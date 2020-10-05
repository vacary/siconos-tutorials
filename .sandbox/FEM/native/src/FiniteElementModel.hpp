/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*! \file FiniteElementLinearTIDS.hpp

 */
#ifndef FINITEELEMENTMODEL_H
#define FINITEELEMENTMODEL_H
#include <vector>
#include <iostream>
#include <map>

#include "SiconosPointers.hpp"                //for SP::
#include "SiconosAlgebraTypeDef.hpp"          //for Index

#include "FemFwd.hpp"
#include "Mesh.hpp"

// a Finite Element node
struct FENode
{
  /* node number */
  size_t _num;

  /* associated Mvertex */
  MVertex * _mVertex;

  /* associated dof number  in the global dof vector*/
  SP::Index _dofIndex;

  /* */

  FENode(size_t num, MVertex *v, SP::Index dofIndex ): _num(num), _mVertex(v), _dofIndex(dofIndex){};


  size_t num()
  {
    return _num;
  }


  
  SP::Index dofIndex(){return _dofIndex;};


  double x()
  {
    return _mVertex->x();
  }
  
  double y()
  {
    return _mVertex->y();
  }
  
  double z()
  {
    return _mVertex->z();
  }
  
  void display()
  {
    std::cout << "     - Fe Node - number: " << _num
              << "               - ndof:" << _dofIndex->size()
              << "               - dofIndex: "<< _dofIndex->front() << ":" << _dofIndex->back()  ;
    std::cout << std::endl;
  };
};

enum FINITE_ELEMENT_TYPE
{
  T3=2,
  T6,
  Q4
};



static const double GP_T3_1_p1[]= {0.333333333333333333, 0.333333333333333333, 0.5 };
static  std::vector<const double* > GaussPointsT3_1 = {GP_T3_1_p1};

static const double GP_T3_2_p1[]= {0.66666666666666667, 0.16666666666666667, 0.1666666666666666 };
static const double GP_T3_2_p2[]= {0.16666666666666667, 0.66666666666666667, 0.1666666666666666 };
static const double GP_T3_2_p3[]= {0.16666666666666667, 0.16666666666666667, 0.1666666666666666 };
static  std::vector<const double* > GaussPointsT3_2 = {GP_T3_2_p1, GP_T3_2_p2, GP_T3_2_p3};





static  std::vector<const double* > GaussPointsEmpty = {};


// a Finite Element
struct FElement
{
  /* element number */
  size_t _num;

  /* Element type */
  FINITE_ELEMENT_TYPE _type;

  /* number of dof by Element  */
  unsigned _ndof;

  /** nodes */
  std::vector<SP::FENode> _nodes;

  /* associated Mesh element */
  MElement * _mElement;

  FElement(FINITE_ELEMENT_TYPE type, unsigned int ndof, MElement *e):
    _num(e->num()), _type(type), _ndof(ndof), _mElement(e){};

  unsigned int ndof()
  {
    return _ndof;
  }

  

  
  int order()
  {
    switch(_type)
    {
    case T3:
      return 1;

      break;
    default:
      RuntimeException::selfThrow("FElement::GaussPoints(). element type not recognized");
    }
    return 0;
  }

  std::vector<const double*> & GaussPoints(int order)
  {
    switch(_type)
    {
    case T3:
      if (order==1)
        return  GaussPointsT3_1;
      else if (order==2)
        return  GaussPointsT3_2;
      break;
    default:
      RuntimeException::selfThrow("FElement::GaussPoints(). element type not recognized");
    }
    return  GaussPointsEmpty;
  }

  std::vector<SP::FENode> & nodes()
  {
    return _nodes;
  }

  void shapeFunctionIso2D(double ksi, double eta, double* N, double * Nksi, double * Neta )
  {
    switch(_type)
    {
    case T3:
    {
      N[0] = 1.0 - ksi - eta;
      N[1] = ksi;
      N[2] = eta;
      Nksi[0] = -1.0;
      Nksi[1] = 1.0;
      Nksi[2] = 0.0;
      Neta[0] = -1.0;
      Neta[1] = 0.0;
      Neta[2] = 1.0;
      break;
    }
    default:
      RuntimeException::selfThrow("FElement::shapeFunctionIso2D(). element type not recognized");
    }
  }

  void display()
  {
    std::cout << " - FElement - number: " << _num
              << "            - type: " << _type
              << "            - ndof: " << _ndof
              << "            - number of nodes: " << _nodes.size() ;
    std::cout << std::endl;
    for (SP::FENode n : _nodes)
    {
      n->display();
    }
  };
};


// a finite element model
class FiniteElementModel
{
protected :

  /** a mesh */
  SP::Mesh _mesh;

  /** nodes */
  std::vector<SP::FENode> _nodes;

  /** elements */
  std::vector<SP::FElement> _elements;

  /** vertex to node map **/
  std::map<MVertex * , SP::FENode > _vertexToNode;

  /** MElement to FElement map **/
  std::map<MElement * , SP::FElement > _mElementTOFElement;

  /** default constructor */
  FiniteElementModel() {};

public:
  FiniteElementModel(SP::Mesh mesh):
    _mesh(mesh){};

  std::vector<SP::FElement > &  elements()
  {
    return _elements;
  }

  std::vector<SP::FENode > &  nodes()
  {
    return _nodes;
  }
  SP::FENode vertexToNode(MVertex * v)
  {
    return _vertexToNode.at(v);
  }
  /* create the FEM model from the mesh and the element type
   * \return the number of dof */
  unsigned int init();

  /* Assembly method for elemetary matrix */
  void AssembleElementaryMatrix(SP::SiconosMatrix M,
                                SimpleMatrix& Me, FElement& fe);

  /** compute Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeMassMatrix(SP::SiconosMatrix, double massDensity);

  /** compute elementary Mass Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe,  double massDensity);

  /** compute Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeStiffnessMatrix(SP::SiconosMatrix, Material& mat );

  /** compute elementary Stiffness Matrix
   * should be computeMass of LagrangianDS ?
   **/
  void computeElementaryStiffnessMatrix(SimpleMatrix& Me, FElement& fe, Material& mat);




  void display(bool brief) const;
};



#endif // FINITEELEMENTMODEL_H
