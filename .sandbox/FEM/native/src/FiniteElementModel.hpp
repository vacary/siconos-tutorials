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

#include "SiconosPointers.hpp"                //for SP::
#include "SiconosAlgebraTypeDef.hpp"          //for Index

#include "FemFwd.hpp"
#include "Mesh.hpp"

// a Finite Element node
struct FENode
{
  /* node number */
  size_t _num;

  /* Element type */

  /* associated Mvertex */
  MVertex * _mVertex;

  /* associated dof number  in the global dof vector*/
  SP::Index _dofIndex;

  /* */

  FENode(size_t num, MVertex *v, SP::Index dofIndex ):_num(num), _mVertex(v), _dofIndex(dofIndex){};

  SP::Index dofIndex(){return _dofIndex;};
  
  void display()
  {
    std::cout << " - Fe Node - number: " << _num ;
    std::cout << std::endl;
  };
};

enum FINITE_ELEMENT_TYPE
{
  T3,
  T6,
  Q4
};


static const double GP_T3_x[] = {0.333333333333333333};
static const double GP_T3_y[] = {0.333333333333333333};
static const double GP_T3_w[] = {0.5};
static int  GP_T3_order = 1;
static  std::vector<const double* > GaussPointsT3 = {GP_T3_x, GP_T3_y, GP_T3_w};

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
  std::vector<FENode *> _nodes;

  /* associated Mesh element */
  MElement * _mElement;

  FElement(FINITE_ELEMENT_TYPE type, unsigned int ndof, MElement *e):
    _num(e->num()), _type(type), _ndof(ndof), _mElement(e){};

  unsigned int ndof()
  {
    return _ndof;
  }

  int GaussIntegrationOrder()
  {
    switch(_type)
    {
    case T3:
      return GP_T3_order;
      
      break;
    default:
      RuntimeException::selfThrow("FElement::GaussPoints(). element type not recognized");
    }
    return 0;
  }
  
  std::vector<const double*> & GaussPoints()
  {
    switch(_type)
    {
    case T3:
      return  GaussPointsT3;
      break;
    default:
      RuntimeException::selfThrow("FElement::GaussPoints(). element type not recognized");
    }
    return  GaussPointsEmpty;
  }
  
  void display()
  {
    std::cout << " - FElement - number: " << _num
              << "            - type: " << _type
              << "            - ndof: " << _ndof
              << "            - number of nodes" << _nodes.size() ;
    std::cout << std::endl;
  };
};


// a finite element model
class FiniteElementModel
{
protected :

  /** a mesh */
  SP::Mesh _mesh;

  /** nodes */
  std::vector<FENode *> _nodes;

  /** elements */
  std::vector<FElement *> _elements;


  /** default constructor */
  FiniteElementModel() {};

public:
  FiniteElementModel(SP::Mesh mesh):
    _mesh(mesh){};

  std::vector<FElement *> & elements()
  {
    return _elements;
  }
  
  /* create the FEM model from the mesh and the element type
   * \return the number of dof */
  unsigned int init();

  void display(bool brief) const;
};



#endif // FINITEELEMENTMODEL_H
