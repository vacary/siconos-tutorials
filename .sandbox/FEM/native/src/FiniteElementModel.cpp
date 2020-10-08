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

#include "FiniteElementModel.hpp"
#include "SimpleMatrix.hpp"
#include "Material.hpp"
#include "SiconosAlgebraProd.hpp"
#include "SimpleMatrixFriends.hpp"


// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "debug.h"

#include <stdio.h>
#include <iostream>
#include <map>

unsigned int FiniteElementModel::init()
{
  DEBUG_BEGIN("FiniteElement::init()\n");


  _nodes.resize(_mesh->vertices().size());

  /* Construction of the elements w.r.t type */
  unsigned int numElement=0;
  unsigned int ndof=0;
  unsigned int dofIdx=0;
  unsigned int dim = _mesh->dim();

  unsigned int num_node=0;

  for (MElement * e : _mesh->elements())
  {
    DEBUG_PRINTF("MElement->num() : %zu\n", e->num() );
    switch(e->type()) // we follow gmsh convention
    {
      case 2: // 3-node triangle.
      {
        /* We should normally ask the user for the type of element we associate with 
         * the type of MElement of gmsh */
        
        SP::FElement fe (new FElement(T3, 6, e));
        _elements.push_back(fe); /* the FE element number is equal to the MElement number */
        _mElementTOFElement[e] = _elements.back();
        for (MVertex * v : e->vertices())
        {

          if (_vertexToNode.find(v) == _vertexToNode.end())
          {
            ndof+=2;
            SP::Index dofIndex (new Index(0));
            dofIndex->push_back(dofIdx++);
            dofIndex->push_back(dofIdx++);

            //DEBUG_PRINTF("  create node num_mode: %lu with dofIndex.size():%lu for vertex num = %lu\n", num_node, dofIndex->size(), v->num() );
            _nodes[num_node].reset(new FENode(num_node, v, dofIndex));
            _vertexToNode[v] = _nodes[num_node];
            num_node++;

          }
          else
          {
            //DEBUG_PRINTF("  node already exists for vertex : %zu \n", v->num() );
            int typeElement = v->elements()[0]->type();
            for (MElement * e : v->elements())
            {
              if (e->type() != typeElement)
              {
                RuntimeException::selfThrow("FiniteElementModel::initDOF(). Element type are consistent for all the elements connected to the vertex.");
                return 0;
              }
            }
          }
          //assert(_nodes[num_node]);
          //DEBUG_EXPR(vertexNode.at(v)->display(););
          fe->nodes().push_back(_vertexToNode.at(v));
        }
        break;
      }
      default:
        RuntimeException::selfThrow("FiniteElementModel::initDOF(). element not recognized");
    }
  }
    // for (MVertex * v : e->vertices())
    // {
    //   if (!_nodes[v->num()])
    //   {
    //     ndof+=2;

    //     SP::Index dofIndex (new Index(2));

    //     dofIndex->push_back(dofIdx++);
    //     dofIndex->push_back(dofIdx++);
    //     _nodes[v->num()].reset(new FENode(v, dofIndex));
    //   }
    //   else
    //   {

    //   }
    //   _elements.back->_nodes.push_back(_nodes[v->num()]);

    // }

  // /* we construct from the mesh :
  //  * - the vector of degree of freedom (ddl)
  //  * - the connectivity table between the nodes and the ddl
  //  */

  // unsigned int num=0;
  // unsigned int ndof=0;
  // unsigned int dofIdx;

  // for (MVertex * v : _mesh->vertices())
  // {
  //   int typeElement = v->elements()[0]->type();
  //   DEBUG_PRINTF("typeElement : %i \n", typeElement);
  //   for (MElement * e : v->elements())
  //   {
  //     if (e->type() != typeElement)
  //     {
  //       RuntimeException::selfThrow("FiniteElementModel::initDOF(). Element type are consistent for all the elements connected to the vertices.");
  //       return 0;
  //     }
  //   }
  //   switch(typeElement) // we follow gmsh convention
  //   {
  //   case 2: // 3-node triangle.
  //   {
  //     ndof+=2;

  //     SP::Index dofIndex (new Index(2));

  //     dofIndex->push_back(dofIdx++);
  //     dofIndex->push_back(dofIdx++);

  //     _nodes.push_back(new FENode(num++, v, dofIndex)); /* the FE node number is equal to the MVertex number */
  //     DEBUG_PRINT("_nodes push back \n");
  //     break;
  //   }
  //   default:
  //     RuntimeException::selfThrow("FiniteElementModel::initDOF(). element not recognized");
  //   }
  //}

  return ndof;
  DEBUG_END("FiniteElement::init()\n");

}
void FiniteElementModel::AssembleElementaryMatrix(SP::SiconosMatrix M,
                                                  SimpleMatrix& Me, FElement& fe)
{
  int node1_cnt =0;



  for(SP::FENode  node1 : fe.nodes())
  {
    Index& dofIndex1 = *node1->dofIndex();
    int node2_cnt =0;
    for(SP::FENode  node2 : fe.nodes())
    {
      Index& dofIndex2 = *node2->dofIndex();
      for (int i = 0; i < dofIndex1.size() ; i++)
      {
        for (int j = 0; j < dofIndex2.size() ; j++)
        {
          //DEBUG_PRINTF("i = %i\t j=%i,  Me.getValue(i,j) = %e\n", i+node1_cnt*dofIndex1.size(), j+node2_cnt*dofIndex2.size(), Me.getValue(i,j));

          M->setValue(dofIndex1[i], dofIndex2[j], Me.getValue(i+node1_cnt*dofIndex1.size(),j+node2_cnt*dofIndex2.size())+ M->getValue(dofIndex1[i], dofIndex2[j]));
        }
      }
      node2_cnt++;
    }
    node1_cnt++;
  }
  //y
  //yDEBUG_EXPR(M->display());
}


void FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )
{
  DEBUG_BEGIN("FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");


  Me.zero();

  int ndof = fe.ndof();
  int order = fe.order();
  std::vector<SP::FENode> & nodes= fe.nodes();
  int nnodes= nodes.size();

  int dim = _mesh->dim();
  double * N =(double *)malloc(nnodes *sizeof(double));

  double * Nksi =(double *)malloc(nnodes *sizeof(double)); // only in 2D
  double * Neta =(double *)malloc(nnodes *sizeof(double));

  double * J =(double *)malloc(dim*dim*sizeof(double));

  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder=2;
  for (const double* gp : fe.GaussPoints(integrationOrder))
  {
    if (_mesh->dim() ==2) //Ugly
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta= gp[0];
      double  gp_ksi= gp[1];
      double  gp_w= gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta );
      // Compute element determinant
      for (int i =0; i<4; i++) J[i]=0.0;
      for (int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Neta[n]*nodes[n]->_mVertex->x();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->y();
      }
      double detJ = J[0]*J[3] - J[1]*J[2];
      // DEBUG_PRINTF("detJ = %e\n", detJ );
      // DEBUG_EXPR(std::cout << "Gauss points : "<< gp_eta << " "  << gp_ksi << " "  << gp_w << " "   << std::endl;);

      double coeff = gp_w * massDensity * detJ;

      /* M += (coeff * Nt N)*/
      for (int i = 0; i < nnodes; i++)
      {
        for (int j= 0; j < nnodes; j++)
        {
          // DEBUG_PRINTF(" N[%i] = %e\t N[%i] = %e\t entry = %e \n", i, N[i], j, N[j], coeff* N[i]*N[j]);
          Me.setValue(i*2,j*2, coeff* N[i]*N[j] + Me.getValue(i*2,j*2));
          Me.setValue(i*2+1,j*2+1, coeff* N[i]*N[j] + Me.getValue(i*2+1,j*2+1));
        }
      }
    }
  }
  free(N);
  free(Nksi);
  free(Neta);
  free(J);

  DEBUG_END("FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");
}



void FiniteElementModel::computeMassMatrix(SP::SiconosMatrix M, double massDensity )
{
  DEBUG_BEGIN("FiniteElementModel::computeMassMatrix(SP::SiconosMatrix M, double massDensity )\n");
  M->zero();

  /* loop over the elements */
  for (SP::FElement fe : elements())
  {
    unsigned int ndofElement  = fe->_ndof;
    SP::SimpleMatrix Me(new SimpleMatrix(ndofElement,ndofElement)); // to be optimized if all the element are similar
    computeElementaryMassMatrix(*Me, *fe, massDensity);
    AssembleElementaryMatrix(M, *Me, *fe);
  }
  DEBUG_END("FiniteElementModel::computeMassMatrix(SP::SiconosMatrix M, double massDensity )\n");
}

void FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke, FElement& fe, Material& mat )
{
  DEBUG_BEGIN("FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");


  Ke.zero();
  // Compute element determinant
  int ndof = fe.ndof();
  int order = fe.order();
  std::vector<SP::FENode> & nodes= fe.nodes();
  int nnodes= nodes.size();

  int dim = _mesh->dim();
  double * N =(double *)malloc(nnodes *sizeof(double));

  double * Nksi =(double *)malloc(nnodes *sizeof(double)); // only in 2D
  double * Neta =(double *)malloc(nnodes *sizeof(double));
  double * Nx =(double *)malloc(nnodes*sizeof(double));
  double * Ny =(double *)malloc(nnodes*sizeof(double));



  double * J =(double *)malloc(dim*dim*sizeof(double));
  double * Jinv =(double *)malloc(dim*dim*sizeof(double));




  /** We perform integration by summing over the gauss points
   * this could be simplified by explicit formulae
   */

  int integrationOrder=1;
  for (const double* gp : fe.GaussPoints(integrationOrder))
  {
    if (_mesh->dim() ==2) //Ugly
    {
      // Compute shape function and derivatives of shape function
      double  gp_eta= gp[0];
      double  gp_ksi= gp[1];
      double  gp_w= gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta );
      // Compute element determinant
      for (int i =0; i<4; i++) J[i]=0.0;
      for (int n =0; n < nnodes; n++)
      {
        // DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        // DEBUG_PRINTF(" x = %e\t y = %e\n", nodes[n]->_mVertex->x(), nodes[n]->_mVertex->y());
        J[0] = J[0] + Nksi[n]*nodes[n]->_mVertex->x();
        J[1] = J[1] + Nksi[n]*nodes[n]->_mVertex->y();
        J[2] = J[2] + Neta[n]*nodes[n]->_mVertex->x();
        J[3] = J[3] + Neta[n]*nodes[n]->_mVertex->y();
      }
      double detJ = J[0]*J[3] - J[1]*J[2];
      // DEBUG_PRINTF("detJ = %e\n", detJ );

      // compute inverse of the Jacobian
      Jinv[0] = J[3]/detJ;
      Jinv[1] =-J[1]/detJ;
      Jinv[2] =-J[2]/detJ;
      Jinv[3] = J[0]/detJ;

      // Compute the derivative w.r.t x and y of the shape function
      for (int n =0; n < nnodes; n++)
      {
        DEBUG_PRINTF(" Nksi[%i] = %e\t Neta[%i] = %e\n", n, Nksi[n], n, Neta[n]);
        Nx[n] = Jinv[0] * Nksi[n] + Jinv[1] * Neta[n];
        Ny[n] = Jinv[2] * Nksi[n] + Jinv[3] * Neta[n];
        DEBUG_PRINTF(" Nx[%i] = %e\t Ny[%i] = %e\n", n, Nx[n], n, Ny[n]);
      }

      // Construct the B matrix (its form is consistent with the choice of the reprensentation of strain)
      SP::SimpleMatrix B(new SimpleMatrix(3,ndof));
      B->zero();
      for (int n =0; n < nnodes; n++)
      {
        B->setValue(0, 2*n,   Nx[n]);
        B->setValue(1, 2*n,   0.0);
        B->setValue(2, 2*n,   Ny[n]);
        B->setValue(0, 2*n+1, 0.0);
        B->setValue(1, 2*n+1, Ny[n]);
        B->setValue(2, 2*n+1, Nx[n]);
      }
 
      //Compute 2D elatic tensor
      SP::SimpleMatrix D(new SimpleMatrix(3,3));
      double E = mat.elasticYoungModulus();
      double nu =  mat.poissonCoefficient();
      double coef = E/(1-nu*nu);
      if (mat.analysisType2D() == PLANE_STRAIN)
      {
        (*D)(0,0) = coef*(1.-nu);
        (*D)(0,1) = coef*nu;
        (*D)(0,2) = 0.0;

        (*D)(1,0) = (*D)(0,1);
        (*D)(1,1) = (*D)(0,0);
        (*D)(1,2) = 0.0;

        (*D)(2,0) = 0.0;
        (*D)(2,1) = 0.0;
        (*D)(2,2) = 0.5*coef*(1.0 - 2* nu);
      }
      else
        RuntimeException::selfThrow("FiniteElementModel::computeElementaryStiffnessMatrix. Other type of analysis not yet implemented");
      // Compte BT D B
      SP::SimpleMatrix DB(new SimpleMatrix(3,ndof));
      prod(*D, *B, *DB, true);
      SP::SimpleMatrix BT(new SimpleMatrix(ndof,3));
      BT->trans(*B);
      SP::SimpleMatrix BTDB(new SimpleMatrix(ndof,ndof));
      prod(*BT, *DB, *BTDB, true);

      double coeff =0.0;

      if (mat.analysisType2D() == PLANE_STRAIN || mat.analysisType2D() == PLANE_STRESS)
      {
        coeff = gp_w  * detJ *  mat.thickness();
      }
      else
        RuntimeException::selfThrow("FiniteElementModel::computeElementaryStiffnessMatrix. Other type of analysis not yet implemented");
      Ke += (coeff * *BTDB) ;

    }

  }
  //yDEBUG_EXPR(Ke.display(););
  free(N);
  free(Nksi);
  free(Neta);
  free(J);
  free(Jinv);


  DEBUG_END("FiniteElementModel::computeElementaryStiffnessMatrix(SimpleMatrix& Ke, FElement& fe, Material& mat  )\n");
}



void FiniteElementModel::computeStiffnessMatrix(SP::SiconosMatrix K, Material& mat  )
{
  DEBUG_BEGIN("FiniteElementModel::computeStiffnessMatrix(SP::SiconosMatrix K, Material& mat )\n");
  K->zero();
  /* loop over the elements */
  for (SP::FElement fe : elements())
  {
    unsigned int ndofElement  = fe->_ndof;
    SP::SimpleMatrix Ke(new SimpleMatrix(ndofElement,ndofElement)); // to be optimized if all the element are similar
    computeElementaryStiffnessMatrix(*Ke, *fe, mat);
    AssembleElementaryMatrix(K, *Ke, *fe);
  }
  //DEBUG_EXPR(K->display(););
  DEBUG_END("FiniteElementModel::computeStiffnessMatrix(SP::SiconosMatrix K, Material& mat )\n");
}

void FiniteElementModel::display(bool brief) const
{
  std::cout << "===== FiniteElementModel display ===== " <<std::endl;
  std::cout << "- numberOfNodes : " << _nodes.size() << std::endl;
  std::cout << "- numberOfElements : " << _elements.size() << std::endl;
  for (SP::FElement f : _elements)
  {
    f->display();
  }
}
