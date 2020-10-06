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


#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "debug.h"

#include <stdio.h>
#include <iostream>


unsigned int FiniteElementModel::init()
{
  DEBUG_BEGIN("FiniteElement::init()\n");


  _nodes.resize(_mesh->vertices().size());

  /* Construction of the elements w.r.t type */
  unsigned int numElement=0;
  unsigned int ndof=0;
  unsigned int dofIdx=0;
  unsigned int dim = _mesh->dim();
  for (MElement * e : _mesh->elements())
  {
    DEBUG_PRINTF("MElement->num() : %zu\n", e->num() );
    switch(e->type()) // we follow gmsh convention
    {
      case 2: // 3-node triangle.
      {
        SP::FElement fe (new FElement(T3, 6, e));
        _elements.push_back(fe); /* the FE element number is equal to the MElement number */
        for (MVertex * v : e->vertices())
        {
          if (!_nodes[v->num()-1])
          {
            ndof+=2;
            SP::Index dofIndex (new Index(0));
            dofIndex->push_back(dofIdx++);
            dofIndex->push_back(dofIdx++);

            DEBUG_PRINTF("create node num: %lu with dofIndex.size():%lu \n", v->num()-1, dofIndex->size() );
            _nodes[v->num()-1].reset(new FENode(v, dofIndex));
          }
          else
          {
            DEBUG_PRINTF(" node already exists v->num(): %zu \n", v->num() );
            int typeElement = v->elements()[0]->type();
            for (MElement * e : v->elements())
            {
              if (e->type() != typeElement)
              {
                RuntimeException::selfThrow("FiniteElementModel::initDOF(). Element type are consistent for all the elements connected to the vertices.");
                return 0;
              }
            }
          }
          assert(_nodes[v->num()-1]);
          fe->nodes().push_back(_nodes[v->num()-1]);
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


void FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )
{
  DEBUG_BEGIN("FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");


  Me.zero();
  // Compute element determinant


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
      double  gp_eta= gp[0];
      double  gp_ksi= gp[1];
      double  gp_w= gp[2];
      fe.shapeFunctionIso2D(gp_eta, gp_ksi, N, Nksi, Neta );
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

      /* M += coeff * Nt N*/
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
  DEBUG_EXPR(Me.display(););
  DEBUG_END("FiniteElementModel::computeElementaryMassMatrix(SimpleMatrix& Me, FElement& fe, double massDensity )\n");
}

void FiniteElementModel::AssembleElementaryMatrix(SP::SiconosMatrix M,
                                                  SimpleMatrix& Me, FElement& fe)
{

  for(SP::FENode  node1 : fe.nodes())
  {
    node1->display();
    Index& dofIndex1 = *node1->dofIndex();
    for(SP::FENode  node2 : fe.nodes())
    {
      Index& dofIndex2 = *node2->dofIndex();
      for (int i = 0; i < dofIndex1.size() ; i++)
      {
        for (int j = 0; j < dofIndex2.size() ; j++)
        {
          M->setValue(dofIndex1[i], dofIndex2[j], Me.getValue(i,j)+ M->getValue(dofIndex1[i], dofIndex2[j]));
        }
      }
    }
  }
  DEBUG_EXPR(M->display());

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


