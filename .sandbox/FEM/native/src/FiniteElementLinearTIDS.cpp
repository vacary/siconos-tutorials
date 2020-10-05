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

#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"



#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "debug.h"

#include <stdio.h>
#include <iostream>

FiniteElementLinearTIDS::FiniteElementLinearTIDS(SP::Mesh mesh, SP::Material material):
_mesh(mesh), _material(material)
{

  _FEModel.reset(new FiniteElementModel(mesh));
  _ndof = _FEModel->init();


  _q0.reset(new SiconosVector(_ndof,0.0));
  _velocity0.reset(new SiconosVector(_ndof,0.0));
  LagrangianDS::_init(_q0,_velocity0);

  computeMassMatrix();
  
}
void FiniteElementLinearTIDS::computeElementaryMassMatrix(SP::SimpleMatrix Me, FElement * fe)
{

  // Compute element determinant
  double  element_det=1.0;
  
  // perform the integration

  std::vector<const double*> gp = fe->GaussPoints();
  
  for (int order =0; order < fe->GaussIntegrationOrder(); order++)
  {
    if (_mesh->dim() ==2) //Ugly
    {
      double  gp_x= gp[0][order];
      double  gp_y= gp[1][order];
      double  gp_w= gp[2][order];
      std::cout << "Gauss points : "<< gp_x << " "  << gp_y << " "  << gp_w << " "   << std::endl;
      //SimpleMatrix N = fe->shapeFunction(gp_x,gp_y);
      double coeff = gp_w * _material->massDensity() * element_det;
      // *Me +=  coef * (N.transpose() * N) ;
    }
  }
}

void FiniteElementLinearTIDS::AssembleElementaryMatrix(SP::SiconosMatrix M,
                                                       SP::SimpleMatrix Me, FElement * fe)
{

}

void FiniteElementLinearTIDS::computeMassMatrix()
{
  
  for (FElement * fe : _FEModel->elements())
  {
    unsigned int ndofElement  = fe->_ndof;
    SP::SimpleMatrix Me(new SimpleMatrix(ndofElement,ndofElement)); // to be optimized
    computeElementaryMassMatrix(Me, fe);
    AssembleElementaryMatrix(_mass, Me, fe);
  }
  
}



void FiniteElementLinearTIDS::display(bool brief) const
{
  std::cout << "===== FiniteElementLinearTIDS display ===== " <<std::endl;
  std::cout << "- ndof : " << _ndof << std::endl;
  _FEModel->display(brief);
}
