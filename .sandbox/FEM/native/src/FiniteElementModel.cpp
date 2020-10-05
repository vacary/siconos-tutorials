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



#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "debug.h"

#include <stdio.h>
#include <iostream>


unsigned int FiniteElementModel::init()
{
  DEBUG_BEGIN("FiniteElement::init()\n");


  /* Construction of the elements w.r.t type */
  unsigned int numElement=0;
  unsigned int dim = _mesh->dim();
  for (MElement * e : _mesh->elements())
  {
    switch(e->type()) // we follow gmsh convention
    {
      case 2: // 3-node triangle.
        _elements.push_back(new FElement(T3, 6, e)); /* the FE element number is equal to the MElement number */
      break;
      default:
        RuntimeException::selfThrow("FiniteElementModel::initDOF(). element not recognized");
    }
  }

  /* we construct from the mesh :
   * - the vector of degree of freedom (ddl)
   * - the connectivity table between the nodes and the ddl
   */

  unsigned int num=0;
  unsigned int ndof=0;
  unsigned int dofIdx;
  
  for (MVertex * v : _mesh->vertices())
  {
    int typeElement = v->elements()[0]->type();
    DEBUG_PRINTF("typeElement : %i \n", typeElement);
    for (MElement * e : v->elements())
    {
      if (e->type() != typeElement)
      {
        RuntimeException::selfThrow("FiniteElementModel::initDOF(). Element type are consistent for all the elements connected to the vertices.");
        return 0;
      }
    }
    switch(typeElement) // we follow gmsh convention
    {
    case 2: // 3-node triangle.
    {
      ndof+=2;

      SP::Index dofIndex (new Index(2));

      dofIndex->push_back(dofIdx++);
      dofIndex->push_back(dofIdx++);

      _nodes.push_back(new FENode(num++, v, dofIndex)); /* the FE node number is equal to the MVertex number */
      DEBUG_PRINT("_nodes push back \n");
      break;
    }
    default:
      RuntimeException::selfThrow("FiniteElementModel::initDOF(). element not recognized");
    }
  }
  
  return ndof;
  DEBUG_END("FiniteElement::init()\n");

}

void FiniteElementModel::display(bool brief) const
{
  std::cout << "===== FiniteElementModel display ===== " <<std::endl;
  std::cout << "- numberOfElements : " << _elements.size() << std::endl;
  for (FElement * f : _elements)
  {
    f->display();
  }
  std::cout << "- numberOfNodes : " << _nodes.size() << std::endl;

}
