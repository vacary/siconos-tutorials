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

#include "MeshUtils.hpp"
#include "Mesh.hpp"
#include "./src/gmsh_io.hpp"
#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "debug.h"

#include <stdio.h>
#include <iostream>


Mesh* createMeshFromGMSH(std::string gmsh_filename)
{

  int *element_node;
  int element_num;
  int element_order;

  int m;
  int node_num;
  double *node_x;

  cout << "\n";
  cout << "  Read data from a file.\n";
//
//  Get the data size.
//
  gmsh_size_read ( gmsh_filename, node_num, m, element_num,
    element_order );
//
//  Print the sizes.
//
  cout << "\n";
  cout << "  Node data read from file \"" << gmsh_filename << "\"\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "  Spatial dimension = " << m << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "  Element order = " << element_order << "\n";

//
//  Allocate memory.
//
  node_x = ( double * ) malloc ( m * node_num * sizeof ( double ) );
  element_node = ( int * )
    malloc ( element_order * element_num * sizeof ( int ) );
  std::vector<MVertex *> vertices;
  //vertices.resize(node_num);
  std::vector<MElement *> elements;
  //elements.resize(elements);
//
//  Get the data.
//
  gmsh_data_read ( gmsh_filename, m, node_num, node_x,
    element_order, element_num, element_node );

  for(int v = 0; v <node_num; v++)
  {
    if (m==2)
    {
      vertices.push_back(new MVertex(v, node_x[0+v*m] , node_x[1+v*m], 0.));
    }
    else if (m==3)
    {
      vertices.push_back(new MVertex(v, node_x[0+v*m] , node_x[1+v*m], node_x[2+v*m]));
    }
  }
  unsigned int element_cnt=0;
  for(int e = 0; e <element_num; e++)
  {
    std::vector<MVertex *> vertices_e ;
    for ( int i = 0; i < element_order; i++ )
    {
      int num_v = element_node[i+e *element_order] - 1;
      vertices_e.push_back(vertices[num_v]);
    }

    elements.push_back(new MElement(element_cnt++, element_order, vertices_e));
  }



//
//  Print some of the data.
//
  r8mat_transpose_print_some ( m, node_num, node_x,
    1, 1, m, 1000, "  Coordinates for first 10 nodes:" );

  i4mat_transpose_print_some ( element_order, element_num, element_node,
    1, 1, element_order, 100, "  Connectivity for first 10 elements:" );
  return new Mesh(m, vertices, elements);

}
