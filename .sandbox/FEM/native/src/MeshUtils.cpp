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
#include <math.h>

#include <fstream>
using namespace std;
#include <stdio.h>
#include <iostream>
#include <sstream>

#include <string>
#include <algorithm>
#include <iterator>

#define DEBUG_STDOUT
#define DEBUG_NOCOLOR
#define DEBUG_MESSAGES
#include "debug.h"


template <class Container>
void split(const std::string& str, Container& cont,
              const std::string& delims = " ")
{
    std::size_t current, previous = 0;
    current = str.find_first_of(delims);
    while (current != std::string::npos) {
        cont.push_back(str.substr(previous, current - previous));
        previous = current + 1;
        current = str.find_first_of(delims, previous);
    }
    cont.push_back(str.substr(previous, current - previous));
}

Mesh* createMeshFromGMSH2(std::string gmsh_filename)
{

  ifstream in(gmsh_filename);

  if(!in) {
    cout << "Cannot open input file.\n";
    return NULL;
  }

  char str[500];



  float gmsh_version;
  int m =3;
  unsigned int number_of_vertices;
  std::vector<MVertex *> vertices;
  unsigned int number_of_elements;
  std::vector<MElement *> elements;
  double min_z= 1e+64, max_z = -1e+64;
  string toto;

  while(in) {
    std::string line;
    std::getline(in, line);
    //if(in) cout << line << endl;
    //std::cout << "line.front()" << line.front() << std::endl;

    // std::vector<std::string> words;
    // split4(line, words);
    // std::copy(words.begin(), words.end(),
    //           std::ostream_iterator<std::string>(std::cout, "\n"));

    if (line.compare("$MeshFormat") == 0 )
    {
      std::getline(in, line);
      std::vector<std::string> words;
      split(line, words);
      // std::copy(words.begin(), words.end(),
      //           std::ostream_iterator<std::string>(std::cout, "\n"));

      stringstream token(words[0]);
      token >> gmsh_version;
      //std::cout << "gmsh_version :" << gmsh_version << endl;
    }

    if (line.compare("$Nodes") == 0 )
    {
      std::getline(in, line);
      stringstream token(line);
      token >> number_of_vertices;
      std::cout << "number_of_vertices : " << number_of_vertices <<  endl;
      while (std::getline(in, line))
      {
        if (line.compare("$EndNodes") == 0 ) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        stringstream token(words[0]);
        int vertex_number ;
        token >> vertex_number;
        double x, y, z ;
        stringstream t_x(words[1]); t_x >> x;
        stringstream t_y(words[2]); t_y >> y;
        stringstream t_z(words[3]); t_z >> z;
        max_z = std::max(max_z,z);
        min_z = std::min(min_z,z);
        vertices.push_back(new MVertex(vertex_number, x, y, z));
        //vertices.back()->display();
      }
    }

    if (fabs(min_z-max_z) < 1e-16)
    {
      m =2;
    }
    if (line.compare("$Elements") == 0 )
    {
      std::getline(in, line);
      stringstream token(line);
      token >> number_of_elements;
      std::cout << "number_of_elements : " << number_of_elements <<  endl;
      while (std::getline(in, line))
      {
        if (line.compare("$EndElements") == 0 ) break;
        std::vector<std::string> words;
        split(line, words);
        // std::copy(words.begin(), words.end(),
        //           std::ostream_iterator<std::string>(std::cout, "\n"));
        stringstream token(words[0]);
        int element_number ;
        token >> element_number;
        stringstream t_type(words[1]);
        int element_type ;
        t_type >> element_type;
        stringstream t_nt(words[2]);
        int number_of_tags ;
        t_nt >> number_of_tags;
        std::vector<int> tags ;
        for (int k =0 ; k < number_of_tags; k++)
        {
          int tag;
          stringstream token(words[k]);
          token >> tag;
          tags.push_back(tag);
        }
        
        std::vector<MVertex *> vertices_e ;
        for (int k = 3+number_of_tags; k < words.size(); k++)
        {
          int node_number;
          stringstream token(words[k]);
          token >> node_number;
          //        std::cout << "node_number" << node_number << std::endl;
          int v;
          for (int k  = 0 ; k <  vertices.size(); k++)
          {
            if (node_number-1+k < vertices.size())
              v  = node_number-1+k;
            else
            {
              v = node_number-1+k -vertices.size();
            }

            if (node_number == vertices[v]->num())
            {
              vertices_e.push_back(vertices[v]); 
              break;
            }
          }
        }
        elements.push_back(new MElement(element_number, element_type, vertices_e));
        //elements.back()->display();
      }
    }




  }

  in.close();

  
  return new Mesh(m, vertices, elements);;
}
