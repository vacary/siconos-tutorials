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

/*! \file Mesh.hpp

 */
#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <iostream>


/** A simple class for material
 */

class Material
{

protected:
  /* serialization hooks */
  //ACCEPT_SERIALIZATION(Mesh);

  /** mass density */
  double _massDensity;

  /** Young Modulus */
  double _ElasticYoungModulus;

  /** Poison coefficient */
  double _poissonCoefficient;

  /** default constructor */
  Material() {};

public:

  /** constructor
   */
  Material(double massDensity, double ElasticYoungModulus,  double poissonCoefficient):
    _massDensity(massDensity),
    _ElasticYoungModulus(ElasticYoungModulus),
    _poissonCoefficient(poissonCoefficient){};


  /** destructor */
  ~Material(){};

  double massDensity()
  {
    return _massDensity;
  }

  
  /** print the data of the Material
   */
  void display(bool brief = true) const;

  //
  //ACCEPT_STD_VISITORS();

};
#endif // MESH_H
