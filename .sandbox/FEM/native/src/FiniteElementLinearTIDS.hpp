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
#ifndef FINITEELEMENTLAGRANGIANTIDS_H
#define FINITEELEMENTLAGRANGIANTIDS_H

#include "LagrangianLinearTIDS.hpp"
#include "Mesh.hpp"
#include "FiniteElementModel.hpp"
//#include "SiconosMatrix.hpp"       // UBLAS_TYPE


#include "FemFwd.hpp"

/** Finite Element discretization of elastic solids that inherits from Lagrangian Linear Systems with time invariant coefficients
 * - \f$M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$
 */

class FiniteElementLinearTIDS : public LagrangianLinearTIDS
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** a mesh */
  SP::Mesh _mesh;

  /** a material */
  SP::Material _material;

  /** a finite element model */
  SP::FiniteElementModel _FEModel;

  /* Storage type for the matrices */
  Siconos::UBLAS_TYPE _storageType;
  
  /** default constructor */
  FiniteElementLinearTIDS():LagrangianLinearTIDS() {};

  
public:

  /** constructor from initial state and all matrix operators.
   * \param mesh the mesh that defined the spatial discretization
   * \param material
   */
  FiniteElementLinearTIDS(SP::Mesh mesh, SP::Material material,  Siconos::UBLAS_TYPE storageType=Siconos::DENSE);

  /** destructor */
  ~FiniteElementLinearTIDS(){};

  void setStorageType(Siconos::UBLAS_TYPE type)
  {
    _storageType=type;
  }

  
  SP::FiniteElementModel FEModel()
  {
    return _FEModel;
  };

  
  
  
  void display(bool brief) const;



  ACCEPT_STD_VISITORS();

};
#endif // FINITEELEMENTLAGRANGIANTIDS_H
