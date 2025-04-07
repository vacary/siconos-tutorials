/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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

#include <SiconosFEM.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <cmath>  // for fabs

#include "FENode.hpp"
#include "FiniteElementModel.hpp"
#include "MeshUtils.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"

using namespace std;
using namespace siconos::mechanics::fem;
using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;
    //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
    //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
    // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
    // string gmsh_filename = "./mesh_data/triangle_reference.msh";
    // string gmsh_filename = "./mesh_data/square_6.msh";
    auto gmsh_filename = "./mesh_data/square_200.msh";
    // string gmsh_filename = "./mesh_data/square_2720.msh";

    auto v0 = std::make_shared<MVertex>(0, 0.0, 0.0, 0.0);
    auto vEnd = std::make_shared<MVertex>(0, 1.0, 0.0, 0.0);
    int nbBeams = 4;
    int dim = 2;
    auto mesh = createBeamMesh(v0, vEnd, nbBeams, 2);
    // mesh->display(false);

    siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

            // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
            // 1/5.); // material for  triangle_felippa.msh
    double density = 7800.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
        density, 210e9, 1 / 3., 1.0);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>>
        materials = {{bulk_material_tag, mat1}};

    auto FEsolid = std::make_shared<siconos::mechanics::fem::FiniteElementLinearTIDS>(
            mesh, materials, siconos::algebra::UblasType::SPARSE);


    return 0;


  } catch (...) {
    std::cerr << "Exception caught in T3.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
