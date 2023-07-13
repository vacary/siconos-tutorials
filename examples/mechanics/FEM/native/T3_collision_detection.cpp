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

using Matrix = siconos::algebra::SimpleMatrix;
using Vector = siconos::algebra::SiconosVector;

class MyContactDetection : public siconos::simulation::InteractionManager {
 protected:
  double _initial_gap{0.};

  std::shared_ptr<Vector> _normal{nullptr};

  std::shared_ptr<Vector> _tangent{nullptr};

  std::shared_ptr<siconos::modeling::NonSmoothLaw> _nslaw{nullptr};
  unsigned int _contact_frame_dimension = 1;

  std::shared_ptr<siconos::mechanics::fem::FiniteElementModel> _femodel{
      nullptr};

  std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> _fesolid{
      nullptr};

  struct contacting_node {
    std::shared_ptr<siconos::mechanics::fem::FENode> _node{nullptr};
    std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
    unsigned int _node_index{0};

    contacting_node(std::shared_ptr<siconos::mechanics::fem::FENode> node,
                    unsigned int node_index)
        : _node(node), _node_index(node_index) {}
  };

  std::list<contacting_node*> _contacting_nodes;

 public:
  MyContactDetection(
      double initial_gap, std::shared_ptr<Vector> normal,
      std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw,
      std::shared_ptr << siconos::mechanics::fem::FiniteElementLinearTIDS >
          fesolid)
      : InteractionManager(),
        _initial_gap(initial_gap),
        _nslaw(nslaw),
        _normal(normal),
        _fesolid(fesolid)

  {
    _femodel = _fesolid->FEModel();

    auto nslaw_with_friction = std::dynamic_pointer_cast < siconos::modeling
        : NewtonImpactFrictionNSL > (nslaw);
    if (nslaw_with_friction) _contact_frame_dimension = 2;

    if (_contact_frame_dimension == 2)  // create tangent vertor
    {
      _tangent = std::make_shared<Vector>(2);
      _tangent->setValue(0, -(*_normal)(1));
      _tangent->setValue(1, (*_normal)(0));
    }
  }
  virtual ~MyContactDetection() noexcept = default;

  // shoud be done with find and a lambda function
  contacting_node* find_contacting_node(unsigned int node_index) {
    for (contacting_node* cnn : _contacting_nodes) {
      if (cnn->_node_index == node_index) {
        // std::cout << "existing interaction" << std::endl;
        return cnn;
      }
    }
    return nullptr;
  }

  /** Called by Simulation after updating positions prior to starting
   * the Newton loop. */
  void updateInteractions(
      std::shared_ptr<siconos::simulation::Simulation> simulation) {
    // std::cout<< "\nCall to updateInteractions in MyContactDetection" <<
    // std::endl;

    Vector& displacement = *(_fesolid->q());
    auto solid = simulation->nonSmoothDynamicalSystem();

    // update the list of contacting node by brute contact detection
    //_contacting_nodes.clear();
    for (auto& n : _femodel->nodes()) {
      // std::cout << "n->num() is: " << n->num() << std::endl;
      if (fabs(n->y()) <= 1e-16 and fabs(n->x()) >= 1e-16) {
        // std::cout << "contact node number : " << n->num() << " " << n->y() <<
        // std::endl;
        contacting_node* cn = new contacting_node(n, (*n->dofIndex())[0]);
        if (!find_contacting_node(cn->_node_index)) {
          // std::cout << "add contact node number : " << n->num() << " is
          // contacting list " <<  std::endl;
          _contacting_nodes.push_back(cn);
        }
      } else {
        contacting_node* cn = find_contacting_node((*n->dofIndex())[0]);
        if (cn) {
          // std::cout << "remove contact node number : " << n->num() << " is
          // contacting list " <<  std::endl;
          _contacting_nodes.remove(cn);
        }
      }
    }

    // update/create Interaction for contacting points
    for (contacting_node* cn : _contacting_nodes) {
      auto n = cn->_node;
      unsigned int node_idx = (*n->dofIndex())[0];
      if (cn->_inter)  // update interaction->relation
      {
        // std::cout << "update interaction" << std::endl;

        std::shared_ptr<Vector> pc2;
        if (_contact_frame_dimension == 2) {
          auto r =
              std::static_pointer_cast<siconos::mechanics::fem::NodeFem2d2DR>(
                  cn->_inter->relation());
          pc2 = r->pc2();
        } else {
          auto r =
              std::static_pointer_cast<siconos::mechanics::fem::NodeFem1d2DR>(
                  cn->_inter->relation());
          pc2 = r->pc2();
        }

        pc2->setValue(0, cn->_node->x() + displacement(node_idx));
        pc2->setValue(1, -_initial_gap);
      } else  // create an interaction and link
      {
        // std::cout << "create interaction" << std::endl;
        auto pc2 = std::make_shared<Vector>(2);
        pc2->setValue(0, cn->_node->x() + displacement(node_idx));
        pc2->setValue(1, -_initial_gap);
        std::shared_ptr<Relation> relation;
        if (_contact_frame_dimension == 2) {
          relation = std::make_shared<siconos::mechanics::fem::NodeFem2d2DR>(
              cn->_node, pc2, _normal, _tangent);
        } else {
          relation = std::make_shared<siconos::mechanics::fem::NodeFem1d2DR>(
              cn->_node, pc2, _normal);
        }
        auto inter =
            std::make_shared<siconos::modeling::Interaction>(_nslaw, relation);
        cn->_inter = inter;
        // link the interaction and the dynamical system
        solid->link(inter, _fesolid);
      }
    }
    // std::cout << "end to updateInteractions in MyContactDetection" <<
    // std::endl;
  }
};

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;
    //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
    //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
    Ly = 1.0;
    // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
    // string gmsh_filename = "./mesh_data/triangle_reference.msh";
    auto gmsh_filename = "./mesh_data/square_200.msh";
    // string gmsh_filename = "./mesh_data/square_6.msh";
    // string gmsh_filename = "./mesh_data/square_2720.msh";

    auto mesh = siconos::mechanics::fem::createMeshFromGMSH2(gmsh_filename);
    // mesh->display(false);

    siconos::mechanics::fem::writeMeshforPython(mesh);

    int bulk_material_tag = 1;
    int boundary_condition_tag = 2;
    int applied_force_tag = 3;

    // std::shared_ptr<Material> mat1 = std::make_shared<Material>(1, 8*36/5.,
    // 1/5.); // material for  triangle_felippa.msh
    double density = 7800.;
    auto mat1 = std::make_shared<siconos::mechanics::fem::Material>(
        density, 210e9, 1 / 3.);
    std::map<unsigned int, std::shared_ptr<siconos::mechanics::fem::Material>>
        materials = {{bulk_material_tag, mat1}};

    auto start = std::chrono::system_clock::now();
    auto FEsolid =
        std::make_shared<siconos::mechanics::fem::FiniteElementLinearTIDS>(
            mesh, materials, siconos::algebra::UblasType::SPARSE);
    auto end = std::chrono::system_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << "Assembly time : " << elapsed << " ms\n";

    auto femodel = FEsolid->FEModel();
    // FEsolid->K()->display();

    /*------------------------------------------------- Applied forces  */

    auto nodal_forces = std::make_shared<Vector>(2);
    (*nodal_forces)(1) = -1e7;
    FEsolid->applyNodalForces(applied_force_tag, nodal_forces);

    /*------------------------------------------------- Boundary Conditions  */
    /* This part should be hidden in a new BC function for a node number
     * and a dof index. */

    auto node_dof_index = std::make_shared<std::vector<int>>(0);
    node_dof_index->push_back(0);
    node_dof_index->push_back(1);

    FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag,
                                              node_dof_index);
    FEsolid->boundaryConditions()->display();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;       // initial computation time
    double T = 1e-02;    // final computation time
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    auto solid =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e = 0.0;
#define WITH_FRICTION
#ifdef WITH_FRICTION
    double mu = 1.0;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(
        e, 0.0, mu, 2);
#else
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
#endif

    double initial_gap = Ly * 5e-4;
    auto displacement = FEsolid->q();
    auto normal = std::make_shared<Vector>(2);
    normal->setValue(0, 0.0);
    normal->setValue(1, 1.0);
    auto collision_detection = std::make_shared<MyContactDetection>(
        initial_gap, normal, nslaw, FEsolid);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
#ifdef WITH_FRICTION
    auto osnspb =
        std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);
#else
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::LCP>();
#endif
    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(solid, t, OSI,
                                                                 osnspb);

    s->insertInteractionManager(collision_detection);

    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    Matrix dataPlot(N + 1, outputSize);

    auto q = FEsolid->q();
    auto v = FEsolid->velocity();
    auto p = FEsolid->p(1);
    // auto lambda = inter->lambda(1);

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(FEsolid->dimension() - 1);
    dataPlot(0, 2) = (*v)(FEsolid->dimension() - 1);
    dataPlot(0, 3) = (*p)(0);
    // dataPlot(0, 4) = (*lambda)(0);

    auto filename =
        siconos::mechanics::fem::prepareWriteDisplacementforPython("T3");
    writeDisplacementforPython(mesh, femodel, q, filename);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      s->computeOneStep();
      // osnspb->display();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(FEsolid->dimension() - 1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension() - 1);
      dataPlot(k, 3) = (*p)(0);

      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);
      // dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      k++;
      siconos::tools::progressBar((double)k / N);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                  .count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.resize(k, outputSize);
    siconos::algebra::io::write("T3_square_200.dat", dataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
#ifdef WITH_FRICTION

    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "T3_square_200_with_friction.ref", eps)) >= eps)
      return 1;
#else
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "T3_square_200.ref", eps)) >= eps)
      return 1;
#endif
    return 0;

  } catch (...) {
    std::cerr << "Exception caught in T3.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
