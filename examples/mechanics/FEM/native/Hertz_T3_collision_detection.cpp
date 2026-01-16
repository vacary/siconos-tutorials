/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

#include <SolverOptions.h>

#include <FENode.hpp>
#include <FemTools.hpp>
#include <FiniteElementLinearTIDS.hpp>
#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <Mesh.hpp>
#include <MeshUtils.hpp>
#include <NodeFem1d2DR.hpp>
#include <NodeFem2d2DR.hpp>
#include <SiconosKernel.hpp>
#include <SiconosMatrix.hpp>
#include <SiconosVector.hpp>
#include <chrono>
#include <cmath>  // for fabs

// class MyContactDetection : public siconos::simulation::InteractionManager {
//  protected:
//   double _initial_gap{0.};
//   ;

//   unsigned int _contact_condition_tag;

//   std::shared_ptr<Vector> _normal{nullptr};

//   std::shared_ptr<Vector> _tangent{nullptr};

//   std::shared_ptr<siconos::modeling::NonSmoothLaw> _nslaw{nullptr};
//   unsigned int _contact_frame_dimension = 1;

//   std::shared_ptr<siconos::mechanics::fem::FiniteElementModel> _femodel{nullptr};

//   std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS> _fesolid;
//   std::list<std::shared_ptr<siconos::mechanics::fem::FENode> > _contacting_node_zone;

//   struct contacting_node {
//     std::shared_ptr<siconos::mechanics::fem::FENode> _node{nullptr};
//     std::shared_ptr<siconos::modeling::Interaction> _inter{nullptr};
//     unsigned int _node_index{0};

//     contacting_node(std::shared_ptr<siconos::mechanics::fem::FENode> node,
//                     unsigned int node_index)
//         : _node(node), _node_index(node_index) {}
//   };

//   std::list<contacting_node*> _contacting_nodes;

//  public:
//   MyContactDetection(double initial_gap, int contact_condition_tag,
//                      std::shared_ptr<siconos::algebra::SiconosVector> normal,
//                      std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw,
//                      std::shared_ptr<siconos::mechanics::fem::FiniteElementLinearTIDS>
//                      fesolid)
//       : InteractionManager(),
//         _initial_gap(initial_gap),
//         _contact_condition_tag(contact_condition_tag),
//         _nslaw(nslaw),
//         _normal(normal),
//         _fesolid(fesolid)

//   {
//     _femodel = _fesolid->FEModel();

//     _contacting_node_zone = *_femodel->contactingNodes(_contact_condition_tag);

//     auto nslaw_with_friction =
//         std::dynamic_pointer_cast<siconos::modeling::NewtonImpactFrictionNSL>(nslaw);

//     if (nslaw_with_friction) _contact_frame_dimension = 2;

//     if (_contact_frame_dimension == 2)  // create tangent vertor
//     {
//       _tangent = std::make_shared<Vector>(2);
//       (*_tangent)(0) = -(*_normal)(1);
//       (*_tangent)(1) = (*_normal)(0);
//     }
//   }
//   virtual ~MyContactDetection() noexcept = default;

//   // shoud be done with find and a lambda function
//   contacting_node* find_contacting_node(unsigned int node_index) {
//     for (contacting_node* cnn : _contacting_nodes) {
//       if (cnn->_node_index == node_index) {
//         // std::cout << "existing interaction" << std::endl;
//         return cnn;
//       }
//     }
//     return nullptr;
//   }

//   /** Called by Simulation after updating positions prior to starting
//    * the Newton loop. */
//   void updateInteractions(std::shared_ptr<siconos::simulation::Simulation> simulation) {
//     // std::cout<< "\nCall to updateInteractions in MyContactDetection" <<
//     // std::endl;

//     Vector& displacement = *(_fesolid->q());
//     auto solid = simulation->nonSmoothDynamicalSystem();

//     // update the list of contacting node by brute contact detection
//     //_contacting_nodes.clear();

//     for (auto n : _contacting_node_zone) {
//       // std::cout << "n->num() is: " << n->num() << std::endl;
//       if (fabs(n->y()) <= 1e-01)  // and fabs(n->x()) >= 1e-05)
//       {
//         // std::cout << "contact node number : " << n->num() << " " << n->y() <<
//         // std::endl;
//         contacting_node* cn = new contacting_node(n, n->global_dof_index()[0]);
//         if (!find_contacting_node(cn->_node_index)) {
//           // std::cout << "add contact node number : " << n->num() << " is
//           // contacting list " <<  std::endl;
//           _contacting_nodes.push_back(cn);
//         }
//       } else {
//         contacting_node* cn = find_contacting_node(n->global_dof_index()[0]);
//         if (cn) {
//           // std::cout << "remove contact node number : " << n->num() << " is
//           // contacting list " <<  std::endl;
//           _contacting_nodes.remove(cn);
//         }
//       }
//     }

//     // update/create Interaction for contacting points
//     for (contacting_node* cn : _contacting_nodes) {
//       auto n = cn->_node;
//       unsigned int node_idx = n->global_dof_index()[0];
//       if (cn->_inter)  // update interaction->relation
//       {
//         // std::cout << "update interaction" << std::endl;

//         std::shared_ptr<Vector> pc2;
//         if (_contact_frame_dimension == 2) {
//           std::shared_ptr<siconos::mechanics::fem::NodeFem2d2DR> r =
//               std::static_pointer_cast<siconos::mechanics::fem::NodeFem2d2DR>(
//                   cn->_inter->relation());
//           pc2 = r->pc2();
//         } else {
//           std::shared_ptr<siconos::mechanics::fem::NodeFem1d2DR> r =
//               std::static_pointer_cast<siconos::mechanics::fem::NodeFem1d2DR>(
//                   cn->_inter->relation());
//           pc2 = r->pc2();
//         }

//         (*pc2)(0) = cn->_node->x() + displacement(node_idx);
//         (*pc2)(1) = -_initial_gap;
//       } else  // create an interaction and link
//       {
//         // std::cout << "create interaction" << std::endl;
//         auto pc2 = std::make_shared<Vector>(2);
//         (*pc2)(0) = cn->_node->x() + displacement(node_idx);
//         (*pc2)(1) = -_initial_gap;
//         std::shared_ptr<siconos::modeling::Relation> relation;
//         if (_contact_frame_dimension == 2) {
//           relation = std::make_shared<siconos::mechanics::fem::NodeFem2d2DR>(
//               cn->_node, pc2, _normal, _tangent);
//         } else {
//           relation =
//               std::make_shared<siconos::mechanics::fem::NodeFem1d2DR>(cn->_node, pc2,
//               _normal);
//         }
//         auto inter = std::make_shared<siconos::modeling::Interaction>(_nslaw, relation);
//         cn->_inter = inter;
//         // link the interaction and the dynamical system
//         solid->link(inter, _fesolid);
//       }
//     }
//     // std::cout << "end to updateInteractions in MyContactDetection" <<
//     // std::endl;
//   }
// };

int main(int argc, char* argv[]) {
  try {
    double Ly = 1.0;
    auto gmsh_filename = "./mesh_data/hertz.msh2";
    // Applied forces
    siconos::algebra::SiconosVector nodal_forces{2};
    nodal_forces << 0., -1e6;

    siconos::mechanics::fem::Tags tags;
    tags[siconos::mechanics::fem::MeshTags::bulk_material] = 6;
    tags[siconos::mechanics::fem::MeshTags::applied_forces] = 5;
    // No BC

    siconos::mechanics::fem::Material mat{7800, 210e6, 1. / 3};

    auto FEsolid = siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
        gmsh_filename, tags, mat, nodal_forces, {});

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;     // initial computation time
    double T = 1e-02;  // final computation time

    auto solid = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */
    double e = 0.0;
    double mu = 1.0;
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(e, 0.0, mu, 2);

    double initial_gap = 0.0;  // Ly*5e-05;
    auto displacement = FEsolid->q();
    auto normal = std::make_shared<siconos::algebra::SiconosVector>(2);
    (*normal)(0) = 0.0;
    (*normal)(1) = 1.0;
    int contact_condition_tag = 4;

    auto contact_condition = [](const siconos::mechanics::fem::FENode& node) {
      return (fabs(node.y()) <= 1e-01);  // and fabs(n->x()) >= 1e-05));
    };

    auto collision_detection = std::make_shared<siconos::mechanics::fem::ContactDetection>(
        initial_gap, normal, nslaw, FEsolid, contact_condition, contact_condition_tag);

    // ------------------
    // --- Simulation ---
    // ------------------
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    // -- (1) OneStepIntegrators --
    auto OSI = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);
    OSI->setGamma(0.0);

    // -- (2) Time discretisation --
    auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
    auto osnspb = std::make_shared<siconos::nonsmooth_formulations::FrictionContact>(2);
    auto options = osnspb->numericsSolverOptions();

    options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    // numerics_set_verbose(2);
    // solver_options_print(options);

    // -- (4) Simulation setup with (1) (2) (3)
    auto s = std::make_shared<siconos::simulation::TimeStepping>(solid, t, OSI, osnspb);

    s->insertInteractionManager(collision_detection);

    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);

    auto q = FEsolid->q_read();
    auto v = FEsolid->velocity_read();
    auto p = FEsolid->p_read(1);
    // auto lambda = inter->lambda(1);

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = q(FEsolid->dimension() - 1);
    dataPlot(0, 2) = v(FEsolid->dimension() - 1);
    dataPlot(0, 3) = p(0);
    // dataPlot(0, 4) = (*lambda)(0);

    auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython(
        "Hertz_T3_collision_detection");
    auto femodel = FEsolid->FEModel();
    auto mesh = femodel->mesh();
    siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);

    // --- Time loop ---
    std::cout << "====> Start computation ... \n";
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    auto start = std::chrono::system_clock::now();
    while (s->hasNextEvent()) {
      s->computeOneStep();
      // siconos::algebra::print(*osnspb);
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = q(FEsolid->dimension() - 1);
      dataPlot(k, 2) = v(FEsolid->dimension() - 1);
      dataPlot(k, 3) = p(0);

      if (k % 1 == 0)
        siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q, filename);

      std::cout << "numerics -- "
                << " iterations: " << options->iparam[SICONOS_IPARAM_ITER_DONE]
                << " precision: " << options->dparam[SICONOS_DPARAM_RESIDU] << std::endl;

      // dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();

      // double y_max = 0.0;

      // auto y_max = (s->y(0, 0))->minCoeff();
      // std::cout << "y_max violation " << std::max(y_max, 0.0) << std::endl;
      k++;
      // progressBar((double)k/N);
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    dataPlot.conservativeResize(k, outputSize);
    siconos::algebra::io::write("Hertz_T3_1612.dat", dataPlot, siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             dataPlot, "Hertz_T3_1612_with_friction.ref", eps)) >= 0.0 &&
        error > eps)
      return 1;
  } catch (...) {
    std::cerr << "Exception caught in Hertz_T3_collision_detection.cpp\n";
    siconos::exception::process();
    return 1;
  }
}
