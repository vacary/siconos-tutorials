/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include <stdio.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <memory>

#include "FiniteElementLinearTIDS.hpp"
#include "MeshUtils.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NodeFem1d2DR.hpp"
#include "NodeFem2d2DR.hpp"
#include "NumericsVerbose.h"
#include "SolverOptions.h"

using namespace std;
using namespace siconos::mechanics::fem::native;

class MyContactDetection : public InteractionManager {
 protected:
  double _initial_gap;

  unsigned int _contact_condition_tag;

  std::shared_ptr<SiconosVector> _normal;

  std::shared_ptr<SiconosVector> _tangent;

  std::shared_ptr<NonSmoothLaw> _nslaw;
  unsigned int _contact_frame_dimension = 1;

  std::shared_ptr<FiniteElementModel> _femodel;

  std::shared_ptr<FiniteElementLinearTIDS> _fesolid;
  std::list<std::shared_ptr<FENode> > _contacting_node_zone;

  struct contacting_node {
    std::shared_ptr<FENode> _node;
    std::shared_ptr<Interaction> _inter;
    unsigned int _node_index;

    contacting_node(std::shared_ptr<FENode> node, unsigned int node_index)
        : _node(node), _node_index(node_index) {}
  };

  std::list<contacting_node*> _contacting_nodes;

 public:
  MyContactDetection(double initial_gap, int contact_condition_tag,
                     std::shared_ptr<SiconosVector> normal,
                     std::shared_ptr<NonSmoothLaw> nslaw,
                     std::shared_ptr<FiniteElementLinearTIDS> fesolid)
      : InteractionManager(),
        _initial_gap(initial_gap),
        _contact_condition_tag(contact_condition_tag),
        _nslaw(nslaw),
        _normal(normal),
        _fesolid(fesolid)

  {
    _femodel = _fesolid->FEModel();

    _contacting_node_zone = *_femodel->contactingNodes(_contact_condition_tag);

    std::shared_ptr<NewtonImpactFrictionNSL> nslaw_with_friction =
        std::dynamic_pointer_cast<NewtonImpactFrictionNSL>(nslaw);
    if (nslaw_with_friction) _contact_frame_dimension = 2;

    if (_contact_frame_dimension == 2)  // create tangent vertor
    {
      _tangent = std::make_shared<SiconosVector>(2);
      _tangent->setValue(0, -(*_normal)(1));
      _tangent->setValue(1, (*_normal)(0));
    }
  }
  virtual ~MyContactDetection() {}

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
  void updateInteractions(std::shared_ptr<Simulation> simulation) {
    // std::cout<< "\nCall to updateInteractions in MyContactDetection" <<
    // std::endl;

    SiconosVector& displacement = *(_fesolid->q());
    std::shared_ptr<NonSmoothDynamicalSystem> solid = simulation->nonSmoothDynamicalSystem();

    // update the list of contacting node by brute contact detection
    //_contacting_nodes.clear();

    for (std::shared_ptr<FENode> n : _contacting_node_zone) {
      // std::cout << "n->num() is: " << n->num() << std::endl;
      if (fabs(n->y()) <= 1e-01)  // and fabs(n->x()) >= 1e-05)
      {
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
      std::shared_ptr<FENode> n = cn->_node;
      unsigned int node_idx = (*n->dofIndex())[0];
      if (cn->_inter)  // update interaction->relation
      {
        // std::cout << "update interaction" << std::endl;

        std::shared_ptr<SiconosVector> pc2;
        if (_contact_frame_dimension == 2) {
          std::shared_ptr<siconos::mechanics::fem::NodeFem2d2DR> r =
              std::static_pointer_cast<siconos::mechanics::fem::NodeFem2d2DR>(
                  cn->_inter->relation());
          pc2 = r->pc2();
        } else {
          std::shared_ptr<siconos::mechanics::fem::NodeFem1d2DR> r =
              std::static_pointer_cast<siconos::mechanics::fem::NodeFem1d2DR>(
                  cn->_inter->relation());
          pc2 = r->pc2();
        }

        pc2->setValue(0, cn->_node->x() + displacement(node_idx));
        pc2->setValue(1, -_initial_gap);
      } else  // create an interaction and link
      {
        // std::cout << "create interaction" << std::endl;
        std::shared_ptr<SiconosVector> pc2 = std::make_shared<SiconosVector>(2);
        pc2->setValue(0, cn->_node->x() + displacement(node_idx));
        pc2->setValue(1, -_initial_gap);
        std::shared_ptr<Relation> relation;
        if (_contact_frame_dimension == 2) {
          relation = std::make_shared<siconos::mechanics::fem::NodeFem2d2DR>(
              cn->_node, pc2, _normal, _tangent);
        } else {
          relation =
              std::make_shared<siconos::mechanics::fem::NodeFem1d2DR>(cn->_node, pc2, _normal);
        }
        std::shared_ptr<Interaction> inter = std::make_shared<Interaction>(_nslaw, relation);
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
  double Ly = 1.0;
  string gmsh_filename = "./mesh_data/hertz.msh2";

  std::shared_ptr<Mesh> mesh(createMeshFromGMSH2(gmsh_filename));
  // mesh->display(false);

  writeMeshforPython(mesh);

  int bulk_material_tag = 6;
  int contact_condition_tag = 4;
  int applied_force_tag = 5;

  double density = 7800.;
  // a very soft materila is used to postprocess large deformations.
  std::shared_ptr<Material> mat1 = std::make_shared<Material>(density, 210e6, 1 / 3.);
  std::map<unsigned int, std::shared_ptr<Material> > materials = {{bulk_material_tag, mat1}};

  try {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::shared_ptr<FiniteElementLinearTIDS> FEsolid =
        std::make_shared<FiniteElementLinearTIDS>(mesh, materials, Siconos::SPARSE);
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << "Assembly time : " << elapsed << " ms" << endl;
    std::cout << " " << std::endl;
    // FEsolid->display(true);

    std::shared_ptr<FiniteElementModel> femodel = FEsolid->FEModel();
    // FEsolid->K()->display();
    // FEsolid->mass()->display();
    //  getchar();

    /*------------------------------------------------- Applied forces  */

    std::shared_ptr<SiconosVector> nodal_forces = std::make_shared<SiconosVector>(2);
    nodal_forces->setZero();
    //(*nodal_forces)(0) = 1e6;
    (*nodal_forces)(1) = -1e6;
    FEsolid->applyNodalForces(applied_force_tag, nodal_forces);

    // FEsolid->fext()->display();
    // getchar();

    // /*------------------------------------------------- Boundary Conditions
    // */
    // /* This part should be hidden in a new BC function for a node number
    //  * and a dof index. */

    // std::shared_ptr<IndexInt> node_dof_index = std::make_shared<IndexInt>(0);
    // node_dof_index->push_back(0);
    // node_dof_index->push_back(1);

    // FEsolid->applyDirichletBoundaryConditions(boundary_condition_tag,
    // node_dof_index); FEsolid->boundaryConditions()->display();

    // -------------
    // --- Model ---
    // -------------
    double t0 = 0;       // initial computation time
    double T = 1e-02;    // final computation time
    double h = 1e-05;    // time step
    double theta = 1.0;  // theta for MoreauJeanOSI integrator

    std::shared_ptr<NonSmoothDynamicalSystem> solid =
        std::make_shared<NonSmoothDynamicalSystem>(t0, T);

    // add the dynamical system in the non smooth dynamical system
    solid->insertDynamicalSystem(FEsolid);

    /*------------------------------------------------- Contact Conditions  */

    double e = 0.0;
#define WITH_FRICTION
#ifdef WITH_FRICTION
    double mu = 1.0;
    std::shared_ptr<NonSmoothLaw> nslaw =
        std::make_shared<NewtonImpactFrictionNSL>(e, 0.0, mu, 2);
#else
    std::shared_ptr<NonSmoothLaw> nslaw = std::make_shared<NewtonImpactNSL>(e);
#endif

    double initial_gap = 0.0;  // Ly*5e-05;
    std::shared_ptr<SiconosVector> displacement = FEsolid->q();
    std::shared_ptr<SiconosVector> normal = std::make_shared<SiconosVector>(2);
    normal->setValue(0, 0.0);
    normal->setValue(1, 1.0);
    std::shared_ptr<MyContactDetection> collision_detection =
        std::make_shared<MyContactDetection>(initial_gap, contact_condition_tag, normal, nslaw,
                                             FEsolid);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    std::shared_ptr<MoreauJeanOSI> OSI = std::make_shared<MoreauJeanOSI>(theta);
    OSI->setIsWSymmetricDefinitePositive(true);
    OSI->setGamma(0.0);

    // -- (2) Time discretisation --
    std::shared_ptr<TimeDiscretisation> t = std::make_shared<TimeDiscretisation>(t0, h);

    // -- (3) one step non smooth problem
#ifdef WITH_FRICTION
    std::shared_ptr<OneStepNSProblem> osnspb = std::make_shared<FrictionContact>(2);
#else
    std::shared_ptr<OneStepNSProblem> osnspb = std::make_shared<LCP>();
#endif
    SolverOptions* options = osnspb->numericsSolverOptions().get();

    options->dparam[SICONOS_DPARAM_TOL] = 1e-10;
    // numerics_set_verbose(2);
    // solver_options_print(options);

    // -- (4) Simulation setup with (1) (2) (3)
    std::shared_ptr<TimeStepping> s = std::make_shared<TimeStepping>(solid, t, OSI, osnspb);

    s->insertInteractionManager(collision_detection);

    // =========================== End of model definition
    // ===========================

    // ================================= Computation
    // =================================

    int N = ceil((T - t0) / h);  // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    SiconosMatrix dataPlot(N + 1, outputSize);

    std::shared_ptr<SiconosVector> q = FEsolid->q();
    std::shared_ptr<SiconosVector> v = FEsolid->velocity();
    std::shared_ptr<SiconosVector> p = FEsolid->p(1);
    // std::shared_ptr<SiconosVector> lambda = inter->lambda(1);

    dataPlot(0, 0) = solid->t0();
    dataPlot(0, 1) = (*q)(FEsolid->dimension() - 1);
    dataPlot(0, 2) = (*v)(FEsolid->dimension() - 1);
    dataPlot(0, 3) = (*p)(0);
    // dataPlot(0, 4) = (*lambda)(0);

    std::string filename = prepareWriteDisplacementforPython("T3");
    writeDisplacementforPython(mesh, femodel, q, filename);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    start = std::chrono::system_clock::now();

    while (s->hasNextEvent()) {
      s->computeOneStep();
      // osnspb->display();
      //  --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      // std::cout << (*q)(0) << std::endl;
      dataPlot(k, 1) = (*q)(FEsolid->dimension() - 1);
      dataPlot(k, 2) = (*v)(FEsolid->dimension() - 1);
      dataPlot(k, 3) = (*p)(0);

      if (k % 1 == 0) writeDisplacementforPython(mesh, femodel, q, filename);

      std::cout << "numerics -- "
                << " iterations: " << options->iparam[SICONOS_IPARAM_ITER_DONE]
                << " precision: " << options->dparam[SICONOS_DPARAM_RESIDU] << std::endl;

      // dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      // std::cout << "y     " ;
      // s->y(0,0)->display();
      // std::cout << "ydot  ";
      // s->y(1,0)->display();
      // std::cout << "lambda";
      // s->lambda(1,0)->display();

      double y_max = 0.0;

      auto y_max = (s->y(0, 0))->minCoeff();
      std::cout << "y_max violation " << std::max(y_max, 0.0) << std::endl;
      k++;
      // progressBar((double)k/N);
    }
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation time : " << elapsed << " ms" << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("Hertz_T3_1612.dat", "ascii", dataPlot, "noDim");
    double error = 0.0, eps = 1e-12;
#ifdef WITH_FRICTION
    if ((error = ioMatrix::compareRefFile(dataPlot, "Hertz_T3_1612_with_friction.ref", eps)) >=
            0.0 &&
        error > eps)
      return 1;
#else
    if ((error = ioMatrix::compareRefFile(dataPlot, "Hertz_T3_1612.ref", eps)) >= 0.0 &&
        error > eps)
      return 1;
#endif

  } catch (...) {
    cerr << "Exception caught in T3.cpp" << endl;
    Siconos::exception::process();
    return 1;
  }
}
