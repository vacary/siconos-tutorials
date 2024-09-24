/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include <SolverOptions.h>

#include <SiconosKernel.hpp>
#include <chrono>
#include <numbers>

using Matrix = siconos::algebra::SiconosMatrix;
using Vector = siconos::algebra::SiconosVector;

using namespace std;

constexpr double GGearth = 9.8100;

//---------------------------------------------------
double LengthBlock = 1.0;                  // Length of the rocking block
double HeightBlock = 0.5;                  // Height of the rocking block
unsigned int Nfreedom = 3;                 // Number of degrees of freedom
unsigned int Ncontact = 2;                 // Number of contacts
double MassBlock = 1.0;                    // Mass of the rocking block
double PosXiniPointA = 0.0;                // Initial coordinate X of the point A
double PosYiniPointA = 0.0;                // Initial coordinate Y of the point A
double AngleThetaIni = numbers::pi / 3.0;  // Initial angle theta of the block
double VelXiniPointA = 0.0;                // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0;                // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;               // Initial angular velocity of the block
double e = 0.9;                            // Restitution coefficient
double TimeInitial = 0.0;                  // Initial time of the simulation
double TimeFinal = 0.58;                   // Final time of the simulation
double _rho = 0.99;                        // used to computer parameters for NewMark Scheme
double StepSize = 0.001;                   // Time step size
unsigned int NpointSave = 1500;            //
unsigned int SizeOutput = 13;              //
unsigned int maxIter = 20000;
bool IsTreatFirstSteps = false;
bool IsHandleVelConstraint = false;
//==========================================================================================================
//                                             Main function
//==========================================================================================================
int main(int argc, char* argv[]) {
  //---------------------------- calculate the computation time
  //--------------------------------------------------
  try {
    //===========================================================================================================
    //                  I: Declare the dynamical systems
    //===========================================================================================================
    // 1. Set the mass matrix
    auto Mass = std::make_shared<Matrix>(Nfreedom, Nfreedom);
    double InertiaBlock;
    InertiaBlock = (MassBlock / 12.0) * ((HeightBlock * HeightBlock) +
                                         (LengthBlock * LengthBlock));  // moment of inertia
    (*Mass)(0, 0) = MassBlock;
    (*Mass)(1, 1) = MassBlock;
    (*Mass)(2, 2) = InertiaBlock;
    // 2. Set the initial position of the block in function of the initial position of the
    // contact point A (left-hand contact)
    auto PosIniBlock = std::make_shared<Vector>(Nfreedom);
    (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) -
                        0.5 * HeightBlock * sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) +
                        0.5 * HeightBlock * cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;
    cout.precision(15);
    cout << "x0: " << (*PosIniBlock)(0) << endl;
    cout << "y0: " << (*PosIniBlock)(1) << endl;
    cout << "theta0: " << (*PosIniBlock)(2) << endl;
    cout.precision(15);
    // (*PosIniBlock)(0) = 0.5;
    // (*PosIniBlock)(1) = 0.5;
    // (*PosIniBlock)(2) = 0.0;

    // 3. Set the initial velocity of the block in function of the initial relative velocity of
    // the contact point A
    auto VelIniBlock = std::make_shared<Vector>(Nfreedom);
    (*VelIniBlock)(0) = VelXiniPointA - (0.5 * LengthBlock * sin(AngleThetaIni) +
                                         0.5 * HeightBlock * cos(AngleThetaIni)) *
                                            RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5 * LengthBlock * cos(AngleThetaIni) -
                                         0.5 * HeightBlock * sin(AngleThetaIni)) *
                                            RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;

    // (*VelIniBlock)(0) = 0.0;
    // (*VelIniBlock)(1) = 0.0;
    // (*VelIniBlock)(2) = 0.0;

    // 4. Instantiate the object of "LagrangianTIDS"
    auto RockingBlock = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
        PosIniBlock, VelIniBlock, Mass);
    // 5. Set the external force
    Vector ForceExtern{Nfreedom};
    ForceExtern.setZero();
    ForceExtern(1) = -MassBlock * GGearth;
    RockingBlock->setConstantFExt(ForceExtern);    std::vector<double> zparams = {LengthBlock, HeightBlock};
    auto zz = std::make_shared<Vector>(zparams);
    RockingBlock->setzPtr(zz);

    //
    //----------------------------- Display variables of the dynamical
    // system---------------------------------------
    cout << "Initial position of the rocking block:\n";
    PosIniBlock->display();
    cout << "Initial velocity of the rocking block:\n";
    VelIniBlock->display();
    cout << "Mass matrix of the rocking block:\n";
    Mass->display();
    cout << "External force applied on the rocking block:" << endl;
    ForceExtern->display();
    //==================================================================================================================
    //              II: Declare the relation et interaction between dynamical systems
    //==================================================================================================================
    //
    /*
    auto H= std::make_shared<Matrix>(1,Nfreedom);
    (*H)(0,1) = 1.0;
    auto E= std::make_shared<Vector>(1);
    (*E)(0) = -0.5*HeightBlock;
    */
    // Impact law
    auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
    // Interaction at contact point 1
    // auto relation1= std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, E));
    auto relation1 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "RockingBlockPlugin:h1", "RockingBlockPlugin:G1", "RockingBlockPlugin:G1dot");
    auto inter1 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation1);
    // Interaction at contact point 2
    // auto relation2= std::make_shared<siconos::modeling::LagrangianLinearTIR>(H, E));
    auto relation2 = std::make_shared<siconos::modeling::LagrangianScleronomousR>(
        "RockingBlockPlugin:h2", "RockingBlockPlugin:G2", "RockingBlockPlugin:G2dot");
    auto inter2 = std::make_shared<siconos::modeling::Interaction>(nslaw, relation2);
    // Interactions for the whole dynamical system
    //================================================================================================================
    //            III. Create the "model" object
    //================================================================================================================
    auto RoBlockModel =
        std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(TimeInitial, TimeFinal);
    RoBlockModel->insertDynamicalSystem(RockingBlock);
    RoBlockModel->link(inter1, RockingBlock);
    RoBlockModel->link(inter2, RockingBlock);
    //================================================================================================================
    //            IV. Create the simulation
    //================================================================================================================
    // 1. Time discretization
    auto TimeDiscret =
        std::make_shared<siconos::simulation::TimeDiscretisation>(TimeInitial, StepSize);
    // 2. Integration solver for one step
    auto NewMarkAlpha =
        std::make_shared<siconos::integrators::NewMarkAlphaOSI>(_rho, IsHandleVelConstraint);
    // 3. Nonsmooth problem
    auto impact = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto position = std::make_shared<siconos::nonsmooth_formulations::LCP>();
    auto acceleration = std::make_shared<siconos::nonsmooth_formulations::LCP>();

    position->numericsSolverOptions()->dparam[0] = 1e-12;
    impact->numericsSolverOptions()->dparam[0] = 1e-12;
    acceleration->numericsSolverOptions()->dparam[0] = 1e-12;

    // 4. Simulation with (1), (2), (3)
    auto EDscheme =
        std::make_shared<siconos::simulation::EventDriven>(RoBlockModel, TimeDiscret);
    EDscheme->insertIntegrator(NewMarkAlpha);
    EDscheme->insertNonSmoothProblem(impact, siconos::simulation::SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_ACC);
    EDscheme->insertNonSmoothProblem(position,
                                     siconos::simulation::SICONOS_OSNSP_ED_SMOOTH_POS);
    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // cout << "Impact law included in the simulation: " << check1 << endl;
    // cout << "LCP at acceleration level included in the simulation: " << check2 << endl;
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization
    // ------------------------------------------------------
    EDscheme->setPrintStat(true);
    auto eventsManager =
        EDscheme->eventsManager();  // ponters point to the "eventsManager" object
    auto PosBlock =
        RockingBlock->q();  // pointer points to the position vector of the rocking block
    auto VelBlock =
        RockingBlock->velocity();  // pointer points to the velocity of the rocking block
    auto AcceBlock =
        RockingBlock->acceleration();  // pointer points to the velocity of the rocking block

    auto indexSet0 = RoBlockModel->topology()->indexSet(0);
    cout << "Size of IndexSet0: " << indexSet0->size() << endl;

    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    //-------------------- Save the output during simulation
    //---------------------------------------------------------
    Matrix DataPlot(NpointSave, SizeOutput);
    //------------- At the initial time
    //-----------------------------------------------------------------------------
    DataPlot(0, 0) = RoBlockModel->t0();
    DataPlot(0, 1) = (*PosBlock)(0);  // Position X
    DataPlot(0, 2) = (*PosBlock)(1);  // Position Y
    DataPlot(0, 3) = (*PosBlock)(2);  // Angle theta
    DataPlot(0, 4) = (*VelBlock)(0);  // Velocity Vx
    DataPlot(0, 5) = (*VelBlock)(1);  // Velocity Vy
    DataPlot(0, 6) = (*VelBlock)(2);  // Angular velocity
    DataPlot(0, 7) = 0.0;             // Gap at first contact
    DataPlot(0, 8) = 0.0;             // Gap at second contact
    DataPlot(0, 9) = 0.0;             // Relative velocity at first contact
    DataPlot(0, 10) = 0.0;            // Relative velocity at second contact
    DataPlot(0, 11) = 0.0;            // Force at first contact
    DataPlot(0, 12) = 0.0;            // Force at second contact
    //----------------------------------- Simulation starts
    //----------------------------------------------------------
    cout << "====> Start computation ... " << endl << endl;
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    double alpha_m, alpha_f, beta, gamma;
    unsigned int k = 1;
    auto start = std::chrono::system_clock::now();
    while (EDscheme->hasNextEvent() && (k < NpointSave)) {
      if (IsTreatFirstSteps) {
        if (k == 1)  // first step
        {
          alpha_m = 0.0;
          alpha_f = 0.0;
          gamma = 1.0 / 2.0 + 1.0 / 10.0;
          beta = 1.0 / 4.0 * (gamma + 1.0 / 2.0) * (gamma + 1.0 / 2.0);
          NewMarkAlpha->setAlpha_m(alpha_m);
          NewMarkAlpha->setAlpha_f(alpha_f);
          NewMarkAlpha->setBeta(beta);
          NewMarkAlpha->setGamma(gamma);
        } else if (k == 2) {
          alpha_m = 0.0;
          alpha_f = -1.0 / 3.0;  // -1/3 <= alpha_f <= 0
          gamma = 1.0 / 2.0 - alpha_f;
          beta = 1.0 / 4.0 * (1.0 - alpha_f) * (1.0 - alpha_f);
          NewMarkAlpha->setAlpha_m(alpha_m);
          NewMarkAlpha->setAlpha_f(alpha_f);
          NewMarkAlpha->setBeta(beta);
          NewMarkAlpha->setGamma(gamma);
        } else {
          NewMarkAlpha->setParametersFromRho_infty(_rho);
        }
      }
      EDscheme->advanceToEvent();  // lead the simulation run from one event to the next
      auto GapCon1 = inter1->y(0);
      auto GapCon2 = inter2->y(0);
      auto VelCon1 = inter1->y(1);
      auto VelCon2 = inter2->y(1);

      auto LambdaCon1 = inter1->lambda(2);
      auto LambdaCon2 = inter2->lambda(2);
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() ==
          siconos::simulation::EventType::NS)  // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events
      //---------------------------
      if (NSEvent) {
        DataPlot(k, 0) = EDscheme->startingTime();  // instant at non-smooth event
        DataPlot(k, 1) = RockingBlock->qMemory().getSiconosVector(0)(0);         // Position X
        DataPlot(k, 2) = RockingBlock->qMemory().getSiconosVector(0)(1);         // Position Y
        DataPlot(k, 3) = RockingBlock->qMemory().getSiconosVector(0)(2);         // Angle theta
        DataPlot(k, 4) = RockingBlock->velocityMemory().getSiconosVector(0)(0);  // Velocity Vx
        DataPlot(k, 5) = RockingBlock->velocityMemory().getSiconosVector(0)(1);  // Velocity Vy
        DataPlot(k, 6) =
            RockingBlock->velocityMemory().getSiconosVector(0)(2);  // Angular velocity
        // EDscheme->update(1);
        k++;
        ++NumberNSEvent;

        NSEvent = false;
        // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events
      //---------------
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = (*PosBlock)(0);     // Position X
      DataPlot(k, 2) = (*PosBlock)(1);     // Position Y
      DataPlot(k, 3) = (*PosBlock)(2);     // Position theta
      DataPlot(k, 4) = (*VelBlock)(0);     // Velocity Vx
      DataPlot(k, 5) = (*VelBlock)(1);     // Velocity Vy
      DataPlot(k, 6) = (*VelBlock)(2);     // Velocity Vtheta
      DataPlot(k, 7) = (*GapCon1)(0);      // Gap at first contact
      DataPlot(k, 8) = (*GapCon2)(0);      // Gap at second contact
      DataPlot(k, 9) = (*VelCon1)(0);      // Relative velocity at first contact
      DataPlot(k, 10) = (*VelCon2)(0);     // Relative velocity at second contact
      DataPlot(k, 11) = (*LambdaCon1)(0);  // Force at first contact
      DataPlot(k, 12) = (*LambdaCon2)(0);  // Force at second contact
      // go to the next time step
      k++;

      // // Display information
      // cout << "********At the end of integation step***************"<< (k - 1) << endl;
      // cout << "Information on Dynamical System\n";
      // cout << "Position: ";
      // PosBlock->display();
      // cout << "Velocity: ";
      // VelBlock->display();
      // cout << "Acceleration: ";
      // AcceBlock->display();
      // cout << "Information on contacts\n";
      // for(std::tie(ui,uiend) = indexSet0->vertices(); ui!=uiend; ++ui)
      //   {
      //     auto inter = indexSet0->bundle(*ui);
      //     cout << "Contact number: " << inter->number() << endl;
      //     cout << "Contact gap: ";
      //     inter->y(0)->display();
      //     cout << "Contact relative velocity: ";
      //     inter->y(1)->display();
      //     cout << "Contact Force: \n";
      //     inter->lambda(2)->display();
      //   }
    };
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
    std::cout << "\nComputation time : " << elapsed << " ms\n";

    // --- Output files ---
    std::cout << "====> Output file writing ...\n";
    DataPlot.resize(k, SizeOutput);
    siconos::algebra::io::write("RockingBlockED_NewMarkAlpha.dat", DataPlot,
                                siconos::algebra::io::ASCII_OUT,
                                siconos::algebra::io::WriteType::nodim);
    double error = 0.0, eps = 1e-12;
    if ((error = siconos::algebra::io::compareRefFile(
             DataPlot, "RockingBlockED_NewMarkAlpha.ref", eps)) > eps)
      return 1;
    return 0;
  }
  //============================== Catch exceptions
  //===================================================================
  catch (...) {
    siconos::exception::process();
    return 1;
  }
}
