// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include <chrono>
#include "SiconosKernel.hpp"
#include <SolverOptions.h>
#include <stdlib.h>
using namespace std;
#define PI 3.141592653589793
#define GGearth  9.8100
//---------------------------------Decalre global variables ---------------------------------------------------
double LengthBlock = 1.0;        // Length of the rocking block
double HeightBlock = 0.5;        // Height of the rocking block
unsigned int Nfreedom = 3;       // Number of degrees of freedom
unsigned int Ncontact = 2;       // Number of contacts
double MassBlock = 1.0;          // Mass of the rocking block
double PosXiniPointA = 0.0;      // Initial coordinate X of the point A
double PosYiniPointA = 0.0;      // Initial coordinate Y of the point A
double AngleThetaIni = PI / 3.0; // Initial angle theta of the block
double VelXiniPointA = 0.0 ;     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0 ;     // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;    // Initial angular velocity of the block
double e = 0.9;       // Restitution coefficient
double TimeInitial = 0.0;        // Initial time of the simulation
double TimeFinal =  0.58;      // Final time of the simulation
double _rho = 0.99;             // used to computer parameters for NewMark Scheme
double StepSize = 0.001;         // Time step size
unsigned int NpointSave = 1500;   //
unsigned int SizeOutput = 13;     //
unsigned int maxIter = 20000;
bool IsTreatFirstSteps = false;
bool IsHandleVelConstraint = false;
//==========================================================================================================
//                                             Main function
//==========================================================================================================
int main(int argc, char* argv[])
{
  //---------------------------- calculate the computation time --------------------------------------------------
  std::chrono::time_point<std::chrono::system_clock> start, end;
  start = std::chrono::system_clock::now();
  try
  {
    //===========================================================================================================
    //                  I: Declare the dynamical systems
    //===========================================================================================================
    //1. Set the mass matrix
    SP::SiconosMatrix Mass(new SimpleMatrix(Nfreedom, Nfreedom));
    double InertiaBlock;
    InertiaBlock = (MassBlock / 12.0) * ((HeightBlock*HeightBlock) + (LengthBlock*LengthBlock)); // moment of inertia
    (*Mass)(0, 0) = MassBlock;
    (*Mass)(1, 1) = MassBlock;
    (*Mass)(2, 2) = InertiaBlock;
    //2. Set the initial position of the block in function of the initial position of the contact point A (left-hand contact)
    SP::SiconosVector PosIniBlock(new SiconosVector(Nfreedom));
    (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;
    cout.precision(15);
    cout << "x0: " << (*PosIniBlock)(0) << endl;
    cout << "y0: " << (*PosIniBlock)(1) << endl;
    cout << "theta0: " << (*PosIniBlock)(2) << endl;
    cout.precision(15);
    cout << "PI: " << PI << endl;
    // (*PosIniBlock)(0) = 0.5;
    // (*PosIniBlock)(1) = 0.5;
    // (*PosIniBlock)(2) = 0.0;

    //3. Set the initial velocity of the block in function of the initial relative velocity of the contact point A
    SP::SiconosVector VelIniBlock(new SiconosVector(Nfreedom));
    (*VelIniBlock)(0) = VelXiniPointA - (0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;

    // (*VelIniBlock)(0) = 0.0;
    // (*VelIniBlock)(1) = 0.0;
    // (*VelIniBlock)(2) = 0.0;

    //4. Instantiate the object of "LagrangianTIDS"
    SP::LagrangianLinearTIDS RockingBlock(new LagrangianLinearTIDS(PosIniBlock, VelIniBlock, Mass));
    //5. Set the external force
    SP::SiconosVector ForceExtern(new SiconosVector(Nfreedom));
    (*ForceExtern)(1) = -MassBlock * GGearth;
    RockingBlock->setFExtPtr(ForceExtern);
    std::vector<double> zparams = {LengthBlock, HeightBlock};
    SP::SiconosVector zz(new SiconosVector(zparams));
    RockingBlock->setzPtr(zz);

    //
    //----------------------------- Display variables of the dynamical system---------------------------------------
    cout << "Initial position of the rocking block:" << endl;
    PosIniBlock->display();
    cout << "Initial velocity of the rocking block:" << endl;
    VelIniBlock->display();
    cout << "Mass matrix of the rocking block:" << endl;
    Mass->display();
    cout << "External force applied on the rocking block:"  << endl;
    ForceExtern->display();
    //==================================================================================================================
    //              II: Declare the relation et interaction between dynamical systems
    //==================================================================================================================
    //
    /*
    SP::SiconosMatrix H(new SimpleMatrix(1,Nfreedom));
    (*H)(0,1) = 1.0;
    SP::SiconosVector E(new SiconosVector(1));
    (*E)(0) = -0.5*HeightBlock;
    */
    // Impact law
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    // Interaction at contact point 1
    //SP::Relation relation1(new LagrangianLinearTIR(H, E));
    SP::Relation relation1(new LagrangianScleronomousR("RockingBlockPlugin:h1", "RockingBlockPlugin:G1", "RockingBlockPlugin:G1dot"));
    SP::Interaction inter1(new Interaction(nslaw, relation1));
    // Interaction at contact point 2
    //SP::Relation relation2(new LagrangianLinearTIR(H, E));
    SP::Relation relation2(new LagrangianScleronomousR("RockingBlockPlugin:h2", "RockingBlockPlugin:G2", "RockingBlockPlugin:G2dot"));
    SP::Interaction inter2(new Interaction(nslaw, relation2));
    // Interactions for the whole dynamical system
    //================================================================================================================
    //            III. Create the "model" object
    //================================================================================================================
    SP::NonSmoothDynamicalSystem RoBlockModel(new NonSmoothDynamicalSystem(TimeInitial, TimeFinal));
    RoBlockModel->insertDynamicalSystem(RockingBlock);
    RoBlockModel->link(inter1, RockingBlock);
    RoBlockModel->link(inter2, RockingBlock);
    //================================================================================================================
    //            IV. Create the simulation
    //================================================================================================================
    //1. Time discretization
    SP::TimeDiscretisation TimeDiscret(new TimeDiscretisation(TimeInitial, StepSize));
    //2. Integration solver for one step
    SP::OneStepIntegrator OSI(new NewMarkAlphaOSI(_rho, IsHandleVelConstraint));
    SP::NewMarkAlphaOSI _NewMarkAlpha = std::static_pointer_cast<NewMarkAlphaOSI>(OSI);
    //3. Nonsmooth problem
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem position(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());

    position->numericsSolverOptions()->dparam[0] = 1e-12;
    impact->numericsSolverOptions()->dparam[0] = 1e-12;
    acceleration->numericsSolverOptions()->dparam[0] = 1e-12;


    //4. Simulation with (1), (2), (3)
    SP::Simulation EDscheme(new EventDriven(RoBlockModel, TimeDiscret));
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    EDscheme->insertNonSmoothProblem(position, SICONOS_OSNSP_ED_SMOOTH_POS);
    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // cout << "Impact law included in the simulation: " << check1 << endl;
    // cout << "LCP at acceleration level included in the simulation: " << check2 << endl;
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization ------------------------------------------------------
    EDscheme->setPrintStat(true);
    SP::EventsManager eventsManager = EDscheme->eventsManager(); // ponters point to the "eventsManager" object
    SP::SiconosVector PosBlock = RockingBlock->q();              // pointer points to the position vector of the rocking block
    SP::SiconosVector VelBlock = RockingBlock->velocity();       // pointer points to the velocity of the rocking block
    SP::SiconosVector AcceBlock = RockingBlock->acceleration();       // pointer points to the velocity of the rocking block

    SP::InteractionsGraph indexSet0 = RoBlockModel->topology()->indexSet(0);
    cout << "Size of IndexSet0: " << indexSet0->size() << endl;


    InteractionsGraph::VIterator ui, uiend;
    //-------------------- Save the output during simulation ---------------------------------------------------------
    SimpleMatrix DataPlot(NpointSave, SizeOutput);
    //------------- At the initial time -----------------------------------------------------------------------------
    DataPlot(0, 0) = RoBlockModel->t0();
    DataPlot(0, 1) = (*PosBlock)(0); // Position X
    DataPlot(0, 2) = (*PosBlock)(1); // Position Y
    DataPlot(0, 3) = (*PosBlock)(2); // Angle theta
    DataPlot(0, 4) = (*VelBlock)(0); // Velocity Vx
    DataPlot(0, 5) = (*VelBlock)(1); // Velocity Vy
    DataPlot(0, 6) = (*VelBlock)(2); // Angular velocity
    DataPlot(0, 7) = 0.0;  // Gap at first contact
    DataPlot(0, 8) = 0.0;  // Gap at second contact
    DataPlot(0, 9) = 0.0;  // Relative velocity at first contact
    DataPlot(0, 10) = 0.0;  // Relative velocity at second contact
    DataPlot(0, 11) = 0.0; // Force at first contact
    DataPlot(0, 12) = 0.0; // Force at second contact
    //----------------------------------- Simulation starts ----------------------------------------------------------
    cout << "====> Start computation ... " << endl << endl;
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    double alpha_m, alpha_f, beta, gamma;
    unsigned int k = 1;
    while(EDscheme->hasNextEvent() && (k < NpointSave))
    {
      if(IsTreatFirstSteps)
      {
        if(k == 1)  // first step
        {
          alpha_m = 0.0;
          alpha_f = 0.0;
          gamma = 1.0/2.0 + 1.0/10.0;
          beta = 1.0/4.0 * (gamma + 1.0/2.0) *(gamma + 1.0/2.0)  ;
          _NewMarkAlpha->setAlpha_m(alpha_m);
          _NewMarkAlpha->setAlpha_f(alpha_f);
          _NewMarkAlpha->setBeta(beta);
          _NewMarkAlpha->setGamma(gamma);
        }
        else if(k == 2)
        {
          alpha_m = 0.0;
          alpha_f = -1.0 / 3.0; // -1/3 <= alpha_f <= 0
          gamma = 1.0/2.0 - alpha_f;
          beta = 1.0/4.0 * (1.0 - alpha_f)*(1.0 - alpha_f);
          _NewMarkAlpha->setAlpha_m(alpha_m);
          _NewMarkAlpha->setAlpha_f(alpha_f);
          _NewMarkAlpha->setBeta(beta);
          _NewMarkAlpha->setGamma(gamma);
        }
        else
        {
          _NewMarkAlpha->setParametersFromRho_infty(_rho);
        }
      }
      EDscheme->advanceToEvent(); // lead the simulation run from one event to the next
      SP::SiconosVector GapCon1 = inter1->y(0);
      SP::SiconosVector GapCon2 = inter2->y(0);
      SP::SiconosVector VelCon1 = inter1->y(1);
      SP::SiconosVector VelCon2 = inter2->y(1);

      SP::SiconosVector LambdaCon1 = inter1->lambda(2);
      SP::SiconosVector LambdaCon2 = inter2->lambda(2);
      //---------- detect the statue of the current event ------------------------------------
      if(eventsManager->nextEvent()->getType() == 2)  // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events ---------------------------
      if(NSEvent)
      {
        DataPlot(k, 0) = EDscheme->startingTime(); // instant at non-smooth event
        DataPlot(k, 1) = RockingBlock->qMemory().getSiconosVector(0)(0);      // Position X
        DataPlot(k, 2) = RockingBlock->qMemory().getSiconosVector(0)(1);      // Position Y
        DataPlot(k, 3) = RockingBlock->qMemory().getSiconosVector(0)(2);      // Angle theta
        DataPlot(k, 4) = RockingBlock->velocityMemory().getSiconosVector(0)(0); // Velocity Vx
        DataPlot(k, 5) = RockingBlock->velocityMemory().getSiconosVector(0)(1); // Velocity Vy
        DataPlot(k, 6) = RockingBlock->velocityMemory().getSiconosVector(0)(2); // Angular velocity
        //EDscheme->update(1);
        k++;
        ++NumberNSEvent;

        NSEvent = false;
        // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events ---------------
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = (*PosBlock)(0); //Position X
      DataPlot(k, 2) = (*PosBlock)(1); //Position Y
      DataPlot(k, 3) = (*PosBlock)(2); // Position theta
      DataPlot(k, 4) = (*VelBlock)(0); // Velocity Vx
      DataPlot(k, 5) = (*VelBlock)(1); // Velocity Vy
      DataPlot(k, 6) = (*VelBlock)(2); // Velocity Vtheta
      DataPlot(k, 7) = (*GapCon1)(0);  // Gap at first contact
      DataPlot(k, 8) = (*GapCon2)(0);  // Gap at second contact
      DataPlot(k, 9) = (*VelCon1)(0);  // Relative velocity at first contact
      DataPlot(k, 10) = (*VelCon2)(0);  // Relative velocity at second contact
      DataPlot(k, 11) = (*LambdaCon1)(0); // Force at first contact
      DataPlot(k, 12) = (*LambdaCon2)(0); // Force at second contact
      // go to the next time step
      k++;

      // // Display information
      // cout << "********At the end of integation step***************"<< (k - 1) << endl;
      // cout << "Information on Dynamical System" << endl;
      // cout << "Position: ";
      // PosBlock->display();
      // cout << "Velocity: ";
      // VelBlock->display();
      // cout << "Acceleration: ";
      // AcceBlock->display();
      // cout << "Information on contacts" << endl;
      // for(std::tie(ui,uiend) = indexSet0->vertices(); ui!=uiend; ++ui)
      //   {
      //     SP::Interaction inter = indexSet0->bundle(*ui);
      //     cout << "Contact number: " << inter->number() << endl;
      //     cout << "Contact gap: ";
      //     inter->y(0)->display();
      //     cout << "Contact relative velocity: ";
      //     inter->y(1)->display();
      //     cout << "Contact Force: " << endl;
      //     inter->lambda(2)->display();
      //   }
    };
    //----------------------- At the end of the simulation --------------------------
    cout << "End of the simulation" << endl;
    cout << "Number of events processed during simulation: " << (k + 1) << endl;
    cout << "Number of non-smooth events: " << NumberNSEvent << endl;
    cout << "====> Output file writing ..." << endl << endl;
    DataPlot.resize(k,SizeOutput);
    ioMatrix::write("RockingBlockED_NewMarkAlpha.dat", "ascii", DataPlot, "noDim");

    Index index(11);
    for(int k =0; k< 11; k++) index.push_back(k);
    // Comparison with a reference file
    double error=0.0, eps=1e-10;
    if((error=ioMatrix::compareRefFile(DataPlot, "RockingBlockED_NewMarkAlpha.ref",
                                       eps, index)) >= 0.0
        && error > eps)
      return 1;
  }
  //============================== Catch exceptions ===================================================================
  catch(SiconosException e)
  {
    cerr << e.report() << endl;
    return 1;
  }
  catch(...)
  {
    cerr << "Exception caught." << endl;
    return 1;
  }
}
