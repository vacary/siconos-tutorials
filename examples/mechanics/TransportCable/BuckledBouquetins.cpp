
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

#include <SiconosKernel.hpp>
#include <TransportCableProfil.h>
#include <TransportCableManager.h>
#include <TransportCableModel.h>
#include <chrono>
#include <ioVector.hpp>
#include <fstream>
#include <string>

int main()

{

  try {

    // Reads model parameters from a json file
    // Mechanical params, geometry, supports positions ...
    std::string modelFile = "buckled_bouquetins.json";
    auto model = std::make_shared<TransportCableModel>(modelFile);

    // Creates the result object: it handles two ropeways (sequence of ropes),
    // the supports (pylons and pulleys), and some extra variables related to
    // contacts.
    auto results = std::make_shared<TransportCableResult>();

    // Creates the profile, object which owns the model and the result and which
    // is responsible for the computation of the initial state/profile using
    // catenary and FEM.
    auto profil = std::make_shared<TransportCableProfil>(*model, *results);

    // -- Applies catenary equations to compute a first profile of the ropeways
    // --
    int nb_nodes = 50;  // Catenary, number of nodes per rope span
    double tol = 1e-10; // Tol. used in Newton-Raphson for catenary equation
    int nmax = 20;      // Newton-Raphson, max number of iterations
    profil->computeInitialProfil(nb_nodes, tol, nmax);

    // Save ropeways variables into json file
    ojson out;
    results->to_json(out, "ropeway");
    // results->to_json(out);

    std::ofstream out_ropes("catenary.json");
    out_ropes << std::setw(4) << out << std::endl;

    // std::ifstream in("bouquetins_ref.json");
    // json reader;
    // in >> reader;
    // auto qref1 = ioVector::readVectorFromJson(reader["rope1"]["q"]);
    // auto qref2 = ioVector::readVectorFromJson(reader["rope2"]["q"]);

    // std::ifstream in2("ropes.json");
    // auto reader2 = nlohmann::json::parse()
    //  std::ifstream in("ropes.json");
    // auto q1 = ioVector::readVectorFromJson(out["rope1"]["q"]);
    // auto q2 = ioVector::readVectorFromJson(out["rope2"]["q"]);

    // std::cout << qref1->norm2() << " " << qref1->size() << std::endl;

    // std::cout << ((*qref1) == (*q1)) << std::endl;

    // -- Fem part --

    nb_nodes = 1400; // FEM number of nodes
    double eps = 0.1;
    tol = 1e-3; // tolerance used to activate constraints
    profil->computeFEM(nb_nodes, eps, tol);

    // ojson out;
    results->to_json(out);
    std::ofstream out_ropes_fem("ropes_fem.json");
    out_ropes_fem << std::setw(4) << out << std::endl;

    auto positions = ioVector::readVectorFromJson(out["q"]);
    std::cout << positions->norm2() << " " << positions->size() << std::endl;
    

    auto manager = std::make_shared<TransportCableManager>(modelFile);
    std::string outFile = "results.json";
    ojson out2;
    json args;
    auto res = manager->computeFEM(args, outFile, out2);
    
    return 0;
  }

  catch (...) {
    Siconos::exception::process();
    return 1;
  }
}
