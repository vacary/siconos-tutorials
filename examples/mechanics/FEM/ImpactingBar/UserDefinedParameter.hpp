namespace user_defined_ref {

// User-defined main parameters
constexpr unsigned int nDof = 100;           // degrees of freedom for the beam
constexpr double t0 = 1e-8;                  // initial computation time
constexpr double T = 0.0015;                 // final computation time
constexpr double h = 2e-6;                   // time step
constexpr double position_init = 0.00005;    // initial position
constexpr double velocity_init = -.1;        // initial velocity
constexpr double epsilon = 0.5;              // 1e-1;
constexpr double theta = 1 / 2.0 + epsilon;  // theta for MoreauJeanOSI integrator
// theta = 1.0;
constexpr double E = 210e9;     // young Modulus
constexpr double S = 0.000314;  //  Beam Section 1 cm  for the diameter
// S=0.1;
constexpr double L = 1.0;       // length of the  beam
constexpr double rho = 7800.0;  // specific mass
// rho=1.0;
// double g = 9.81;  // Gravity
constexpr double g = 0.0;
}  // namespace user_defined



namespace user_defined {

// User-defined main parameters
constexpr unsigned int nDof = 10;            // degrees of freedom for the beam
constexpr double t0 = 1e-8;                  // initial computation time
constexpr double T = 0.0015;                 // final computation time
constexpr double h = 1e-7;                   // time step
constexpr double position_init = 0.00005;    // initial position
constexpr double velocity_init = -.1;        // initial velocity
constexpr double epsilon = 0.0;              // 1e-1;
constexpr double theta = 1 / 2.0 + epsilon;  // theta for MoreauJeanOSI integrator
// theta = 1.0;
constexpr double E = 210e9;     // young Modulus
constexpr double S = 0.000314;  //  Beam Section 1 cm  for the diameter
// S=0.1;
constexpr double L = 1.0;       // length of the  beam
constexpr double rho = 7800.0;  // specific mass
// rho=1.0;
// constexpr double g = 9.81; // Gravity
constexpr double g = 0.0;
}  // namespace user_defined



