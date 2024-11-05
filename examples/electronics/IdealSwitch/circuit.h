#ifndef CIRCUIT_H
#define CIRCUIT_H

#define CLSC_CIRCUIT

// #define SICONOS_DEBUG

#ifdef CLSC_CIRCUIT
/*
 *An example from Paolo Maffezzoni nov 2006.
 *
 *
 */

namespace user_defined {

constexpr double sR = 1;
constexpr double sR1s = 0.001;
constexpr double sR1d = 0.001;
constexpr double sR2 = 1000;
constexpr double sL = 200e-6;
constexpr double sT = 20e-6;
constexpr double sStep = 0.1e-6;
constexpr double sE_plus = 7.5;
constexpr double sE_moins = -2.5;
constexpr double sTf = 10 * sT;
constexpr int sNSLawSize = 9;
constexpr int sN = 4;
constexpr int sM = 5;
constexpr double sAmpli = 100.0;

#else
/*    V1     /      V2
 * ___|_____/  _____|________
 * |         ^               |
 * |         |               |
 * |         V1              |
 * |                         |
 *_____                     ---
 *|u=  |                       C
 *|S(t)|                    ___
 *-----                      |
 * |                         |
 * |                         |
 * |                         |
 * |__________________________
 *
 *
 * x=V2
 * L=(V1,i,L3,L4)^t=(L1,L2,L3,L4)^t=
 *
 * C dx/dt = 0.x +r
 * r = g(L)=L2
 *
 *           |-L1+s(t)
 * Y=h(x,L)= |L1-X+(L4+R1)L2
 *           |R2-L4-R1
 *           |-L1+L3
 *
 *
 *
 */
constexpr int sNSLawSize = 4;
constexpr int sN = 2;
constexpr int sM = 2;

constexpr double sR1 = 1;
constexpr double sR2 = 1000;
constexpr double sC = 1e-2;
constexpr double sW = 1000;
constexpr double sTf = 0.07;
constexpr double sStep = 0.00001;

#endif
}

#endif
