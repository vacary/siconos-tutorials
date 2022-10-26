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
#ifdef _WIN32
#define SICONOS_EXPORT extern "C" __declspec(dllexport)
#else
#define SICONOS_EXPORT extern "C"
#endif
#include <stdio.h>

extern "C" double OneCubicFextFunction(double time)
{
  double t_start;
  double res = 10000.*time;
  if (time < 0.05)
  {
    res = 0.25e+00;
  }
  else if ((time >= 0.05) and (time < 0.1)) 
  {
    res = 0.5e+00;
  }
  else if ((time >= 0.1) and (time < 0.2)) 
  {
    res = 1e+00;
  }
  else if ((time >= 0.2) and (time < 0.21)) 
  {
    res = 1.5e+00;
  }
  else if ((time >= 0.21) and (time < 0.33)) 
  {
    res = -1e+00;
  }
  else
    res = 1.5e+00;
  return 2*res;
}


SICONOS_EXPORT void OneCubicFExt(double time, unsigned int sizeOfq, double *fExt, unsigned int sizeZ, double* z)
{
  for(unsigned int i = 0; i < sizeOfq; i++)
    fExt[i] = 0.0;
  fExt[0] =  OneCubicFextFunction(time);
}
