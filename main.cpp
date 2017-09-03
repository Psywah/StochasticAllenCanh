// Copyright (C) 2005-2007 Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Anders Logg 2011
#include <iostream>
#include <fstream>
#include <dolfin.h>
#include <math.h>
#include "AllenCahn2D.h"

#define EPS 0.01
#define DT 0.0001
#define SAVE_DT 1*DT
#define END_T 5*DT
#define SIGMA 0.1

#define N 256

using namespace dolfin;

// Initial conditions
class InitialConditions : public Expression
{
public:
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double d = sqrt((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5)) - 0.25;
   // values[0] = tanh(d/sqrt(2)/EPS);
   values[0] = tanh(d/(sqrt(2)*EPS));
  }
};


//class InitialConditions : public Expression
//{
//public:
//  InitialConditions() { dolfin::seed(200); }
//  void eval(Array<double>& values, const Array<double>& x) const
//  {
//    values[0] = 2.0*dolfin::rand() - 1.0; 
//
//  }
//};


/**
 * @brief main function
 *
 * @param argc
 * @param argv[]
 *
 * @return 
 */
int main()
{
  //init(argc, argv);

  // Mesh
  //RectangleMesh mesh(-1, -1, 1, 1, 256, 256);
//   Point p0(-1, -1);
//   Point p1(1, 1);
//   RectangleMesh mesh(p0, p1, N, N);
  auto mesh = std::make_shared<UnitSquareMesh>(32, 32);   
  //UnitSquareMesh mesh(32,32);
  // Time stepping and model parameters
  // Constant epsilon(EPS);
  // Constant dt(DT);
  
  
  auto epsilon = std::make_shared<Constant>(EPS);
  auto sigma = std::make_shared<Constant>(SIGMA);
  //auto dw = std::make_shared<Constant>(dolfin::rand() - .5);
  auto dw = std::make_shared<Constant>(0.0);
  auto dt = std::make_shared<Constant>(DT); 

  
  double t = 0.0;
  double T = END_T;

  // Create user-defined function space 
  // AllenCahn2D::FunctionSpace V(mesh);
  // Function phi(V);
  // Function phi0(V);
  auto V = std::make_shared<AllenCahn2D::FunctionSpace>(mesh); 
  auto phi = std::make_shared<Function>(V);
  auto phi0 = std::make_shared<Function>(V);
  

  // Initial value 
  InitialConditions u_initial;
  phi0->interpolate(u_initial);
  // Bilinear form 
  AllenCahn2D::BilinearForm a(V, V);
  a.dt = dt;
  AllenCahn2D::LinearForm L(V); 
  L.phi0 = phi0; L.eps = epsilon; L.dt = dt; L.dw=dw; L.sigma= sigma;

  
  // Save initial condition to file
  File file("./allen_cahn.pvd");
  //  file << phi0;
  file << *phi0;
  std::ofstream sol_phi;
  sol_phi.open("./phi.txt",std::ios_base::app);
  for (int i=0; i<phi0->vector()->size(); ++i) {
      sol_phi << (*phi0->vector())[i] << std::endl;}
  // Solve
  int n_save = 1;
  while (t < T-DOLFIN_EPS)
  {
      // Update for next time step
      solve(a == L, *phi);
      // update time 
      t += DT;
      *phi0->vector() = *phi->vector();
      std::cout << "t = " << t << std::endl;
      // if (t > n_save*SAVE_DT) 
      if ((t - n_save*SAVE_DT<1e-5)&&(n_save*SAVE_DT-t<1e-5)){
          //      file << phi0;
          file << *phi0;
          n_save ++;
          std::cout << "n_save = " << n_save  << std::endl;
          for (int i=0; i<phi0->vector()->size(); ++i) {
              sol_phi  << (*phi0->vector())[i] << std::endl;
          }
      }
  }
  sol_phi.close();
  //  file << phi0;
  file << *phi0;
  // Plot solution
  plot(phi0);

  //plot(mesh);
  interactive();
  
  return 0;
}
