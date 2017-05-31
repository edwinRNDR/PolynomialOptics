/* Polynomial Optics
 * (C) 2012, Matthias Hullin <hullin@cs.ubc.ca>
 * (University of British Columbia)
 *
 * All rights reserved. If you end up using this code 
 * or ideas contained therein, please cite our EGSR 2012 paper. 
 * More info:
 * 
 * http://www.cs.ubc.ca/labs/imager/tr/2012/PolynomialOptics/
 * 
 * ===========================
 * Example_BasicArithmetic.cpp 
 * ===========================
 * Example application introducing the Polynomial Optics library 
 * with its classes PolyTerm, TruncPoly, TruncPolySystem,
 * and the library of optical elements and materials.
 */
 
#include <TruncPoly/TruncPolySystem.hh>

#include <OpticalElements/OpticalMaterial.hh>
#include <OpticalElements/Spherical5.hh>
#include <OpticalElements/Cylindrical5.hh>
#include <OpticalElements/Propagation5.hh>
#include <OpticalElements/TwoPlane5.hh>

#include <OpticalElements/FindFocus.hh>

#include <iostream>

using namespace std;

int main() {

  cout << endl 
       << "PART 1: Simple polynomial operations" << endl 
       << "====================================" << endl << endl;
  
  Poly2f sum = Term2f(1.f, 1, 0) + Term2f(1.f, 0, 1); // x0' = x0 + x1
  Poly1f halfsquare = Term1f(0.5f, 2); // x0' = 0.5 * x0^2
  
  cout << "sum = " << sum << endl << endl; // prints "1 * x1 + 1 * x0"

  // Concatenate polynomials:
  
  Poly2f halfsquare_of_sum = sum >> halfsquare;
  cout << "halfsquare_of_sum = " << halfsquare_of_sum << endl << endl; // prints "0.5 * x1^2 + 1 * x0 * x1 + 0.5 * x0^2"

  // Derivative with respect to x0
  Poly2f derivative = halfsquare_of_sum.get_derivative(0);
  
  float input[2] = {2.f, 1.f};
  cout << "Derivative w.r.t. x0, evaluated at {2.f,1.f}: " 
       << derivative << " = " << derivative.evaluate(input) << endl << endl;
  
  // Take the 5th power of (x0 + x1 + 1):  
  Poly2f sum5 = (Term2f(1.f, 1, 0) + Term2f(1.f, 0, 1) + 1) ^ 5;
  
  cout << "(x0 + x1 + 1)^5 = " << sum5 << endl << endl; // prints a lengthy polynomial
  
  cout << "After truncation to degree 3: " << (sum5 % 3) << endl << endl;


  cout << "PART 2: Optical systems" << endl 
       << "=======================" << endl << endl;

  // Retrieve glass from database by name
  OpticalMaterial glass1("BK7", true);
  
  // Retrieve glass from database by refractive index and Abbe number
  OpticalMaterial glass2(1.6, 47, true);
  
  // Construct a simple composite lens from two materials
  
  const int degree = 3; // Truncate behind cubic terms
  const float lambda = 550; 
  
  const float d0 = 5000;
  const float R1 = 65.22;
  const float d1 = 9.60;
  const float R2 = -62.03;
  const float d2 = 4.20;
  const float R3 = -1240.67;

  // Construct primitive optical elements and combine them to a system:
  System44f lens =
  	two_plane_5(d0, degree) // Two-plane parametrization of rays for distance d0 between planes
  	>> refract_spherical_5(R1, 1.f, glass1.get_index(lambda), degree)
    >> propagate_5(d1, degree)
    >> refract_spherical_5(R2, glass1.get_index(lambda), glass2.get_index(lambda), degree)
    >> propagate_5(d2, degree)
    >> refract_spherical_5(R3, glass2.get_index(lambda), 1.f, degree);

  cout << endl << lens << endl << endl;
  lens.print_stats(cout);
  
  // Obtain back focal length
  float bfl = find_focus_X(lens);
  cout << "Back focal length: " << bfl << endl;
  
  // Integrate focus into total optical system
  lens = lens >> propagate_5(bfl);
  
  // Output linear terms of resulting optical system. 
  // Note that x0', x1' don't depend on x2, x3, 
  // the coordinates in the entrance pupil plane.
  cout << "Linear system for lens: " << (lens%1) << endl << endl;



  cout << "PART 3: Fancy stuff (no output, see source)" << endl 
       << "===========================================" << endl << endl;

  System44f spectral_lens[3];
  
  float lambdas[3] = {450.f, 550.f ,650.f};
  
  for (int i = 0; i < 3; ++i) {
    
    float lambda = lambdas[i];

    spectral_lens[i] =
      two_plane_5(d0, degree) // Two-plane parametrization of rays for distance d0 between planes
      >> refract_spherical_5(R1, 1.f, glass1.get_index(lambda), degree)
      >> propagate_5(d1, degree)
      >> refract_spherical_5(R2, glass1.get_index(lambda), glass2.get_index(lambda), degree)
      >> propagate_5(d2, degree)
      >> refract_spherical_5(R3, glass2.get_index(lambda), 1.f, degree);
  }

  // Quadratic interpolation into a new system. 
  // (Experimental implementation; requires same term ordering in all three systems)
  System54f querped_spectral_lens 
    = spectral_lens[0].querp_with(spectral_lens[1],spectral_lens[2],lambdas[0],lambdas[1],lambdas[2]);


  // Get Jacobian coefficients for a certain ray:
  float jacobian[20]; // 4x5 matrix
  float spectral_ray[5] = {1.0, 2.0, 0.1, 0.1, 500.0}; // rx, ry, dx, dy, lambda
  querped_spectral_lens.get_jacobian(jacobian, spectral_ray);

  // In order to get Jacobian in analytical form, use (TruncPoly|TruncPolySystem)::get_derivative(int wrt);
  
}
