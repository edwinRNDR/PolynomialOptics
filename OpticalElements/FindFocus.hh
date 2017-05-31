/* Polynomial Optics
 * (C) 2012, Matthias Hullin <hullin@cs.ubc.ca>
 * (University of British Columbia)
 *
 * Feel free to do what you want with this code, as long as this 
 * header remains intact. Needless to say, we appreciate proper 
 * attribution, be it through citation of our EGSR 2012 paper, movie 
 * credits, or otherwise.
 *
 * http://www.cs.ubc.ca/labs/imager/tr/2012/PolynomialOptics/
 * 
 * ============================
 * OpticalElements/FindFocus.hh
 * ============================
 * Find the back focal length of a given system from 1st-order terms
 * (matrix optics). That is, we want to find the distance at which
 * film location (rx', ry') becomes independent of ray direction (dx,
 * dy). 
 */

#ifndef findfocus_hh
#define findfocus_hh

#include "../TruncPoly/TruncPolySystem.hh"

template <typename scalar>
inline scalar find_focus_X(TruncPolySystem<scalar, 4, 4> system) {

  // Find terms that are linear in dx and nothing else
  scalar b = system[0].get_coeff(0,0,1,0);
  scalar d = system[2].get_coeff(0,0,1,0);

  return -b / d;
}

template <typename scalar>
inline scalar find_focus_Y(TruncPolySystem<scalar, 4, 4> system) {

  // Find terms that are linear in dy and nothing else
  scalar b = system[1].get_coeff(0,0,0,1);
  scalar d = system[3].get_coeff(0,0,0,1);

  return -b / d;
}

template <typename scalar>
inline scalar get_magnification_X(TruncPolySystem<scalar, 4, 4> system) {

  return system[0].get_coeff(1,0,0,0);
}

template <typename scalar>
inline scalar get_magnification_Y(TruncPolySystem<scalar, 4, 4> system) {

  return system[1].get_coeff(0,1,0,0);
}

template <typename scalar>
inline scalar get_stability_param_X(TruncPolySystem<scalar, 4, 4> system) {
  scalar a = system[0].get_coeff(1,0,0,0);
  scalar b = system[0].get_coeff(0,0,1,0);
  scalar c = system[2].get_coeff(1,0,0,0);
  scalar d = system[2].get_coeff(0,0,1,0);

  return (a+d)/2;
}



#endif
