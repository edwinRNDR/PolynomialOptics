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
 * ===============================
 * OpticalElements/Propagation5.hh
 * ===============================
 * Free-space propagation up to degree 5. 
 * Defines function propagate_5.
 */


#ifndef propagate5_hh
#define propagate5_hh

#include "../TruncPoly/TruncPolySystem.hh"

// These functions do the actual work; C++ / Transform4f wrappers below
inline void propagate_5(PT4fData **polynomials, int *numTerms, float l);

#ifdef __cplusplus // C++ versions return Transform4f

#ifndef _make_tpt4f_
#define _make_tpt4f_

inline const Transform4f make_TPT4f(PT4fData **polynomials, int *num_terms, int trunc = 5) {
  Transform4f pt;
  pt.trunc_degree = trunc;
  for (int i = 0; i < 4; ++i) {
    Poly4f poly(num_terms[i],polynomials[i],sizeof(PT4fData));
    pt[i] = poly; 
    pt[i].trunc_degree = trunc;
  }
  pt.consolidate_terms();
  return pt;
}

inline const Poly4f make_TP4f(PT4fData *polynomial, int num_terms, int trunc = 5) {
  Poly4f poly(num_terms,polynomial,sizeof(PT4fData));
  poly.trunc_degree = trunc;
  poly.consolidate_terms();
  return poly;
}

#endif

inline const Transform4f propagate_5(float dist, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  propagate_5(element,num_terms,dist);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}


#endif

inline void propagate_5(PT4fData **polynomials, int *num_terms, float l) {
  // Maple output


  // How does x' depend?
  num_terms[0] = 7;

  polynomials[0][0] = PT4fData( 1,
				1,0,0,0);  // x coefficient
  polynomials[0][1] = PT4fData( l,
				0,0,1,0);  // dx coefficient
  polynomials[0][2] = PT4fData( l/2.f,
				0,0,3,0);  // dx^3 coefficient
  polynomials[0][3] = PT4fData( l/2.f,
				0,0,1,2);  // dx * dy^2 coefficient
  polynomials[0][4] = PT4fData( 3*l/8.f,
				0,0,1,4);  // dx * dy^4 coefficient
  polynomials[0][5] = PT4fData( 3*l/4.f,
				0,0,3,2);  // dx^3 * dy^2 coefficient
  polynomials[0][6] = PT4fData( 3*l/8.f,
				0,0,5,0);  // dx^5 coefficient


  // How does y' depend?
  num_terms[1] = 7;

  polynomials[1][0] = PT4fData( 1,
				0,1,0,0);  // y coefficient
  polynomials[1][1] = PT4fData( l,
				0,0,0,1);  // dy coefficient
  polynomials[1][2] = PT4fData( l/2.f,
				0,0,0,3);  // dy^3 coefficient
  polynomials[1][3] = PT4fData( l/2.f,
				0,0,2,1);  // dx^2 * dy coefficient
  polynomials[1][4] = PT4fData( 3*l/8.f,
				0,0,4,1);  // dx^4 * dy coefficient
  polynomials[1][5] = PT4fData( 3*l/4.f,
				0,0,2,3);  // dx^2 * dy^3 coefficient
  polynomials[1][6] = PT4fData( 3*l/8.f,
				0,0,0,5);  // dy^5 coefficient

  
  // How does dx' depend?
  num_terms[2] = 1;

  polynomials[2][0] = PT4fData( 1,
				0,0,1,0);  // dx remains unchanged


  // How does dy' depend?
  num_terms[3] = 1;

  polynomials[3][0] = PT4fData( 1,
				0,0,0,1);  // dy remains unchanged
}

#endif
