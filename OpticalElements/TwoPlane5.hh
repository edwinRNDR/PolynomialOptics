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
 * OpticalElements/TwoPlane5.hh
 * ============================
 * Two-plane ray parametrization. Defines function two_plane_5 
 * as 4in/4out system of degree 5, mapping world and entrance pupil 
 * positions (xw yw xp yp) to ray (xp yp dx dy).
 */


#ifndef TwoPlane_hh
#define TwoPlane_hh

#include "../TruncPoly/TruncPolySystem.hh"
// These functions do the actual work; C++ / Transform4f wrappers below
inline void two_plane_5(PT4fData **polynomials, int *num_terms, float dz);


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


inline const Transform4f two_plane_5(float dz, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  two_plane_5(element,num_terms,dz);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

#endif


inline void two_plane_5(PT4fData **polynomials, int *num_terms, float dz) {
  // Aperture coordinates are the new ray position
  num_terms[0] = 1;
  polynomials[0][0] = PT4fData(1,0,0,1,0);
  num_terms[1] = 1;
  polynomials[1][0] = PT4fData(1,0,0,0,1);

  // Obtain direction components a la (x_ap-x_world) / dist:

#define dz3 (dz*dz*dz)
#define dz5 (dz*dz*dz*dz*dz)

  num_terms[2] = 0;
  polynomials[2][num_terms[2]++] = PT4fData(1/dz,
					    0,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/dz,
					    1,0,0,0);

  polynomials[2][num_terms[2]++] = PT4fData(-0.5/dz3  ,0,0,1,2);
  polynomials[2][num_terms[2]++] = PT4fData(-0.5/dz3  ,0,0,3,0);
  polynomials[2][num_terms[2]++] = PT4fData(1/dz3     ,0,1,1,1);
  polynomials[2][num_terms[2]++] = PT4fData(-0.5/dz3  ,0,2,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5/dz3   ,1,0,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(1.5/dz3   ,1,0,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/dz3    ,1,1,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.5/dz3   ,1,2,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1.5/dz3  ,2,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5/dz3   ,3,0,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.375/dz5 ,0,0,1,4);
  polynomials[2][num_terms[2]++] = PT4fData(0.75/dz5  ,0,0,3,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.375/dz5 ,0,0,5,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1.5/dz5  ,0,1,1,3);
  polynomials[2][num_terms[2]++] = PT4fData(-1.5/dz5  ,0,1,3,1);
  polynomials[2][num_terms[2]++] = PT4fData(2.25/dz5  ,0,2,1,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.75/dz5  ,0,2,3,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1.5/dz5  ,0,3,1,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.375/dz5 ,0,4,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-0.375/dz5,1,0,0,4);
  polynomials[2][num_terms[2]++] = PT4fData(-2.25/dz5 ,1,0,2,2);
  polynomials[2][num_terms[2]++] = PT4fData(-1.875/dz5,1,0,4,0);
  polynomials[2][num_terms[2]++] = PT4fData(1.5/dz5   ,1,1,0,3);
  polynomials[2][num_terms[2]++] = PT4fData(4.5/dz5   ,1,1,2,1);
  polynomials[2][num_terms[2]++] = PT4fData(-2.25/dz5 ,1,2,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(-2.25/dz5 ,1,2,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(1.5/dz5   ,1,3,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(-0.375/dz5,1,4,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(2.25/dz5  ,2,0,1,2);
  polynomials[2][num_terms[2]++] = PT4fData(3.75/dz5  ,2,0,3,0);
  polynomials[2][num_terms[2]++] = PT4fData(-4.5/dz5  ,2,1,1,1);
  polynomials[2][num_terms[2]++] = PT4fData(2.25/dz5  ,2,2,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-0.75/dz5 ,3,0,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(-3.75/dz5 ,3,0,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(1.5/dz5   ,3,1,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(-0.75/dz5 ,3,2,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(1.875/dz5 ,4,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-0.375/dz5,5,0,0,0);


  num_terms[3] = 0;
  polynomials[3][num_terms[3]++] = PT4fData(1/dz,
					    0,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/dz,
					    0,1,0,0);

  polynomials[3][num_terms[3]++] = PT4fData(-0.5/dz3  ,0,0,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(-0.5/dz3  ,0,0,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(1/dz3     ,1,0,1,1);
  polynomials[3][num_terms[3]++] = PT4fData(-0.5/dz3  ,2,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.5/dz3   ,0,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(1.5/dz3   ,0,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/dz3    ,1,1,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.5/dz3   ,2,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1.5/dz3  ,0,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.5/dz3   ,0,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.375/dz5 ,0,0,4,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.75/dz5  ,0,0,2,3);
  polynomials[3][num_terms[3]++] = PT4fData(0.375/dz5 ,0,0,0,5);
  polynomials[3][num_terms[3]++] = PT4fData(-1.5/dz5  ,1,0,3,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1.5/dz5  ,1,0,1,3);
  polynomials[3][num_terms[3]++] = PT4fData(2.25/dz5  ,2,0,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.75/dz5  ,2,0,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(-1.5/dz5  ,3,0,1,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.375/dz5 ,4,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-0.375/dz5,0,1,4,0);
  polynomials[3][num_terms[3]++] = PT4fData(-2.25/dz5 ,0,1,2,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1.875/dz5,0,1,0,4);
  polynomials[3][num_terms[3]++] = PT4fData(1.5/dz5   ,1,1,3,0);
  polynomials[3][num_terms[3]++] = PT4fData(4.5/dz5   ,1,1,1,2);
  polynomials[3][num_terms[3]++] = PT4fData(-2.25/dz5 ,2,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-2.25/dz5 ,2,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(1.5/dz5   ,3,1,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(-0.375/dz5,4,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(2.25/dz5  ,0,2,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(3.75/dz5  ,0,2,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(-4.5/dz5  ,1,2,1,1);
  polynomials[3][num_terms[3]++] = PT4fData(2.25/dz5  ,2,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-0.75/dz5 ,0,3,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-3.75/dz5 ,0,3,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(1.5/dz5   ,1,3,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(-0.75/dz5 ,2,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(1.875/dz5 ,0,4,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-0.375/dz5,0,5,0,0);

}

#endif
