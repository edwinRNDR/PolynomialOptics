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
 * =============================
 * OpticalElements/Spherical5.hh
 * =============================
 * Provides polynomial system for spherical refraction and reflection
 * up to degree 5. Defines functions refract_spherical_5, 
 * reflect_spherical_5, and cos_angle_spherical_5.
 */


#ifndef spherical5_hh
#define spherical5_hh

#include "../TruncPoly/TruncPolySystem.hh"

// These functions do the actual work; C++ / Transform4f wrappers below
inline void refract_spherical_5(PT4fData **polynomials, int *num_terms, float R, float n1, float n2);
inline void reflect_spherical_5(PT4fData **polynomials, int *num_terms, float R);
inline void cos_angle_spherical_5(PT4fData *polynomial, int &num_terms, float R);


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


inline const Transform4f refract_spherical_5(float R, float n1, float n2, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  refract_spherical_5(element,num_terms,R,n1,n2);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

inline const Transform4f reflect_spherical_5(float R, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  reflect_spherical_5(element,num_terms,R);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

inline const Poly4f cos_angle_spherical_5(float R, int trunc = 5) {
  
  int num_terms;
  PT4fData *element = new PT4fData[1024];
  cos_angle_spherical_5(element,num_terms,R);
  Poly4f result = make_TP4f(element, num_terms, trunc);
  delete[] element;
  return result;
}

#endif


inline void refract_spherical_5(PT4fData **polynomials, int *num_terms, float R, float n1, float n2) {

  float R2 = R*R;
  float R3 = R2*R;
  float R4 = R2*R2;
  float R5 = R4*R;
  float n12 = n1*n1;
  float n13 = n12*n1;
  float n14 = n13*n1;
  float n22 = n2*n2;
  float n23 = n22*n2;
  float n24 = n23*n2;
  // Maple output     // cf. Spherical5.mw
  
  // How does x' depend?

  num_terms[0] = 0;

  // Linear (paraxial) term:
  float x1_xCoeff    = 1;
  polynomials[0][num_terms[0]++] = PT4fData( x1_xCoeff,
					    1,0,0,0);  // x coefficient

  // Degree-3 terms:
  float x1_x3Coeff   = 0.5*(n2-n1)/(R2*n2);
  float x1_xy2Coeff  = x1_x3Coeff;
  float x1_x2dxCoeff = 0.5*(n2-n1)/(R*n2);
  float x1_y2dxCoeff = x1_x2dxCoeff;
  polynomials[0][num_terms[0]++] = PT4fData( x1_x3Coeff,
					    3,0,0,0);  // x^3 coefficient
  polynomials[0][num_terms[0]++] = PT4fData( x1_xy2Coeff,
					    1,2,0,0);  // x * y^2 coefficient
  polynomials[0][num_terms[0]++] = PT4fData( x1_x2dxCoeff,
					    2,0,1,0);  // x^2 * dx coefficient
  polynomials[0][num_terms[0]++] = PT4fData( x1_y2dxCoeff,
					    0,2,1,0);  // y^2 * dx coefficient

  // Degree-5 terms:
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   0,2,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   0,2,3,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
					   0,3,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(-0.125*(-3*n23+5*n1*n22-4*n12*n2+2*n13)/(n23*R3),
					   0,4,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*n1*(-n22+n12)/(n23*R2),
					   1,2,0,2);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
					   1,2,2,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
					   1,3,0,1);
  polynomials[0][num_terms[0]++] = PT4fData(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   1,4,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   2,0,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   2,0,3,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
					   2,1,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-5*n23-6*n12*n2+4*n13+7*n1*n22)/(n23*R3),
					   2,2,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*n1*(-n22+n12)/(n23*R2),
					   3,0,0,2);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
					   3,0,2,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
					   3,1,0,1);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   3,2,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.125*(-7*n23+9*n1*n22-8*n12*n2+6*n13)/(n23*R3),
					   4,0,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   5,0,0,0);

  
  // How does y' depend?
  // Axial symmetry - just copy stuff from x'
  
  num_terms[1] = 5;

  // Linear term:
  float y1_yCoeff    = x1_xCoeff;
  polynomials[1][0] = PT4fData( y1_yCoeff,
				0,1,0,0);  // y coefficient

  // Degree-3 terms:
  float y1_y3Coeff   = x1_x3Coeff;
  float y1_x2yCoeff  = x1_x3Coeff;
  float y1_y2dyCoeff = x1_x2dxCoeff;
  float y1_x2dyCoeff = x1_x2dxCoeff;
  polynomials[1][1] = PT4fData( y1_y3Coeff,
				0,3,0,0);  // y^3 coefficient
  polynomials[1][2] = PT4fData( y1_x2yCoeff,
				2,1,0,0);  // x^2 * y coefficient
  polynomials[1][3] = PT4fData( y1_y2dyCoeff,
				0,2,0,1);  // y^2 * dy coefficient
  polynomials[1][4] = PT4fData( y1_x2dyCoeff,
				2,0,0,1);  // x^2 * dy coefficient

  // Degree-5 Terms:
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   2,0,2,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   2,0,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
					   3,0,1,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(-3*n23+5*n1*n22-4*n12*n2+2*n13)/(n23*R3),
					   4,0,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*n1*(-n22+n12)/(n23*R2),
					   2,1,2,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
					   2,1,0,2);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
					   3,1,1,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   4,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   0,2,2,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-n23+n13)/(n23*R),
					   0,2,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(-n23+n1*n22-n12*n2+n13)/(n23*R2),
					   1,2,1,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-5*n23-6*n12*n2+4*n13+7*n1*n22)/(n23*R3),
					   2,2,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*n1*(-n22+n12)/(n23*R2),
					   0,3,2,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-2*n23+n1*n22-2*n12*n2+3*n13)/(n23*R2),
					   0,3,0,2);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(-n2+n1)*(n12+n22)/(n23*R3),
					   1,3,1,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   2,3,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(-7*n23+9*n1*n22-8*n12*n2+6*n13)/(n23*R3),
					   0,4,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(-n2+n1)*(-2*n1*n2+3*n22+2*n12)/(n23*R4),
					   0,5,0,0);
  
  // How does dx' depend?
  
  num_terms[2] = 9;

  // Linear terms (paraxial):
  float dx1_dxCoeff   = n1/n2;
  float dx1_xCoeff    = (n1-n2)/(R*n2);
  polynomials[2][0] = PT4fData( dx1_dxCoeff,
				0,0,1,0);  // dx coefficient
  polynomials[2][1] = PT4fData( dx1_xCoeff,
				1,0,0,0);  // x coefficient

  // Degree-3 terms:
  float dx1_x3Coeff   = 0.5*n1*(n1-n2)/(n22*R3);
  float dx1_x2dxCoeff = 0.5*(n1-n2)*(2*n1+n2)/(n22*R2);
  float dx1_xy2Coeff  = dx1_x3Coeff;
  float dx1_xydyCoeff = n1*(n1-n2)/(n22*R2);
  float dx1_y2dxCoeff = 0.5*(n1-n2)/(R2*n2);
  float dx1_xdx2Coeff = 0.5*n1*(n1-n2)/(R*n22);
  float dx1_xdy2Coeff = dx1_xdx2Coeff;
  polynomials[2][2] = PT4fData( dx1_x3Coeff,
				3,0,0,0);  // x^3 coefficient
  polynomials[2][3] = PT4fData( dx1_x2dxCoeff,
				2,0,1,0);  // x^2 * dx coefficient
  polynomials[2][4] = PT4fData( dx1_xy2Coeff,
				1,2,0,0);  // x * y^2 coefficient
  polynomials[2][5] = PT4fData( dx1_xydyCoeff,
				1,1,0,1);  // x * y * dy coefficient
  polynomials[2][6] = PT4fData( dx1_y2dxCoeff,
				0,2,1,0);  // y^2 * dx coefficient
  polynomials[2][7] = PT4fData( dx1_xdx2Coeff,
				1,0,2,0);  // x * dx^2 coefficient
  polynomials[2][8] = PT4fData( dx1_xdy2Coeff,
				1,0,0,2);  // x * dy^2 coefficient


  // Degree-5 Terms:

  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n22+n12)/(n22*R2),
					   0,2,1,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n22+n12)/(n22*R2),
					   0,2,3,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*(-n22+n12)/(n22*R3),
					   0,3,1,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*(-n2+n1)*(2*n1+n2)/(n22*R4),
					   0,4,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*n1*(-n23+n13)/(R*n24),
					   1,0,0,4);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*n1*(-n23+n13)/(R*n24),
					   1,0,2,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*n1*(-n23+n13)/(R*n24),
					   1,0,4,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*n12*(-n22+n12)/(R2*n24),
					   1,1,0,3);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*n12*(-n22+n12)/(R2*n24),
					   1,1,2,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*n1*(-n23+3*n13-2*n1*n22)/(R3*n24),
					   1,2,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n1*n23+2*n12*n22+n14-2*n24)/(R3*n24),
					   1,2,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*n1*(-n23+n13)/(R4*n24),
					   1,3,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*n1*(-n23+n13)/(R5*n24),
					   1,4,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
					   2,0,1,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
					   2,0,3,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*(-n12*n22+2*n14-n24)/(R3*n24),
					   2,1,1,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-3*n1*n23+2*n12*n22+2*n14-n24)/(R4*n24),
					   2,2,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*n1*(-n23+n13)/(R3*n24),
					   3,0,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*(-n1*n23+3*n14-2*n24)/(R3*n24),
					   3,0,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.5*n1*(-n23+n13)/(R4*n24),
					   3,1,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(0.25*n1*(-n23+n13)/(R5*n24),
					   3,2,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*(-5*n1*n23+4*n14+2*n12*n22-n24)/(R4*n24),
					   4,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(0.125*n1*(-n23+n13)/(R5*n24),
					   5,0,0,0);


  // How does dy' depend?


  num_terms[3] = 9;

  // Linear terms (paraxial):
  float dy1_dyCoeff   = n1/n2;
  float dy1_yCoeff    = (n1-n2)/(R*n2);
  polynomials[3][0] = PT4fData( dy1_dyCoeff,
				0,0,0,1);  // dy coefficient
  polynomials[3][1] = PT4fData( dy1_yCoeff,
				0,1,0,0);  // y coefficient

  // Degree-3 terms:
  float dy1_y3Coeff   = 0.5*n1*(n1-n2)/(n22*R3);
  float dy1_y2dyCoeff = 0.5*(n1-n2)*(2*n1+n2)/(n22*R2);
  float dy1_x2yCoeff  = dy1_y3Coeff;
  float dy1_xydxCoeff = n1*(n1-n2)/(n22*R2);
  float dy1_x2dyCoeff = 0.5*(n1-n2)/(R2*n2);
  float dy1_ydy2Coeff = 0.5*n1*(n1-n2)/(R*n22);
  float dy1_ydx2Coeff = dy1_ydy2Coeff;
  polynomials[3][2] = PT4fData( dy1_y3Coeff,
				0,3,0,0);  // y^3 coefficient
  polynomials[3][3] = PT4fData( dy1_y2dyCoeff,
				0,2,0,1);  // y^2 * dy coefficient
  polynomials[3][4] = PT4fData( dy1_x2yCoeff,
				2,1,0,0);  // x^2 * y coefficient
  polynomials[3][5] = PT4fData( dy1_xydxCoeff,
				1,1,1,0);  // x * y * dx coefficient
  polynomials[3][6] = PT4fData( dy1_x2dyCoeff,
				2,0,0,1);  // x^2 * dy coefficient
  polynomials[3][7] = PT4fData( dy1_ydy2Coeff,
				0,1,0,2);  // y * dy^2 coefficient
  polynomials[3][8] = PT4fData( dy1_ydx2Coeff,
				0,1,2,0);  // y * dx^2 coefficient

  // Degree-5 Terms:

  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n22+n12)/(n22*R2),
					   2,0,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n22+n12)/(n22*R2),
					   2,0,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*(-n22+n12)/(n22*R3),
					   3,0,1,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*(-n2+n1)*(2*n1+n2)/(n22*R4),
					   4,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(-n23+n13)/(R*n24),
					   0,1,4,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(-n23+n13)/(R*n24),
					   0,1,2,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(-n23+n13)/(R*n24),
					   0,1,0,4);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n12*(-n22+n12)/(R2*n24),
					   1,1,3,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n12*(-n22+n12)/(R2*n24),
					   1,1,1,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(-n23+3*n13-2*n1*n22)/(R3*n24),
					   2,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n1*n23+2*n12*n22+n14-2*n24)/(R3*n24),
					   2,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n1*(-n23+n13)/(R4*n24),
					   3,1,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(-n23+n13)/(R5*n24),
					   4,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
					   0,2,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n12*n22+2*n14-n24)/(R2*n24),
					   0,2,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*(-n12*n22+2*n14-n24)/(R3*n24),
					   1,2,1,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-3*n1*n23+2*n12*n22+2*n14-n24)/(R4*n24),
					   2,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(-n23+n13)/(R3*n24),
					   0,3,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n1*n23+3*n14-2*n24)/(R3*n24),
					   0,3,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n1*(-n23+n13)/(R4*n24),
					   1,3,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(-n23+n13)/(R5*n24),
					   2,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*(-5*n1*n23+4*n14+2*n12*n22-n24)/(R4*n24),
					   0,4,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(-n23+n13)/(R5*n24),
					   0,5,0,0);
}

inline void reflect_spherical_5(PT4fData **polynomials, int *num_terms, float R) {

  float R2 = R*R;
  float R3 = R2*R;
  float R4 = R2*R2;
  float R5 = R4*R;

  // Maple output     // cf. Spherical5.mw


  // x dependency:
  num_terms[0] = 0;

  // 1 linear term, 4x degree-3 terms, 16x degree-5 terms:
  polynomials[0][num_terms[0]++] = PT4fData(1,
					   1,0,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/R,
					   0,2,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/R2,
					   1,2,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/R,
					   2,0,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/R2,
					   3,0,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(0.5/R,
					   0,2,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(0.5/R,
					   0,2,3,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R2,
					   0,3,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(7/(4*R3),
					   0,4,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R2,
					   1,2,2,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R3,
					   1,3,0,1);
  polynomials[0][num_terms[0]++] = PT4fData(7/(4*R4),
					   1,4,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/(2*R),
					   2,0,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(1/(2*R),
					   2,0,3,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R2,
					   2,1,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(11/(2*R3),
					   2,2,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R2,
					   3,0,2,0);
  polynomials[0][num_terms[0]++] = PT4fData(2/R3,
					   3,1,0,1);
  polynomials[0][num_terms[0]++] = PT4fData(7/(2*R4),
					   3,2,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(15/(4*R3),
					   4,0,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(7/(4*R4),
					   5,0,0,0);

  // y dependency:
  num_terms[1] = 0;

  // 1 linear term, 4x degree-3 terms, 16x degree-5 terms:
  polynomials[1][num_terms[1]++] = PT4fData(1,
					   0,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(1/R,
					   2,0,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(1/R2,
					   2,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(1/R,
					   0,2,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(1/R2,
					   0,3,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(0.5/R,
					   2,0,2,1);
  polynomials[1][num_terms[1]++] = PT4fData(0.5/R,
					   2,0,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(2/R2,
					   3,0,1,1);
  polynomials[1][num_terms[1]++] = PT4fData(7/(4*R3),
					   4,0,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(2/R2,
					   2,1,0,2);
  polynomials[1][num_terms[1]++] = PT4fData(2/R3,
					   3,1,1,0);
  polynomials[1][num_terms[1]++] = PT4fData(7/(4*R4),
					   4,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(1/(2*R),
					   0,2,2,1);
  polynomials[1][num_terms[1]++] = PT4fData(1/(2*R),
					   0,2,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(2/R2,
					   1,2,1,1);
  polynomials[1][num_terms[1]++] = PT4fData(11/(2*R3),
					   2,2,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(2/R2,
					   0,3,0,2);
  polynomials[1][num_terms[1]++] = PT4fData(2/R3,
					   1,3,1,0);
  polynomials[1][num_terms[1]++] = PT4fData(7/(2*R4),
					   2,3,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(15/(4*R3),
					   0,4,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(7/(4*R4),
					   0,5,0,0);

  // x' dependency:
  num_terms[2] = 0;

  polynomials[2][num_terms[2]++] = PT4fData(1,
					   0,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(2/R,
					   1,0,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(1/R2,
					   0,2,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R,
					   1,0,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R,
					   1,0,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(-2/R2,
					   1,1,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R3,
					   1,2,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R2,
					   2,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R3,
					   3,0,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(4*R4),
					   0,4,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(4*R),
					   1,0,0,4);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R),
					   1,0,2,2);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(4*R),
					   1,0,4,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R3),
					   1,2,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R3),
					   1,2,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R4,
					   1,3,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(4*R5),
					   1,4,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(-3/(2*R4),
					   2,2,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R3),
					   3,0,0,2);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R3),
					   3,0,2,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/R4,
					   3,1,0,1);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(2*R5),
					   3,2,0,0);
  polynomials[2][num_terms[2]++] = PT4fData(-5/(4*R4),
					   4,0,1,0);
  polynomials[2][num_terms[2]++] = PT4fData(-1/(4*R5),
					   5,0,0,0);

  // y' dependency:
  num_terms[3] = 0;

  polynomials[3][num_terms[3]++] = PT4fData(1,
					   0,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(2/R,
					   0,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(1/R2,
					   2,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R,
					   0,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R,
					   0,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(-2/R2,
					   1,1,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R3,
					   2,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R2,
					   0,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R3,
					   0,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R4),
					   4,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R),
					   0,1,4,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R),
					   0,1,2,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R),
					   0,1,0,4);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R3),
					   2,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R3),
					   2,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R4,
					   3,1,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R5),
					   4,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-3/(2*R4),
					   2,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R3),
					   0,3,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R3),
					   0,3,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R4,
					   1,3,1,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R5),
					   2,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-5/(4*R4),
					   0,4,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R5),
					   0,5,0,0);
}

inline void cos_angle_spherical_5(PT4fData *polynomial, int &num_terms, float R) {  

  float R2 = R*R;
  float R3 = R2*R;
  float R4 = R3*R;
  
  num_terms = 0;
  
  polynomial[num_terms++] = PT4fData(1         ,0,0,0,0);

  polynomial[num_terms++] = PT4fData(-0.5      ,0,0,0,2);
  polynomial[num_terms++] = PT4fData(-0.5      ,0,0,2,0);
  polynomial[num_terms++] = PT4fData(-1/R      ,0,1,0,1);
  polynomial[num_terms++] = PT4fData(-0.5/R2   ,0,2,0,0);
  polynomial[num_terms++] = PT4fData(-1/R      ,1,0,1,0);
  polynomial[num_terms++] = PT4fData(-0.5/R2   ,2,0,0,0);

  polynomial[num_terms++] = PT4fData(-0.125    ,0,0,0,4);
  polynomial[num_terms++] = PT4fData(-0.125    ,0,0,4,0);

  polynomial[num_terms++] = PT4fData(-0.25     ,0,0,2,2);

  polynomial[num_terms++] = PT4fData(-0.25/R4  ,2,2,0,0);

  polynomial[num_terms++] = PT4fData(-0.25/R2  ,0,2,0,2);
  polynomial[num_terms++] = PT4fData(-0.25/R2  ,0,2,2,0);
  polynomial[num_terms++] = PT4fData(-0.25/R2  ,2,0,0,2);
  polynomial[num_terms++] = PT4fData(-0.25/R2  ,2,0,2,0);
 
  polynomial[num_terms++] = PT4fData(-0.5/R3   ,0,3,0,1);
  polynomial[num_terms++] = PT4fData(-0.5/R3   ,3,0,1,0);

  polynomial[num_terms++] = PT4fData(-0.125/R4 ,0,4,0,0);
  polynomial[num_terms++] = PT4fData(-0.125/R4 ,4,0,0,0);

  polynomial[num_terms++] = PT4fData(-0.5/R3   ,1,2,1,0);
  polynomial[num_terms++] = PT4fData(-0.5/R3   ,2,1,0,1);
}









#endif
