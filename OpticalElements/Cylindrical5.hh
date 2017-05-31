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
 * OpticalElements/Cylindrical5.hh
 * ===============================
 * Provides polynomial system for cylindrical refraction and reflection
 * up to degree 5. Defines functions refract_cylindrical_{x|y}_5 
 * and reflect_cylindrical_{x|y}_5.
 */


#ifndef cylindrical5_hh
#define cylindrical5_hh

#include "../TruncPoly/TruncPolySystem.hh"

// These functions do the actual work; C++ / Transform4f wrappers below
inline void refract_cylindrical_x_5(PT4fData **polynomials, int *num_terms, float R, float n1, float n2);
inline void reflect_cylindrical_x_5(PT4fData **polynomials, int *num_terms, float R);
inline void refract_cylindrical_y_5(PT4fData **polynomials, int *num_terms, float R, float n1, float n2);
inline void reflect_cylindrical_y_5(PT4fData **polynomials, int *num_terms, float R);




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

inline const Transform4f refract_cylindrical_x_5(float R, float n1, float n2, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  refract_cylindrical_x_5(element,num_terms,R, n1, n2);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

inline const Transform4f refract_cylindrical_y_5(float R, float n1, float n2, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  refract_cylindrical_x_5(element,num_terms,R, n1, n2);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

inline const Transform4f reflect_cylindrical_x_5(float R, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  reflect_cylindrical_x_5(element,num_terms,R);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

inline const Transform4f reflect_cylindrical_y_5(float R, int trunc = 5) {
  
  int num_terms[4];
  PT4fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[1024];
  }
  reflect_cylindrical_x_5(element,num_terms,R);
  Transform4f result = make_TPT4f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}


#endif



inline void refract_cylindrical_x_5(PT4fData **polynomials, int *num_terms, float R, float n1, float n2){
  float R2 = R*R;
  float R3 = R2*R;
  float R4 = R2*R2;
  float R5 = R3*R2;
  float n22 = n2*n2;
  float n23 = n22*n2;
  float n24 = n23*n2;
  float n12 = n1*n1;
  float n13 = n12*n1;
  float n14 = n13*n1;

  // Maple output    

  // Dependencies of x:
  num_terms[0] = 0;
  polynomials[0][num_terms[0]++] = PT4fData(1,
					   1,0,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(-n2+n1)/(R*n2),
					   0,2,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(n13-n23)/(R*n23),
					   0,2,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(-0.25*(n13-n23)/(R*n23),
					   0,2,3,0);
  polynomials[0][num_terms[0]++] = PT4fData(-0.5*(n13-n2*n12+n1*n22-n23)/(R2*n23),
					   0,3,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(-0.125*(-R2*n23+n1*n22*R2-4*n12*R2*n2+2*n13*R2+2*n1*n22*R2)/(R3*n23*(R2)),
					   0,4,1,0);

  // Dependencies of y:
  num_terms[1] = 0;
  polynomials[1][num_terms[1]++] = PT4fData(1,
					   0,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(n1-n2)/(n2*R),
					   0,2,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.5*(n1-n2)/(n2*R*R),
					   0,3,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(n13-n23)/(n23*R),
					   0,2,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(n13-n23)/(n23*R),
					   0,2,2,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*(-2*n23*R-2*n12*R*n2-n1*R*n22+2*n1*R*n22+3*n13*R)/(R2*n23*R),
					   0,3,0,2);
  polynomials[1][num_terms[1]++] = PT4fData(-0.25*n1*(n12-n22)/(n23*R*R),
					   0,3,2,0);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(-8*n12*R2*n2+6*n13*R2+4*n1*R2*n22+(4*n22*R*R)*n1-4*n23*R*R-R2*n23-2*R2*n23+n1*R2*n22)/(R3*n23*(R2)),
					   0,4,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(-0.125*(n1-n2)*(-2*n1*n2*R2+2*n12*R2+2*n22*R2+n22*R2)/(R3*n23*(R3)),
					   0,5,0,0);

  // Dependencies of x':
  num_terms[2] = 0;
  polynomials[2][num_terms[2]++] = PT4fData(n1/n2, 
					   0,0,1,0);

  // Dependencies of y':
  num_terms[3] = 0;
  polynomials[3][num_terms[3]++] = PT4fData(n1/n2,
					   0,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData((n1-n2)/(n2*R),
					   0,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n1*(n1-n2)/(n22*R),
					   0,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n1*(n1-n2)/(n22*R),
					   0,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*(n1-n2)*(2*n1+n2)/(n22*R2),
					   0,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.5*n1*(n1-n2)/(n22*R3),
					   0,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(n13-n23)/(R*n24),
					   0,1,0,4);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(n13-n23)/(R*n24),
					   0,1,2,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(n13-n23)/(R*n24),
					   0,1,4,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(2*n14-n12*n22-n24)/(R2*n24),
					   0,2,0,3);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(2*n14-n12*n22-n24)/(R2*n24),
					   0,2,2,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*(-n1*n23+3*n14-2*n24)/(R3*n24),
					   0,3,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(0.25*n1*(n23+n13-2*n1*n22)/(R3*n24),
					   0,3,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*(4*n14-5*n1*n23+2*n12*n22-n24)/(R4*n24),
					   0,4,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(0.125*n1*(n13-n23)/(R5*n24),
					   0,5,0,0);
}


inline void reflect_cylindrical_x_5(PT4fData **polynomials, int *num_terms, float R){
  float R2 = R*R;
  float R3 = R2*R;

  // Maple output

  // Dependencies of x:
  num_terms[0] = 0;
  polynomials[0][num_terms[0]++] = PT4fData(1,
					   1,0,0,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/R,
					   0,2,1,0);
  polynomials[0][num_terms[0]++] = PT4fData(1/(2*R),
					   0,2,1,2);
  polynomials[0][num_terms[0]++] = PT4fData(1/(2*R),
					   0,2,3,0);
  polynomials[0][num_terms[0]++] = PT4fData((R+R)/(R2*R),
					   0,3,1,1);
  polynomials[0][num_terms[0]++] = PT4fData(0.25*(R2+4*R2)/(R3*(R2)),
					   0,4,1,0);

  // Dependencies of y:
  num_terms[1] = 0;
  polynomials[1][num_terms[1]++] = PT4fData(1,
					   0,1,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(1/R,
					   0,2,0,1);
  polynomials[1][num_terms[1]++] = PT4fData(1/R2,
					   0,3,0,0);
  polynomials[1][num_terms[1]++] = PT4fData(1/(2*R),
					   0,2,0,3);
  polynomials[1][num_terms[1]++] = PT4fData(1/(2*R),
					   0,2,2,1);
  polynomials[1][num_terms[1]++] = PT4fData((R+R)/(R2*R),
					   0,3,0,2);
  polynomials[1][num_terms[1]++] = PT4fData((1/4)*(R2+10*R2+4*R*R)/(R3*(R2)),
					   0,4,0,1);
  polynomials[1][num_terms[1]++] = PT4fData((1/4)*(6*R2+R2)/(R3*(R3)),
					   0,5,0,0);

  // Dependencies of x':
  num_terms[2] = 0;
  polynomials[2][num_terms[2]++] = PT4fData(1,
					   0,0,1,0);

  // Dependencies of y':
  num_terms[3] = 0;
  polynomials[3][num_terms[3]++] = PT4fData(1,
					   0,0,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(2/R,
					   0,1,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R,
					   0,1,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R,
					   0,1,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R2,
					   0,2,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/R3,
					   0,3,0,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R),
					   0,1,0,4);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R),
					   0,1,2,2);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R),
					   0,1,4,0);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(2*R3),
					   0,3,0,2);
  polynomials[3][num_terms[3]++] = PT4fData(1/(2*R3),
					   0,3,2,0);
  polynomials[3][num_terms[3]++] = PT4fData(-5/(4*R2*R2),
					   0,4,0,1);
  polynomials[3][num_terms[3]++] = PT4fData(-1/(4*R3*R2),
					   0,5,0,0);
}


inline void refract_cylindrical_y_5(PT4fData **polynomials, int *num_terms_rotated, float R, float n1, float n2){
  PT4fData *element[4];
  int num_terms[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[64];
  }

  refract_cylindrical_x_5(element, num_terms, R, n1, n2);

  // juggle X and Y coordinates
  for (int i = 0; i < 4; ++i) {
    num_terms_rotated[i] = num_terms[i ^ 0x01];
    for (int t = 0; t < num_terms_rotated[i]; ++t) {\
      polynomials[i][t].coeff = element[i ^ 0x01][t].coeff;
      polynomials[i][t].e0    = element[i ^ 0x01][t].e1;
      polynomials[i][t].e1    = element[i ^ 0x01][t].e0;
      polynomials[i][t].e2    = element[i ^ 0x01][t].e3;
      polynomials[i][t].e3    = element[i ^ 0x01][t].e2;
    }
  }
  for (int i = 0; i < 4; ++i) {
    delete element[i];
  }
}


inline void reflect_cylindrical_y_5(PT4fData **polynomials, int *num_terms_rotated, float R){
  PT4fData *element[4];
  int num_terms[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT4fData[64];
  }

  reflect_cylindrical_x_5(element, num_terms, R);

  // juggle X and Y coordinates
  for (int i = 0; i < 4; ++i) {
    num_terms_rotated[i] = num_terms[i ^ 0x01];
    for (int t = 0; t < num_terms_rotated[i]; ++t) {\
      polynomials[i][t].coeff = element[i ^ 0x01][t].coeff;
      polynomials[i][t].e0    = element[i ^ 0x01][t].e1;
      polynomials[i][t].e1    = element[i ^ 0x01][t].e0;
      polynomials[i][t].e2    = element[i ^ 0x01][t].e3;
      polynomials[i][t].e3    = element[i ^ 0x01][t].e2;
    }
  }
  for (int i = 0; i < 4; ++i) {
    delete element[i];
  }
}


#endif
