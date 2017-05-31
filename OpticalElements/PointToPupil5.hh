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
 * ================================
 * OpticalElements/PointToPupil5.hh
 * ================================
 * Point-to-pupil mapping. Defines function point_to_pupil_5 
 * as 2in/4out system, mapping pupil position (xy) to ray (xydxdy).
 */


#ifndef p2p5_hh
#define p2p5_hh

#include "../TruncPoly/TruncPolySystem.hh"

// These functions do the actual work; C++ / Transform4f wrappers below
inline void point_to_pupil_5(PT2fData **polynomials, int *num_terms, float x0, float y0, float z0);



#ifdef __cplusplus // C++ versions return Transform4f

#ifndef _make_system24f_
#define _make_system24f_

inline const System24f make_System24f(PT2fData **polynomials, int *num_terms, int trunc = 5) {
  System24f pt;
  pt.trunc_degree = trunc;
  for (int i = 0; i < 4; ++i) {
    Poly2f poly(num_terms[i],polynomials[i],sizeof(PT2fData));
    pt[i] = poly; 
    pt[i].trunc_degree = trunc;
  }
  pt.consolidate_terms();
  return pt;
}

#endif 

inline const System24f point_to_pupil_5(float x0, float y0, float z0, int trunc = 5) {
  
  int num_terms[4];
  PT2fData *element[4];
  for (int i = 0; i < 4; ++i) {
    element[i] = new PT2fData[64];
  }
  point_to_aperture_5(element,num_terms,x0,y0,z0);
  System24f result = make_System24f(element, num_terms, trunc);
  for (int i = 0; i < 4; ++i) {
    delete[] element[i];
  }
  return result;
}

#endif // __cplusplus



inline void point_to_pupil_5(PT2fData **polynomials, int *num_terms, float x0, float y0, float z0) {

  // How does x' depend?
  num_terms[0] = 1;

  polynomials[0][0] = PT2fData( 1,
				1,0);  // x coefficient


  // How does y' depend?
  num_terms[1] = 1;

  polynomials[1][0] = PT2fData( 1,
				0,1);  // y coefficient

  
  float x02 = x0*x0;
  float y02 = y0*y0;
  float z02 = z0*z0;
  float x04 = x02*x02;
  float y04 = y02*y02;
  float z04 = z02*z02;
  float x06 = x02*x04;
  float y06 = y02*y04;
  float z06 = z02*z04;

  float R2 = x02+y02+z02;
  float R = sqrt(R2);
  float R3 = R2 * R;
  float R5 = R3 * R2;
  float R7 = R5 * R2;
  float R9 = R7 * R2;
  float R11 = R9 * R2;

  // How does dx' depend?  dx' = (x-x0)/R
  num_terms[2] = 0;

  polynomials[2][num_terms[2]++] = PT2fData( -x0/R,
					    0,0); 
  polynomials[2][num_terms[2]++] = PT2fData(-x0*y0/R3,
					   0,1);
  polynomials[2][num_terms[2]++] = PT2fData((y02+z02)/R3,
					   1,0);
  polynomials[2][num_terms[2]++] = PT2fData(0.5*x0*(-2*y02+x02+z02)/R5,
					   0,2);
  polynomials[2][num_terms[2]++] = PT2fData(y0*(-2*x02+y02+z02)/R5,
					   1,1);
  polynomials[2][num_terms[2]++] = PT2fData(1.5*x0*(y02+z02)/R5,
					   2,0);
  polynomials[2][num_terms[2]++] = PT2fData(0.5*x0*y0*(-2*y02+3*x02+3*z02)/R7,
					   0,3);
  polynomials[2][num_terms[2]++] = PT2fData(-0.5*(11*x02*y02-2*y04-y02*z02-2*x04-x02*z02+z04)/R7,
					   1,2);
  polynomials[2][num_terms[2]++] = PT2fData(1.5*x0*y0*(-2*x02+3*y02+3*z02)/R7,
					   2,1);
  polynomials[2][num_terms[2]++] = PT2fData(-0.5*(-4*x02*y02-4*x02*z02+y04+2*y02*z02+z04)/R7,
					   3,0);
  polynomials[2][num_terms[2]++] = PT2fData(-0.125*x0*(8*y04-24*x02*y02-24*y02*z02+3*x04+6*x02*z02+3*z04)/R9,
					   0,4);
  polynomials[2][num_terms[2]++] = PT2fData(-0.5*y0*(21*x02*y02-2*y04+y02*z02-12*x04-9*x02*z02+3*z04)/R9,
					   1,3);
  polynomials[2][num_terms[2]++] = PT2fData(-0.75*x0*(21*x02*y02-12*y04-9*y02*z02-2*x04+x02*z02+3*z04)/R9,
					   2,2);
  polynomials[2][num_terms[2]++] = PT2fData(-0.5*y0*(8*x04-24*x02*y02-24*x02*z02+3*y04+6*y02*z02+3*z04)/R9,
					   3,1);
  polynomials[2][num_terms[2]++] = PT2fData(-0.625*x0*(-4*x02*y02-4*x02*z02+3*y04+6*y02*z02+3*z04)/R9,
					   4,0);
  polynomials[2][num_terms[2]++] = PT2fData(-0.125*x0*y0*(8*y04-40*x02*y02-40*y02*z02+15*x04+30*x02*z02+15*z04)/R11,
					   0,5);
  polynomials[2][num_terms[2]++] = PT2fData(0.125*(-16*y04*z02+159*y02*x04-21*y02*z04-21*x04*z02-6*x02*z04-136*x02*y04+138*x02*y02*z02+8*y06-12*x06+3*z06)/R11,
					   1,4);
  polynomials[2][num_terms[2]++] = PT2fData(-3.75*x0*y0*(13*x02*y02-4*y04-y02*z02-4*x04-x02*z02+3*z04)/R11,
					   2,3);
  polynomials[2][num_terms[2]++] = PT2fData(0.25*(-21*y04*z02-136*y02*x04-6*y02*z04-16*x04*z02-21*x02*z04+159*x02*y04+138*x02*y02*z02-12*y06+8*x06+3*z06)/R11,
					   3,2);
  polynomials[2][num_terms[2]++] = PT2fData(-0.625*x0*y0*(8*x04-40*x02*y02-40*x02*z02+15*y04+30*y02*z02+15*z04)/R11,
					   4,1);
  polynomials[2][num_terms[2]++] = PT2fData(0.375*(8*y02*x04+8*x04*z02-12*x02*y04-24*x02*y02*z02-12*x02*z04+y06+3*y04*z02+3*y02*z04+z06)/R11,
					   5,0);



  // How does dy' depend?   dy' = (y-y0)/R
  num_terms[3] = 0;

  polynomials[3][num_terms[3]++] = PT2fData(-y0/R,
					   0,0); 
  polynomials[3][num_terms[3]++] = PT2fData(-y0*x0/R3,
					   1,0);
  polynomials[3][num_terms[3]++] = PT2fData((x02+z02)/R3,
					   0,1);
  polynomials[3][num_terms[3]++] = PT2fData(0.5*y0*(-2*x02+y02+z02)/R5,
					   2,0);
  polynomials[3][num_terms[3]++] = PT2fData(x0*(-2*y02+x02+z02)/R5,
					   1,1);
  polynomials[3][num_terms[3]++] = PT2fData(1.5*y0*(x02+z02)/R5,
					   0,2);
  polynomials[3][num_terms[3]++] = PT2fData(0.5*y0*x0*(-2*x02+3*y02+3*z02)/R7,
					   3,0);
  polynomials[3][num_terms[3]++] = PT2fData(-0.5*(11*y02*x02-2*x04-x02*z02-2*y04-y02*z02+z04)/R7,
					   2,1);
  polynomials[3][num_terms[3]++] = PT2fData(1.5*y0*x0*(-2*y02+3*x02+3*z02)/R7,
					   1,2);
  polynomials[3][num_terms[3]++] = PT2fData(-0.5*(-4*y02*x02-4*y02*z02+x04+2*x02*z02+z04)/R7,
					   0,3);
  polynomials[3][num_terms[3]++] = PT2fData(-0.125*y0*(8*x04-24*y02*x02-24*x02*z02+3*y04+6*y02*z02+3*z04)/R9,
					   4,0);
  polynomials[3][num_terms[3]++] = PT2fData(-0.5*x0*(21*y02*x02-2*x04+x02*z02-12*y04-9*y02*z02+3*z04)/R9,
					   3,1);
  polynomials[3][num_terms[3]++] = PT2fData(-0.75*y0*(21*y02*x02-12*x04-9*x02*z02-2*y04+y02*z02+3*z04)/R9,
					   2,2);
  polynomials[3][num_terms[3]++] = PT2fData(-0.5*x0*(8*y04-24*y02*x02-24*y02*z02+3*x04+6*x02*z02+3*z04)/R9,
					   1,3);
  polynomials[3][num_terms[3]++] = PT2fData(-0.625*y0*(-4*y02*x02-4*y02*z02+3*x04+6*x02*z02+3*z04)/R9,
					   0,4);
  polynomials[3][num_terms[3]++] = PT2fData(-0.125*y0*x0*(8*x04-40*y02*x02-40*x02*z02+15*y04+30*y02*z02+15*z04)/R11,
					   5,0);
  polynomials[3][num_terms[3]++] = PT2fData(0.125*(-16*x04*z02+159*x02*y04-21*x02*z04-21*y04*z02-6*y02*z04-136*y02*x04+138*y02*x02*z02+8*x06-12*y06+3*z06)/R11,
					   4,1);
  polynomials[3][num_terms[3]++] = PT2fData(-3.75*y0*x0*(13*y02*x02-4*x04-x02*z02-4*y04-y02*z02+3*z04)/R11,
					   3,2);
  polynomials[3][num_terms[3]++] = PT2fData(0.25*(-21*x04*z02-136*x02*y04-6*x02*z04-16*y04*z02-21*y02*z04+159*y02*x04+138*y02*x02*z02-12*x06+8*y06+3*z06)/R11,
					   2,3);
  polynomials[3][num_terms[3]++] = PT2fData(-0.625*y0*x0*(8*y04-40*y02*x02-40*y02*z02+15*x04+30*x02*z02+15*z04)/R11,
					   1,4);
  polynomials[3][num_terms[3]++] = PT2fData(0.375*(8*x02*y04+8*y04*z02-12*y02*x04-24*y02*x02*z02-12*y02*z04+x06+3*x04*z02+3*x02*z04+z06)/R11,
					   0,5);
};



#endif
