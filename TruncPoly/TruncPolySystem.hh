/* Polynomial Optics
 * (C) 2012, Matthias Hullin <hullin@cs.ubc.ca>
 * (University of British Columbia)
 * (C) 2012, Johannes Hanika <hanatos@gmail.com>
 * (Weta Digital)
 *
 * Feel free to do what you want with this code, as long as this 
 * header remains intact. Needless to say, we appreciate proper 
 * attribution, be it through citation of our EGSR 2012 paper, movie 
 * credits, or otherwise.
 *
 * More info:
 * 
 * http://www.cs.ubc.ca/labs/imager/tr/2012/PolynomialOptics/
 * 
 * ============================
 * TruncPoly/TruncPolySystem.hh
 * ============================
 * Defines classes polyTerm, truncPoly, truncPolySystem (C++ only)
 * as well as data container structs.
 */

#ifndef TruncPolySystem_hh
#define TruncPolySystem_hh

// Data type for powers. x^255 should be enough for all reasonable purposes. 
// Gives nice 8-byte packing for 4-variable polynomial terms.
#define echar unsigned char    

#ifdef __cplusplus // Templated class definitions requires C++

#include <vector>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdio>

// Abbreviations to reduce visual clutter of templates
#define _tS template <typename scalar>
#define _S <scalar>
#define _tO template <typename other_scalar>
#define _O <other_scalar>
#define _tSN template <typename scalar, int num_vars>
#define _SN <scalar, num_vars>
#define _ON <other_scalar, num_vars>
#define _tSIO template <typename scalar, int num_vars_in, int num_vars_out>
#define _SIO <scalar, num_vars_in, num_vars_out>
#define _OIO <other_scalar, num_vars_in, num_vars_out>
#define _O1O <other_scalar, (num_vars_in-1), num_vars_out>

#define _sc scalar
#define _nv num_vars
#define _nvi num_vars_in
#define _nvo num_vars_out

using namespace std;


_tSN  class PolyTerm;
_tSN  class TruncPoly;
_tSIO class TruncPolySystem;


_tS
inline scalar intpow(const scalar base, const echar exponent) {
  if (exponent == 0) return scalar(1);
  else if (exponent == 1)  return base;
  else if (exponent == 2)  return base*base;
  else if (exponent % 2 == 1)  return base * intpow(base,exponent/2) * intpow(base,exponent/2); 
  else return intpow(base,exponent/2) * intpow(base,exponent/2);
}


// A single term of a polynomial. 
// Example: alpha * x1^3 * x2 * x4^2
// coefficient = alpha; exponents = {0, 3, 1, 0, 4}

_tSN
class PolyTerm {

public:
  scalar coefficient;    // Scalar factor
  echar exponents[num_vars];    // Powers of each term -- 0..255 is plenty for all purposes
    
  // Adds up all exponents to get total degree of term. 
  // Example: alpha * x1^3 * x2 * x4^2 => degree = 8
  int degree() const;
    
  // Constructor sets everything to zero
  PolyTerm();

  // Popular constructors
  PolyTerm(scalar coeff);                                       // constant
  PolyTerm(scalar coeff, echar ex0);                              // univariate
  PolyTerm(scalar coeff, echar ex0, echar ex1);                     // bivariate
  PolyTerm(scalar coeff, echar ex0, echar ex1, echar ex2);            // trivariate
  PolyTerm(scalar coeff, echar ex0, echar ex1, echar ex2, echar ex3);   // quadrivariate

    


  /**************     POLYTERM ARITHMETIC     ***************/

  // Assignment operator, also works for different scalar type
  _tO
  PolyTerm _SN & operator= (const PolyTerm _ON &rhs);

  // Assignment operator
  PolyTerm _SN & operator= (const scalar rhs);

  // Multiplication-assignment operator: adds exponents, multiplies coefficients
  PolyTerm  _SN & operator*= (const PolyTerm _SN &rhs);

  // Addition operator: term + term = polynomial
  const TruncPoly _SN operator+ (const PolyTerm _SN &rhs) const;

  // Addition operator: term - term = polynomial
  const TruncPoly _SN operator- (const PolyTerm _SN &rhs) const;

  // Multiplication operator: adds exponents, multiplies coefficients
  const PolyTerm  _SN  operator* (const PolyTerm _SN &rhs) const;

  // Relation operator (needed for sorting). Sorting order is first to last exponent
  // (x0 * x1^4 < x0^2 * x1^3)
  const bool  operator< (const PolyTerm _SN &rhs) const;

  // Check exponents for equality (NOT coefficients)
  const bool  operator== (const PolyTerm _SN &rhs) const;
    
  // Same for inequality (NOT coefficients)
  inline bool  operator!= (const PolyTerm _SN &rhs) const {
    return !(*this == rhs); // just invert equality!
  }

  // Access exponents by [] brackets:
  echar& operator[] (int var);       
  int  operator[] (int var) const; 

};

  
// A polynomial

_tSN
class TruncPoly {

public:
    
  int trunc_degree;
  vector<PolyTerm _SN > terms;
  bool consolidated;
  
  // Return count of terms
  int getNumTerms() {return terms.size();} 

  // Default constructor: empty polynomial
  TruncPoly() {trunc_degree = -1; consolidated = true;} // do not truncate by default

  // Construct from PolyTerm
  inline TruncPoly(PolyTerm _SN term) {
    trunc_degree = -1;
    terms.push_back(term);
    consolidated = true;
  }
  
  // Construct from TruncPoly
  _tO
  inline TruncPoly(TruncPoly _ON rhs) {
    operator=(rhs);
  }

  // Construct from scalar
  inline TruncPoly(scalar constant) {
    trunc_degree = -1;
    terms.push_back(PolyTerm _SN(constant));
    consolidated = true;
  }


  /********     RAW DATA DUMPING (for easy interfacing to C-like languages such as GLSL, CUDA)     *******/

  // Construct by raw data dump
  inline TruncPoly(int num_terms, void* data, int elementSize);

  // Dump in raw data
  inline void assign(int num_terms, void* data, int elementSize) {TruncPoly(num_terms, data, elementSize);}

  // Dump out raw data
  inline void read_data(void* data, int elementSize);


  /**************     MISC ASSIGNMENT FUNCTIONS     ***************/

  // Assign TruncPoly to TruncPoly, also works for other scalar types
  _tO
  TruncPoly _SN &  operator= (const TruncPoly _ON &rhs);

  // Assignment PolyTerm to TruncPoly (a single term is also a polynomial)
  TruncPoly _SN &  operator= (const PolyTerm _SN &rhs);


  /**************     POLYNOMIAL ARITHMETIC     ***************/

  // Add to polynomial
  TruncPoly  _SN & operator+= (const TruncPoly _SN &rhs);

  // Subtract from polynomial
  TruncPoly  _SN & operator-= (const TruncPoly _SN &rhs);

  // Multiply with polynomial
  TruncPoly  _SN & operator*= (const TruncPoly _SN &rhs);

  // Add two polynomials
  const TruncPoly  _SN operator+ (const TruncPoly _SN &rhs) const;

  // Subtract two polynomials
  const TruncPoly  _SN operator- (const TruncPoly _SN &rhs) const;

  // Multiply two polynomials
  const TruncPoly  _SN operator* (const TruncPoly _SN &rhs) const;

  // nth power of a polynomial
  inline const TruncPoly _SN operator^ (const echar power) const {
    //assert(power >= 0); // always true, because unsigned
    TruncPoly _SN result;
    if (power == 0) result=TruncPoly _SN ((scalar)1);
    else if (power == 1)  result=*this;
    else if (power == 2)  result=*this * *this;
    else if (power % 2 == 1)  result= *this * (*this^(power/2)) * (*this^(power/2)); 
    else result = (*this^(power/2)) * (*this^(power/2)); 
    return result;
  }

  // Linearly interpolate with other poly, returning new poly with num_vars+1 input vars
  const TruncPoly <scalar, num_vars+1> lerp_with(  const TruncPoly _SN &poly2, const scalar position1, const scalar position2) const;

  // Quadratically interpolate with two other polys, returning new
  // poly with num_vars+1 input vars
  const TruncPoly <scalar, num_vars+1> querp_with(const TruncPoly _SN &poly2, const TruncPoly _SN &poly3, const scalar pos1, const scalar pos2, const scalar pos3) const;

  // Bake one input variable as constant into coefficients
  const TruncPoly <scalar, num_vars-1> bake_input_variable(int varIdx, scalar value) const;

  // Evaluate the poly: Map input (x0,...) to output scalar
  scalar evaluate(const scalar x[num_vars]);
    
  // Compute derivative with respect to var. wrt at location x.
  scalar get_derivative(int wrt, scalar x[num_vars]);

  // Derive system analytically wrt ("with respect to") variable wrt.
  const TruncPoly _SN  get_derivative(int wrt) const;

  // Swap columns
  void swap_columns(int col1, int col2) {
    assert(col1 >= 0 && col1 < num_vars && col2 >= 0 && col2 < num_vars);
    for (int i = 0; i < terms.size(); ++i) {
      echar tmp = terms[i].exponents[col1];
      terms[i].exponents[col1] = terms[i].exponents[col2];
      terms[i].exponents[col2] = tmp;
    }
  }


  // Truncates poly degree, sorts terms, and consolidates duplicates 
  // Example: 4 x0^3 + 5 x0^3 + 2 x0^2 => 9 x0^3 + 2 x0^2
  void consolidate_terms();//int degree=-1);    


  // Builds the concatenation of two polynomials P1 and P2 (P2 must be in 1 variable):
  // x'' = P12(x) = P2(P1(x));
  const TruncPoly _SN plug_into(const TruncPoly <scalar, 1> &P2) const;

  // Concatenate polynomials in operator syntax
  const TruncPoly _SN operator>>(const TruncPoly <scalar, 1> &P2) const {
    return plug_into(P2);
  }

  // Cut polynomial to desired degree
  void cut(int degree) {
    trunc_degree = degree;
    consolidate_terms();
  }

  TruncPoly  _SN & operator%= (int degree) {
    cut(degree); return *this;
  }

  const TruncPoly  _SN operator% (int degree) const {
    TruncPoly  _SN result = *this;
    result.cut(degree);
    return result;
  }

  // Needs to be friends with PolyTerm+PolyTerm => TruncPoly addition
  template <typename s, int nv> 
  friend  TruncPoly<s, nv> operator+(const PolyTerm<s,nv>& lhs, const PolyTerm<s,nv>& );

  // Search term with desired exponents, and return corresponding coefficient
  scalar get_coeff(const echar *e0) const;
  scalar get_coeff(echar e0) const;
  scalar get_coeff(echar e0, echar e1) const;
  scalar get_coeff(echar e0, echar e1, echar e2) const;
  scalar get_coeff(echar e0, echar e1, echar e2, echar e3) const;
  scalar get_coeff(echar e0, echar e1, echar e2, echar e3, echar e4) const;
  scalar get_coeff(echar e0, echar e1, echar e2, echar e3, echar e4, echar e5) const;

  // Output some trivia about the polynomial: 
  // #in variables, trunc degree, highest degree term, #terms
  void print_stats(std::ostream& os) const;  

  // Access terms by [] brackets:
  PolyTerm _SN & operator[] (int idx);       
  PolyTerm _SN  operator[] (int idx) const; 

};


// Polynomial system: rectangular equation system of the form
//  x0'         = P1        (x0, ... , x{nvIn-1})
//  ...           ...
//  x{nvOut-1}' = P{nvOut-1}(x0, ... , x{nvIn-1})
//
//  For systems P1 and P2 applied in sequence, x'' = P2(x') = P2(P1(x)), 
//  the combined system P12 can be obtained as P1.plug_into(P2);
//  (If P1 and P2 contain only terms of degree 1, this is equivalent to a matrix product P2P1.)

_tSIO
class TruncPolySystem {

public:
  TruncPoly <scalar, num_vars_in> equations[num_vars_out];
  int trunc_degree;

  TruncPolySystem(int degree=-1){trunc_degree=degree;}

  // Construct from other scalar tps
  _tO
  inline TruncPolySystem(TruncPolySystem _OIO rhs) {
    operator=(rhs);
  }


  // Assign TruncPoly to TruncPolySystem: only defined for num_vars_out==1
  _tO
  TruncPolySystem _SIO &  operator= (const TruncPoly <other_scalar, num_vars_in> &rhs) {
    assert(num_vars_out==1);
    equations[0] = rhs;
    trunc_degree = rhs.trunc_degree;

    return *this;
  }

  // Assign TruncPolySystem to TruncPolySystem
  _tO
  TruncPolySystem _SIO &  operator= (const TruncPolySystem _OIO &rhs) {
    if (this == (TruncPolySystem _SIO *) &rhs)      // Self-assignment?
      return *this;        // Yes, so skip assignment, and just return *this.
    for (int i = 0; i < num_vars_out; ++i)
      equations[i] = rhs.equations[i];
    trunc_degree = rhs.trunc_degree;
    return *this;
  }


  // Access equations by [] brackets:
  TruncPoly <scalar, num_vars_in>& operator[] (int i);       
  TruncPoly <scalar, num_vars_in>  operator[] (int i) const; 




  // Builds the concatenation of two polynomial systems P1 and P2:
  // x'' = P12(x) = P2(P1(x));
  template <int num_vars_result>
  const TruncPolySystem<scalar, num_vars_in, num_vars_result> plug_into(const TruncPolySystem<scalar, num_vars_out, num_vars_result> &P2) const;

  // Concatenate systems in operator syntax
  template <int num_vars_result>
  const TruncPolySystem<scalar, num_vars_in, num_vars_result> operator>>(const TruncPolySystem<scalar, num_vars_out, num_vars_result> &P2) const {
    TruncPolySystem<scalar, num_vars_in, num_vars_result> result;
    result = (*this).plug_into(P2);
    return result;
  }

  // Builds the concatenation of a polynomial system P1 and a polynomial P2:
  // x'' = P12(x) = P2(P1(x));
  const TruncPoly <scalar, num_vars_in> plug_into(const TruncPoly <scalar, num_vars_out> &P2) const {
    TruncPolySystem <scalar, num_vars_out, 1> proxy2;
    proxy2.equations[0] = P2;  proxy2.trunc_degree = P2.trunc_degree;
    TruncPolySystem <scalar, num_vars_in, 1> proxyOut = plug_into(proxy2);
    return proxyOut[0];

  }

  // Concatenate system and polynomial in operator syntax
  const TruncPoly <scalar, num_vars_in> operator>>(const TruncPoly <scalar, num_vars_out> &P2) const {
    return plug_into(P2);
  }

  // Linearly interpolate with other system, returning new system with num_vars_in+1 input vars
  const TruncPolySystem <scalar, num_vars_in+1, num_vars_out> lerp_with(  const TruncPolySystem _SIO &system2, const scalar position1, const scalar position2) const;

  // Quadratically interpolate with two other systems, returning new
  // system with num_vars_in+1 input vars
  const TruncPolySystem <scalar, num_vars_in+1, num_vars_out> querp_with(const TruncPolySystem _SIO &system2, const TruncPolySystem _SIO &system3, const scalar pos1, const scalar pos2, const scalar pos3) const;

  // Bake one input variable as constant into coefficients
  const TruncPolySystem <scalar, num_vars_in-1, num_vars_out> bake_input_variable(int varIdx, scalar value) const;

  // Add a new input variable that is passed through to output without doing anything else
  const TruncPolySystem <scalar, num_vars_in+1, num_vars_out+1> add_passthrough_variable() const;

  // Drop a line from the equation system
  const TruncPolySystem <scalar, num_vars_in, num_vars_out-1> drop_equation(int eqIdx) const;

  // Evaluate the system: Map input (x0,...) to output (x1,...)
  void evaluate(const scalar x0[num_vars_in], scalar x1[num_vars_out], bool print = false);

  // Evaluate the system: Map input (x0,...) to output (x1,...)
  void evaluate_double(const scalar x0[num_vars_in], scalar x1[num_vars_out], bool print = false);

  // Evaluate the system and return in same variable
  void evaluate(scalar x0[num_vars_in]);

  // Evaluate the i-th equation of the system
  scalar evaluate_line(const scalar x0[num_vars_in], int i);

  // Derive system analytically wrt ("with respect to") variable wrt.
  const TruncPolySystem _SIO  get_derivative(int wrt) const;

  // Computes the system's Jacobian at location x, and stores it in coeffs (row-major indexing).
  void get_jacobian(scalar coeffs[num_vars_in*num_vars_out], scalar x[num_vars_in]);

  // Swap rows
  void swap_rows(int row1, int row2) {
    assert (row1 >= 0 && row1 < num_vars_out && row2 >= 0 && row2 < num_vars_out);
    TruncPoly <scalar, num_vars_in> temp = equations[row1];
    equations[row1] = equations[row2];
    equations[row2] = temp;
  }

  // Swap columns
  void swap_columns(int col1, int col2) {
    for (int i = 0; i < num_vars_out; ++i)
      equations[i].swap_columns(col1,col2);  
  }

  // Limits poly degree, sorts terms, and consolidates duplicates 
  // Example: 4 x0^3 + 5 x0^3 + 2 x0^2 => 9 x0^3 + 2 x0^2
  inline void consolidate_terms() {
    for (int i = 0; i < num_vars_out; ++i) {
      equations[i].trunc_degree = trunc_degree;
      equations[i].consolidate_terms();
    }
  } 

  // Cut polynomial system to desired degree
  void cut(int degree) {
    trunc_degree = degree;
    consolidate_terms();
  }

  TruncPolySystem  _SIO & operator%= (int degree) {
    cut(degree);
    return *this;
  }

  const TruncPolySystem  _SIO operator% (int degree) const {
    TruncPolySystem _SIO result = *this;
    result.cut(degree);
    return result;
  }

  // Add to this polynomial systems
  TruncPolySystem   _SIO & operator+= (const TruncPolySystem  _SIO &rhs);

  // Add two polynomial systems
  const TruncPolySystem   _SIO operator+ (const TruncPolySystem  _SIO &rhs) const;

  // Output some trivia about the system: 
  // #in variables, #out variables, trunc degree, highest degree term, #terms

  void print_stats(std::ostream& os) const;

};



/**************     FORMATTED OUTPUT     ***************/

// Print PolyTerm as text
_tSN
std::ostream& operator<<(std::ostream& os, const PolyTerm _SN& term);

// Print TruncPoly as text (polynomial)
_tSN
std::ostream& operator<<(std::ostream& os, const TruncPoly _SN& poly);

// Print TruncPolySystem as text (equation system)
_tSIO
std::ostream& operator<<(std::ostream& os, const TruncPolySystem _SIO& tpt);


/**************************************************************************/
/*                         IMPLEMENTATION PART                            */
/**************************************************************************/



/**************     MIXED-CLASS ARITHMETIC     ***************/




// Addition of polynomial and term (in general) gives a polynomial
template <typename s, int nv>
TruncPoly<s, nv> operator+(const TruncPoly<s,nv>& lhs, const PolyTerm<s,nv>& rhs){

  TruncPoly <s, nv> poly = lhs;
  poly.terms.push_back(rhs);
  poly.consolidate_terms();
      
  return poly;
}

// Addition of term and polynomial (in general) gives a polynomial
template <typename s, int nv>
TruncPoly<s, nv> operator+(const PolyTerm<s,nv>& lhs, const TruncPoly<s,nv>& rhs){

  return rhs+lhs;
}

// Left-multiplication of PolyTerm with a scalar
template <typename s, int nv>
PolyTerm <s,nv> operator*(const s lhs, const PolyTerm <s,nv> & rhs) {
  PolyTerm <s,nv> result = rhs;
  result.coefficient *= lhs;
  return result;
}

// Left-multiplication of TruncPoly with a scalar
template <typename s, int nv>
TruncPoly <s,nv> operator*(const s lhs, const TruncPoly <s,nv> & rhs) {

  if (lhs == 0) { TruncPoly <s,nv> zero; return zero;} // multiplication by zero: return empty polynomial
  TruncPoly <s,nv> poly = rhs;
  for (unsigned int i = 0; i < poly.terms.size(); ++i) {
    poly.terms[i].coefficient *= lhs;
  }
  return poly;
}









/**************  PolyTerm CLASS FUNCTIONS  ****************/


_tSN
int PolyTerm _SN::degree() const {
  int d = 0;
  for (int i = 0; i < num_vars; ++i) {
    d += exponents[i];
  }
  return d;
}    

_tSN
inline PolyTerm _SN::PolyTerm() {
  coefficient = 0;
  for (int i = 0; i < num_vars; ++i)
    exponents[i] = 0;
}    

  // Constructor from scalar
  _tSN
inline PolyTerm _SN::PolyTerm(scalar constant) {
    coefficient = constant;
    for (int i = 0; i < num_vars; ++i)
      exponents[i] = 0;
  }

  // Construct a term with 1 independent variable
  _tSN
inline PolyTerm _SN::PolyTerm(scalar coeff, echar ex0) {
    assert(num_vars == 1);
    coefficient = coeff;
    exponents[0] = ex0; 
  }

  // Construct a term with 2 independent variables
  _tSN
inline PolyTerm _SN::PolyTerm(scalar coeff, echar ex0, echar ex1) {
    assert(num_vars == 2);
    coefficient = coeff;
    exponents[0] = ex0; 
    exponents[1] = ex1; 
  }

  // Construct a term with 3 independent variables
  _tSN
inline PolyTerm _SN::PolyTerm(scalar coeff, echar ex0, echar ex1, echar ex2) {
    assert(num_vars == 3);
    coefficient = coeff;
    exponents[0] = ex0; 
    exponents[1] = ex1; 
    exponents[2] = ex2; 
  }

  // Construct a term with 4 independent variables
  _tSN
inline PolyTerm _SN::PolyTerm(scalar coeff, echar ex0, echar ex1, echar ex2, echar ex3) {
    assert(num_vars == 4);
    coefficient = coeff;
    exponents[0] = ex0; 
    exponents[1] = ex1; 
    exponents[2] = ex2; 
    exponents[3] = ex3;
  }


  _tSN
  _tO
  PolyTerm  _SN & PolyTerm _SN::operator= (const PolyTerm _ON &rhs) {
  if (this == (PolyTerm _SN *) &rhs)      // Same object?
    return *this;        // Yes, so skip assignment, and just return *this.
  coefficient = rhs.coefficient;
  for (int i = 0; i < num_vars; ++i)
    exponents[i] = rhs.exponents[i];
  return *this;
}

_tSN
PolyTerm  _SN & PolyTerm _SN::operator= (const scalar rhs) {
  coefficient = rhs;
  for (int i = 0; i < num_vars; ++i)
    exponents[i] = 0;
  return *this;
}

_tSN
const TruncPoly  _SN PolyTerm _SN::operator+ (const PolyTerm _SN &rhs) const {
  TruncPoly _SN result;
  result = result + *this + rhs;
  return result;
} 

_tSN
const TruncPoly  _SN PolyTerm _SN::operator- (const PolyTerm _SN &rhs) const {
  TruncPoly _SN result;
  result = result + *this - rhs;
  return result;
}  

_tSN
PolyTerm  _SN & PolyTerm _SN::operator*= (const PolyTerm _SN &rhs) {
  coefficient *= rhs.coefficient;
  for (int i = 0; i < num_vars; ++i) {
    exponents[i] += rhs.exponents[i];
  }
  return *this;
} 

_tSN
const PolyTerm  _SN PolyTerm _SN::operator* (const PolyTerm _SN &rhs) const {
  PolyTerm _SN result = *this;
  result *= rhs;           
  return result;           
} 

_tSN
const bool PolyTerm _SN::operator< (const PolyTerm _SN &rhs) const {
  for (int i = 0; i < num_vars; ++i) { 
    // The first differing exponent decides
    if (exponents[i] != rhs.exponents[i]) 
      return (exponents[i] < rhs.exponents[i]);
  }
  return false; // else equal, hence "not less than"
}

_tSN
const bool PolyTerm _SN::operator== (const PolyTerm _SN &rhs) const {
  
  for (int i = 0; i < num_vars; ++i) { 
    // If one of the exponents differs, not equal
    if (exponents[i] != rhs.exponents[i])
      return false;
  }
  
  return true; // else equal
}

  // Bracket operator for exponent access:
_tSN
inline echar& PolyTerm _SN::operator[] (int var) {
  assert(var>=0 && var<num_vars);
  return exponents[var];
}
 
_tSN
inline int PolyTerm _SN::operator[] (int var) const {
  assert(var>=0 && var<num_vars);
  int result;
  return exponents[var];
} 




/**************  TruncPoly CLASS FUNCTIONS  ****************/

_tSN
inline TruncPoly _SN::TruncPoly(int num_terms, void* data, int elementSize) {
  assert(elementSize==sizeof(PolyTerm _SN)); // make sure the data is of the right format!
  trunc_degree = -1;
  terms.resize(num_terms);
  for (int i = 0; i < num_terms; ++i) { // avoid memcpy
    terms[i] = ((PolyTerm _SN*) data) [i];
  }
  consolidated = false;
}

_tSN
inline void TruncPoly _SN::read_data(void* data, int elementSize) {
  assert(elementSize==sizeof(PolyTerm _SN)); // make sure the data is of the right format!
  for (int i = 0; i < terms.size(); ++i) { // avoid memcpy
    ((PolyTerm _SN*) data) [i] = terms[i];
  }
}

_tSN
_tO
TruncPoly _SN & TruncPoly _SN::operator= (const TruncPoly _ON &rhs) {
  
if (this == (TruncPoly _SN *)&rhs)      // Same object?
    return *this;        // Yes, so skip assignment, and just return *this.
  trunc_degree = rhs.trunc_degree;
  terms.resize(rhs.terms.size());
  for (int i = 0; i < (int)terms.size(); ++i) terms[i] = rhs.terms[i];
  consolidated = rhs.consolidated;
  return *this;
}

_tSN
TruncPoly _SN & TruncPoly _SN::operator= (const PolyTerm _SN &rhs) {
  trunc_degree = -1;
  terms.clear();
  terms.push_back(rhs);
  consolidated = true;
  return *this;
}

_tSN
TruncPoly  _SN & TruncPoly _SN::operator+= (const TruncPoly _SN &rhs) {
  // Be conservative about truncation. User can always re-truncate later
  trunc_degree = max(trunc_degree, rhs.trunc_degree); 
  int lhsSize = terms.size();
  terms.resize(lhsSize + rhs.terms.size());

  // Simply add rhs's terms:
  int i = lhsSize;

  for (int j = 0; j < (int)rhs.terms.size(); ++j, ++i) 
    terms[i] = rhs.terms[j];

  consolidate_terms();  
  return *this;
}

_tSN
const TruncPoly  _SN TruncPoly _SN::operator+ (const TruncPoly _SN &rhs) const {
  TruncPoly _SN result = *this;
  result += rhs;           
  return result;           
}

_tSN
TruncPoly  _SN & TruncPoly _SN::operator-= (const TruncPoly _SN &rhs) {
  *this += ((scalar)(-1) * rhs);
  return *this;
}

_tSN
const TruncPoly  _SN TruncPoly _SN::operator- (const TruncPoly _SN &rhs) const {
  TruncPoly _SN result = *this;
  result -= rhs;           
  return result;           
}

_tSN
TruncPoly  _SN & TruncPoly _SN::operator*= (const TruncPoly _SN &rhs) {
  TruncPoly _SN result;
  result.terms.resize(terms.size() * rhs.terms.size());
  // Be conservative about truncation. User can always re-truncate later
  result.trunc_degree = max(trunc_degree, rhs.trunc_degree); 
  int t = 0;
  for (int i = 0; i < (int)terms.size(); ++i)
    for (int j = 0; j < (int)rhs.terms.size(); ++j) {
      if (trunc_degree < 0 || (terms[i].degree() + rhs.terms[j].degree()) <= trunc_degree) {
	PolyTerm _SN product = terms[i] * rhs.terms[j];
	result.terms[t++] = product;
      }
    }

  result.consolidate_terms();
  terms = result.terms;
  return *this;
  
}

_tSN
const TruncPoly  _SN TruncPoly _SN::operator* (const TruncPoly _SN &rhs) const {
  TruncPoly _SN result = *this;
  result *= rhs;           
  return result;           
}

_tSN
void TruncPoly _SN::consolidate_terms() {

  int i; // index input side
  int j; // index output side

  // Pass 0: Sort terms according to their exponents
  sort(terms.begin(), terms.end());

  // Pass 1: Only keep terms with tolerable degree and nonzero coefficient
  vector<PolyTerm _SN> newTerms;
  j = 0;
  newTerms.resize(terms.size());
  for (i = 0; i < (int)terms.size(); ++i) {
    if ((trunc_degree < 0 || terms[i].degree() <= trunc_degree) && terms[i].coefficient!=(scalar)0) { 
      newTerms[j++] = terms[i];
    }
  }
  newTerms.resize(j);

  // Pass 2: Add up terms with identical exponents
  if (newTerms.size() > 0) { 
    terms[0] = newTerms[0];
    j = 1;
    for (i = 1; i < (int)newTerms.size(); ++i) {
      if (newTerms[i-1] < newTerms[i]) {
	terms[j++] = newTerms[i];
      } else {
	terms[j-1].coefficient += newTerms[i].coefficient;
      } 
    }
    terms.resize(j);
  } else terms.resize(0);

  consolidated = true;

}


  // Linearly interpolate with other poly, returning new poly with num_vars+1 input vars
_tSN
const TruncPoly <scalar, num_vars+1> TruncPoly _SN::lerp_with(const TruncPoly _SN &poly2, const scalar position1, const scalar position2) const {

  TruncPoly _SN poly1 = *this;

  TruncPoly<scalar, num_vars+1> result;

  // Reserve space for terms from input systems
  int num_terms = poly1.terms.size() + poly2.terms.size();
  result.terms.resize(2 * num_terms);
  int j = 0;
  // Copy terms from poly 1, and add linear dependency
  for (int k = 0; k < (int)poly1.terms.size(); ++k) {
    scalar c1 = poly1.terms[k].coefficient;
    result.terms[j].coefficient = c1 * (1-(position1/(position2-position1)));
    for (int ie = 0; ie < num_vars; ++ie)
      result.terms[j].exponents[ie] = poly1.terms[k].exponents[ie];
    result.terms[j].exponents[num_vars] = 0;
    
    // Copy term
    result.terms[num_terms + j] = result.terms[j];
    result.terms[num_terms + j].coefficient = - c1 / (position2-position1);
    result.terms[num_terms + j].exponents[num_vars] = 1;
    ++j;
  }
  
  
  
  // Copy terms from poly 2
  for (int k = 0; k < (int)poly2.terms.size(); ++k) {
    scalar c2 = poly2.terms[k].coefficient;
    result.terms[j].coefficient = c2 * position1 / (position2-position1);
    for (int ie = 0; ie < num_vars; ++ie)
      result.terms[j].exponents[ie] = poly2.terms[k].exponents[ie];
    result.terms[j].exponents[num_vars] = 0;
    
    // Copy term
    result.terms[num_terms + j] = result.terms[j];
    result.terms[num_terms + j].coefficient = c2 / (position2-position1);
    result.terms[num_terms + j].exponents[num_vars] = 1;
    ++j;
  }

  // Truncate resulting degree to max_degree+1 to preserve additional linear dependency
  if (poly1.trunc_degree >= 0 && poly2.trunc_degree >= 0)
    result %= max(poly1.trunc_degree,poly2.trunc_degree) + 1;
  else result %= -1;
  return result;

};

  // Quadratically interpolate with two other polys, returning new
  // poly with num_vars+1 input vars
_tSN
const TruncPoly <scalar, num_vars+1> TruncPoly _SN::querp_with(const TruncPoly _SN &poly2, const TruncPoly _SN &poly3, const scalar pos1, const scalar pos2, const scalar pos3) const {

};

  // Bake one input variable as constant into coefficients
_tSN
const TruncPoly <scalar, num_vars-1> TruncPoly _SN::bake_input_variable(int varIdx, scalar value) const {
};











// Evaluate one line (equation) of the system: Map input (x0,...) to output x_i
_tSN
scalar TruncPoly _SN::evaluate(const scalar x[num_vars]){
  scalar result = 0; 

  
  for (int j = terms.size() - 1; j >= 0; --j) {
    PolyTerm<scalar,num_vars> term = terms[j];
    scalar term_value = term.coefficient;
    for (int k = 0; k < num_vars; ++k)
      if (term.exponents[k]) 
	term_value *= intpow(x[k],term.exponents[k]);
    result += term_value;
  }
  return result;
}


_tSN
scalar TruncPoly _SN::get_derivative(int wrt, scalar x[num_vars]) {
  scalar der = 0;
  for (int i = 0; i < (int)terms.size(); ++i) {
    
    scalar term = terms[i].coefficient;

    for (int v = 0; v < num_vars; ++v) { // go through variables
      
      if (v == wrt) {
	echar e = terms[i].exponents[v];
	if (e) {
	  term *= e * intpow(x[v], e-1);
	} else term = 0;
      } else {
	echar e = terms[i].exponents[v];
	term *= intpow(x[v], e);
      }
    }
    
    der += term;
  }
  return der;
}

_tSN
const TruncPoly _SN  TruncPoly _SN::get_derivative(int wrt) const {
  assert(wrt >= 0 && wrt < num_vars);
  
  TruncPoly _SN result;
  
  for (int i = 0; i < (int)terms.size(); ++i) {
    
    PolyTerm _SN term;
    echar e = terms[i].exponents[wrt];

    if (e) {
      term.coefficient = terms[i].coefficient * (scalar)e;
      for (int j = 0; j < num_vars; ++j)
	term.exponents[j] = terms[i].exponents[j];
      --term.exponents[wrt];
      if (term.coefficient!=(scalar)0) 
	result.terms.push_back(term);
    }
  }
  
  //result.consolidate_terms();
  return result;

}

  // Search term with desired exponents, and return corresponding coefficient
_tSN
scalar TruncPoly _SN::get_coeff(const echar *e) const {
  bool found = false;
  int i = 0;
  for (i = 0; i < (int)terms.size() && !found; ++i) {
    found = true;
    for (int j = 0; j < num_vars; ++j)
      if (terms[i].exponents[j] != e[j]) found = false;
  }

  if (!found) return 0;
  return terms[i-1].coefficient;
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0) const {
  assert(num_vars == 1);
  const echar exponents[1] = {e0};
  return get_coeff(exponents);
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0, echar e1) const {
  assert(num_vars == 2);
  const echar exponents[2] = {e0,e1};
  return get_coeff(exponents);
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0, echar e1, echar e2) const {
  assert(num_vars == 3);
  const echar exponents[3] = {e0,e1,e2};
  return get_coeff(exponents);
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0, echar e1, echar e2, echar e3) const {
  assert(num_vars == 4);
  const echar exponents[4] = {e0,e1,e2,e3};
  return get_coeff(exponents);
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0, echar e1, echar e2, echar e3, echar e4) const {
  assert(num_vars == 5);
  const echar exponents[5] = {e0,e1,e2,e3,e4};
  return get_coeff(exponents);
}

_tSN
scalar TruncPoly _SN::get_coeff(echar e0, echar e1, echar e2, echar e3, echar e4, echar e5) const {
  assert(num_vars == 6);
  const echar exponents[6] = {e0,e1,e2,e3,e4,e5};
  return get_coeff(exponents);
}



_tSN
void TruncPoly _SN::print_stats(std::ostream& os) const {
  int num_terms = terms.size();
  int highest_degree = 0;
  for (int j = 0; j < num_terms; ++j) {
    int term_degree = terms[j].degree();
    if (term_degree > highest_degree) 
      highest_degree = term_degree;
    }

  os << "Polynomial (" 
     << num_vars << "in, num_terms=" 
     << num_terms << ", highest_degree=" 
     << highest_degree <<", trunc=";
  if (trunc_degree >= 0) os << trunc_degree;
  else os << "NOTRUNC";
  os << ")" << endl;
}



  // Builds the concatenation of two polynomials P1 and P2 (P2 must be in 1 variable)
  // x'' = P12(x) = P2(P1(x));
_tSN
const TruncPoly _SN
TruncPoly _SN::plug_into(const TruncPoly<scalar,1> &P2) const {

  TruncPoly _SN P12; // start with empty polynomial 
  
  // Be conservative about truncation. User can always re-truncate later
  P12.trunc_degree = max(trunc_degree, P12.trunc_degree); 

  for (int j = 0; j < (int)P2.terms.size(); ++j) { // then, term by term
    TruncPoly _SN P2Term(P2.terms[j].coefficient); // keep the coefficient
    P2Term.trunc_degree = P12.trunc_degree;
     
    P2Term *= (*this)^(P2.terms[j].exponents[0]);
    P12 += P2Term;
    P12.consolidate_terms();
  }
  return P12;
}

_tSN
inline
PolyTerm _SN & TruncPoly _SN::operator[] (int idx) {
  assert(idx>=0 && idx<terms.size());
  return terms[idx];
}
 
_tSN
inline
PolyTerm _SN TruncPoly _SN::operator[] (int idx) const {
  assert(idx>=0 && idx<terms.size());
  return terms[idx];
} 


/**************  TruncPolySystem CLASS FUNCTIONS  ****************/

// Constructs system by linear blending of two systems with one less input variable
_tSIO
const TruncPolySystem<scalar, num_vars_in+1, num_vars_out>
TruncPolySystem _SIO::lerp_with( const TruncPolySystem _SIO &system2, 
				 const scalar position1,
				 const scalar position2) const {
  TruncPolySystem _SIO system1 = *this;

  TruncPolySystem<scalar, num_vars_in+1, num_vars_out> result;
  
  int trunc_degree = -1;

  // Lerping a system means to lerp each of its equations:
  for (int i = 0; i < num_vars_out; ++i) {
    result[i] = system1[i].lerp_with(system2[i], position1, position2);
    trunc_degree = max(trunc_degree, result[i].trunc_degree);
  }

  result %= trunc_degree;
  return result;
}


// Constructs system by quadratic blending ("querping") of three systems
_tSIO
const TruncPolySystem<scalar, num_vars_in+1, num_vars_out>
TruncPolySystem _SIO::querp_with( const TruncPolySystem _SIO &system2, 
				  const TruncPolySystem _SIO &system3, 
				  const scalar pos1,
				  const scalar pos2, 
				  const scalar pos3) const {

  TruncPolySystem<scalar, num_vars_in+1, num_vars_out> result;

  bool TermsSame = true;
  for (int i = 0; i < num_vars_out; ++i) {
	// Reserve space for terms from input systems
    int num_terms = system2[i].terms.size();
    result[i].terms.resize(3 * num_terms);
    for (int k = 0; k < (int)system2[i].terms.size(); ++k) {
      for (int e = 0; e < num_vars_in; ++e) 
	if (equations[i].terms[k].exponents[e] != system2[i].terms[k].exponents[e]
	    || equations[i].terms[k].exponents[e] != system3[i].terms[k].exponents[e])
	  TermsSame = false;

      double x1 = pos2;
      double y1 = system2[i].terms[k].coefficient;
      double x2 = pos1;
      double y2 = equations[i].get_coeff(equations[i].terms[k].exponents);
      double x3 = pos3;
      double y3 = system3[i].get_coeff(equations[i].terms[k].exponents);


      double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
      // Determine quadratic (A), linear (B) and constant (C) coefficients
      double A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
      double B = (x3 * x3 * (y1 - y2) 
		  + x2 * x2 * (y3 - y1) 
		  + x1 * x1 * (y2 - y3)
		  ) / denom;
      double C = (
		  x2 * x3 * (x2 - x3) * y1 
		  + x3 * x1 * (x3 - x1) * y2 
		  + x1 * x2 * (x1 - x2) * y3
		  ) / denom;

      // Make three new terms (parabola in the new variable)
      result[i].terms[k].coefficient = A;
      for (int e = 0; e < num_vars_in; ++e) {
	result[i].terms[k].exponents[e] = system2[i].terms[k].exponents[e];
      }
      result[i].terms[k].exponents[num_vars_in] = 2;
      
      result[i].terms[k+num_terms].coefficient = B;
      for (int e = 0; e < num_vars_in; ++e) {
	result[i].terms[k+num_terms].exponents[e] = system2[i].terms[k].exponents[e];
      }
      result[i].terms[k+num_terms].exponents[num_vars_in] = 1;

      result[i].terms[k+2*num_terms].coefficient = C;
      for (int e = 0; e < num_vars_in; ++e) {
	result[i].terms[k+2*num_terms].exponents[e] = system2[i].terms[k].exponents[e];
      }
      result[i].terms[k+2*num_terms].exponents[num_vars_in] = 0;
    }
  }

  if (trunc_degree >= 0 && system3.trunc_degree >= 0 && system2.trunc_degree >= 0)
    result %= max(max(trunc_degree,system3.trunc_degree),system2.trunc_degree) + 2;
  else result %= -1;

  if (!TermsSame) {
    printf("[TruncPolySystem::querp_with] Warning: Querping is experimental and only works reliably when the terms of all three systems are in the same order (which not the case here!)\n");
  }

  return result;

}



_tSIO
const 
TruncPolySystem <scalar, num_vars_in-1, num_vars_out> TruncPolySystem _SIO::bake_input_variable(int varIdx, scalar value) const {
  TruncPolySystem <scalar, num_vars_in-1, num_vars_out> result;
  assert (varIdx >= 0 && varIdx < num_vars_in);

  // Copy terms over, removing "column" varIdx
  for (int i = 0; i < num_vars_out; ++i) {
    result[i].terms.resize(equations[i].terms.size());
    for (int k = 0; k < (int)equations[i].terms.size(); ++k) {
      result[i].terms[k].coefficient = equations[i].terms[k].coefficient 
	* intpow(value, equations[i].terms[k].exponents[varIdx]);
      int j = 0;
      for (int ie = 0; ie < num_vars_in; ++ie) {
	if (ie != varIdx) {
	  result[i].terms[k].exponents[j++] = equations[i].terms[k].exponents[ie];
	}
      }
    }
    result[i].consolidate_terms();
  }
  result %= trunc_degree;
  return result;
}

_tSIO
const 
TruncPolySystem <scalar, num_vars_in+1, num_vars_out+1> TruncPolySystem _SIO::add_passthrough_variable() const {
  TruncPolySystem <scalar, num_vars_in+1, num_vars_out+1> result;
  // Copy old system over
  for (int i = 0; i < num_vars_out; ++i) {
    result[i].terms.resize(equations[i].terms.size());
    for (int k = 0; k < (int)equations[i].terms.size(); ++k) {
      result[i].terms[k].coefficient = equations[i].terms[k].coefficient;
      for (int ie = 0; ie < num_vars_in; ++ie) {
	result[i].terms[k].exponents[ie] = equations[i].terms[k].exponents[ie];
      }
      result[i].terms[k].exponents[num_vars_in] = 0;
    }
  }
  // Add diagonal entry for new variable:
  result[num_vars_out] = PolyTerm<scalar,num_vars_in+1>(1);
  result[num_vars_out].terms[0].exponents[num_vars_in] = 1;
  result %= trunc_degree;
  return result;
}

_tSIO
const 
TruncPolySystem <scalar, num_vars_in, num_vars_out-1> TruncPolySystem _SIO::drop_equation(int eqIdx) const {
  TruncPolySystem <scalar, num_vars_in, num_vars_out-1> result;
  assert (eqIdx > 0 && eqIdx < num_vars_out);
  // Copy terms over, removing "column" varIdx
  int j = 0;

  for (int i = 0; i < num_vars_out; ++i) {
	  if (i != eqIdx) {
		result[j] = equations[i]; 
		++j;
	  }
  }
  result %= trunc_degree;
  return result;
}


_tSIO
inline
TruncPoly<scalar,num_vars_in>& TruncPolySystem _SIO::operator[] (int row) {
  assert(row>=0 && row<num_vars_out);
  return equations[row];
}
 
_tSIO
inline
TruncPoly<scalar,num_vars_in> TruncPolySystem _SIO::operator[] (int row) const {
  assert(row>=0 && row<num_vars_out);
  return equations[row];
} 


  // Builds the concatenation of two polynomial systems P1
  // x'' = P12(x) = P2(P1(x));
_tSIO
template <int num_vars_result>
const TruncPolySystem<scalar,num_vars_in,num_vars_result> 
TruncPolySystem _SIO::plug_into(const TruncPolySystem<scalar,num_vars_out,num_vars_result> &P2) const {

  TruncPolySystem<scalar,num_vars_in,num_vars_result> P12;

  P12.trunc_degree = max(trunc_degree, P2.trunc_degree); // Be conservative about truncation. User can always re-truncate later

  for (int i = 0; i < num_vars_result; ++i) { // Row by row through the equation system
    
    // Now, unroll the polynomial by plugging in the corresponding rows of P1:

    TruncPoly <scalar, num_vars_in> P12i; // start with empty polynomial 
    P12i.trunc_degree = trunc_degree;

    for (int j = 0; j < (int)P2[i].terms.size(); ++j) { // then, term by term

      TruncPoly <scalar, num_vars_in> P2Term(P2[i].terms[j].coefficient); // keep the coefficient
      
      for (int k = 0; k < num_vars_out; ++k) { // but replace all variables with their polynomials as found in P1
	TruncPoly <scalar, num_vars_in> temp = equations[k] ^ P2[i].terms[j].exponents[k];
	P2Term *= temp;
      }
      
      P12i += P2Term;
    } 

    P12.equations[i] = P12i;
    P12.equations[i].trunc_degree = trunc_degree;
    P12.equations[i].consolidate_terms();
 
  }
  return P12;
}




// Evaluate the system: Map input (x0,...) to output (x1,...)
_tSIO
void TruncPolySystem _SIO::evaluate(const scalar x0[num_vars_in], scalar x1[num_vars_out], bool print){
  for (int i = 0; i < num_vars_out; ++i) {
    scalar result = 0;
    for (int j = (int)equations[i].terms.size() - 1; j >= 0; --j) {
      PolyTerm<scalar,num_vars_in> *term = &equations[i].terms[j];
      scalar term_value = term->coefficient;
      for (int k = 0; k < num_vars_in; ++k) 
	term_value *= intpow(x0[k],term->exponents[k]);
      
      result += term_value;
    }
    x1[i] = result;
  }
}

// Evaluate the system with intermediate steps taken in double precision: Map input (x0,...) to output (x1,...)
_tSIO
void TruncPolySystem _SIO::evaluate_double(const scalar x0[num_vars_in], scalar x1[num_vars_out], bool print){
  for (int i = 0; i < num_vars_out; ++i) {
    // Accumulate in double temp to avoid precision issues
    double result = 0;
    for (int j = (int)equations[i].terms.size() - 1; j >= 0; --j) {
      PolyTerm<scalar,num_vars_in> *term = &equations[i].terms[j];
      double term_value = term->coefficient;
      for (int k = 0; k < num_vars_in; ++k) {
	double base = x0[k];
	term_value *= intpow(base,term->exponents[k]);
      }
      
      result += (double)term_value;
    }
    x1[i] = result;
  }
}

// Evaluate the system "in place" (return in same variable)
_tSIO
void TruncPolySystem _SIO::evaluate(scalar x0[num_vars_in]) {
  assert(num_vars_in == num_vars_out);
  scalar x1[num_vars_in];
  evaluate(x0,x1);
  for (int i = 0; i < num_vars_out; ++i)
    x0[i] = x1[i];
}

// Evaluate one line (equation) of the system: Map input (x0,...) to output x_i
_tSIO
scalar TruncPolySystem _SIO::evaluate_line(const scalar x0[num_vars_in], int i){
  assert(i >= 0 && i < num_vars_out);
  return evaluate(equations[i], x0);
}

_tSIO
const TruncPolySystem _SIO TruncPolySystem _SIO::get_derivative(int wrt) const {
  assert(wrt >= 0 && wrt < num_vars_in);
  TruncPolySystem _SIO result;
  for (int i = 0; i < num_vars_out; ++i) {
    result.equations[i] = equations[i].get_derivative(wrt);
  }
  return result;
}

_tSIO
void TruncPolySystem _SIO::get_jacobian(scalar coeffs[num_vars_in*num_vars_out], scalar x[num_vars_in]) {

  for (int i = 0; i < num_vars_out; ++i) {
    for (int j = 0; j < num_vars_in; ++j) {

      coeffs[j+i*num_vars_in] = equations[i].get_derivative(j, x);

    }
  }

}

_tSIO
TruncPolySystem   _SIO & TruncPolySystem _SIO::operator+= 
(const TruncPolySystem  _SIO &rhs) {
  trunc_degree = max(trunc_degree, rhs.trunc_degree); // Be conservative about truncation. User can always re-truncate later  
  for (int i = 0; i < num_vars_out; ++i)
    equations[i] += rhs.equations[i];

  consolidate_terms();
  return *this;
}

_tSIO
const TruncPolySystem   _SIO TruncPolySystem _SIO::operator+ 
(const TruncPolySystem  _SIO &rhs) const {
      
  TruncPolySystem  _SIO result = *this;
  result += rhs;
  return result;
}


_tSIO
void TruncPolySystem _SIO::print_stats(std::ostream& os) const {
  int num_terms[num_vars_out];
  int highest_degree = 0;
  for (int i = 0; i < num_vars_out; ++i) {
    num_terms[i] = equations[i].terms.size();
    for (int j = 0; j < num_terms[i]; ++j) {
      int term_degree = equations[i].terms[j].degree();
      if (term_degree > highest_degree) 
	highest_degree = term_degree;
    }
  }
  os << "Polynomial system (" 
     << num_vars_in << "in, " 
     << num_vars_out << "out, num_terms=";

  for (int i = 0; i < num_vars_out; ++i)
    os << num_terms[i] << ((i==num_vars_out-1)?"":";");
  os << ", highest_degree=" 
     << highest_degree << ", trunc=";
  if (trunc_degree >= 0) os << trunc_degree;
  else os << "NOTRUNC";
  os << ")" << endl;
      
}















/* Non-member functions for I/O and cross-class arithmetic*/


_tSN
void print(std::ostream& os, const PolyTerm _SN term, string *var_names = NULL) {

  os << setprecision(8) << term.coefficient;

  if (var_names) {
    for(int i=0; i<num_vars; i++) {
      echar e = term.exponents[i];
      string vn = var_names[i];
      if (e)
	os << " " <<vn;
      if (e>1)
	os << "^" << int(e);
    }
    
  } else {
    for(int i=0; i<num_vars; i++) {
      echar e = term.exponents[i];
      if (e)
	os << " * x" << i;
      if (e>1)
	os << "^" << int(e);
    }

  }
}

_tSN
std::ostream& operator<<(std::ostream& os, const PolyTerm _SN& term) {

  print(os, term);
  return os;
}

_tSN
void print(std::ostream& os, const TruncPoly _SN poly, string *var_names = NULL, bool linebreak = false) {

  int N = poly.terms.size();

  if (N) {
    for(int i = 0; i < N - 1; i++) {
      print(os, poly.terms[i], var_names);
      if (linebreak) os << endl;
      os << " + ";
    }
    
    print(os, poly.terms[N-1], var_names);
   
  } else os << "[Null]";

}

_tSN
std::ostream& operator<<(std::ostream& os, const TruncPoly _SN& poly) {

  print(os, poly, NULL, false);      
  return os;
}

_tSIO
void print(std::ostream& os, const TruncPolySystem _SIO& tpt, string *var_names = NULL, bool linebreak = false)
{
  for(int i = 0; i < num_vars_out; i++) {
    if (var_names) { 
      string vn = var_names[i];
      os << vn << "' = ";
    } else {
      os << "x" << i << "' = ";
    }

    print(os,tpt.equations[i],var_names,linebreak);
    os << endl;
  }
}


_tSIO
std::ostream& operator<<(std::ostream& os, const TruncPolySystem _SIO& tpt) {

  print(os, tpt);
  return os;
}


// Abbreviated class names

typedef PolyTerm<float,1> Term1f;
typedef PolyTerm<float,2> Term2f;
typedef PolyTerm<float,3> Term3f;
typedef PolyTerm<float,4> Term4f;
typedef PolyTerm<float,5> Term5f;
typedef PolyTerm<float,6> Term6f;
typedef PolyTerm<double,1> Term1d;
typedef PolyTerm<double,2> Term2d;
typedef PolyTerm<double,3> Term3d;
typedef PolyTerm<double,4> Term4d;
typedef PolyTerm<double,5> Term5d;
typedef PolyTerm<double,6> Term6d;
typedef PolyTerm<long double,1> Term1l;
typedef PolyTerm<long double,2> Term2l;
typedef PolyTerm<long double,3> Term3l;
typedef PolyTerm<long double,4> Term4l;
typedef PolyTerm<long double,5> Term5l;
typedef PolyTerm<long double,6> Term6l;


typedef TruncPoly<float,1> Poly1f;
typedef TruncPoly<float,2> Poly2f;
typedef TruncPoly<float,3> Poly3f;
typedef TruncPoly<float,4> Poly4f;
typedef TruncPoly<float,5> Poly5f;
typedef TruncPoly<float,6> Poly6f;
typedef TruncPoly<double,1> Poly1d;
typedef TruncPoly<double,2> Poly2d;
typedef TruncPoly<double,3> Poly3d;
typedef TruncPoly<double,4> Poly4d;
typedef TruncPoly<double,5> Poly5d;
typedef TruncPoly<double,6> Poly6d;
typedef TruncPoly<long double,1> Poly1l;
typedef TruncPoly<long double,2> Poly2l;
typedef TruncPoly<long double,3> Poly3l;
typedef TruncPoly<long double,4> Poly4l;
typedef TruncPoly<long double,5> Poly5l;
typedef TruncPoly<long double,6> Poly6l;

// Systems:
typedef TruncPolySystem<float,1,1> System11f;
typedef TruncPolySystem<float,1,2> System12f;
typedef TruncPolySystem<float,1,3> System13f;
typedef TruncPolySystem<float,1,4> System14f;
typedef TruncPolySystem<float,1,5> System15f;
typedef TruncPolySystem<float,1,6> System16f;
typedef TruncPolySystem<float,2,1> System21f;
typedef TruncPolySystem<float,2,2> System22f;
typedef TruncPolySystem<float,2,3> System23f;
typedef TruncPolySystem<float,2,4> System24f;
typedef TruncPolySystem<float,2,5> System25f;
typedef TruncPolySystem<float,2,6> System26f;
typedef TruncPolySystem<float,3,1> System31f;
typedef TruncPolySystem<float,3,2> System32f;
typedef TruncPolySystem<float,3,3> System33f;
typedef TruncPolySystem<float,3,4> System34f;
typedef TruncPolySystem<float,3,5> System35f;
typedef TruncPolySystem<float,3,6> System36f;
typedef TruncPolySystem<float,4,1> System41f;
typedef TruncPolySystem<float,4,2> System42f;
typedef TruncPolySystem<float,4,3> System43f;
typedef TruncPolySystem<float,4,4> System44f;
typedef TruncPolySystem<float,4,5> System45f;
typedef TruncPolySystem<float,4,6> System46f;
typedef TruncPolySystem<float,5,1> System51f;
typedef TruncPolySystem<float,5,2> System52f;
typedef TruncPolySystem<float,5,3> System53f;
typedef TruncPolySystem<float,5,4> System54f;
typedef TruncPolySystem<float,5,5> System55f;
typedef TruncPolySystem<float,5,6> System56f;
typedef TruncPolySystem<float,6,1> System61f;
typedef TruncPolySystem<float,6,2> System62f;
typedef TruncPolySystem<float,6,3> System63f;
typedef TruncPolySystem<float,6,4> System64f;
typedef TruncPolySystem<float,6,5> System65f;
typedef TruncPolySystem<float,6,6> System66f;
typedef TruncPolySystem<double,1,1> System11d;
typedef TruncPolySystem<double,1,2> System12d;
typedef TruncPolySystem<double,1,3> System13d;
typedef TruncPolySystem<double,1,4> System14d;
typedef TruncPolySystem<double,1,5> System15d;
typedef TruncPolySystem<double,1,6> System16d;
typedef TruncPolySystem<double,2,1> System21d;
typedef TruncPolySystem<double,2,2> System22d;
typedef TruncPolySystem<double,2,3> System23d;
typedef TruncPolySystem<double,2,4> System24d;
typedef TruncPolySystem<double,2,5> System25d;
typedef TruncPolySystem<double,2,6> System26d;
typedef TruncPolySystem<double,3,1> System31d;
typedef TruncPolySystem<double,3,2> System32d;
typedef TruncPolySystem<double,3,3> System33d;
typedef TruncPolySystem<double,3,4> System34d;
typedef TruncPolySystem<double,3,5> System35d;
typedef TruncPolySystem<double,3,6> System36d;
typedef TruncPolySystem<double,4,1> System41d;
typedef TruncPolySystem<double,4,2> System42d;
typedef TruncPolySystem<double,4,3> System43d;
typedef TruncPolySystem<double,4,4> System44d;
typedef TruncPolySystem<double,4,5> System45d;
typedef TruncPolySystem<double,4,6> System46d;
typedef TruncPolySystem<double,5,1> System51d;
typedef TruncPolySystem<double,5,2> System52d;
typedef TruncPolySystem<double,5,3> System53d;
typedef TruncPolySystem<double,5,4> System54d;
typedef TruncPolySystem<double,5,5> System55d;
typedef TruncPolySystem<double,5,6> System56d;
typedef TruncPolySystem<double,6,1> System61d;
typedef TruncPolySystem<double,6,2> System62d;
typedef TruncPolySystem<double,6,3> System63d;
typedef TruncPolySystem<double,6,4> System64d;
typedef TruncPolySystem<double,6,5> System65d;
typedef TruncPolySystem<double,6,6> System66d;
typedef TruncPolySystem<long double,1,1> System11l;
typedef TruncPolySystem<long double,1,2> System12l;
typedef TruncPolySystem<long double,1,3> System13l;
typedef TruncPolySystem<long double,1,4> System14l;
typedef TruncPolySystem<long double,1,5> System15l;
typedef TruncPolySystem<long double,1,6> System16l;
typedef TruncPolySystem<long double,2,1> System21l;
typedef TruncPolySystem<long double,2,2> System22l;
typedef TruncPolySystem<long double,2,3> System23l;
typedef TruncPolySystem<long double,2,4> System24l;
typedef TruncPolySystem<long double,2,5> System25l;
typedef TruncPolySystem<long double,2,6> System26l;
typedef TruncPolySystem<long double,3,1> System31l;
typedef TruncPolySystem<long double,3,2> System32l;
typedef TruncPolySystem<long double,3,3> System33l;
typedef TruncPolySystem<long double,3,4> System34l;
typedef TruncPolySystem<long double,3,5> System35l;
typedef TruncPolySystem<long double,3,6> System36l;
typedef TruncPolySystem<long double,4,1> System41l;
typedef TruncPolySystem<long double,4,2> System42l;
typedef TruncPolySystem<long double,4,3> System43l;
typedef TruncPolySystem<long double,4,4> System44l;
typedef TruncPolySystem<long double,4,5> System45l;
typedef TruncPolySystem<long double,4,6> System46l;
typedef TruncPolySystem<long double,5,1> System51l;
typedef TruncPolySystem<long double,5,2> System52l;
typedef TruncPolySystem<long double,5,3> System53l;
typedef TruncPolySystem<long double,5,4> System54l;
typedef TruncPolySystem<long double,5,5> System55l;
typedef TruncPolySystem<long double,5,6> System56l;
typedef TruncPolySystem<long double,6,1> System61l;
typedef TruncPolySystem<long double,6,2> System62l;
typedef TruncPolySystem<long double,6,3> System63l;
typedef TruncPolySystem<long double,6,4> System64l;
typedef TruncPolySystem<long double,6,5> System65l;
typedef TruncPolySystem<long double,6,6> System66l;

// Transforms: square version of system

typedef TruncPolySystem<float,1,1> Transform1f;
typedef TruncPolySystem<float,2,2> Transform2f;
typedef TruncPolySystem<float,3,3> Transform3f;
typedef TruncPolySystem<float,4,4> Transform4f;
typedef TruncPolySystem<float,5,5> Transform5f;
typedef TruncPolySystem<float,6,6> Transform6f;
typedef TruncPolySystem<double,1,1> Transform1d;
typedef TruncPolySystem<double,2,2> Transform2d;
typedef TruncPolySystem<double,3,3> Transform3d;
typedef TruncPolySystem<double,4,4> Transform4d;
typedef TruncPolySystem<double,5,5> Transform5d;
typedef TruncPolySystem<double,6,6> Transform6d;
typedef TruncPolySystem<long double,1,1> Transform1l;
typedef TruncPolySystem<long double,2,2> Transform2l;
typedef TruncPolySystem<long double,3,3> Transform3l;
typedef TruncPolySystem<long double,4,4> Transform4l;
typedef TruncPolySystem<long double,5,5> Transform5l;
typedef TruncPolySystem<long double,6,6> Transform6l;



#endif // #ifdef cplusplus


// Additional data structure definitions, might come in handy at some point

typedef struct pterm1f {
  float coeff;
  echar e0;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm1f(): coeff(0),e0(0){};
  inline pterm1f(float c, echar e0): coeff(c),e0(e0){};
#endif
} PT1fData;

typedef struct pterm2f {
  float coeff;
  echar e0, e1;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm2f(): coeff(0),e0(0),e1(0){};
  inline pterm2f(float c, echar e0, echar e1): coeff(c),e0(e0),e1(e1){};
#endif
} PT2fData;

typedef struct pterm3f {
  float coeff;
  echar e0, e1, e2;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm3f(): coeff(0),e0(0),e1(0),e2(0){};
  inline pterm3f(float c, echar e0, echar e1, echar e2): coeff(c),e0(e0),e1(e1),e2(e2){};
#endif
} PT3fData;

typedef struct pterm4f {
  float coeff;
  echar e0, e1, e2, e3;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm4f(): coeff(0),e0(0),e1(0),e2(0),e3(0){};
  inline pterm4f(float c, echar e0, echar e1, echar e2, echar e3): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3){};
#endif
} PT4fData;

typedef struct pterm5f {
  float coeff;
  echar e0, e1, e2, e3, e4;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm5f(): coeff(0),e0(0),e1(0),e2(0),e3(0),e4(0){};
  inline pterm5f(float c, echar e0, echar e1, echar e2, echar e3, echar e4): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3),e4(e4){};
#endif
} PT5fData;

typedef struct pterm6f {
  float coeff;
  echar e0, e1, e2, e3, e4, e5;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm6f(): coeff(0),e0(0),e1(0),e2(0),e3(0),e4(0),e5(0){};
  inline pterm6f(float c, echar e0, echar e1, echar e2, echar e3, echar e4, echar e5): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3),e4(e4),e5(e5){};
#endif
} PT6fData;





typedef struct pterm1d {
  double coeff;
  echar e0;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm1d(): coeff(0),e0(0){};
  inline pterm1d(double c, echar e0): coeff(c),e0(e0){};
#endif
} PT1dData;

typedef struct pterm2d {
  double coeff;
  echar e0, e1;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm2d(): coeff(0),e0(0),e1(0){};
  inline pterm2d(double c, echar e0, echar e1): coeff(c),e0(e0),e1(e1){};
#endif
} PT2dData;

typedef struct pterm3d {
  double coeff;
  echar e0, e1, e2;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm3d(): coeff(0),e0(0),e1(0),e2(0){};
  inline pterm3d(double c, echar e0, echar e1, echar e2): coeff(c),e0(e0),e1(e1),e2(e2){};
#endif
} PT3dData;

typedef struct pterm4d {
  double coeff;
  echar e0, e1, e2, e3;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm4d(): coeff(0),e0(0),e1(0),e2(0),e3(0){};
  inline pterm4d(double c, echar e0, echar e1, echar e2, echar e3): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3){};
#endif
} PT4dData;

typedef struct pterm5d {
  double coeff;
  echar e0, e1, e2, e3, e4;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm5d(): coeff(0),e0(0),e1(0),e2(0),e3(0),e4(0){};
  inline pterm5d(double c, echar e0, echar e1, echar e2, echar e3, echar e4): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3),e4(e4){};
#endif
} PT5dData;

typedef struct pterm6d {
  double coeff;
  echar e0, e1, e2, e3, e4, e5;
#ifdef __cplusplus     // add a constructor, if C++
  inline pterm6d(): coeff(0),e0(0),e1(0),e2(0),e3(0),e4(0),e5(0){};
  inline pterm6d(double c, echar e0, echar e1, echar e2, echar e3, echar e4, echar e5): coeff(c),e0(e0),e1(e1),e2(e2),e3(e3),e4(e4),e5(e5){};
#endif
} PT6dData;






#endif
