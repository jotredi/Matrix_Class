// Matrix class
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

template <class X>
class matrix{
public:
  //Constructors
  matrix(unsigned, unsigned, X = 0.0);
  matrix(unsigned);
  matrix(void);
  //Destructor
  ~matrix();

  //Resize
  void resize(unsigned, unsigned, X = 0.0);

  //Rows & Cols
  void setRow(unsigned, std::vector<X>);
  void setCol(unsigned, std::vector<X>);
  std::vector<X> getRow(const unsigned &);
  std::vector<X> getCol(const unsigned &);
  X sumRow(const unsigned &);
  X sumCol(const unsigned &);

  //Matrix Operations
  matrix<X> operator+(matrix<X>);
  matrix<X> operator-(matrix<X>);
  matrix<X> operator*(matrix<X>);
  matrix<X> operator&(matrix<X>);
  matrix<X> T();

  //Scalar Operations
  matrix<X> operator+(X);
  matrix<X> operator-(X);
  matrix<X> operator*(X);
  matrix<X> operator/(X);
  matrix<X> operator^(X);

  //Other methods
  matrix<X> cofactor(unsigned, unsigned);
  matrix<double> inverse();
  matrix<X> adjoint();
  X det();
  X sum();
  X max();
  X min();

  //Access methods
  X& operator()(const unsigned &, const unsigned &);
  unsigned getRows() const;
  unsigned getCols() const;
  void print() const;

private:
  unsigned rows, cols;
  std::vector<std::vector<X> > M;
};

#include "Matrix.cpp"

#endif
