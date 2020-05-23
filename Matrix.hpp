/*********************************************************************
Matrix Class
============
Copyright (C) 2020  Jose Trescastro Diaz

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*********************************************************************/
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <limits>
#include <fstream>
#include <iterator>

template <class X>
class matrix{
public:
  //Constructors
  matrix(unsigned, unsigned, X = 0.0);
  matrix(unsigned, unsigned, std::string);
  matrix(unsigned);
  matrix(void);
  //Destructor
  ~matrix(void);

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

  //Math methods
  X sum();
  X max();
  X min();
  double mean();
  matrix<double> ln();

  //Other methods
  matrix<X> cofactor(unsigned, unsigned);
  matrix<double> inverse();
  matrix<X> adjoint();
  X det();

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
