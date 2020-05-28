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
  matrix(unsigned, unsigned, const X& = 0.0);
  matrix(unsigned, unsigned, const std::string);
  matrix(unsigned);
  matrix(void);
  //Destructor
  virtual ~matrix(void);

  //Change
  void resize(unsigned, unsigned, const X& = 0.0);
  void insertRow(unsigned, const X& = 0.0);
  void insertCol(unsigned, const X& = 0.0);

  //Rows & Cols
  void setRow(unsigned, const std::vector<X> &);
  void setCol(unsigned, const std::vector<X> &);
  std::vector<X> getRow(const unsigned &) const;
  std::vector<X> getCol(const unsigned &) const;
  X sumRow(const unsigned &) const;
  X sumCol(const unsigned &) const;

  //Matrix Operations
  matrix<X> operator+(const matrix<X> &) const;
  matrix<X> operator-(const matrix<X> &) const;
  matrix<X> operator*(const matrix<X> &) const;
  matrix<X> operator&(const matrix<X> &) const;

  void operator+=(const matrix<X> &);
  void operator-=(const matrix<X> &);

  //Scalar Operations
  matrix<X> operator+(const X&) const;
  matrix<X> operator-(const X&) const;
  matrix<X> operator*(const X&) const;
  matrix<X> operator/(const X&) const;
  matrix<X> operator^(const X&) const;

  void operator+=(const X&);
  void operator-=(const X&);
  void operator*=(const X&);
  void operator/=(const X&);
  void operator^=(const X&);

  //Math methods
  X sum() const;
  X max() const;
  X min() const;
  double mean() const;
  matrix<double> ln() const;
  matrix<double> sqr() const;

  //Other methods
  matrix<X> cofactor(const unsigned &, const unsigned &) const;
  matrix<double> inverse() const;
  matrix<X> adjoint() const;
  matrix<X> T() const;
  X det() const;

  //Access methods
  X& operator()(const unsigned &, const unsigned & = 0);
  const X& operator()(const unsigned &, const unsigned & = 0) const;
  unsigned getRows() const;
  unsigned getCols() const;
  void print() const;

private:
  unsigned rows, cols;
  std::vector<std::vector<X> > M;
};

#include "Matrix.cpp"

#endif
