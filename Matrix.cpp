/****************************************
Matrix Class
============
Copyright (C) 2020  Jose Trescastro Diaz
****************************************/

#ifndef MATRIX_CPP
#define MATRIX_CPP

#include "Matrix.hpp"

//=================== CONSTRUCTORS ===================//
//Usual constructor
template <class X>
matrix<X>::matrix(unsigned r, unsigned c, const X& initial){
  rows = r;
  cols = c;
  M.resize(r);
  for (unsigned i=0; i<M.size(); ++i)  M[i].resize(c, initial);
}

//Read matrix from file
template <class X>
matrix<X>::matrix(unsigned r, unsigned c, const std::string file){
  // Initialize matrix dimmensions
  rows = r;
  cols = c;
  M.resize(r);
  for (unsigned i=0; i<M.size(); ++i)  M[i].resize(c);

  // Read file
  std::ifstream data_file(file);
  std::istream_iterator<X> start(data_file), end;
  std::vector<X> data(start, end);

  // Check dimmensions of file
  if(data.size() != rows*cols){
    throw std::logic_error( "Dimmensions are not correct" );
  }

  // Load data
  unsigned pos = 0;
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      M[i][j] = data[pos++];
    }
  }
}

//Identity matrix of size n
template <class X>
matrix<X>::matrix(unsigned n){
  rows = n;
  cols = n;
  M.resize(n);
  for (unsigned i=0; i<M.size(); ++i){
    M[i].resize(n,0.0);
    M[i][i]=1;
  }
}

//Empty Matrix
template <class X>
matrix<X>::matrix(void){
  rows=0;
  cols=0;
}

//Destructor
template <class X>
matrix<X>::~matrix(void){}


//==================== CHANGE METHODS ===================//
template <class X>
void matrix<X>::resize(unsigned r, unsigned c, const X& initial){
  rows = r;
  cols = c;
  M.resize(r);
  for (unsigned i=0; i<M.size(); ++i)  M[i].resize(c, initial);
}

//Insert row
template <class X>
void matrix<X>::insertRow(unsigned r, const X& value){
  M.insert(M.begin() + r, std::vector<X> (cols, value));
  rows += 1;
}

//Insert column
template <class X>
void matrix<X>::insertCol(unsigned c, const X& value){
  for(unsigned i=0; i<rows; ++i){
    M[i].insert(M[i].begin() + c, value);
  }
  cols += 1;
}


//=================== ROWS & COLS ===================//
//Set a row
template <class X>
void matrix<X>::setRow(unsigned i, const std::vector<X> &v){
  if(cols==v.size())  M[i]=v;
  else  throw std::invalid_argument( "setRow: Invalid row size" );
}

//Set a column
template <class X>
void matrix<X>::setCol(unsigned j, const std::vector<X> &v){
  if(rows==v.size()){
    for(unsigned i=0; i<rows; ++i)  M[i][j]=v[i];
  }
  else  throw std::invalid_argument( "setCol: Invalid column size" );
}

//Get a row
template <class X>
std::vector<X> matrix<X>::getRow(const unsigned &i) const{
  std::vector<X> v = M[i];
  return v;
}

//Get a column
template <class X>
std::vector<X> matrix<X>::getCol(const unsigned &j) const{
  std::vector<X> v(rows);
  for(unsigned i=0; i<rows; ++i) v[i]=M[i][j];
  return v;
}

//Sum of elements from a row
template <class X>
X matrix<X>::sumRow(const unsigned &i) const{
  X sum = 0;
  for (unsigned c=0; c<cols; ++c) sum += M[i][c];
  return sum;
}

//Sum of elements from a column
template <class X>
X matrix<X>::sumCol(const unsigned &j) const{
  X sum = 0;
  for (unsigned r=0; r<rows; ++r) sum += M[r][j];
  return sum;
}


//=============== MATRIX OPERATIONS ===============//
//Addition of two matrices
template <class X>
matrix<X> matrix<X>::operator+(const matrix<X> &B) const{
  if((rows != B.getRows()) || (cols != B.getCols())){
    throw std::logic_error( "Matrix dimensions don't match" );
  }
  matrix<X> sum(rows, cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      sum(i,j) = M[i][j] + B(i,j);
    }
  }
  return sum;
}

//Subtraction of two matrices
template <class X>
matrix<X> matrix<X>::operator-(const matrix<X> &B) const{
  if((rows != B.getRows()) || (cols != B.getCols())){
    throw std::logic_error( "Matrix dimensions don't match" );
  }
  matrix<X> diff(rows, cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      diff(i,j) = M[i][j] - B(i,j);
    }
  }
  return diff;
}

//Matrix multiplication
template <class X>
matrix<X> matrix<X>::operator*(const matrix<X> &B) const{
  //Number of rows and cols for each matrix
  unsigned r1, c1, r2, c2;
  r1 = rows;
  c1 = cols;
  r2 = B.getRows();
  c2 = B.getCols();

  if (c1 == r2){
    //We can perform matrix multiplication
    X sum;
    matrix<X> R(r1,c2);

    for (unsigned i = 0; i<r1; ++i){
      for (unsigned j = 0; j<c2; ++j){
        sum = 0;
        for (unsigned k = 0; k<c1; ++k)  sum += M[i][k]*B(k,j);
        R(i,j) = sum;
      }
    }
    return R;
  }
  else{
    //Matrix dimensions are not the same so we return an error message
    throw std::logic_error( "Matrix dimensions don't match" );
  }
}

//Element-wise multiplication
template <class X>
matrix<X> matrix<X>::operator&(const matrix<X> &B) const{
  if((rows != B.getRows()) || (cols != B.getCols())){
    throw std::logic_error( "Matrix dimensions don't match" );
  }
  matrix<X> R(rows, cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      R(i,j) = M[i][j] * B(i,j);
    }
  }
  return R;
}


//Cumulative addition
template <class X>
void matrix<X>::operator+=(const matrix<X> &B){
  if((rows != B.getRows()) || (cols != B.getCols())){
    throw std::logic_error( "Matrix dimensions don't match" );
  }
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      M[i][j] += B(i,j);
    }
  }
}

//Cumulative subtraction
template <class X>
void matrix<X>::operator-=(const matrix<X> &B){
  if((rows != B.getRows()) || (cols != B.getCols())){
    throw std::logic_error( "Matrix dimensions don't match" );
  }
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      M[i][j] -= B(i,j);
    }
  }
}


//=============== SCALAR OPERATIONS ===============//
//Scalar Addition
template <class X>
matrix<X> matrix<X>::operator+(const X& s) const{
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]+s;
  }
  return R;
}

//Scalar Subtraction
template <class X>
matrix<X> matrix<X>::operator-(const X& s) const{
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]-s;
  }
  return R;
}

//Scalar Multiplication
template <class X>
matrix<X> matrix<X>::operator*(const X& s) const{
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]*s;
  }
  return R;
}

//Scalar Division
template <class X>
matrix<X> matrix<X>::operator/(const X& s) const{
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]/s;
  }
  return R;
}

//Power
template <class X>
matrix<X> matrix<X>::operator^(const X& p) const{
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = pow(M[i][j],p);
  }
  return R;
}


//Cumulative Scalar Addition
template <class X>
void matrix<X>::operator+=(const X& s){
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) M[i][j] += s;
  }
}

//Cumulative Scalar Subtraction
template <class X>
void matrix<X>::operator-=(const X& s){
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) M[i][j] -= s;
  }
}

//Cumulative Scalar Multiplication
template <class X>
void matrix<X>::operator*=(const X& s){
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) M[i][j] *= s;
  }
}

//Cumulative Scalar Division
template <class X>
void matrix<X>::operator/=(const X& s){
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) M[i][j] /= s;
  }
}

//Cumulative Power
template <class X>
void matrix<X>::operator^=(const X& p){
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) M[i][j] = pow(M[i][j],p);
  }
}


//=============== MATH METHODS ===============//
//Sum of all elements
template <class X>
X matrix<X>::sum() const{
  X sum = 0;
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) sum += M[i][j];
  }
  return sum;
}

//Max value
template <class X>
X matrix<X>::max() const{
  X max = - std::numeric_limits<X>::infinity();
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      if (M[i][j] > max)  max = M[i][j];
    }
  }
  return max;
}

//Min value
template <class X>
X matrix<X>::min() const{
  X min = std::numeric_limits<X>::infinity();
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      if (M[i][j] < min)  min = M[i][j];
    }
  }
  return min;
}

// Mean
template <class X>
double matrix<X>::mean() const{return (double) sum()/(rows*cols);}

// Natural Log
template <class X>
matrix<double> matrix<X>::ln() const{
  matrix<double> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      R(i,j) = log(M[i][j]);
    }
  }
  return R;
}

// Square Root
template <class X>
matrix<double> matrix<X>::sqr() const{
  matrix<double> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      R(i,j) = sqrt(M[i][j]);
    }
  }
  return R;
}


//=============== OTHER METHODS ===============//
//Cofactor of a Matrix
template <class X>
matrix<X> matrix<X>::cofactor(const unsigned &r, const unsigned &c) const{
  if(rows != cols){
    throw std::logic_error( "Matrix must be a square matrix" );
  }
  // New matrix
  matrix<X> C(rows-1, cols-1);
  // New indexes
  unsigned i = 0, j = 0;
  for(unsigned row=0; row<rows; ++row){
    for(unsigned col=0; col<cols; ++col){
      if ((row != r) && (col != c)){
        C(i,j++) = M[row][col];
        // When row filled
        if(j == cols-1){
          j = 0;
          i++;
        }
      }
    }
  }
  return C;
}

//Inverse
template <class X>
matrix<double> matrix<X>::inverse() const{
  // Check if determinat is 0 (Singular matrix)
  X d = det();
  if (d == 0){
    throw std::logic_error( "Inverse doesn't exist" );
  }
  // Inverse matrix
  matrix<double> inv(rows,cols);
  // Adjoint
  matrix<X> adj = adjoint();
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      inv(i,j) = (double) adj(i,j)/d;
    }
  }
  return inv;
}

//Adjoint
template <class X>
matrix<X> matrix<X>::adjoint() const{
  if(rows != cols){
    throw std::logic_error( "Matrix must be a square matrix" );
  }
  // Adjoint matrix
  matrix<X> A(rows, cols);
  int sign = 1;

  if(rows==1){
    A(0,0) = 1;
    return A;
  }
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      sign = ((i+j)%2==0)? 1: -1;
      A(j,i) = sign * cofactor(i,j).det();
    }
  }
  return A;
}

//Matrix transpose
template <class X>
matrix<X> matrix<X>::T() const{
  matrix<X> transposed(cols,rows);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) transposed(j,i)=M[i][j];
  }
  return transposed;
}

//Determinant of a Matrix
template <class X>
X matrix<X>::det() const{
  if(rows != cols){
    throw std::logic_error( "Matrix must be a square matrix" );
  }
  if(rows==1) return M[0][0];
  X d = 0; // determinat
  int sign = 1;
  // For each element of first row
  for (unsigned j = 0; j<cols; ++j){
    d += sign * M[0][j] * cofactor(0,j).det();
    sign = -sign; // alternate sign
  }
  return d;
}


//=============== ACCESS METHODS ===============//
//Matrix access in the form M(i,j)
template <class X>
X& matrix<X>::operator()(const unsigned &i, const unsigned &j){return M[i][j];}

//Matrix access in the form M(i,j) for constant
template <class X>
const X& matrix<X>::operator()(const unsigned &i, const unsigned &j) const{
  return M[i][j];
}

//Get number of rows and columns
template <class X>
unsigned matrix<X>::getRows() const{return rows;}
template <class X>
unsigned matrix<X>::getCols() const{return cols;}

//Print Matrix
template <class X>
void matrix<X>::print() const{
  for (unsigned i=0; i<rows; ++i){
    for (unsigned j=0; j<cols; ++j)  std::cout << M[i][j] << " ";
    std::cout << std::endl;
  }
}

#endif
