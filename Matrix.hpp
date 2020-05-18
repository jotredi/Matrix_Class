// Matrix class

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
  X& operator()(const unsigned &, const unsigned &);
  unsigned getRows() const;
  unsigned getCols() const;
  void print() const;
  X sum();
  X max();
  X sumRow(const unsigned &);
  X sumCol(const unsigned &);

private:
  unsigned rows, cols;
  std::vector<std::vector<X> > M;
};

//Usual constructor
template <class X>
matrix<X>::matrix(unsigned r, unsigned c, X initial){
  rows = r;
  cols = c;
  M.resize(r);
  for (unsigned i=0; i<M.size(); ++i)  M[i].resize(c, initial);
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

//Resize
template <class X>
void matrix<X>::resize(unsigned r, unsigned c, X initial){
  rows = r;
  cols = c;
  M.resize(r);
  for (unsigned i=0; i<M.size(); ++i)  M[i].resize(c, initial);
}

//Set a row
template <class X>
void matrix<X>::setRow(unsigned i, std::vector<X> v){
  if(cols==v.size())  M[i]=v;
  else  throw std::invalid_argument( "setRow: Invalid row size" );
}
//Set a column
template <class X>
void matrix<X>::setCol(unsigned j, std::vector<X> v){
  if(rows==v.size()){
    for(unsigned i=0; i<rows; ++i)  M[i][j]=v[i];
  }
  else  throw std::invalid_argument( "setCol: Invalid column size" );
}

//Get a row
template <class X>
std::vector<X> matrix<X>::getRow(const unsigned &i){
  std::vector<X> v = M[i];
  return v;
}
//Get a column
template <class X>
std::vector<X> matrix<X>::getCol(const unsigned &j){
  std::vector<X> v(rows);
  for(unsigned i=0; i<rows; ++i) v[i]=M[i][j];
  return v;
}

//Addition of two matrices
template <class X>
matrix<X> matrix<X>::operator+(matrix<X> B){
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

//Substraction of two matrices
template <class X>
matrix<X> matrix<X>::operator-(matrix<X> B){
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
matrix<X> matrix<X>::operator*(matrix<X> B){
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
    //Matrix dimenssions are not the same so we return an error message
    throw std::logic_error( "Matrix dimensions don't match" );
  }
}

//Element-wise multiplication
template <class X>
matrix<X> matrix<X>::operator&(matrix<X> B){
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

//Matrix transpose
template <class X>
matrix<X> matrix<X>::T(){
  matrix<X> transposed(cols,rows);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) transposed(j,i)=M[i][j];
  }
  return transposed;
}

//Scalar Addition
template <class X>
matrix<X> matrix<X>::operator+(X s){
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]+s;
  }
  return R;
}

//Scalar Substraction
template <class X>
matrix<X> matrix<X>::operator-(X s){
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]-s;
  }
  return R;
}

//Scalar Multiplication
template <class X>
matrix<X> matrix<X>::operator*(X s){
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]*s;
  }
  return R;
}

//Scalar Division
template <class X>
matrix<X> matrix<X>::operator/(X s){
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = M[i][j]/s;
  }
  return R;
}

//Power
template <class X>
matrix<X> matrix<X>::operator^(X p){
  matrix<X> R(rows,cols);
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) R(i,j) = pow(M[i][j],p);
  }
  return R;
}

//Sum of all elements
template <class X>
X matrix<X>::sum(){
  X sum = 0;
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j) sum += M[i][j];
  }
  return sum;
}

//Max value
template <class X>
X matrix<X>::max(){
  X max = - std::numeric_limits<X>::infinity();
  for(unsigned i=0; i<rows; ++i){
    for(unsigned j=0; j<cols; ++j){
      if (M[i][j] > max)  max = M[i][j];
    }
  }
  return max;
}

//Sum of elements from a row
template <class X>
X matrix<X>::sumRow(const unsigned &i){
  X sum = 0;
  for (unsigned c=0; c<cols; ++c) sum += M[i][c];
  return sum;
}
//Sum of elements from a column
template <class X>
X matrix<X>::sumCol(const unsigned &j){
  X sum = 0;
  for (unsigned r=0; r<rows; ++r) sum += M[r][j];
  return sum;
}

//Matrix access in the form M(i,j)
template <class X>
X& matrix<X>::operator()(const unsigned &i, const unsigned &j){return M[i][j];}

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
