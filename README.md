# Matrix_Class
This is a very useful, extensible and reusable Matrix class with the most common matrix/vector operations and methods for C++.

### Constructors
* `matrix(r, c, initial)`: create a matrix with `r` rows, `c` columns & `initial` as initial value.
* `matrix(n)`: create an identity matrix of size `n`.
* `matrix()`: empty matrix with 0 rows & 0 cols.

Example:
```c++
// Creation of 3x3 matrix with initial value of 1 and type int:
matrix<int> A (3,3,1);

// Identity matrix of size 3:
matrix<int> I (3);

// Creation of a column vector of size 3 and type double:
matrix<double> v (3,1);
// If initial value is not given, it is initialized with 0's.
```

### Access Methods
* `M (i,j)`: access the ith, jth element.
* `getRows()`: return number of rows.
* `getCols()`: return number of columns.
* `print()`: print matrix.

### Rows & Cols Methods
* `resize(r, c, initial)`: resize the matrix with `r` rows, `c` columns & `initial` value.
* `setRow(i, v)`: set the ith row with vector `v`.
* `setCol(j, v)`: set the jth column with vector `v`.
* `getRow(i)`: get ith row.
* `getCol(j)`: get jth column.
* `sumRow(i)`: sum of elements from ith row.
* `sumCol(j)`: sum of elements from jth column.

### Matrix Operations
* `A + B`: matrix sum.
* `A - B`: matrix subtraction.
* `A * B`: matrix multiplication.
* `A & B`: element-wise matrix multiplication.
* `A.T()`: transpose of matrix A.

### Scalar Operation
* `A + s`: scalar sum.
* `A - s`: scalar subtraction.
* `A * s`: scalar multiplication.
* `A / B`: scalar division.
* `A ^ p`: raise each element of A to the power of `p`.

### Other methods
* `inverse()`: return inverse of the matrix (if exists).
* `cofactor(i,j)`: return matrix formed without the ith row & jth column.
* `adjoint()`: return the adjoint matrix (transpose of the cofactor matrix).
* `det()`: determinant of the matrix.
* `sum()`: return the sum of all elements.
* `max()`: return the maximum value.
* `min()`: return the minimum value.

### Reference
* https://medium.com/@furkanicus/how-to-create-a-matrix-class-using-c-3641f37809c7
* https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File/
* https://www.geeksforgeeks.org/adjoint-inverse-matrix/?ref=lbp
