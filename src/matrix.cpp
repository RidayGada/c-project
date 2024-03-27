#include <iostream>   //includes the basic header files
#include <cmath>
#include <climits>
#include <algorithm>
#include <exception>
#include <stdexcept>

#include "matrix.hpp"  //includes the header file matrix.hpp

using namespace simple_matrix;  

// BEGIN CLASS

#define index(i, j) ((n_ * (i)) + (j))
//a macro which calculates the linear index of an element in a 2D array given its row and column indices.
#define check_size(mat) if (!((mat).m_ == m_ && (mat).n_ == n_)) throw bad_size();
//a macro  to check if the dimensions of two matrices are compatible for operations
#define for_ij(m, n) for (uint i = 0; i < (m); i++) for (uint j = 0; j < (n); j++)
//a macro to iterate over the matrix indices

matrix::matrix(uint rows, uint cols) //constructor that takes 2 arguments
		: m_{rows}, n_{cols}, buf_{(m_==0||n_==0) ? nullptr : new double[m_*n_]} { //if eitherrows or cols is zero it sets buf_ to nullptr
	std::fill(buf_, buf_ + (m_*n_), 0); //uses fill to initialize all elements of matrix to zero
}

matrix::matrix(uint rows, uint cols, std::initializer_list<double> list)//constructor that take 3 arguments
		: m_{rows}, n_{cols}, buf_{(m_==0||n_==0) ? nullptr : new double[m_*n_]} {
	if (list.size() != (m_*n_)) //if size of initializer list and matrix is not same it throws exception
		throw std::invalid_argument{"List must be of same size as matrix"};
	std::copy(list.begin(), list.end(), buf_); //copies list elements to buf_
}

matrix::matrix(uint rows, uint cols, const double *values) //constructor with 3 arg, 3rd is pointer to an array of double
		: m_{rows}, n_{cols}, buf_{(m_==0||n_==0) ? nullptr : new double[m_*n_]} {
	std::copy((const double *)values, (const double *)(values + (m_*n_)), buf_);
} //copies array elements to buffer

matrix::matrix(uint rows, uint cols, std::initializer_list<int> list)//constructor with 3 arg, 3rd is initializer list of int datatype
		: m_{rows}, n_{cols}, buf_{(m_==0||n_==0) ? nullptr : (new double[m_*n_])} {
	if (list.size() != (m_*n_)) //throws exception if size of initializer list and matrix is unequal
		throw std::invalid_argument{"List must be of same size as matrix"};
	std::copy(list.begin(), list.end(), buf_); //copies elements of list to buf_
}

matrix::matrix(uint rows, uint cols, const int *values) //constr. with 3 args, 3rd is pointer to array of int
		: m_{rows}, n_{cols}, buf_{(m_==0||n_==0) ? nullptr : new double[m_*n_]} {
	std::copy((const int *)values, (const int *)(values + (m_*n_)), buf_);
} //copies array elements to buffer

//these constructors ensure proper memory allocation and initialization of matrix elements.

matrix::matrix(const matrix& mat) //copy constructor
		: m_{mat.m_}, n_{mat.n_}, buf_{(m_==0||n_==0) ? nullptr : new double[m_*n_]} {
	std::copy(mat.buf_, mat.buf_ + (m_*n_), buf_);
}

matrix::matrix(matrix&& mat) //move constructor
		: m_{0}, n_{0}, buf_{nullptr} {
	std::swap(m_, mat.m_);
	std::swap(n_, mat.n_);
	std::swap(buf_, mat.buf_);
}

matrix::matrix(const std::string& matstr) //constructor takes std::string mastr as argument
		: m_{0}, n_{0}, buf_{nullptr} { //initialize m and n to 0 and buff to nullptr
	std::string::const_iterator start = matstr.cbegin();
	std::string::const_iterator end;
	*this = parse_matrix(start, end); //assigns the result of *this to parse_matrix
}

matrix::~matrix() { //destructor
	delete[] buf_; //deallocates memory for buf_ using delete[] operator
}

uint matrix::m() const {
	return m_; //function that returns the no. of rows of matrix
}

uint matrix::n() const {
	return n_; //function that returns the no. of column of matrix
}

double matrix::get(uint i, uint j) const { //function that returns the element at that particular row i and column j
	if (i >= m_ || j >= n_) 
		throw std::out_of_range("Term isn't within matrix"); //throws exception out of range if i or j are out of range
	return buf_[index(i, j)];
}

void matrix::set(uint i, uint j, double value) { //function that sets value of element at row i and col j
	if (i >= m_ || j >= n_)
		throw std::out_of_range("Term isn't within matrix"); //throws exception out of range if i or j are out of range
	buf_[index(i, j)] = value;
}

bool matrix::is_empty() const { //funtion that returns true if matrix is empty 
	return m_ == 0 || n_ == 0;  //or else returns false
}

bool matrix::is_square() const { //function that returns true if matrix is square 
	return m_ == n_;			 //or else returns false
}

bool matrix::is_diagonal() const { //function that returns true if matrix is a diagonal matrix
	if (!is_square())
		return false;
	return is_upper_triangular() && is_lower_triangular(); //checks both if its upper and lower triangular matrix
}

//function that returns true if matrix is upper triangular
bool matrix::is_upper_triangular() const {
	if (!is_square())
		return false;
	for (uint i = 1; i < m_; i++) //all elements below diagonal must be equal to 0
	for (uint j = 0; j < i; j++)
		if (!EQUAL(get(i, j), 0))
			return false;
	return true;
}

//function that returns true if matrix is lower triangular
bool matrix::is_lower_triangular() const {
	if (!is_square())
		return false;
	for (uint j = 1; j < n_; j++) //all elements above diagonal must be equal to 0
	for (uint i = 0; i < j; i++)
		if (!EQUAL(get(i, j), 0))
			return false;
	return true;
}

bool matrix::is_invertible() const {
	return det() != 0; //function returns true if matrix is invertible i.e det is not equal to 0
}

matrix matrix::get_row(uint i) const { // function returns a row vector of the matrix at index i
	matrix a(1, n_); //creates a row matrix 
	for (uint j = 0; j < n_; j++) //assigns values to the row matrix
		a(0, j) = (*this)(i, j);
	return a; //returns row matrix
}

matrix matrix::get_col(uint j) const { //function returns a column vector of the matrix at index j
	matrix a(m_, 1);  //creates a column matrix
	for (uint i = 0; i < m_; i++) //assigns values to the column matrix
		a(i, 0) = (*this)(i, j);
	return a; //returns col matrix
}

void matrix::set_row(uint i, const matrix& row) { //function sets the values of a row in the matrix
	if (row.m_ != 1 || row.n_ != n_)
		throw bad_size(); //if the dimensions don't match, it throws a bad_size exception
	for (int j = 0; j < n_; j++)   //assigns the corresponding value from the 
		(*this)(i, j) = row(0, j); //input row matrix to row of the current matrix
} 

void matrix::set_col(uint j, const matrix& col) { //function sets the values of a col in the matrix
	if (col.n_ != 1 || col.m_ != m_)
		throw bad_size(); //if the dimensions don't match, it throws a bad_size exception
	for (uint i = 0; i < m_; i++)  //assigns the corresponding value from the
		(*this)(i, j) = col(i, 0); //input col matrix to col of the current matrix
}

double matrix::trace() const { //function to calculate trace of the matrix
	if (!is_square())
		throw not_square(); //first checks if matrix is square or not
	double trace = 0;
	for (int i = 0; i < m_; i++)
		trace += (*this)(i, i); //calculates trace i.e the sum of diaginal elements
	return trace; //retruns the value of trace calculated.
}

double matrix::det() const { //function to calculate determinant of matrix
	if (!is_square())
		throw not_square(); //first checks if matrix is square or not
	if (m_ < 2) //checks matrix size if its less than 2 
		throw bad_size(); //throws bad_seize exception if true
	if (m_ == 2) // Might as well speed things up
		return ((*this)(0, 0) * (*this)(1, 1)) - ((*this)(0, 1) * (*this)(1, 0));
	//for 2x2 matrix calculates deteminant directly 
	uint p[n_], v[n_];
	// iterates through all possible permutations of indices to calculate the 
	//determinant using the permutation sign and the product of matrix elements
	permutation_init(n_, p, v);
	double detsum = 0, prod;
	int sgn = 1;
	do {
		prod = 1;
		for (uint i = 0; i < n_; i++)
			prod *= (*this)(i, p[i]);
		prod *= sgn;
		sgn = -sgn;
		detsum += prod;
	} while (permutation_permute(n_, p, v));
	return detsum; //returns the deteminant value
}

matrix matrix::transpose() const { //function that returns transpose of matrix
	matrix mat(n_, m_); //creates new matrix mat with swapped dimensions
	for_ij(m_, n_) //iterates through each element of original matrix and assigns value
		mat(j, i) = (*this)(i, j);  //to the newly created mat matrix
	return mat; //returns transpose of matrix
}

matrix matrix::adj() const { //function returns the adjugate of the matrix
	matrix mat = cofactor_matrix().transpose(); //calculates cofactor matrix using cofactor_matrix() and then transposes it
	return mat; //returns adjugate matrix
}

//function returns a submatrix by removing the i-th row and j-th column
matrix matrix::submatrix(uint i, uint j) const {
	matrix mat(m_ - 1, n_ - 1); //new matrix if dimension reduced by 1
	uint is = 0, js;
	for (uint ip = 0; ip < m_; ip++) 
	{   //iterates through each element of the original matrix, skipping the i-th row and j-th column,
		if (ip == i)                   //and assigns it to the corresponding position in the submatrix
			continue;
		js = 0;
		for (uint jp = 0; jp < n_; jp++) {
			if (jp == j)
				continue;
			mat(is, js) = (*this)(ip, jp);
			++js;
		}
		++is;
	}
	return mat; //returns submatrix
}

//function calculates the determinant of the submatrix obtained by removing row i and column j
double matrix::minordet(uint i, uint j) const {
	if (!is_square())
		throw not_square(); //throws not_square() if matrix is not square
	if (m_ < 3)
		throw bad_size(); //throws a bad_size() if matrix size is less than 3
	return submatrix(i, j).det(); //returns the determinant value of the submatrix
}

//function calculates the cofactor of the row i and column j
double matrix::cofactor(uint i, uint j) const {
	double deter = minordet(i, j);         //determines the sign of the cofactor
	return ((i + j) % 2) ? -deter : deter; //and returns the calculated cofactor value
}

matrix matrix::minor_matrix() const { //function returns the minor matrix of the square matrix
	if (!is_square())
		throw not_square(); //checks if square or not
	matrix mat(m_, m_); //creates new matrix mat
	for (uint i = 0; i < m_; i++) //iterates through each element of the original matrix and calculates
	for (uint j = 0; j < m_; j++) //the corresponding minor determinant using the minordet function
		mat(i, j) = minordet(i, j);
	return mat; //returns the minor matrix
}

matrix matrix::cofactor_matrix() const { //function returns the cofactor matrix of the square matrix
	matrix mat = minor_matrix(); //calculates minor of matrix mat
	int alt = 1;
	for (uint k = 0; k < (m_ * n_); k++) {
		mat.buf_[k] *= alt; //iterates through each element of the 
		alt = -alt;         //minor matrix and applies the alternating sign pattern to each element
	}
	return mat; //returns the cofactor matrix
}

matrix matrix::invert() const { //function returns the inverse of the square matrix
	if (!is_square())
		throw not_square(); //checks if its square or not
	double deter = det();
	if (deter == 0)
		throw not_invertible(); //throws exception if not invertible
	matrix mat;
	if (m_ == 2) { //for a 2x2 matrix, it calculates the inverse directly and returns it
		mat = matrix(2, 2);
		mat(0, 0) = (*this)(1, 1);
		mat(1, 1) = (*this)(0, 0);
		mat(0, 1) = -(*this)(0, 1);
		mat(1, 0) = -(*this)(1, 0);
	} else
		mat = adj(); //for greater sized matrices it calculate adjugate 
	mat /= deter;    //and divides it by deteminant to calculate inverse matrix
	return mat;      //returns inverse matrix
}

// Reimplement?
//function solves a system of linear equations represented by the matrix equation Ax = b
matrix matrix::solve(const matrix& ans) const {
	if (!is_square())
		throw not_square(); //checks if matrix is not square
	if (ans.m_ != m_ || ans.n_ != 1)
		throw bad_size(); //checks if the dimensions of the solution vector ans match with the matrix
	double deter = det(); //calculates the determinant of the matrix
	if (deter == 0)
		throw not_solvable(); //throws not solvable if determinant is equal to 0
	//solving the system of equations
	matrix res(m_, 1);
	matrix mat;
	for (uint j = 0; j < n_; j++) {
		mat = *this;
		mat.set_col(j, ans);
		res(j, 0) = mat.det() / deter;
	}
	return res; //returns the result matrix containing the solutions
}

void matrix::swap(matrix& other) { //function swaps the contents of the current
	std::swap(m_, other.m_);	   //matrix with another matrix passed as a reference
	std::swap(n_, other.n_);
	std::swap(buf_, other.buf_); //also swaps the internal buffer buf_ which stores the matrix elements
}

//these overloaded function call operators allow accessing 
//the elements of the matrix using the (i, j) syntax
double& matrix::operator()(uint i, uint j) {
	if (i >= m_ || j >= n_) 
		throw std::out_of_range("Term isn't within matrix"); //throws exception if out of range
	return buf_[index(i, j)];
}

double matrix::operator()(uint i, uint j) const {
	if (i >= m_ || j >= n_)
		throw std::out_of_range("Term isn't within matrix"); //throws exception if out of range
	return buf_[index(i, j)];
}

//these overloaded assignment operators assign one matrix to another
matrix& matrix::operator=(const matrix& a) { //for assignment from another matrix
	matrix tmp(a);
	this->swap(tmp);
	return *this;
}

matrix& matrix::operator=(matrix&& a) { //for move assignment from an rvalue matrix
	this->swap(a);
	return *this;
}

matrix matrix::operator-() { //this overloaded unary minus operator negates all elements of the matrix
	matrix a(m_, n_); //new matrix a
	for_ij(m_, n_) buf_[index(i, j)] *= -1;
	return a; //returns the newly created matrix a
}

matrix& matrix::operator+=(const matrix& a) { //overloaded compound assignment operator which adds
	check_size(a);							  //matrix to another matrix
	for_ij(m_, n_) buf_[index(i, j)] += a.buf_[index(i, j)];
	return *this;
}

matrix& matrix::operator-=(const matrix& a) { //overloaded compound assignment operator which subtracts
	check_size(a);							  //matrix to another matrix
	for_ij(m_, n_) buf_[index(i, j)] -= a.buf_[index(i, j)];
	return *this;
}

matrix& matrix::operator*=(const matrix& a) { //overloaded compound assignment operator
	matrix mat = a * *this;					  //to multiply matrix by another matrix
	*this = mat;
	return *this;
}

matrix& matrix::operator*=(double a) {     //overloaded compound assignment operator
	for_ij(m_, n_) buf_[index(i, j)] *= a; //to multiply matrix by a scalar value
	return *this;
}

matrix& matrix::operator/=(double a) {     //overloaded compound assignment operator
	for_ij(m_, n_) buf_[index(i, j)] /= a; //to divide matrix by a scalar value
	return *this;
}

bool matrix::operator==(const matrix& a) { //overloaded comparison operator 
	if (a.m_ != m_ || a.n_ != n_) //checks if the dimensions of both matrices are the same ot not
		return false;
	for (uint k = 0; k < (m_ * n_); k++) //iterates through all elements of the matrices and
		if (!EQUAL(a.buf_[k], buf_[k]))  //checks if they are equal using the EQUAL macro
			return false; //returns false if not equal
	return true; //returns true if equal
}

bool matrix::operator!=(const matrix& a) { 
	return !operator==(a);//inequality operator simply negates the result of the equality operator
}

//here we undefine the previously declared macros to avoid name clashes
#undef index 
#undef check_size
#undef for_ij

// END CLASS

namespace simple_matrix {

	matrix identity_matrix(uint m) { //function to genearte indentity matrix of mxm size
		matrix mat(m, m);
		for (uint k = 0; k < m; k++)
			mat(k, k) = 1; //diagonal elements equal to 1
		return mat;
	}

//function overloads the output stream insertion operator (<<) to allow streaming 
//of matrix objects to output streams, such as std::cout
	std::ostream& operator<<(std::ostream& out, const matrix& a) { //takes 2 args
		out << '[';
		for (uint i = 0; i < a.m(); ++i) { //iterates through each row and column of matrix 
			if (i)
				out << ";  ";
			for (uint j = 0; j < a.n(); ++j) {
				if (j)
					out << ',';
				out << ' ' << a(i, j);
			}
		}
		out << ' ' << ']';
		return out; //basicaly gives us human readable form of matrix
	}

//function overloads the input stream extraction operator (>>) to allow 
//reading matrix objects from input streams, such as std::cin
	std::istream& operator>>(std::istream& in, matrix& a) {
		std::istream_iterator<char> start(in);
		std::istream_iterator<char> end;
		a = parse_matrix(start, end);
		//calls a function parse_matrix with the iterators as arguments
		//to parse the matrix from the input stream
		return in;
	}

	matrix operator+(const matrix& a, const matrix& b) {
		matrix mat(a); //overloaded arithmetic operator to add two matrix objects
		mat += b;
		return mat;
	}

	matrix operator-(const matrix& a, const matrix& b) {
		matrix mat(a); //overloaded arithmetic operator to subtract two matrix objects
		mat -= b;
		return mat;
	}

	matrix operator*(const matrix& a, double b) {
		matrix mat(a); //overloaded arithmetic operator to multiply matrix by scalar value
		mat *= b;
		return mat;
	}

	matrix operator*(double a, const matrix& b) {
		return b * a; 	//overloaded arithmetic operator to multiply matrix by scalar value
	}

	//overloaded arithmetic operator allows multiplication of two matrices
	matrix operator*(const matrix& a, const matrix& b) {
		if (a.n() != b.m())
			throw bad_size();
		matrix mat(a.m(), b.n());
		for (uint i = 0; i < a.m(); i++) //multiplies each column to each row
		for (uint j = 0; j < b.n(); j++) {
			double sum = 0;
			for (uint k = 0; k < a.n(); k++)
				sum += a(i, k) * b(k, j);
			mat(i, j) = sum;
		}
		return mat; //returns the resultant matrix
	}

	//overloaded arithmetic operator allows division of a matrix by a scalar value of double datatype
	matrix operator/(const matrix& a, double b) {
		matrix mat(a);
		mat /= b;
		return mat; //returns resultant matrix
	}

}
