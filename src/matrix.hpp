#ifndef SIMPLE_MATRIX_MATRIX_H_
#define SIMPLE_MATRIX_MATRIX_H_

#include <iostream>  //these 3 are basic c++ header files included
#include <iterator>
#include <vector>

// Error used for comparisons
#ifndef EPSILON
#define EPSILON 0.0000000001  //defines a macro which is of float
#endif
#ifndef EQUAL
#define EQUAL(a, b) (abs((a) - (b)) < EPSILON)  //another macro to basically check equality of two no.s such that
#endif											//their difference less than epsilon macro deifned above

namespace simple_matrix {  //defined a namespace here for various user defined exceptions
	struct bad_size : public std::exception {  //user defined exception for matrix size
		virtual const char *what() const throw() {  //if matrix size is not correct
			return "matrix/matricies not compatible sizes";  //returns that size is not compatible
		}
	};
	struct not_square : public bad_size {  //user defined exception is matrix is not square
		const char *what() const throw() {
			return "matrix must be square";  //returns that matrix should be a square matrix
		}
	};
	struct not_invertible : public std::exception {  //user defined exception for matrix inverse
		const char *what() const throw() {  //if matrix is not invertible
			return "matrix is not invertible";  //returns invertibility of the matrix
		}
	};
	struct not_solvable : public std::exception {  //user defined exception for unsolvable systems 
		const char *what() const throw() {
			return "System is not solvable";  //returns system is unsolvable
		}
	};

	typedef unsigned int uint;  //defined another way to define an unsigned int as 'uint'

	class matrix {  //created a class named matrix 
	private:  //for private members of the class;
		// m_ and n_ store the number of rows
		// and columns respectively.
		uint m_, n_;  //2 unsigned integers 
		// buf_ stores each term in the matrix
		double *buf_;  //a pointer declaration for a double variable
		// Permutation functions, used for determintant calculation
		static inline uint permutation_position(uint n, uint i, uint *p);  //decl. of an inline function which is of returntype unsigned int
		static void permutation_init(uint n, uint *p, uint *v);  //a void return type function declaration
		static int permutation_permute(uint n, uint *p, uint *v);  // a int returntype func decl.
	public:
		typedef double* iterator;  //these are two aliases for iterationg over the matrix elements
		typedef const double* const_iterator;

		// Basic
		// =====

		// Default constructor, creates empty matrix
		matrix()  //default constructor (has same name as that of class)
			: m_{0}, n_{0}, buf_{nullptr} {}  //initializing member variables to default values

		// Creates an empty matrix
		matrix(uint rows, uint cols); //constructor for creating matrix with specified rows and columns

		// Creates a matrix and fills horizontally with
		// values from values
		matrix(uint rows, uint cols, std::initializer_list<double>); //constructors for filling matrix
		matrix(uint rows, uint cols, const double *values);

		// Creates a matrix and fills horizontally with
		// integer values from values
		matrix(uint rows, uint cols, std::initializer_list<int>);   //constructors for filling matrix
		matrix(uint rows, uint cols, const int *values);

		// Copy constructor
		matrix(const matrix&);
		
		// Move constructor
		matrix(matrix&&);  //prevents unnecessary copying of data 

		// String initialiser
		matrix(const std::string&);  //constructor that initializes a matrix with string representation

		~matrix(); //destructor

		// Pretty prints the matrix
		std::string pretty() const; //prints matrix in such a way, we can read

		// Get number of rows
		uint m() const;

		// Get number of columns
		uint n() const;

		// Gets the M_ij element
		double get(uint i, uint j) const;  //gets the values of matrix

		// Sets the M_ij element
		void set(uint i, uint j, double value);  //sets matrix with values got using get

		iterator begin() {
			return &buf_[0]; //for mutable access to beginning of matrix
		}

		iterator end() {
			return &buf_[m_ * n_];  //for mutable access to ending of matrix
		}

		const_iterator cbegin() const {
			return &buf_[0];  //for immutable access to beginning of matrix
		}

		const_iterator cend() const {
			return &buf_[m_ * n_];  //for immutable access to ending of matrix
		}

		// Returns true if either m or n are zero, useful
		// for returning an undefined value
		bool is_empty() const;

		// Returns true if matrix is a square matrix
		bool is_square() const;

		// Returns true if matrix is a diagonal matrix
		bool is_diagonal() const;

		// Returns true if matrix is an upper triangular matrix
		bool is_upper_triangular() const;

		// Returns true if matrix is a lower triangular matrix
		bool is_lower_triangular() const;

		// Returns true if matrix is invertible
		// that is, M.det() != 0 , determinant of matrix should be unequal to zero
		bool is_invertible() const;

		// Gets the ith row as a row vector
		matrix get_row(uint i) const;

		// Gets the jth col as a column vector
		matrix get_col(uint j) const;

		// Replaces ith row with a row vector
		void set_row(uint i, const matrix& row);

		// Replaces jth col with a column vector
		void set_col(uint j, const matrix& col);

		// Advanced
		// ========

		// Calculates trace
		double trace() const;  //calculates sum ofdiagonal elements of the matrix

		// Calculates determinant of the matrix
		double det() const;  

		// Transposes matrix
		matrix transpose() const;

		// Adjugates matrix , which can be used to fing the inverse of matrix
		matrix adj() const;

		// Creates a sub matrix, excluding row i and col j
		matrix submatrix(uint i, uint j) const;

		// Calculates minor M_ij
		double minordet(uint i, uint j) const;

		// Calculates cofactor C_ij
		double cofactor(uint i, uint j) const;

		// Creates a matrix of minors
		matrix minor_matrix() const;

		// Creates a matrix of cofactors
		matrix cofactor_matrix() const;

		// Inverts a square matrix
		matrix invert() const;

		// Solves a system of equations in a square matrix
		matrix solve(const matrix& ans) const;

		void swap(matrix& other);  //swaps the contents of two matrices

		// Operator Overloads
		// ==================

		// Element access operators
		double& operator()(uint i, uint j);
		double operator()(uint i, uint j) const; //these operator overloads are for individual access to the matrix elements


		//These are operator overloads for various matrix operations such as assignment,
		//unary negation, add, subtract, multiply, divide, equality, and inequality.
		matrix& operator=(const matrix&); // Copy assignment
		matrix& operator=(matrix&&); // Move assignment
		matrix operator-();
		matrix& operator+=(const matrix&);
		matrix& operator-=(const matrix&);
		matrix& operator*=(const matrix&);
		matrix& operator*=(double);
		matrix& operator/=(double);
		bool operator==(const matrix&);
		bool operator!=(const matrix&);

		friend std::ostream& operator<<(std::ostream& out, const matrix& a);  //friend function that has access to private members of the matrix class
	};

	// 0x0 matrix
	const matrix EMPTY_MATRIX;

	// Creates an identity matrix with the given size
	matrix identity_matrix(uint m);


	//These lines declare standalone operator overloads for input/output streaming
	// (<< and >>) and arithmetic operations (+, -, *, /) involving matrices and scalar values.
	std::ostream& operator<<(std::ostream&, const matrix&);
	std::istream& operator>>(std::istream&, matrix&);
	matrix operator+(const matrix&, const matrix&);
	matrix operator-(const matrix&, const matrix&);
	matrix operator*(const matrix&, double);
	matrix operator*(double, const matrix&);
	matrix operator*(const matrix&, const matrix&);
	matrix operator/(const matrix&, double);

	// STL style parser

	template <class InputIterator>
	simple_matrix::matrix parse_matrix(InputIterator& start, const InputIterator& end) {
		uint m = 0;
		uint n = 0;
		uint cn = 0;
		char ch;

		enum {VOID, CAPTURE} STATE = VOID;

		std::vector<char> inp;
		std::vector<double> mat;

		InputIterator& i = start;
		bool cancel = false;
		while (!(i == end || cancel)) {
			switch (STATE) {
			case VOID:
				switch (*i) {
				case '[':
				case '(':
					STATE = CAPTURE;
					break;
				case ' ':
				case '\n':
				case '\t':
				case '\v':
				case '\f':
				case '\r':
					break;
				default:
					// Unexpected value
					throw std::invalid_argument{"Unexpected character encountered: '" + *i + '"'};
				}
				break;
			case CAPTURE:
				switch (*i) {
				case ',':
					inp.push_back('\0');
					mat.push_back(atof(inp.data()));
					inp.clear();
					cn++;
					break;
				case ']':
					inp.push_back('\0');
					mat.push_back(atof(inp.data()));
					inp.clear();
					cn++;
					m++;
					if (!(n == 0 || cn == n)) {
						throw std::invalid_argument{"Column lengths must be consistent"};
					}
					n = cn;
					cn = 0;
					STATE = VOID;
					cancel = true;
					break;
				case ';':
					inp.push_back('\0');
					mat.push_back(atof(inp.data()));
					inp.clear();
					cn++;
					m++;
					if (!(n == 0 || cn == n)) {
						throw std::invalid_argument{"Column lengths must be consistent"};
					}
					n = cn;
					cn = 0;
					break;
				default:
					inp.push_back(*i);
				}
				break;
			}
			if (!(i == end || cancel))
				++i;
		}
		if (STATE == CAPTURE) {
			throw std::invalid_argument{"Premature end"};
		}
		matrix M{m, n};
		for (uint i = 0; i < m; i++)
		for (uint j = 0; j < n; j++) {
			M.set(i, j, mat[j + (i * n)]);
		}
		return M;
	}
}

#endif // SIMPLE_MATRIX_MATRIX_H_
