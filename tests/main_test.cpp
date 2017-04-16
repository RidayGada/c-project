#include <iostream>

#include "matrix.hpp"

using namespace std;
using namespace simple_matrix;

int test_determinant();
int test_minor();
int test_cofactor();
int test_inverse();
int test_parse();

int main() {
	cout << "Testing determinants..." << flush;;
	if (!test_determinant()) {
		cout << " FAILED" << endl;
		return -1;
	}
	cout << " Passed!" << endl;

	cout << "Testing minors..." << flush;;
	if (!test_minor()) {
		cout << " FAILED" << endl;
		return -1;
	}
	cout << " Passed!" << endl;

	cout << "Testing cofactors..." << flush;;
	if (!test_cofactor()) {
		cout << " FAILED" << endl;
		return -1;
	}
	cout << " Passed!" << endl;

	cout << "Testing inverses..." << flush;;
	if (!test_inverse()) {
		cout << " FAILED" << endl;
		return -1;
	}
	cout << " Passed!" << endl;

	cout << "Testing parser..." << flush;;
	if (!test_parse()) {
		cout << " FAILED" << endl;
		return -1;
	}
	cout << " Passed!" << endl;
}
