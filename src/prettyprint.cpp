//Includes the basic header files

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>

#include "matrix.hpp"  //includes the header file matrix.hpp

using namespace simple_matrix;

std::string matrix::pretty() const {  //method that returns a readable form of matrix
	if (is_empty()) {
		return std::string("Empty");  //returns the string 'Empty' if matrix is empty
	}
	std::ostringstream osstr, main;  //declaration of two string streams
	std::vector<std::string> strs;  //declaration of a vector of strings
	uint mwidth = 0;  //an unsigned int initialized to 0
	for (uint j = 0; j < n_; ++j) {
		for (uint i = 0; i < m_; ++i) { //loops through each element
			double term = (*this)(i, j);
			if (EQUAL(term,0))  //checks if the element is approximately equal to 0 using macro EQUAL
				term = 0;
			osstr << term;  //appends element to osstr
			std::string str = osstr.str();
			osstr.str(std::string());
			strs.push_back(str); //stores resultant string in the vector
			uint nl = str.size();
			mwidth = nl > mwidth ? nl : mwidth;  //to store max width of element so far found
		}
	}
	uint midwidth = (mwidth * n_) + ((n_ + 1) << 1); //for width of border line between upper and lower parts of matrix
	main << "\u250c\u2500" << std::setw(midwidth - 2) << "" << "\u2500\u2510";
	for (uint i = 0; i < m_; i++) { //Constructs the middle part of the matrix, row by row, using Unicode characters
		main << std::endl;
		main << "\u2502";
		for (uint j = 0; j < n_; j++) {
			std::string& str = strs[i + (j * m_)];
			int wlen = (str.size() + mwidth + 1) >> 1; //calculates width of each cell
			main << "  ";
			main << std::setw(wlen); 
			main << str;
			main << std::setw(mwidth - wlen) << ""; //calculates max width found so far
		}
		main << "  \u2502";

	}
	main << std::endl << "\u2514\u2500" << std::setw(midwidth - 2) << "" << "\u2500\u2518";
	//Constructs the bottom border line of the matrix using Unicode characters. Finally, 
	//it returns the constructed string representing the pretty-printed matrix.
	return main.str();
}
