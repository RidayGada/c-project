/*
This is based upon the Steinhaus-Johnson-Trotter algorithm
https://en.wikipedia.org/wiki/Steinhaus–Johnson–Trotter_algorithm
This is used to iterate through every permutation (p) of the set S_(n).

The algorithm is slightly modified to use zero indexing instead, so that
{0, 1, 2} \in S_3. This makes more sense as we are using the results of the
permutations as indexes for the matrix.
*/

#include "matrix.hpp"
#define SWAP(a, b) tmp=(a); (a)=(b); (b)=tmp  //a macro defined to swap two numbers using temp variable

namespace simple_matrix {
	//This function basically calculates the position of element i in array p of size n
	inline uint matrix::permutation_position(uint n, uint i, uint *p) { //an inline function inside class matrix
		for (uint j = 0; j < n; j++)  //a for loop to iterate through the array p
			if (p[j] == i)  //find the position of element i
				return j;  //returns the index of array in which that particular element i is stored
		return 0;  //returns 0 if not found
	}

	//initialized array p(permutation) and v(direction) with default values
	void matrix::permutation_init(uint n, uint *p, uint *v) {  //a function defined in class Matrix
		for (uint i = 0; i < n; i++) {
			p[i] = i;  //sets every element of p as its index value
			v[i] = 0;  //sets every element of v to 0
		}
	}

	//defines a function in class Matrix
	int matrix::permutation_permute(uint n, uint *p, uint *v) {
		uint i = n;
		while (i) { //while loop that iterates from n-1 to 0
			--i;
			uint pi = permutation_position(n, i, p);  //calculates position element i.
			// Determine if mobile that can move left or right based on its position and direction
			if (pi == 0 && !v[pi])
				continue;
			if (pi == n - 1 && v[pi])
				continue;
			if (!v[pi] && p[pi - 1] > i)
				continue;
			if (v[pi] && p[pi + 1] > i)
				continue;
			// Is largest mobile int
			uint tmp;
			uint sp = v[pi] ? pi + 1 : pi - 1;
			SWAP(p[pi], p[sp]); //swaps the element with its adjacent element in the direction
			SWAP(v[pi], v[sp]); //it an move and updates direction accordingly

			// Flip directions of larger ints
			//After swapping, this loop flips the directions of elements
			//larger than the current element i in the permutation array p
			for (uint j = i + 1; j < n; j++) {
				pi = permutation_position(n, j, p);
				v[pi] = !v[pi];
			}
			return 1;  //indicates success that valis permutation is found and generated
		}
		return 0;  //indicates the end of permutation sequence
	}
}
