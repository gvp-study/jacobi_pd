#include "jacobi_pd.hpp"
#include <iostream>
using namespace jacobi_pd;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{

// ...
    int n = 3;       // Matrix size
    double **M;      // A symmetric n x n matrix you want to diagonalize
    double *evals;   // Store the eigenvalues here.
    double **evecs;  // Store the eigenvectors here.
// Allocate space for M, evals, and evecs (omitted)...
    M = new double*[n];
    evecs = new double*[n];
    evals = new double[n];
    for(int i = 0; i < n; i++)
    {
	M[i] = new double[n];
	evecs[i] = new double[n];
    }
    M[0][0] = 2.0; M[0][1] = 1.0; M[0][2] = 1.0;
    M[1][0] = 1.0; M[1][1] = 2.0; M[1][2] =-1.0;  //Note: The matrix
    M[2][0] = 1.0; M[2][1] =-1.0; M[2][2] = 2.0;  //must be symmetric.

// Now create an instance of Jacobi ("eigen_calc").  This will allocate space
// for storing intermediate calculations.  Once created, it can be reused
// multiple times without incurring the cost of allocating memory on the heap.

    Jacobi<double, double*, double**> eigen_calc(n);

// Note:
// If the matrix you plan to diagonalize (M) is read-only, use this instead:
//   Jacobi<double, double*, double**, double const*const*> eigen_calc(n);
// If you prefer using C++ vectors over C-style pointers, this works also:
//   Jacobi<double, vector<double>&, vector<vector<double>>&,
//          const vector<vector<Scalar>>&>  eigen_calc(n);

// Now, calculate the eigenvalues and eigenvectors of M

    eigen_calc.Diagonalize(M, evals, evecs);  //(successful if return value is != 0)

    std::cout << "eigenvalues:  ";
    for (int i=0; i < n; i++)
	cout << evals[i] << " ";
    cout << endl;
    for (int i=0; i < n; i++) {
	cout << "eigenvector" <<i+1<< ": ";
	for (int j=0; j < n; j++)
	    cout << evecs[i][j] << " ";
	cout << endl;
    }
}
