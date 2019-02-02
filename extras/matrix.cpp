#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;

int main(){
	int n, u;
	cout << "Give me a number of columns" << endl;
	cin >> n;
	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n,n);
	cout << "Give me a MU" << endl;
	double mu;
    cin >> mu;
	
	for (int i = 0; i < n; i++){
		m(i,i) = -2 - mu;
		m(0,n-1) = 1;
		m(n-1,0) = 1;	
		if (i>0) m(i,i-1) = 1;
		if (i<n-1) m(i,i+1) = 1;
	}

	subMu(&m, mu);
	Eigen::SelfAdjointEigenSolver<MatrixXd> eigenSolver(m,false);
	
	cout << "	Matrix:" << '\n';
	cout << m << '\n';
	cout << "The eigenvalues of the"<< n <<"x"<< n << "matrix of ones are:"<<endl;

	double* eigenValues = (double *)eigenSolver.eigenvalues().data();
	
	sort(eigenValues, eigenValues + n);
	ofstream myfile;
	myfile.open("eigenvalues.txt")
    for (size_t i = 0; i != n; ++i) {
        cout << eigenValues[i] << endl;
        myfile eigenValues[i] << endl;
	}
	myfile.close();
}
