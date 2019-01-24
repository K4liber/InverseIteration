#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>
#include <fstream>

using namespace std;
using Eigen::MatrixXd;
using namespace Eigen;

/* *** Substract mu value from diagonal of the m Matrix *** */
void subMu(Eigen::MatrixXd *m, double mu){
	int n = sqrt(m->size());
	for (int i=0; i<n; i++){
		m->operator()(i,i) = m->operator()(i,i)-mu;
	}
}

/* *** Push radnom floats from 0 to 1 to the vector v elements *** */
void random_vector(Eigen::VectorXd *v){
	for (int i=0; i<v->size(); i++){
		v->operator()(i) = float(rand())/float(RAND_MAX);
	}
}

/* *** Count eigenvalue of matrix A using inverse iteration algorithm *** */
double inverse_iteration (Eigen::MatrixXd A, double epsilon){
	int n = sqrt(A.size());
	VectorXd x0(n);
	VectorXd x(n);
	random_vector(&x0);
	x0 = x0/x0.norm();
	double res = 1.0;
	int i=0;
	while ((res>epsilon) && (i<1000)){
		i++;
		x = A.colPivHouseholderQr().solve(x0);
		x = x/x.norm();
		res = (x-x0).norm();
		cout<<"itaration: "<<i<<" residual: "<<res<<endl;
		x0 = x;
	}
	VectorXd s = A*x;
	return s(0)/x(0);
}

/* *** Save matrix of class Eigen::MatrixXd to mtx format *** */
void mtx(Eigen::MatrixXd tab, int n){
	ofstream file;
	file.open("matrix.mtx");
	int nozeros=0;
	for (int i = 0; i<n; ++i){
		for(int j = 0; j<n; ++j){
			if( tab(i,j) != 0){
				 nozeros++;
				 file << i+1 <<" "<< j+1 <<" "<< tab(i,j) << '\n';
			 }
		}
	}
	file << n << " " << n << " " << nozeros;	
	file.close();	
}

/* 
	*** Compile: g++ -std=c++11 matrix.cpp -o matrix.x *** 
  	*** Run: ./matrix.x
*/
int main(){
	int n, u;
	cout << "Give me a number of columns" << endl;
	cin >> n;
	Eigen::MatrixXd m = Eigen::MatrixXd::Zero(n,n);
	cout << "Give me a MU" << endl;
	double sjfa;
    cin >> sjfa;
	
	for (int i=0; i<n; i++){
		m(i,i)=-2;
		m(0,n-1)=1;
	    m(n-1,0)=1;	
	    if (i>0) m(i,i-1)=1;
	    if (i<n-1) m(i,i+1)=1;
	 }

	subMu(&m,sjfa);
	Eigen::EigenSolver<MatrixXd> es (m,false);
	
	cout << "	Matrix:" << '\n';
	cout << m << '\n';
	cout << "The eigenvalues of the"<< n <<"x"<< n << "matrix of ones are:"<<endl;

	cout<< es.eigenvalues() << endl;
	double* v=(double *)es.eigenvalues().data();
	double ev[n];
	for(int i=0;i<n;i++){
		ev[i]=v[2*i];
	}
	
	sort(ev, ev + n);
	ofstream myfile;
    myfile.open ("eigensolves.txt");
    for (size_t i = 0; i != n; ++i){
        cout << ev[i] << " "<<endl;
        myfile << ev[i] << '\n';
	}
        
	myfile.close();
	mtx(m, n);
	double value = inverse_iteration(m, 0.0001);
	cout<<"Wartosc wlasna: "<<value<<endl;
	
}
