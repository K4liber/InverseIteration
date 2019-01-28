#include <iostream>
using namespace std;

class LinearSolver {

    public:

    LinearSolver(double** A, double *b, int N) {
        this->A = A;
        this->b = b;
        this->N = N;
    }

    double** reduceSystem(double E, double* b, int size) {
        // Step 0 - create matrix filled with zeros
        double** A = new double*[size];

        for (int i = 0; i < size; ++i) {
            A[i] = new double[size];

            for (int j = 0; j < size; ++j)
                A[i][j] = 0.0;
        }


        // 1 Step - create almost triangular matrix
        double divider = 2.0;
        A[0][0] = E/divider;
        A[0][1] = 1.0/divider;
        A[0][size-1] = 1.0/divider;
        b[0] = b[0]/divider;

        for (int i = 1; i < size-1; ++i) {
            A[i][i] = A[i-1][i] - A[i-1][i-1]*A[i-1][i-1];
            A[i][i+1] = -1.0*A[i-1][i-1];
            A[i][size-1] = 1.0;
            b[i] = b[i-1]-1.0*b[i]*A[i-1][i-1];
            if (i%6 == 0) {
                A[i][i] = A[i][i]/divider;
                A[i][i+1] = A[i][i+1]/divider;
                A[i][size-1] = A[i][size-1]/divider;
                b[i] = b[i]/divider;
            }
        }

        A[size-1][0] = 1.0-1.0*A[size-2][size-2];
        A[size-1][size-1] = A[size-2][size-1]-E*A[size-2][size-2];
        b[size-1] = b[size-2]-b[size-1]*A[size-2][size-2];

        printVector(b, size);

        // Step 2
        return A;
    }

    private:

    void printVector(double *vector, int N) {
        cout<<"(";
        for (int i = 0; i < N-1; i++) {
            cout<<vector[i]<<", ";
        }
        cout<<vector[N-1]<<")"<<endl;
    }

    double** A;
    double* b;
    int N;

};

void printMatrix(double **matrix, int N) {
    for (int i = 0; i < N ; i++) {
        for (int j = 0; j < N-1; j++) {
            cout<<matrix[i][j]<<", ";
        }
        cout<<matrix[i][N-1]<<endl;
    }
}

void printVector(double *vector, int N) {
    cout<<"(";
    for (int i = 0; i < N-1; i++) {
        cout<<vector[i]<<", ";
    }
    cout<<vector[N-1]<<")"<<endl;
}

/* 
	*** Compile: g++ -std=c++11 linearSolver.cpp -o linearSolver.x *** 
  	*** Run: ./linearSolver.x
*/
int main()  
{ 
    int N = 5;
    double** A = new double*[N];
    double* b = new double[N];

    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
        b[i] = 1;

        for (int j = 0; j < N; ++j)
            A[i][j] = i+j+1;
    }

    LinearSolver *ls = new LinearSolver(A, b, N);
    double** reduced = ls->reduceSystem(2, b, N);
    printMatrix(reduced, N);

    return 0; 
}