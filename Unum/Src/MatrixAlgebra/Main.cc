#include "Matrix.hh"
#include "mmio.hh"
#include <iostream>
#include <boost/algorithm/string.hpp> 
#include <cmath>
#include <sstream>
using namespace std;

void toMatrix(string filename, Matrix<bool> &A) {
    ifstream matrixfile(filename);
    string line;
    int s=0;
    vector<string> result(5);

    while(matrixfile.is_open()) {
        getline(matrixfile, line);
        if (matrixfile.eof()) {
            matrixfile.close();
        } else {
            boost::split(result, line, boost::is_any_of(","));
            s=max(s, max(stoi(result[0]), stoi(result[1])));
        }
    }
    s++;

    A.setSize(s, s);

    while(matrixfile.is_open()) {
        getline(matrixfile, line);
        if (matrixfile.eof()) {
            matrixfile.close();
        } else {
            boost::split(result, line, boost::is_any_of(","));
            A.m[stoi(result[0])][stoi(result[1])] = 1;
        }
    }

    for (int i=0;i<s;i++) A.m[i][i]=1;

}

void convert(const Matrix<uint16_t> &in, Matrix<bool> &out) {
    out.setSize(in.nRows(), in.nCols());
    for (int i=0; i<in.nRows(); i++) {
        for (int j=0; j<in.nCols(); j++) {
            if (in.m[i][j]!=0)
                out.m[i][(int) (in.m[i][j]-1)]=1;
        }
    }
}

//input matrix filename
int main(int argc, char*argv[]) {
    Matrix<bool> L;
    //Matrix<double> L, U;
    
    toMatrix("A30_size1.csv", L);
    //L.printMatrix();
    //A.LUdecomposition(L, U);
    //A.printMatrix();
    Matrix<uint16_t> adjList;
    Matrix<bool> nonzeros;
    L.findIndependentOps<1>(adjList);
    
    //L.printMatrix();
    convert(adjList, nonzeros);
    

    for (int i=0;i<L.nRows();i++) {
        for (int j=0;j<L.nCols();j++) {
            if ((L.m[i][j]!=0 && nonzeros.m[i][j]==0 || L.m[i][j]==0 && nonzeros.m[i][j]!=0)) cout << i << ", " << j << endl;
        }
    }

    //Matrix<bool> A(25356, 25356);

    //toMatrix("A30_size1.csv", A);


    //A.findIndependentOps<1>();
    
    return 0;
}