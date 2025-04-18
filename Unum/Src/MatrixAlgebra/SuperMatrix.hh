#ifndef __SUPERMATRIX_HH_
#define __SUPERMATRIX_HH_
// #include <omp.h>
#include <queue>
#include "Matrix.hh"

template <class T>
class SuperMatrix {
    using vec = vector<T>;
    public:
        vector<Matrix<T>*> val;
        vector<int> col_ptr;
        vector<int> row_idx;
        
        int nnz;
        int n;

        //Convert non-blocked matrix into superMatrix.
        //Each item will be a 1x1 matrix. 
        SuperMatrix(const Matrix<T> &A) {
            nnz=A.nnz();
            n=A.nCols();
            
            val.resize(nnz);
            row_idx.resize(nnz);
            col_ptr.resize(n+1);

            int val_idx=0;
            //Matrix<T> entry(1, 1);
            for (int j=0;j<n;j++) {
                col_ptr[j]=val_idx;
                for (int i=0;i<n;i++) {  
                    if (A.m[i][j] != 0) {
                        val[val_idx]=new Matrix<T>(1,1); //Matrix<T> (1, 1); //Will this copy data?
                        val[val_idx]->m[0][0]=A.m[i][j];
                        row_idx[val_idx]=i;
                        val_idx++;
                    }
                }
            }
            col_ptr[n]=nnz;
        }

        void add(vec &x, T scale, const vec &y) const {
            for (int i=0;i<x.size();i++) x[i]+=scale*y[i];
        }

        void print() const {
            puts("val: ");
            for (Matrix<T>* value: val) { (*value).print(); }
            puts("row_idx: ");
            printVector(row_idx);
            puts("col_ptr: ");
            printVector(col_ptr);
        }

        void solveLowerTriangularSerial(vector<vec> &x, const vector<vec> &y) const {
            vector<vec> left_sum(n);
            vector<vec> b(y);
            for (int i=0;i<n;i++) { left_sum[i].resize(y[i].size()); x[i].resize(y[i].size()); }
            
            for (int i=0;i<n;i++) {
                add(b[i], -1, left_sum[i]);
                val[col_ptr[i]]->solveLowerTriangularSystem(x[i], b[i]);
                for (int j=col_ptr[i]+1;j<col_ptr[i+1];j++) {
                    add(left_sum[row_idx[j]], 1, matVec(*(val[j]), x[i]));
                }
            }
        }

        void solveLowerTriangularRandom(vector<vec> &x, const vector<vec> &y) const {
            vector<int> in_degree(n);
            for (int i=0;i<nnz;i++) in_degree[row_idx[i]]++;

            vector<vec> left_sum(n);
            vector<vec> b(y);
            for (int i=0;i<n;i++) { left_sum[i].resize(y[i].size()); x[i].resize(y[i].size()); }
            
            vector<int> r(n);
            queue<int> q;
            for (int i=0;i<n;i++) r[i]=i;
            random_shuffle(r.begin(), r.end());
            for (int x : r) q.push(x);

            int i;
            while (!q.empty()) {
                cout << q.size() << endl;
                i=q.front();
                q.pop();
                if (in_degree[i]>1) { q.push(i); continue; } 
                add(b[i], -1, left_sum[i]);         
                val[col_ptr[i]]->solveLowerTriangularSystem(x[i], b[i]);
                for (int j=col_ptr[i]+1;j<col_ptr[i+1];j++) { //column broadcast. 
                    add(left_sum[row_idx[j]], 1, matVec(*(val[j]), x[i])); //send message x[i] to process holding this block. 
                    in_degree[row_idx[j]]--;
                }
            }
            
        }
        
    
};

#endif
