#include  <fstream>
#include "helpers.hh"

//Get random mpf value in [0, 1]
mpf_class rand(gmp_randstate_t r) {
    mpf_class m;
    mpf_urandomb(m.get_mpf_t(), r, 128);  
    return m;
}

//Get vector of random values in [0, itemSize]
void getVector(vector<mpf_class> &v, int itemSize) {
    gmp_randstate_t r;
    gmp_randinit_mt(r);
    int n = v.size();
    for (int i = 0; i < n; i++) v[i] = (rand(r)-.5)*2*itemSize;
}

void getVectorComplex(vector<complex<mpf_class> > &v, int itemSize) {
    gmp_randstate_t r;
    gmp_randinit_mt(r);
    int n = v.size();
    for (int i = 0; i < n; i++) v[i] = complex<mpf_class>((rand(r)-.5)*2*itemSize, 0);
}
/*
mpf_class mse(vector<complex<mpf_class> > modified, vector<complex<mpf_class> > original) {
    mpf_class sum(0), pow;
    int n = modified.size();

    if (original.size() != n) { fprintf(stderr, "Invalid dimensions."); return 0; }
    
    for (int i=0;i<n;i++) {
        m_pow(pow, abs(modified[i]-original[i]), 2); 
        sum += pow;
    }

    sum /= n;
    return sum;
}

mpf_class mse(vector<vector<complex<mpf_class> > > modified, vector<vector<complex<mpf_class> > > original) {
    mpf_class sum(0), pow;
    int m = modified.size();
    int n = modified[0].size();

    if (original.size() != m || original[0].size() != n) { fprintf(stderr, "Invalid dimensions."); return 0; }
    
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            m_pow(pow, abs(modified[i][j]-original[i][j]), 2); 
            sum += pow;
        }
    }

    sum /= m*n;
    return sum;
}
 */

mpf_class mse(vector<vector<mpf_class> > modified, vector<vector<mpf_class> > original) {
    mpf_class sum(0), pow;
    int m = modified.size();
    int n = modified[0].size();

    if (original.size() != m || original[0].size() != n) { fprintf(stderr, "Invalid dimensions."); return 0; }
    
    for (int i=0;i<m;i++) {
        for (int j=0;j<n;j++) {
            m_pow(pow, modified[i][j]-original[i][j], 2); 
            sum += pow;
        }
    }

    sum /= m*n;
    return sum;
}


//Set matrix to be diagonal with a given condition number.
//In this case condition number = max(cond, n);
//void setConditionNumber(Matrix<mpf_class> &M, mpf_class cond) {
//    if (!(M.isSquare())) {fprintf(stderr, "Please enter square matrix."); return; }
//    gmp_randstate_t r;
//    gmp_randinit_mt(r);
//
//    int n = M.nCols();
//    M.m[0][0]=1;
//    for (int i=1;i<n-1;i++) M.m[i][i]=rand(r)+i;
//    M.m[n-1][n-1] = cond;
//}

int closestPowerTwo(int n) { return pow(2, round(log2(n))); }
double closestPowerTwo(double d) { return pow(2.0, round(log2(d))); }

double closestPowerFour(double d) { 
    double a = log2(d)/log2(4);
    return pow(4.0, round(a)); 
}

void cast(float &a, half       b) { a = float(b); }
void cast(float &a, bfloat16   b) { a = b.f;      }
void cast(float &a, Posit32gmp b) { a = b.toDouble(); }
void cast(half &a,       float b) { a = b; }
void cast(bfloat16 &a,   float b) { a = b; }
void cast(Posit32gmp &a, float b) { a = b; }

void cast(double &a, half       b) { a = float(b); }
void cast(double &a, bfloat16   b) { a = b.f;      }
void cast(double &a, Posit32gmp b) { a = b.toDouble(); }
void cast(double &a, double     b) { a = b; }

void cast(half &a,       double b) { a = b; }
void cast(bfloat16 &a,   double b) { a = b; }
void cast(Posit32gmp &a, double b) { a = b; }

void cast(double &a,     mpf_class b) { a=b.get_d(); }
void cast(float &a,      mpf_class b) { a=b.get_d(); }
void cast(bfloat16 &a,   mpf_class b) { a=b;         }
void cast(half &a,       mpf_class b) { a=b.get_d(); }
void cast(Posit32gmp &a, mpf_class b) { a.set(b);    }

void cast(mpf_class &a, double     b) { a=b;        }
void cast(mpf_class &a, float      b) { a=b;        }
void cast(mpf_class &a, bfloat16   b) { a=b.f;      }   
void cast(mpf_class &a, half       b) { a=float(b); }
void cast(mpf_class &a, Posit32gmp b) { b.get(a);   }

void downcast(vector<double    > &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i] = in[i].get_d(); }
void downcast(vector<float     > &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i] = in[i].get_d(); }
void downcast(vector<bfloat16  > &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i] = in[i];         }
void downcast(vector<half      > &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i] = in[i].get_d(); }
void downcast(vector<Posit32gmp> &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i].set(in[i]);      }
void downcast(vector<mpf_class > &out, vector<mpf_class> in) { for(int i=0;i<in.size();i++) out[i] = in[i];         }

void upcast(vector<mpf_class> &out, vector<double    > in){ for(int i=0;i<in.size();i++) out[i] = in[i];       }
void upcast(vector<mpf_class> &out, vector<float     > in){ for(int i=0;i<in.size();i++) out[i] = in[i];       }
void upcast(vector<mpf_class> &out, vector<bfloat16  > in){ for(int i=0;i<in.size();i++) out[i] = in[i].f;     }
void upcast(vector<mpf_class> &out, vector<half      > in){ for(int i=0;i<in.size();i++) out[i] = half(in[i]); }
void upcast(vector<mpf_class> &out, vector<Posit32gmp> in){ for(int i=0;i<in.size();i++) in[i].get(out[i]);    }
void upcast(vector<mpf_class> &out, vector<mpf_class > in){ for(int i=0;i<in.size();i++) out[i] = in[i];       }

/* 
Record mtx file into posit traffic counter, taking symmetry, complex entries, and
'structure only' matrices into account. 
*/
/*
bool recordToPosit(const char *filename){
    // from MMIO example
    int ret_code;
    FILE *f;
    int M, N, nz,n;
    MM_typecode matcode;
    
    if ((f = fopen(filename, "r")) == NULL)
        return 0;
    
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return 0;
    }

    if (mm_is_pattern(matcode)) return 0;
    
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        return 0;

    Posit32 p;

    //true if we must record each entry twice.
    bool twice = 0;
    if (mm_is_symmetric(matcode) || mm_is_skew(matcode) || mm_is_hermitian(matcode)) twice=1;
    
    char* real    = (char*) calloc(200, sizeof(char));
    char* complex = (char*) calloc(200, sizeof(char));
    for (n=0; n<nz; n++)
    {
        int i,j;
        if (mm_is_complex(matcode)) {
            fscanf(f, "%d %d %s %s\n", &i, &j, real, complex);
            p = string(real); 
            p = string(complex); 
            if (twice && i != j) {
                p = string(real); 
                p = string(complex); 
            }
        } else {
            fscanf(f, "%d %d %s\n", &i, &j, real);
            p = string(real); 
            if (twice && i != j) {
                p = string(real); 
            }
        }
        memset(real,    0, sizeof(char)*200);
        memset(complex, 0, sizeof(char)*200);
    }
    free(real);
    free(complex);
    
    if (f !=stdin) fclose(f);

    return 1;
}*/

/*
Takes all mtx files from directory and writes numerical advantage into second argument file. 
*/
/*
void writeAllMatrixInfo(string directory, string filename) {
    Posit32::initializeCounter();
    auto iter = directory_iterator(directory);
    for (auto e : iter) {
        cout << e.path().u8string().c_str() << endl;
        recordToPosit(e.path().u8string().c_str());
    }
    ofstream file;
    file.open(filename, ofstream::trunc);
    Posit32::writeAdvantage(file);
    file.close();
}*/
