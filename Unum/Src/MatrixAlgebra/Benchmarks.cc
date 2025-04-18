
#include <vector>
#include <stdio.h>
#include <iostream>
#include <CImg.h>
#include <fenv.h>
#include <dirent.h>
#include <string.h>
#include <boost/filesystem.hpp>
#include "FFT.hh"
#include "Quire.hh"

using namespace std;
using namespace cimg_library;
using half_float::half;

//Plot residual after each iteration of CG. Compare different numerical representations.
void CGTest(Matrix<mpf_class> M, string matrixname="unknown matrix", \
    string identifier="", bool quire=0, bool general=0, bool plot=0, bool traffic=0, \
    bool plotTraffic=0, bool clean=0, double relativeTolerance=1e-5, bool scale=0) 
{
    if (traffic && plotTraffic) { fprintf(stderr, "Invalid Arguments."); return; }

    string infoFilename = "plots/CG" + identifier;
    string plotfile, trafficfile, trafficPlotFile;       
    if (plot)        plotfile        = "plots/" + (matrixname + identifier) + ".csv"; 
    if (traffic)     trafficfile     = "plots/" + (matrixname + identifier) + "Traffic.csv";
    if (plotTraffic) trafficPlotFile = "plots/" + (matrixname + identifier) + "TrafficPlot.csv";
    ofstream infoFile;
    
    puts("\nRunning CG benchmark...");
    if (!(M.isSymmetric())) {fprintf(stderr, "Please input symmetric matrix for CG test."); return;}
    
    int n = M.nCols();

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);

    if (scale) {
        M.scaleNorm(M, bM);
        //M=M.rescaleSymmetricOne(3);
        //M.diagScaleAvg(3, M, bM);
        //D.scaleNorm(D, bD);
        //F.scaleNorm(F, bF);
        //P.scaleNorm(P, bP);
    }

    Matrix<double    > D;
    Matrix<float     > F;
    Matrix<Posit32gmp> P;

    D.set(M);
    F.set(M);
    P.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);
    vector<Posit32gmp> xP(n);

    vector<double    > bD(n);
    vector<float     > bF(n);
    vector<Posit32gmp> bP(n);

    downcast(bD, bM);
    downcast(bF, bM);
    downcast(bP, bM);

    double tolerance = M.vectorNorm(bM).get_d()*relativeTolerance;
    cout << "Aiming for relative error of " << relativeTolerance << endl;

    Posit32::clearCounter();
    int d, f, p, b;
    //d = D.conjugateGradientSolver(tolerance, D,  bD, xD, plotfile, "", clean); //leaving double for reference. 
    if (quire) {
        //f = conjugateGradientSolverQ(tolerance, F,  bF, xF, plotfile, "", clean);
        p = conjugateGradientSolverQ(tolerance, P,  bP, xP, plotfile, trafficPlotFile, clean);
    } else {
        f = F.conjugateGradientSolver(tolerance, F, bF, xF, plotfile, "", clean);
        p = P.conjugateGradientSolver(tolerance, P, bP, xP, plotfile, trafficPlotFile, clean);
    }

    if (traffic) {
        ofstream file(trafficfile, ofstream::trunc);
        Posit32::writeAdvantage(file);
        file.close();
    }

    puts("\nNumber of iterations.");
    cout << "Scale?: " << scale << endl;
    cout << "ES: " << ES << endl;
    cout << "Quire?: " << quire << endl;
    cout << "Max/min entry: " << M.getMax() << " : " << M.getMin() << "Range: " << M.getMax()/M.getMin() << endl;
    cout << '\t' << "double iterations = "  << d << endl;
    cout << '\t' << "float iterations = "   << f << endl;
    cout << '\t' << "posit iterations = "   << p << endl;
    cout << '\t' << "average advantage = "  << Posit32::distillAdvantage() << endl;

    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << matrixname << endl;
        infoFile << "relative error tolerance: " << relativeTolerance << endl;
        infoFile << "Scale?: " << scale << endl;
        infoFile << "ES: " << ES << endl;
        infoFile << "Quire?: " << quire << endl;
        infoFile << "Max/min entry: " << M.getMax() << " : " << M.getMin() << "Range: " << M.getMax()/M.getMin() << endl;
        infoFile << '\t' << "double iterations = " << d << endl;
        infoFile << '\t' << "float iterations = "  << f << endl;
        infoFile << '\t' << "posit iterations = "  << p << endl;
        infoFile << '\t' << "average advantage = "  << Posit32::distillAdvantage() << endl;
        infoFile.close();
    }

    //Computed Solutions.
    vector<mpf_class> bmD(n), bmF(n), bmP(n), bmB(n);
    if (quire) {
        upcast(bmD,matVec(D, xD)); 
        vector<float     > fv(n);
        vector<Posit32gmp> fp(n);
        matVecQ(fv, F, xF);
        matVecQ(fp, P, xP);

        upcast(bmF,fv);
        upcast(bmP,fp);
    } else {
        upcast(bmD,matVec(D, xD)); 
        upcast(bmF,matVec(F, xF));
        upcast(bmP,matVec(P, xP));
    }

    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));
    mpf_class rP = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmP));
		
    puts("\nFinal Residuals:");
    cout << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    cout << '\t' << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;

    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << endl;
        infoFile << '\t' << "double relative residual = " << rD/M.vectorNorm(bM) << endl;
        infoFile << '\t' << "float relative residual = "  << rF/M.vectorNorm(bM) << endl;
        infoFile << '\t' << "posit relative residual = "  << rP/M.vectorNorm(bM) << endl;
        infoFile << endl;
        infoFile.close();
    }
    Posit32::writeAdvantage(cout);
    cout << Posit32::distillAdvantage() << endl;
}

//Run trisolve and report final residual. 
void trisolveTest(Matrix<mpf_class> M, string matrixname="unknown matrix", string identifier="", \
    bool cholesky=0, bool general=0, bool traffic=0, bool scale=0) {
    if (!(M.isSquare())) {fprintf(stderr, "Please input square matrix for Tri-Solve test."); return;}

    puts("\nRunning direct solve benchmark...");
    
    int n = M.nCols();

    string infoFilename, trafficFilename;
    if (general) infoFilename= "plots/TriSolve" + identifier;
    if (traffic) trafficFilename = "plots/" + matrixname + identifier + "directTraffic.csv";
    ofstream infoFile; 

    Matrix<double    > D;
    Matrix<float     > F;
    Matrix<Posit32gmp> P;

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);
    
    D.set(M);
    F.set(M);
    P.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);
    vector<Posit32gmp> xP(n);

    vector<double> bD(n);
    vector<float> bF(n);
    vector<Posit32gmp> bP(n);

    downcast(bD, bM);
    downcast(bF, bM);
    downcast(bP, bM);

    if (scale) {
        M.diagScaleAvg(3, M, bM);
        D.diagScaleAvg(3, D, bD);
        F.diagScaleAvg(3, F, bF);
        P.diagScaleAvg(3, P, bP);
    }

    Posit32::clearCounter();

    if (cholesky) {
        D.symmetricTriSolve(xD, bD);
        F.symmetricTriSolve(xF, bF);
        P.symmetricTriSolve(xP, bP);
    } else {
        D.triSolve(xD, bD);
        F.triSolve(xF, bF);
        P.triSolve(xP, bP);
    }

    if (traffic) {
        ofstream file(trafficFilename, ofstream::trunc);
        Posit32::writeAdvantage(file);
        file.close();
    }

    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << matrixname << endl;
        infoFile << "Scale?: " << scale << endl;
        infoFile << "ES: " << ES << endl;
        infoFile << "Cholesky?: " << cholesky << endl;
        infoFile << "Average posit advantage: " << Posit32::distillAdvantage() << endl;
        infoFile.close();
    }

    vector<mpf_class> bmD(n), bmF(n), bmP(n), bmB(n);
    upcast(bmD,matVec(D, xD));
    upcast(bmF,matVec(F, xF));
    upcast(bmP,matVec(P, xP));

    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));
    mpf_class rP = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmP));
		
    puts("\nFinal Residuals:");
    cout << "Scale?: " << scale << endl;
    cout << "ES: " << ES << endl;
    cout << "Cholesky?: " << cholesky << endl;
    cout << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    cout << '\t' << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;
    cout  << log10(rF.get_d()/rP.get_d()) << endl;

    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << '\t' << "double relative residual = " << rD/M.vectorNorm(bM) << endl;
        infoFile << '\t' << "float relative residual = "  << rF/M.vectorNorm(bM) << endl;
        infoFile << '\t' << "posit relative residual = "  << rP/M.vectorNorm(bM) << endl;
        infoFile << log10(rF.get_d()/rP.get_d()) << endl;
        infoFile.close();
    }

    cout << Posit32::distillAdvantage() << endl;
}

/*
void convolveTest() {
    puts("\nRunning convolution benchmark. ");

    int n = 1<<10;

    vector<complex<mpf_class>> inputM(n);
    getVectorComplex(inputM, 1);
    
    vector<complex<float     > > inputF(n);
    vector<complex<Posit32gmp> > inputP(n);
    vector<complex<half      > > inputH(n);
    vector<complex<bfloat16  > > inputB(n);
    
    downcastComplex(inputF, inputM);
    downcastComplex(inputP, inputM);
    downcastComplex(inputH, inputM);
    downcastComplex(inputB, inputM);

    vector<complex<float     > > cF(n);
    vector<complex<Posit32gmp> > cP(n);
    vector<complex<half      > > cH(n);
    vector<complex<bfloat16  > > cB(n);
    vector<complex<mpf_class > > cM(n);

    Posit32::clearCounter();

    convolveFFT(cM, inputM, inputM);
    convolveFFT(cF, inputF, inputF);
    convolveFFT(cP, inputP, inputP);
    convolveFFT(cH, inputH, inputH);
    convolveFFT(cB, inputB, inputB);

    vector<complex<mpf_class> > cfM(n);
    vector<complex<mpf_class> > cpM(n);
    vector<complex<mpf_class> > chM(n);
    vector<complex<mpf_class> > cbM(n);

    upcastComplex(cfM, cF);
    upcastComplex(cpM, cP);
    upcastComplex(chM, cH);
    upcastComplex(cbM, cB);
     
    puts("Errors (mse)...\n");

    cout << "float: "    << mse(cfM, cM) << endl;
    cout << "Posit: "    << mse(cpM, cM) << endl;
    cout << "half: "     << mse(chM, cM) << endl;
    cout << "bfloat16: " << mse(cbM, cM) << endl;

    Posit32::writeAdvantage(cout);
}
*/

/*
void fftImage(string imgFilename, string infofile, bool general=0) {
    CImg<float> img(imgFilename.c_str());
    
    int n = closestPowerTwo(max(img.width(), img.height()));
    cout << n << endl;

    Matrix<complex<mpf_class > > imgM(n, n);

    int startx = (n-img.width())/2;
    int stopx  = img.width() + startx;
    int starty = (n-img.height())/2;
    int stopy  = img.height() + starty;

    double red, green, blue;
    for (int i=starty;i<stopy;i++) {
        for (int j=startx;j<stopx;j++) {
            red   = .3*(*img.data(j-startx, i-starty, 0, 0));
            green = .59*(*img.data(j-startx, i-starty, 0, 1));
            blue  = .11*(*img.data(j-startx, i-starty, 0, 2));
            imgM.m[i][j] = complex<mpf_class>(red+green+blue, 0);
        }
    }

    Matrix<complex<bfloat16  > > imgB(n,n);
    Matrix<complex<half      > > imgH(n,n);
    Matrix<complex<Posit32gmp> > imgP(n,n);
    
    for (int i=0;i<n;i++) downcastComplex(imgB.m[i], imgM.m[i]);
    for (int i=0;i<n;i++) downcastComplex(imgH.m[i], imgM.m[i]);
    for (int i=0;i<n;i++) downcastComplex(imgP.m[i], imgM.m[i]);
    
    fft2D(imgB, imgB);
    ifft2D(imgB, imgB);

    //fft2D(imgH, imgH);
    //fft2D(imgH, imgH);
    
    fft2D(imgP, imgP);
    ifft2D(imgP, imgP);

    Matrix<complex<mpf_class> > imgBm(n, n), imgHm(n, n), imgPm(n, n);
    for (int i=0; i<n; i++) upcastComplex(imgBm.m[i], imgB.m[i]);
    for (int i=0; i<n; i++) upcastComplex(imgHm.m[i], imgH.m[i]);
    for (int i=0; i<n; i++) upcastComplex(imgPm.m[i], imgP.m[i]);

    
    //vector<pair<string, Matrix<complex<mpf_class> > > > files(1);
    //files[0] = pair<string, Matrix<complex<mpf_class> > >(string("testhalf") + ".ppm", imgBm);
    //files[0] = pair<string, Matrix<complex<mpf_class> > >(string("testhalf") + ".ppm", imgHm);
    //files[2] = pair<string, Matrix<complex<mpf_class> > >(output + "Posit"  + ".ppm", imgPm);

    /*
    ofstream imgfile;
    for (auto p : files) {
        imgfile.open(p.first.c_str(), ofstream::trunc);
        imgfile << "P3" << endl;
        imgfile << n << " " << n << endl;
        imgfile << "255" << endl;
        for (int i=0;i<n;i++) {
            for (int j=0;j<n;j++) {
                mpf_class value = p.second.m[i][j].real();
                imgfile << ((int) value.get_d()) << " " << ((int) value.get_d()) << \
                    " " << ((int) value.get_d()) << endl;
            }
        }
        imgfile.close();
    }
    */
    /*
    if (general) {
        ofstream info(infofile, ofstream::app);
        info << imgFilename << endl;
        info << '\t' << "error bfloat: " << mse(imgBm.m, imgM.m) << endl;
        info << '\t' << "error half: " << mse(imgHm.m, imgM.m) << endl;
        info << '\t' << "error posit: " << mse(imgPm.m, imgM.m) << endl;
        info << '\t' << "average posit advantage" << Posit32::distillAdvantage() << endl;
    }

    cout << "error bfloat: " << mse(imgBm.m, imgM.m) << endl;
    cout << "error half: " << mse(imgHm.m, imgM.m) << endl;
    cout << "error posit: " << mse(imgPm.m, imgM.m) << endl;
    cout << '\t' << "average posit advantage" << Posit32::distillAdvantage() << endl;

    //for (auto p : files) system((string("xdg-open ") + p.first).c_str());
}
*/

/*
void correlateImage(string imgFilename, string kernelFilename, string output) {
    CImg<unsigned char> img   (imgFilename.c_str());
    CImg<unsigned char> kernel(kernelFilename.c_str());
    
    int n = closestPowerTwo(max(img.width(), img.height()));
    Matrix<complex<half> > imgM(n, n);
    Matrix<complex<half> > kernelM(n, n);
    imgM.setZero();
    kernelM.setZero();

    int startx = (n-img.width()) / 2;
    int stopx  = img.width() + startx;
    int starty = (n-img.height()) / 2;
    int stopy  = img.height() + starty;

    for (int i=starty;i<stopy;i++) {
        for (int j=startx;j<stopx;j++) {
            imgM.m[i][j]    = complex<half>(*img.data(j-startx, i-starty),    0);
            kernelM.m[i][j] = complex<half>(*kernel.data(j-startx, i-starty), 0); 
        }
    }
    
    Matrix<complex<half> > correlation(n,n);
    correlation.setZero();
    correlate2DFFT(correlation, imgM, kernelM);
    
    ofstream imgfile(output.c_str(), ofstream::trunc);
    imgfile << "P3" << endl;
    imgfile << n << " " << n << endl;
    imgfile << "255" << endl;

    double max=0;
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            complex<half> entry = correlation.m[i][j];
            double value = abs(complex<double>(float(entry.real()), float(entry.imag())));
            if (value > max) max = value;
        }
    }
    
    for (int i=0;i<n;i++) {
        for (int j=0;j<n;j++) {
            complex<half> entry = correlation.m[i][j];
            double value = abs(complex<double>(float(entry.real()), float(entry.imag())));
            imgfile << ((int) value) << " " << ((int) value) << " " << ((int) value) << endl;
        }
    }

    imgfile.close();
    system((string("xdg-open ") + output).c_str());
}
*/


//type 0:=half, 1:=bfloat, 2:=posit
void mixedSolve(Matrix<mpf_class> M, int type, string matrixname="unknown matrix", string identifier="", \
    bool general=0, bool traffic=0) {
    double machineEpsilon=1e-16;
    double tolerance;

    int n = M.nCols();

    string infoFilename(""), trafficfile("");
    if (general) infoFilename= "plots/" + identifier + "squeezeMixed.txt";
    if (traffic) trafficfile = "plots/" + matrixname + identifier + "mixedTraffic.csv";
    ofstream infoFile; 

    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << matrixname << endl;
    }

    Matrix<double> D;
    D.set(M);
    //D = D.rescaleSymmetric(type);
    
    vector<double> x = vector<double>(n, 1/sqrt(n));
    vector<double> b = matVec(D, x);

    tolerance = n*machineEpsilon;
    
    int iterations;
    if (type == 0) {
        try {
            if (general) infoFile << '\t' << "largest/smallest entry: " << \
                D.getMax() << ", " << D.getMin() << ". Range: " << D.getMax()/D.getMin() << endl;
            
            pair<int, double> p = D.symmetricTriSolveMixed<bfloat16>(x, b, tolerance, infoFilename);
            double r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));

            if (general) {
                infoFile << '\t' << "bfloat16 relative residual = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << p.first << " iterations" << endl;
            }

            cout << '\t' << "R*R backward error: " << p.second << endl;
            cout << '\t' << "in " << p.first << " iterations" << endl;  
            cout << '\t' << "bfloat16 relative residual = " << r/D.vectorNorm(b) << endl;
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running bfloat16" << endl;
            infoFile.close();
        }
    }

    if (type == 1) {
        try {
            if (general) infoFile << '\t' << "largest/smallest entry: " << \
                D.getMax() << ", " << D.getMin() << ". Range: " << D.getMax()/D.getMin() << endl;

            pair<int, double> p = D.symmetricTriSolveMixed<half>(x, b, tolerance, infoFilename);
            
            double r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));

            if (general) {
                infoFile << '\t' << "half relative residual = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << p.first << " iterations" << endl;
            }

            cout << '\t' << "R*R backward error: " << p.second << endl;
            cout << '\t' << "half relative residual = " << r/D.vectorNorm(b) << endl;
            cout << '\t' << "in " << p.first << " iterations" << endl; 
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running half precision float" << endl;
            infoFile.close();
        }
    }

    if (type == 2) {
        try {
            if (general) infoFile << '\t' << "largest/smallest entry: " << \
                D.getMax() << ", " << D.getMin() << ". Range: " << D.getMax()/D.getMin() << endl;

            pair<int, double> p = D.symmetricTriSolveMixed<Posit32gmp>(x, b, tolerance, infoFilename);
            double r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));
            
            if (traffic) {
                ofstream file(trafficfile, ofstream::trunc);
                Posit32::writeAdvantage(file);
                file.close();
            }
            
            if (general) {
                infoFile << '\t' << "posit16 relative error = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << p.first << " iterations" << endl;
            }

            cout << '\t' << "R*R backward error: " << p.second << endl;
            cout << '\t' << "posit16 relative error = " << r/D.vectorNorm(b) << endl;
            cout << '\t' << "in " << p.first << " iterations" << endl;
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running Posit." << endl;
            infoFile.close();
        }
    }

    infoFile.close();
}

//type 0:=half, 1:=bfloat, 2:=posit
void triplePrecisionSolve(Matrix<mpf_class> M, int type, string matrixname="unknown matrix", string identifier="", \
    bool general=0, bool traffic=0) {
    double machineEpsilon=5.96e-8;
    double tolerance;

    int n = M.nCols();

    string infoFilename= "plots/" + identifier + "squeezeTriplePrecision.txt";
    string trafficfile = "plots/" + matrixname + identifier + "squeezeTriplePrecisionTraffic.csv";
    ofstream infoFile; 

    
    if (general) {
        infoFile.open(infoFilename, ofstream::app);
        infoFile << matrixname << endl;
    }

    Matrix<float> D;
    D.set(M);
    D = D.rescaleSymmetric(type);
    
    tolerance = n*machineEpsilon;
    
    vector<float> x = vector<float>(n, 1/sqrt(n));
    vector<float> b = matVec(D, x);
    
    int iterations;
    if (type == 0) {
        try {
            iterations = D.symTriSolveThreePrecision<bfloat16, double>(x, b, tolerance);
            printVector(x);
            
            float r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));

            cout << '\t' << "bfloat16 relative error = " << r/D.vectorNorm(b) << endl;
            cout << '\t' << "in " << iterations << " iterations" << endl;
            if (general) {
                infoFile << '\t' << "bfloat16 relative error = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << iterations << " iterations" << endl;
            }
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running bfloat16" << endl;
            infoFile.close();
        }
    }

    if (type == 1) {
        try {
            iterations = D.symTriSolveThreePrecision<half, double>(x, b, tolerance);
            printVector(x);

            float r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));
            cout << '\t' << "half relative error = " << r/D.vectorNorm(b) << endl;
            cout << '\t' << "in " << iterations << " iterations" << endl;
            if (general) {
                infoFile << '\t' << "half relative error = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << iterations << " iterations" << endl;
            }
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running half precision float" << endl;
            infoFile.close();
        }
    }

    if (type == 2) {
        try {
            Posit32::clearCounter();
            iterations = D.symTriSolveThreePrecision<Posit32gmp, double>(x, b, tolerance);
            printVector(x);
            
            float r = D.vectorNorm(D.vectorCombination(1.0, b, -1.0, matVec(D, x)));

            if (traffic) {
                ofstream file(trafficfile, ofstream::trunc);
                Posit32::writeAdvantage(file);
                file.close();
            }
            cout << '\t' << "posit16 relative error = " << r/D.vectorNorm(b) << endl;
            cout << '\t' << "in " << iterations << " iterations" << endl;
            
            if (general) {
                infoFile << '\t' << "posit16 relative error = " << r/D.vectorNorm(b) << endl;
                infoFile << '\t' << "in " << iterations << " iterations" << endl;
                infoFile << '\t' << "average posit advantage: " << \
                    Posit32::distillAdvantage() << endl << endl;
            }
        } catch (...) {
            if (general) infoFile << '\t' << "Exception while running Posit." << endl;
            infoFile.close();
        }
    }

    infoFile.close();
}

void loadAll(string filename) {
    ofstream file(filename, ofstream::app);
    string path = string("allMatrices");
    Matrix<Posit32gmp> A;
    map<int, double> m;
    for (boost::filesystem::directory_entry& entry : boost::filesystem::directory_iterator(path.c_str())){
        cout << entry.path() << endl;
        if (A.recordMatrix(entry.path().string().c_str())) {
            map<int, double> distr = Posit32::distribution();
            for (pair<int, double> item : distr) {
                m[item.first] = m.count(item.first) ? m[item.first]+item.second : item.second;
            }
            Posit32::clearCounter();
        }
    } 

    double total=0;
    for (auto item : m) total += item.second;
    for (auto item : m) m[item.first] /= total;
    for (auto item : m) file << item.first << ",";
    file << endl;
    for (auto item : m) file << item.second << ",";
}

/*
Enter command line arguments as: 
    test ID:
        0->CG (matrixname, other arguments...)
        1->Tri-Solve full precision (matrixname, other arguments...)
        2->Tri-Solve mixed precision (matrixname, tolerance, type, other arguments...)
        3->fft image (image filename, file to write information, other arguments...)
*/
int main(int argc, char*argv[]) {

	if (argc < 2) { fprintf(stderr, "Please enter test id as first argument"); return 0; }
    int test = stoi(argv[1]);
    
    Matrix<mpf_class> systemM;
    string matrixname;
    if (test==0 || test==1 || test==2 || test==4) {
        for (int i = 12; i<strlen(argv[2])-4; i++) matrixname.push_back(argv[2][i]);
        systemM.loadMPF(argv[2]);
    }

    vector<string> args(argc);
    for (int i=0;i<argc;i++) args[i]=argv[i]; 
    string identifier="";
    mpf_class tolerance;
    bool quire, general, plot, traffic, plotTraffic, clean, cholesky, scale;
    int type;
    
    string imgFile; 
    string infofile;

    Posit32::initializeCounter();
    switch(test) 
    {
        case 0:
            tolerance   = args[3];
            quire       = find(begin(args), end(args), "quire") != end(args);
            general     = find(begin(args), end(args), "general") != end(args);
            plot        = find(begin(args), end(args), "plot") != end(args);
            traffic     = find(begin(args), end(args), "traffic") != end(args);
            plotTraffic = find(begin(args), end(args), "plotTraffic") != end(args);
            clean       = find(begin(args), end(args), "clean") != end(args);
            scale       = find(begin(args), end(args), "scale") != end(args);
            
            if (quire) identifier += "quire";
            if (clean) identifier += "Clean";
            //for (vector<mpf_class> row : systemM.m) {
            //    int max = 0;
            //    for (int i=0;i<row.size();i++)
            //        if (row[i]>row[max]) max =i;
            //    cout << max << endl;
            //}
            CGTest(systemM, matrixname, identifier, quire, general, \
                plot, traffic, plotTraffic, clean, tolerance.get_d(), scale);
            break;
        case 1:
            identifier = "After";
            general  = find(begin(args), end(args), "general") != end(args);
            traffic  = find(begin(args), end(args), "traffic") != end(args);
            cholesky = find(begin(args), end(args), "cholesky") != end(args);
            scale    = find(begin(args), end(args), "scale") != end(args); 
            trisolveTest(systemM, matrixname, identifier, cholesky, general, traffic, scale);
            break;
        case 2:
            type = stoi(args[3]);
            general = find(begin(args), end(args), "general") != end(args);
            traffic = find(begin(args), end(args), "traffic") != end(args);
            mixedSolve(systemM, type, matrixname, identifier, general, traffic);
            break;
        case 3:
            imgFile = args[2];
            infofile = args[3];
            general = find(begin(args), end(args), "general") != end(args);
           // fftImage(imgFile, infofile, general);
            break;
        case 4:
            type = stoi(args[3]);
            general = find(begin(args), end(args), "general") != end(args);
            traffic = find(begin(args), end(args), "traffic") != end(args);
            triplePrecisionSolve(systemM, type, matrixname, identifier, general, traffic);
            break;
        case 5:
            loadAll("precisionDistribution.txt");
    }
    return 1;
}
