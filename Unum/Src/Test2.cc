
#include <stdio.h>
#include <string.h>
#include <gmpxx.h>
#include <iostream>
#include <fstream>
#include <set>
#include <bitset>
#include "Posit32.hh"
#include "half.hpp"
#include "bfloat16.hh"
#include "gmp_helpers.hh"
#include "SoftPosit/source/include/softposit_cpp.h"
#include "MatrixAlgebra/helpers.hh"

using namespace std;
using half_float::half;

//For verifying that Posit32(int) and Posit32(double) are consistent with other initalization methods. 
int constructorsTest() {
    cout << "Checking that different intialization methods are consistent..." << endl;
    mpf_class a("0"), b("1"), c("-2");
    Posit32 ma, mb, mc;

    ma.set(a); //initialization using mpf. 
    mb.set(b);
    mc.set(c);
    cout << "ma, mb, mc = " << ma << ", " << mb << ", " << mc << endl;
    Posit32 ia(0), ib(1), ic(-2); //constructing from int
    cout << "ia, ib, ic = " << ia << ", " << ib << ", " << ic << endl;
    Posit32 da(0.0), db(1.0), dc(-2.0); //constructing from double
    cout << "da, db, dc = " << da << ", " << db << ", " << dc << endl;
    if (!(ma == ia && ma == da && ia == da)
        || !(mb == ib && mb == db && ib == db)
        || !(mc == ic && mc == dc && ic == dc)) {
        cout << "Failed" << endl << endl; 
        return 0;
    } 
    else cout << "Success" << endl << endl;
    return 1;
}

int operatorsTest1() {
    cout << "Testing +=, -=... " << endl;

    bool success = 1;
    Posit32gmp a(0), b(1), c(-2);
    a += b;
    if (!(a == b)) success = 0; 
    cout << a << endl;
    a -= b;
    if (!(a.d == 0)) success = 0;
    a += c;
    if (!(a == c)) success = 0;
    a -= c;
    if (!(a.d == 0)) success = 0;

    Posit32gmp d(94.5786); //using double constructor. 
    a += d;
    if (!(a == d)) success = 0;
    a -= d;
    if (!(a.d == 0)) success = 0;
    
    if (success) cout << "+= and -= are good..." << endl;
    
    cout << "Testing *= and /=... ";
    Posit32gmp one(1);
    Posit32gmp aa(0), bb(1), cc(-2);
    aa *= bb; 
    if (!(aa.d == 0)) success = 0;
    aa /= bb;
    if (!(aa.d == 0)) success = 0;
    bb *= cc;
    if (!(bb == cc)) success = 0;
    bb /= cc;
    if (!(bb == one)) success = 0;

    Posit32gmp dd(94.5786);
    bb *= dd;
    if (!(bb == dd)) success = 0;
    bb /= dd;
    if (!(bb == one)) success = 0;
    
    if (success) cout << "*= and /= are good" << endl << endl;
    else cout << "Failed" << endl << endl;

    return success;
}

int operatorsTest2() {
    cout << "checking >, >=, <, <=, == operators..." << endl;
    Posit32 a(0), b(1), c(-2), d(-10);
    bool success = 1;
    if (a > b || !(a < b) || c > a || !(c < a)) success = 0; //mixture of negative and positive
    if (a >= b || !(a <= b) || c >= a || !(c <= a)) success = 0; //inclusive
    if (a == b || b == a || c == b || c == a) success = 0;
    if (c < d || d > c || c <= d || d >= c) success = 0; //both negative. 
    //if (success) cout << "Success" << endl << endl; 
    //else cout << "Failed" << endl << endl;
    cout << (a > b) << endl;
    cout << !(a < b) << endl;
    cout << (c > a) << endl;
    cout << !(c < a) << endl;
    return success;
}

int operatorsTest3() {
    cout << "Checking unary minus..." << endl;
    Posit32gmp b(1.5), c(-1.5);
    if (!(-b == c)) {
        Posit32gmp d = -b;
        cout << d.d << ' ' << c.d;
        cout << "Failed" << endl << endl;
        return 0;
    }
    cout << "Success" << endl << endl;
    return 1;
}

int sqrtTest() {
    cout << "Checking sqrt function" << endl;
    Posit32gmp a(4), aSq(2), b(-4), c(7.3);
    bool success = 1;
    if (!(sqrt(a) == aSq)) success = 0;
    if (success) cout << "Worked for 4" << endl << endl;
    else cout << "Failed" << endl << endl;
    return success;
}

//For making sure that we can use posits and ints/doubles interchangably in code. 
int interoperabilityTest() {
    cout << "Testing compatibility with other number types..." << endl;
    Posit32gmp a(1), b(2), c(-3);
    bool success = 1;
    if (!(a == 1) || !(b==2) || !(c == -3)) success = 0;
    if (!(a == 1.0) || !(b==2.0) || !(c == -3.0)) success = 0;
    if (!(-(a + b) == -3)) success = 0;
    if (!(a < 2) || !(b > 1) || !(c < 0)) success = 0;
    if (!(a <= 2) || !(b >= 1) || !(c <= 0)) success = 0;
    if (!(a <= 2.0) || !(b >= 1.0) || !(c <= 0.0)) success = 0;

    //Testing Arithmetic
    if (!(a + 1 == 2)) success = 0;
    if (!(a - 1 == 0)) success = 0;
    if (!(a * 2 == 2)) success = 0;
    if (!(a / 2 == .5)) success = 0;
    if (success) cout << "Passed all tests" << endl;
    else cout << "Failed at least one test" << endl;
    return success;
}

bool mpf_sub_test() {
    mpf_class a("-29.73");
    mpf_class b("-30.00");
    m_abs_diff(a, a, b);
    cout << a.get_d() << endl;
    mpf_class c("-30.29");
    m_abs_diff(c, b, c);
    cout << c.get_d() << endl;
    cout << mpf_lt(c, a) << endl;
    return 1;
}

void testHalfPrecision() {
    Posit32gmp p1(string(".188273456")), p2(string(".01986")), p3(string("-42.5"));
    cout << p1 << ", " << p2 << ", " << p3 << endl;
    bitset<32> bits = p3.d;
    cout << bits << endl;

    Posit32gmp p4 = p3 + p1;
    cout << p4 << endl;
    bits = p4.d;
    cout << bits << endl;

    half h1(0.000146606), h2(.01986), h3(-42.5);
    cout << h1 << ", " << h2 << ", " << h3 << endl;
    half h4 = h3 + h1;
    cout << h4 << endl;
}

void testbFloat() {
    puts("\nTesting bfloat\n");

    for (int i=0;i<10;++i) {
        float f = ((float) rand()) / ((float) RAND_MAX);
        bfloat16 b = f;
        cout << f << " -> " << b << endl;
        uint32_t lowerbits, upperbits;
        memcpy(&lowerbits, &f, sizeof(f));
        memcpy(&upperbits, &f, sizeof(f));
        
        lowerbits = ((lowerbits >> 16) << 16); //truncate f
        upperbits = ((upperbits >> 16) << 16) + (1 << 16); //truncate and move to next representable.

        memcpy(&f, &upperbits, sizeof(f));
        cout << "upper " << f << endl;
        memcpy(&f, &lowerbits, sizeof(f));
        cout << "lower " << f << endl;
    }

    cout << endl;
}

//Testing Round to nearest tie to even for bfloat16.
int testRounding() {
    uint32_t max = (1 << 16);
    float f1, f2;
    bfloat16 b1, b2;
    Posit32gmp p1, p2;
    mpf_class mb1, mb2;
    mpf_class mp1, mp2;

    unsigned long successesPosit(0), successesbfloat(0);
    unsigned long posit_exact(0), bfloat_exact(0);
    int stride = 64; 
    for (uint32_t i=0;i<max;i+=stride) {
        for (uint32_t j=0;j<max;j+=stride) {
            p1.set(i); p2.set(j);
            if (p1.isInfinity() || p2.isInfinity()) {
                successesPosit++; successesbfloat++; continue;
            }
            //------------------
            uint32_t i2 = i << 16; //encoding bfloats. (bfloat.f uses 16 msb, so must bit shift.)
            uint32_t j2 = j << 16;
            memcpy(&f1, &i2, sizeof(i2));
            memcpy(&f2, &j2, sizeof(j2));
            b1 = bfloat16(f1);
            b2 = bfloat16(f2);
            if (isnan(b1.f) || isnan(b2.f) || isinf(b1.f) || isinf(b2.f)) {
                successesPosit += 1; successesbfloat += 1; continue;
            }
            //------------------
            p1.get(mp1); //convert posit representation to mpf.
            p2.get(mp2);
            mb1 = b1.f; //convert bfloat16 representation to mpf.
            mb2 = b2.f;
            mpf_class p_result_precise = mp1 + mp2; //loss free result. (Posit)
            mpf_class b_result_precise = mb1 + mb2; //loss free result. (bfloat16)
            //------------------
            bfloat16 bfloat_result = b1 + b2; //lossy results
            Posit32gmp posit_result = p1 + p2;
            if (isnan(bfloat_result.f) || isinf(bfloat_result.f) || posit_result.isInfinity()) {
                successesPosit += 1; successesbfloat += 1; continue;
            }
            
            mpf_class b_result_mpf, p_result_mpf; //lossy result in mpf format. 
            b_result_mpf = bfloat_result.f;
            posit_result.get(p_result_mpf);
            //------------------
            //Now we need to find the two adjacent numbers for posit and bfloat and store them in mpf. . . 
            //First, for posit. . .
            Posit32gmp upperP(posit_result), lowerP(posit_result);
            upperP.d += 1; lowerP.d -= 1;
            if (upperP.isInfinity() || lowerP.isInfinity()) {successesPosit += 1; successesbfloat += 1; continue;}

            mpf_class mpf_upper_posit, mpf_lower_posit;
            upperP.get(mpf_upper_posit); lowerP.get(mpf_lower_posit);
            
            //Not sure if this is necessary.
            if (mpf_lt(mpf_upper_posit, mpf_lower_posit)) swap(mpf_upper_posit, mpf_lower_posit);
            //------------------
            //Now we need to do same thing for bfloat. Note that bfloat.f uses 16 most signifcant bits. 
            uint32_t bits_bfloat_result; 
            memcpy(&bits_bfloat_result, &bfloat_result.f, sizeof(bits_bfloat_result));
            uint32_t bitsUpper = bits_bfloat_result + (1 << 16); uint32_t bitsLower = bits_bfloat_result - (1 << 16);
            float upperf, lowerf;
            memcpy(&upperf, &bitsUpper, sizeof(upperf));
            memcpy(&lowerf, &bitsLower, sizeof(lowerf));
            if (isnan(upperf) || isinf(upperf) || isnan(lowerf) || isinf(lowerf)) {
                successesPosit += 1; successesbfloat += 1; continue;
            }
            mpf_class mpf_upper_bfloat16 = upperf; mpf_class mpf_lower_bfloat16 = lowerf;
            
            //If these are flipped, then they are both negative. So just flip them. 
            if (mpf_lt(mpf_upper_bfloat16, mpf_lower_bfloat16)) swap(mpf_upper_bfloat16, mpf_lower_bfloat16);
            //------------------
            
            mpf_class diffUpper, diffLower;
            //Now we have mpf representations of adjacent numbers of result. So time to check that rounding 
            //worked correctly.
            //For posits . . .
            if (mpf_lt(p_result_precise, mpf_lower_posit) || 
                mpf_lt(mpf_upper_posit, p_result_precise)) { //then we're way off. 
                cout << "Proper posit result is not in between two closest representables." << endl;
                return 0;
            }
            if (mpf_lt(p_result_mpf, p_result_precise)) { //then we should have rounded down. 
                m_abs_diff(diffUpper, p_result_precise, mpf_upper_posit); //diffUpper := distance from perfect result to rounded up result 
                m_abs_diff(diffLower, p_result_precise, p_result_mpf); //p_result_mpf is our actual result (rounded down theoretically)
                if (mpf_lt(diffLower, diffUpper)) successesPosit += 1;
                else if (!mpf_lt(diffUpper, diffLower)) successesPosit += (posit_result.d % 2 == 0);
                else { cout << mp1 << " + " << mp2 << " failed." << endl; } 
            }
            else if (mpf_lt(p_result_precise, p_result_mpf)) { //then we should have rounded up. 
                m_abs_diff(diffUpper, p_result_precise, p_result_mpf); //p_result_mpf is our actual result (rounded up theoretically)
                m_abs_diff(diffLower, p_result_precise, mpf_lower_posit); //diffUpper := distance from perfect result to rounded down result
                if (mpf_lt(diffUpper, diffLower)) successesPosit += 1;
                else if (!mpf_lt(diffLower, diffUpper)) {
                    successesPosit += (posit_result.d % 2 == 0);
                }
                else { cout << mp1 << " + " << mp2 << " failed." << endl; } 
            } 
            else {
                posit_exact += 1;
                successesPosit += 1; //Then we should have represented the value exactly. 
            }
            //For bfloat . . . (exactly same idea as above.)
            
            if (mpf_lt(b_result_precise, mpf_lower_bfloat16) || 
                mpf_lt(mpf_upper_bfloat16, b_result_precise)) { //then we're way off. 
                cout << "Proper bfloat result is not in between two closest representables." << endl;
                return 0;
            } 
            else if (mpf_lt(b_result_mpf, b_result_precise)) { //then we should have rounded down. 
                m_abs_diff(diffUpper, b_result_precise, mpf_upper_bfloat16); //diffUpper := distance from perfect result to rounded up result 
                m_abs_diff(diffLower, b_result_precise, b_result_mpf); //p_result_mpf is our actual result (rounded down theoretically)
                if (mpf_lt(diffLower, diffUpper)) successesbfloat += 1;
                else if (!mpf_lt(diffUpper, diffLower)) successesbfloat += ((bits_bfloat_result >> 16) % 2 == 0);
                else { cout << mb1 << " + " << mb2 << " failed." << endl; }
            }
            else if (mpf_lt(b_result_precise, b_result_mpf)) { //then we should have rounded up. 
                m_abs_diff(diffUpper, b_result_precise, b_result_mpf); //p_result_mpf is our actual result (rounded up theoretically)
                m_abs_diff(diffLower, b_result_precise, mpf_lower_bfloat16); //diffUpper := distance from perfect result to rounded down result
                if (mpf_lt(diffUpper, diffLower)) successesbfloat += 1;
                else if (!mpf_lt(diffLower, diffUpper)) successesbfloat += ((bits_bfloat_result >> 16) % 2 == 0); 
                else { cout << mb1 << " + " << mb2 << " failed." << endl; }  
            } 
            else {
                bfloat_exact += 1;
                successesbfloat += 1; //Then we should have represented the value exactly. 
            }

            /*cout << endl << "Posit: " << endl;
            cout << "\t" << mp1 << "+" << mp2 << "=" << posit_result << " : " << "one less=" << mpf_lower_posit
                << " one above=" << mpf_upper_posit << endl;
            cout << "\t" << "Perfect result: " << p_result_precise << endl;

            cout << endl << "bfloat: " << endl;
            cout << "\t" << mb1 << "+" << mb2 << "=" << bfloat_result << " : " << "one less=" << mpf_lower_bfloat16
                << " one above=" << mpf_upper_bfloat16<< endl;
            cout << "\t" << "Perfect result: " << b_result_precise << endl;*/

        }
    }
    cout << "Posits had " << successesPosit << " successes." << endl;
    cout << "bfloats had " << successesbfloat << " successes." << endl;
    cout << "Posit exact representations " << posit_exact << endl;
    cout << "Float exact representations " << bfloat_exact << endl << endl; 
    return 1;    
}

/*
void testSin() {
    Posit32gmp p(string("0.00561523")), p2(string("3.8147e-06"));
    Posit32gmp r = p*p2;

    posit16 sp, sp2;
    sp.value = p.d;
    sp2.value = p2.d;

    cout  << sp << " " << sp2 << endl;
    posit16 rp = sp*sp2;

    mpf_class m1, m2, rmp;
    p.get(m1);
    p2.get(m2);

    cout << m1 << " " << m2 << endl;
    rmp = m1*m2;

    cout << rmp << endl;

    mpf_class rm, rsm;
    r.get(rm);
    Posit32gmp t;
    t.d = rp.value;
    t.get(rsm);

    cout << rm << " " << rsm << endl;

    mpf_class gmean = rm*rsm;
    m_sqrt(gmean, gmean);

    mpf_class amean = (rm+rsm)/2;

    cout << rmp - gmean << endl;

    cout << r << " " << rp << endl;

    bitset<32> bits = r.d;
    bitset<32> bits2 = rp.value;
    cout << bits << " " << bits2 << endl;

    cout << gmean << endl;
}*/


void testArithmetic() {
    uint64_t max = ((uint64_t) 1) << ((uint64_t) SIZE);

    uint stride=8;

    Posit32gmp a, b, ab; 
    posit16    as, bs, abs;
    //mpf_class  amean, gmean;

    for (uint64_t i=0;i<max;i+=stride) {
        a.d = i;
        as.value = a.d;
        if (a.isInfinity() && as.isNaR()) continue;
        
        
        for (uint64_t j=0;j<max;j+=stride) {
            b.d = j;
            bs.value = b.d;

            ab  = a+b;
            abs = as+bs;

            if (ab.d != abs.value) {

                cout << "gmp" << a << " + " << b << " = " << ab << endl;
                cout << "soft" << as << " + " << bs << " = " << abs << endl << endl;

            } 

            ab  = a-b;
            abs = as-bs;

            if (ab.d != abs.value) {

                cout << "gmp" << a << " - " << b << " = " << ab << endl;
                cout << "soft" << as << " - " << bs << " = " << abs << endl << endl;

            } 

            ab  = a*b;
            abs = as*bs;

            if (ab.d != abs.value) {

                cout << "gmp" << a << " * " << b << " = " << ab << endl;
                cout << "soft" << as << " * " << bs << " = " << abs << endl << endl;
            } 

            ab  = a/b;
            abs = as/bs;
            
            //abs.value
            if (ab.d != abs.value) {
                //mpf_class m1, m2, rm;
                //a.get(m1);
                //b.get(m2);
                //rm = m1*m2;

                cout << "gmp" << a << " / " << b << " = " << ab << endl;
                cout << "soft" << as << " / " << bs << " = " << abs << endl << endl;

                //bitset<32> bits  =  ab.d;
                //bitset<32> bits2 = abs.value;
                //cout << bits << " " << bits2 << endl;
            } 
        }
        
    }
}

void testRound() {
    Posit32gmp a=-440;
    Posit32gmp b=125;
    cout << (a-b) << endl;
}

int main(int argc, char*argv[]){
    /*
    Posit32 p(1,6);
    ofstream file;
    if(argc > 1) file.open(argv[1], ofstream::trunc); 
        
    //a : input1, b : input2, result : a + b. 
    Posit32gmp posit_a(p), posit_b(p), posit_result(p); // Posits
    mpf_class mpf_a, mpf_b, mpf_result; // Maximum precisions
    float float_a, float_b, float_result; // Single precisions
    double double_a, double_b, double_result; // Double precisions
    
    unsigned int maxbits=((1<<(p.size)) - 1);
    int i,j;

    //storing all representations in set to later determine if final outputs from positGMP
    //arithmetic are either equal to a representable number (no loss) or are equal to the closest
    //representable number rounded down. 
    
    //comparator to ensure ordering in set is correct. 
    struct mpf_comp { bool operator () (mpf_class a, mpf_class b) {return mpf_gt(a, b);} };  
    set<mpf_class, mpf_comp> representables;
    for (int i=0; i<maxbits; i++) {
        p.set(i);
        if (p.isInfinity()) continue;
        mpf_class representation;
        p.get(representation); //store representation in infinite precision. 
        representables.insert(representation); // add representation to set. storing "wheel values".
    }

    for(i=0; i<maxbits; i++){
        
        posit_a.set(i); //posit_a <-- i
        posit_a.get(mpf_a); //posit_a --> mpf_a 
        double_a = posit_a.getValueDouble(); // double_a <-- posit_a as double
        float_a = (float) double_a;
        if(posit_a.isInfinity()) continue;
        
        for(j=0; j<maxbits; j++){
            posit_b.set(j); //posit_b <-- j
            posit_b.get(mpf_b); //posit_b --> mpf_b 
            double_b = posit_b.getValueDouble(); // double_b <-- posit_b as double
            float_b = (float) double_b;
            if(posit_b.isInfinity()) continue;
            
            //Now performing arithmetic.
            posit_result = posit_a + posit_b; // cast to mpf, add mpfs, cast back to posit_result. 
            mpf_result = mpf_a + mpf_b; // perfect precision reference. 
            double_result = double_a + double_b; // double
            float_result = float_a + float_b; // single
            // printf("da+db==>DR  %le+%le==>%le\n",da,db,dr);
            {
                // using mpf to store float, double, posit errors.
                mpf_class mpf_double_result, mpf_float_result, mpf_posit_result; 
                mpf_class posit_diff, double_diff, float_diff;
                
                // converting results of arithmetic above to mpf. 
                mpf_double_result = double_result; 
                mpf_float_result = float_result;
                posit_result.get(mpf_posit_result);
                bitset<32> bits(posit_result.d);
                
                m_diff(posit_diff, mpf_posit_result, mpf_result); // error for posit precision
                m_diff(double_diff, mpf_double_result, mpf_result); // error for double precision
                m_diff(float_diff, mpf_float_result, mpf_result); // error for single precision
                bool success = mpf_gte(mpf_posit_result, (*representables.lower_bound(mpf_result)))
                    && mpf_lte(mpf_posit_result, (*representables.upper_bound(mpf_result)));
                
                //std::cout << mpf_a << "+" << mpf_b << "=" << mpf_posit_result << "\t" //posit addition
                //    << mpf_a << "+" << mpf_b << "=" << mpf_result << "\t" //vs. mpf addition.
                //    << "posit error: " << posit_diff << " success: " << success << "\t"
                //    << "double result: " << mpf_double_result << ", double error: " << double_diff << "\t" 
                //    << "float result: " << mpf_float_result << ", float error: " << float_diff  
                //    << "bits: " << bits << endl << endl; 

                //std::cout << "Lower bound representation is: " << (*representables.lower_bound(mpf_result)) << endl;
                //Updating input file.
                
                if (file && !(file.fail())){
                    if(i==0 && j==0){ 
                        file << "posit a," << "posit b," << "posit a+b," 
                            << "mpf a," << "mpf b," << "mpf a+b,"
                            << "posit error," << "did it work?,"
                            << "double result," << "double error,"
                            << "float result," << "float error," << "bits," << "rounding" << endl;
                    }
                    file << mpf_a << "," << mpf_b << "," << mpf_posit_result << ","
                        << mpf_a << "," << mpf_b << "," << mpf_result << ","
                        << posit_diff << "," << success << ","
                        << mpf_double_result << "," << double_diff << ","
                        << mpf_float_result << "," << float_diff << "," << bits
                        << endl;
                }
                
            }
        }
    }
    file.close();
    */

    //cout << "Running unit tests..." << endl;
    //Posit32gmp p = 2.5;
    //cout << p.getUseed() << endl;
    //constructorsTest();
    //operatorsTest1();
    //operatorsTest2();
    //operatorsTest3();
    //sqrtTest();
    //interoperabilityTest();
    //testHalfPrecision();
    //testbFloat();
    ////testRounding();
    ////testSin();
    //testArithmetic();
    //testRound();
    ofstream file;
    file.open("MatrixAlgebra/positprecision.csv", ofstream::trunc);
    int fESbits=8;
    for (double value=pow(10, -12);value<=pow(10, 12);value*=1.5) {     
        Posit32gmp p=value;
        float f=value;
        int advantage =  fESbits - (1 + ES + p.getRegimeNbits());
        if (!p.isInfinity())
            file << value << "," <<  pow(10, log10(pow(2,-(24+advantage))))*value << endl; //pow(10, log10(pow(2,-(53))))*value << endl;
        else
            file << value << "," << "-" << endl;
    } 
    file.close();
    
    /*
        have encode and decode functions as members
        encodePostFromFloat32, encodePositFromMPF()
        and complementary decode functions
        decodePositToFloat32 decodePositToFloat64(), decodePositToMPF().
     
        don"t forget to do descriptions of PathForward work
     */
    
    return 1;
}
