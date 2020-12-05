#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <stdio.h>
#include "gmp_helpers.hh"
#define _USE_MATH_DEFINES
using namespace std;

// GMP docs recommend you do NOT use the EQ operation as it is poorly defined (intervals are better)

// lets prototype all of the helpers for GetM using overloading
void m_abs(mpf_class &r,mpf_class &a) {mpf_abs(r.get_mpf_t(),a.get_mpf_t());}
int m_sgn(mpf_class &a) {
    return mpf_sgn(a.get_mpf_t());
}

void m_recip(mpf_class &r){
    mpf_class t=r;
    mpf_ui_div(r.get_mpf_t(),(unsigned long)1,t.get_mpf_t());
}
void m_recip(mpf_class &r,mpf_class &op){
    mpf_ui_div(r.get_mpf_t(),(unsigned long)1,op.get_mpf_t());
}
void m_div(mpf_class &r,long num,mpf_class &denom){
    mpf_class t=r;
    if(num<0) {fprintf(stderr,"Error: negative numerator not allowed in GMP.  numerator=%lu",num); abort();}
    mpf_ui_div(r.get_mpf_t(),(unsigned long)num,denom.get_mpf_t());
}
void m_div(mpf_class &r, unsigned long a, unsigned long b){
    mpf_class t=b;
    mpf_ui_div(r.get_mpf_t(),a,t.get_mpf_t());
}
void m_pow(mpf_class &r,mpf_class op,long p){
    if(p<0) { // generate recip and multiply by positive p
        mpf_class iop,ir;
        unsigned long ip=-p;
        m_recip(iop,op);
        m_pow(r,iop,ip); // multiply recip by positive power

    }
    else {
        mpf_pow_ui(r.get_mpf_t(),op.get_mpf_t(),(unsigned long)p); // convert to an unsigned long
    }
    // std::cout << "m_pow()::Result = " << r << "\tfor input base=" << op << "^" << p << "\n";
}
void m_pow(mpf_class &r,long op,long p){  // accepts integers (for power of 2 or power of N)
    mpf_class t=op;
    m_pow(r,t,p);
}
// mpf_reldiff (mpf_t rop, const mpf_t op1, const mpf_t op2);
void m_diff(mpf_class &r,const mpf_class op1,const mpf_class op2){
    mpf_reldiff(GetM(r), GetM(op1), GetM(op2));
}

void m_sqrt(mpf_class &r, const mpf_class op) {
    mpf_sqrt(GetM(r), GetM(op));
}

//computes the absolute distance between op1 and op2.
void m_abs_diff(mpf_class &r, const mpf_class op1, const mpf_class op2) {
    mpf_sub(GetM(r), GetM(op1), GetM(op2));
    m_abs(r, r);
}

//r := sin(op). Use taylor expansion of n terms. 
void m_sin(mpf_class &r, mpf_class op, int n) {
    mpf_class pi(PI);
    int p = 1;
    mpf_class term;
    unsigned long denom=1;
    mpf_class sum(0), sign(1);
    mpf_class nhp = pi/2; 
    if (mpf_lt(op, nhp)) { op *= -1; op -= pi; }
    
    //Keep in range [-pi/2, pi/2] for best accuracy.
    mpf_class halfPeriods = ceil(op/pi);
    int opSign = ((int) halfPeriods.get_d()) % 2 == 0 ? 1 : -1; 
    op -= halfPeriods*pi;
    mpf_class op2 = op + pi;
    mpf_class distZeroLower, distZeroUpper;
    m_abs_diff(distZeroLower, op, 0);
    m_abs_diff(distZeroUpper, op2, 0);

    if (mpf_lt(distZeroUpper, distZeroLower)) {op = op2; opSign = -opSign;}
    
    op *= opSign;
    
    for(int i=0;i<n;i++) {
        m_pow(term, op, p);
        term /= denom;
        term *= sign;
        sign = -sign;
        denom *= (p+1)*(p+2);
        sum += term;
        p += 2;
    }
    r = sum;
}

//r := cos(op). Use taylor expansion of n terms. 
void m_cos(mpf_class &r, mpf_class op, int n) {
    mpf_class pi(PI);
    int p = 2;
    mpf_class term;
    unsigned long denom=2;
    mpf_class sum(1), sign(-1);
    mpf_class nhp = pi/2; 
    if (mpf_lt(op, nhp)) op *= -1;
    
    //Keep in range [-pi/2, pi/2] for best accuracy.
    mpf_class halfPeriods = ceil(op/pi);
    int opSign = ((int) halfPeriods.get_d()) % 2 == 0 ? 1 : -1;
    op -= halfPeriods*pi;
    mpf_class op2 = op + pi;
    mpf_class distZeroLower, distZeroUpper;
    m_abs_diff(distZeroLower, op, 0);
    m_abs_diff(distZeroUpper, op2, 0);

    if (mpf_lt(distZeroUpper, distZeroLower)) {op = op2; opSign = -opSign;}
    
    for(int i=0;i<n;i++) {
        m_pow(term, op, p);
        term /= denom;
        term *= sign;
        sign = -sign;
        denom *= (p+1)*(p+2);
        sum += term;
        p += 2;
    }
    r = opSign*sum;
}

mpf_class sin(mpf_class r) {
    m_sin(r, r, 20);
    return r;
}

mpf_class cos(mpf_class r) {
    m_cos(r, r, 20);
    return r;
}

mpf_class log2(mpf_class r) {
    return log2(r.get_d());
}

mpf_class pow(long b, mpf_class p) {
    return pow(b, p.get_d());
}

// function for shifting in bits to the posit (shift in from right so we populate the MSB first)
void AppendBitTo(int32_t &r,bool b) {r<<=1; if(b) r|=1;}
void AppendBitTo(unsigned long &r,bool b) {r<<=1; if(b) r|=1;}
