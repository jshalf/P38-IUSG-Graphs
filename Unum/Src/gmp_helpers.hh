#ifndef __GMP_HELPERS_HH_
#define __GMP_HELPERS_HH_

#include <sys/types.h>
#include <gmpxx.h>
#include <cmath>

// convenience functions
#define GetM(x) (x.get_mpf_t())

#define mpf_lt(a,b)  ((mpf_cmp(a.get_mpf_t(),b.get_mpf_t()) <  0)?1:0)
#define mpf_lte(a,b) ((mpf_cmp(a.get_mpf_t(),b.get_mpf_t()) <= 0)?1:0)
#define mpf_gt(a,b)  ((mpf_cmp(a.get_mpf_t(),b.get_mpf_t()) == 1)?1:0)
#define mpf_gte(a,b) ((mpf_cmp(a.get_mpf_t(),b.get_mpf_t()) >= 0)?1:0)
#define PI "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798\
214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881"

// GMP docs recommend you do NOT use the EQ operation as it is poorly defined (intervals are better)

// lets prototype all of the helpers for GetM using overloading
void m_abs(mpf_class &r,mpf_class &a);
int m_sgn(mpf_class &a);

void m_recip(mpf_class &r);
void m_recip(mpf_class &r,mpf_class &op);
void m_div(mpf_class &r,long num,mpf_class &denom);
void m_div(mpf_class &r, unsigned long a, unsigned long b);
void m_pow(mpf_class &r,mpf_class op,long p);
void m_pow(mpf_class &r,long op,long p);
void m_sqrt(mpf_class &r, const mpf_class op);
void m_abs_diff(mpf_class &r, const mpf_class op1, const mpf_class op2);
void m_sin(mpf_class &r, const mpf_class op, int n);
void m_cos(mpf_class &r, const mpf_class op, int n);

// Measures the relative difference between two GMP numbers (how much error)
// mpf_reldiff (mpf_t rop, const mpf_t op1, const mpf_t op2);
void m_diff(mpf_class &r,const mpf_class op1,const mpf_class op2);
// can create sub-classes that accept comparisons of various kinds of posits and doubles

// function for shifting in bits to the posit (shift in from right so we populate the MSB first)
void AppendBitTo(int32_t &r,bool b);
void AppendBitTo(unsigned long &r,bool b);

mpf_class sin(mpf_class r);
mpf_class cos(mpf_class r);
mpf_class log2(mpf_class r);
mpf_class pow(long b, mpf_class p);

#endif
