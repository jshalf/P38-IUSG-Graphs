
#include <math.h>
#include <gmpxx.h>
#include <iostream>
#include "bfloat16.hh"
using namespace std;

bfloat16::bfloat16()         {this->f = 0;}
bfloat16::bfloat16(string s) {this->f = round(stof(s));}
bfloat16::bfloat16(int n)    {this->f = round((float) n);}
bfloat16::bfloat16(float f)  {this->f = round(f);}
bfloat16::bfloat16(double d) {
    if (((float) d) < d)
        this->f = round(d, 1);
    else  
        this->f = round(d, 0);
}

bfloat16::bfloat16(mpf_class m) {
    if (mpf_class((float) m.get_d()) < m)
        this->f = round(m.get_d(), 1);
    else  
        this->f = round(m.get_d(), 0);
}

//Round float to bfloat representation. Round to nearest tie to even.
float bfloat16::round(float in, bool sticky) {
    int trunc = 16; //the number of bits to truncate from 32 bit float. Should be 16. But its nice to tinker with. 
    if (isnan(in) || isinf(in)) { fprintf(stderr, "bfloat overflow.\n"); throw exception(); }

    f = abs(in);
    
    uint32_t fb, lower, upper;
    memcpy(&fb, &f, sizeof(f));
    if (prevent_overflow) {
        if (fb <= min_bfloat16) return getMin(); 
        if (fb >= max_bfloat16) return getMax();
    }
    
    lower = ((fb >> trunc) << trunc); //truncate f
    upper = ((fb >> trunc) << trunc) + (1 << trunc); //truncate and move to next representable.

    mpf_class lowerMPF, upperMPF, mean;
    
    float upperF, lowerF;
    memcpy(&upperF, &upper, sizeof(f)); //convert bits back to float.
    memcpy(&lowerF, &lower, sizeof(f));
    if (isinf(upperF) || isnan(upperF)) {return f==lowerF ? lowerF : upperF; }

    upperMPF=upperF;
    lowerMPF=lowerF;
    
    mean = (upperMPF+lowerMPF)/2;

    float out;
    if      (f > mean) out = upperF;
    else if (f < mean) out = lowerF;
    else               out = (((lower << (31-trunc)) == 0) || sticky) ? lowerF : upperF; //tie to even.
    

    if (in < 0) return -out;
    else        return  out;
    
    //float upperF, lowerF;
    //memcpy(&upperF, &upper, sizeof(f)); 
    //memcpy(&lowerF, &lower, sizeof(f)); 

    /*
    bool lsb, guard;
    lsb    = (fb << ((32-trunc)-1)) >> 31;
    guard  = (fb << ((32-trunc)  )) >> 31;
    sticky = sticky || ((fb << ((32-trunc)+1)) != 31);

    float out;
    if   (guard && (sticky || lsb)) out = upperF;
    else                            out = lowerF;
    */
}

float bfloat16::getMin() { uint32_t min = min_bfloat16; float f1; memcpy(&f1, &min, sizeof(f1)); return f1; }
float bfloat16::getMax() { uint32_t max = max_bfloat16; float f1; memcpy(&f1, &max, sizeof(f1)); return f1; }

bfloat16 operator+(bfloat16 a, bfloat16 b) {return bfloat16(a.f + b.f);}
bfloat16 operator-(bfloat16 a, bfloat16 b) {return bfloat16(a.f - b.f);}
bfloat16 operator-(bfloat16 a            ) {return bfloat16(0   - a.f);}
bfloat16 operator*(bfloat16 a, bfloat16 b) {return bfloat16(a.f * b.f);}
bfloat16 operator/(bfloat16 a, bfloat16 b) {return bfloat16(a.f / b.f);}

bool operator<  (bfloat16 a, bfloat16 b) {return a.f < b.f; }
bool operator>  (bfloat16 a, bfloat16 b) {return a.f > b.f; }
bool operator<= (bfloat16 a, bfloat16 b) {return a.f <= b.f;}
bool operator>= (bfloat16 a, bfloat16 b) {return a.f >= b.f;}
bool operator== (bfloat16 a, bfloat16 b) {return a.f == b.f;}
bool operator!= (bfloat16 a, bfloat16 b) {return a.f != b.f;}

bfloat16& bfloat16::operator+= (const bfloat16& b) {this->f = ((*this) + b).f; return *this;}
bfloat16& bfloat16::operator-= (const bfloat16& b) {this->f = ((*this) - b).f; return *this;}
bfloat16& bfloat16::operator*= (const bfloat16& b) {this->f = ((*this) * b).f; return *this;}
bfloat16& bfloat16::operator/= (const bfloat16& b) {this->f = ((*this) / b).f; return *this;}
bfloat16& bfloat16::operator=  (const bfloat16& b) {this->f = b.f;             return *this;}

ostream& operator<< (ostream& os, bfloat16 a) {return os << a.f;}

bfloat16 sqrt(bfloat16 a) { if (a < 0) puts("sqrt of negative."); return bfloat16(sqrt(a.f));}
bfloat16 sin(bfloat16 a)  {return bfloat16(sin(a.f)); }
bfloat16 cos(bfloat16 a)  {return bfloat16(cos(a.f)); }
bfloat16 abs (bfloat16 a) { 
    if (a < 0) return -a;
    else       return  a; 
}
