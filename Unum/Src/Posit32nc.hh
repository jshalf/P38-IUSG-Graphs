#ifndef __POSIT32NC_HH_
#define __POSIT32NC_HH_

#include "PositBase.hh"

/************************************
 Class: Posit32nc
 Extends: PositBase
 Description: There is an ambiguity in the Posit representation
    as described in the Mathematica notebook.  The regime bits could
    be shown as two's-complement representation or in a no-complement
    representation.  This variation of the Posit32 class is to understand
    what the ramifications are of leaving the format as a natural integer
    rather than as a twos-complement representation.  It is unclear what
    the purpose is of the twos-complement representation as there is rarely
    any addition or subtraction for the regime (so seems to be little point
    for having the twos-comp representation when we could just read off the
    sign bit and be done with it).
 See Also:  Posit32 (the version implemented using the representation that
    is in the mathematic notebook)
 ************************************/
class Posit32nc : public PositBase {
    // int32_t twosComp(int32_t v);
public:
    int32_t d; // Storage for the POSIT bit representation
    // public:
    Posit32nc(size_t _es,size_t maxsize=32);
    Posit32nc(const Posit32nc &a):PositBase(a){this->d=a.d;}
    // Posit32nc(const Posit32nc &other);
    // Posit32(size_t _es):PositBase(sizeof(int32_t)*8,_es),d(0){}
    virtual bool getSignBit(){
        register uint32_t msb=1L<<(size-1);
        return (this->d & msb)?1:0;
    }
    inline virtual void set(int32_t v) { d=v; }
    // virtual void set(int32_t v){ d=(int32_t)v; }
    inline bool regimePolarity(){ // returns 1 if positive and 0 if negative?
        uint32_t mask=1L<<(size-1);
        mask>>=1;
        return polarity(mask & this->d);
    }
    virtual inline void setToInfinity(){
        register uint32_t mask=1L<<(size-1);
        this->d=mask;
    }
    virtual inline bool isInfinity(){
        register uint32_t mask=1L<<(size-1);
        if(this->d == mask) { return 1;}
        else return 0;
    }
    virtual inline void setToZero(){this->d=0;}
    virtual inline bool isZero(){ if(this->d==0) return 1; else return 0; }
    virtual unsigned int getRegimeNbits();
    virtual long getRegimeBits();
    virtual unsigned int getExponentNbits();
    virtual unsigned int getExponentBits();
    virtual unsigned int getExponentValue(); // get the actual value
    virtual unsigned int getFractionNbits();
    virtual unsigned int getFractionBits();
    virtual unsigned int getMaxFractionBits();
    virtual int getUseed();
    virtual void printBinary();
    void printInfoCompact();
    // not yet filled out
    friend Posit32nc operator+ (Posit32nc a,Posit32nc b); // should we do ref to a or copy to stack?
    friend Posit32nc operator- (Posit32nc a,Posit32nc b);
    friend Posit32nc operator* (Posit32nc a,Posit32nc b);
    friend Posit32nc operator/ (Posit32nc a,Posit32nc b);
    Posit32nc &operator= (const Posit32nc &a);
    
    // Encode/Decode
    virtual void set(char *s);
    void set(mpf_class a); // set using exact version
    void get(mpf_class &a); // transcode Posit to an MPF object
    void get(float &a);
    // void get(double &a);
    // void get(int32_t &a); // get to nearest integer?
};



#endif

