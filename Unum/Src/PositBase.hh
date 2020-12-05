#ifndef __POSITBASE_HH_
#define __POSITBASE_HH_

#include <sys/types.h>
#include <gmpxx.h>
#include "gmp_helpers.hh"

// need bit-manip
// extract bits range
// mask bits

#define debug(x,y) std::cout<<x; printBin(y); putchar('\n');

// using namespace std;
void printBin(int32_t v);
int64_t twosComp(int64_t v, size_t nbits);

int32_t twosComp(int32_t v,size_t nbits);

class PositBase { // a base class for Posits of various sizes (might eliminate and
public:
    // globals
    size_t size;  // perhaps change name to nbits to match with Posits reference documentation
    size_t es; // 256
    size_t useed;
    //void *data;
    // local
    inline bool polarity(int a) const { if(a) return 1; else return 0;} // ret 1 positive and 0 if negative
    // should change name to polarity
    inline size_t ipow(size_t base, size_t exp){
        size_t result = 1;
        while (exp) {
            if (exp & 1)
                result *= base;
            exp >>= 1;
            base *= base;
        }
        return result;
    }
    PositBase(size_t _es,size_t _size);
    PositBase(const PositBase &a):size(a.size),es(a.es),useed(a.useed){}
    void CountBitsInPosit(int &r, int &e, int &f) const;
    virtual void setToInfinity()=0;
    virtual void setToZero()=0;
    // min posit? or positable check?  can only happen in infinite precision library
    virtual void set(char *s)=0;
    void get(double &d);
    virtual void printBinary(void)=0;
    void printInfo();
    virtual void set(int32_t v)=0;
    // virtual void set(int32_t v)=0;
    virtual bool getSignBit() const =0;
    virtual unsigned int getRegimeNbits() const =0; // need to count bits after sign (must be virtual).. this is always pos
    virtual long getRegimeBits() const =0; // this might be negative
    long getRegimeValue() const { return getRegimeBits();}
    inline virtual unsigned int getES() const {return this->es;}
    inline virtual unsigned int getExponentMaxNbits() const {return this->es;}
    virtual unsigned int getExponentNbits() const =0;
    virtual unsigned int getExponentBits() const =0;
    virtual unsigned int getExponentValue() const =0;
    virtual unsigned int getFractionNbits() const =0;
    // { return this->size-1-this->getRegimeNbits()-this->getExponentNbits(); }
    virtual unsigned int getFractionBits() const =0;
    virtual unsigned int getMaxFractionBits() const =0;
    virtual bool isZero() const =0;
    virtual bool isInfinity() const =0;
    double getMax();
    double getMin();
    double getValueDouble();
};

#endif

