#ifndef __POSIT32_HH_
#define __POSIT32_HH_

#include "PositBase.hh"
#include <string>
#include <map>
using namespace std;
#define SIZE 32
#define ES 2
#define FAST 0 //Super fast downcasting.
#define USES 0 //Collect traffic data? (overhead)

/************************************
 Class: Posit32
 Extends: PositBase
 Description: This class
 ************************************/
class Posit32 : public PositBase {
    // int32_t twosComp(int32_t v);
public:
	static map<pair<int32_t, int32_t>, long> uses;
    static bool counterInitialized;
    int32_t d;  // Storage for the POSIT bit representation
    // public:
    // This function allows us to directly compare the numerical representation
    // of the posit in the Positnc (non-complement) representation against this
    // version, which always represents the number in twos-complement if it is a negative number.
    // In theory, if I encode a number in Posit32 and Posit32nc, and then undo the twos-complement
    // represenation for the Posit32, it should be the identical bit pattern to the Posit32nc
    int32_t undoTwosComp(){ // debug function
        if(this->getSignBit()) { // must twos comp the lower bits of the result and OR back in the sign bit if it is nuked
            int32_t tc = twosComp(this->d,this->size-1);
            int32_t sgn=1L<<(this->size-1);
            return (tc | sgn); // put the sign bit back in
        }
        else return this->d;
    }
    Posit32();
    Posit32(size_t _es,size_t maxsize=32);
    Posit32(const Posit32 &a):d(a.d),PositBase(a){}
    Posit32(const int n);
    Posit32(const double d);
    Posit32(string s);
    
    // Posit32(size_t _es):PositBase(sizeof(int32_t)*8,_es),d(0){}
    virtual bool getSignBit() const;
    inline virtual void set(int32_t v) { d=v; }
    // virtual void set(int32_t v){ d=(int32_t)v; }
    inline bool regimePolarity() const{ // returns 1 if positive and 0 if negative?
        bool p;
        uint32_t mask=1L<<(size-1);
        mask>>=1;
        p=polarity(mask & this->d);
        // but if sign is negative, the invert the polarity because 2s complement
        if(this->getSignBit()) { // reverse polarity
            return p?0:1;
        }
        else return p; // return normal polarity (no inversion)
    }
    virtual inline void setToInfinity() {
        register uint32_t mask=1L<<(size-1);
        this->d=mask;
    }
    virtual inline bool isInfinity() const {
        register uint32_t mask=1L<<(size-1);
        if(this->d == mask) { return 1;}
        else return 0;
    }
    
    virtual inline void setToZero(){this->d=0;}
    virtual inline bool isZero() const { if(this->d==0) return 1; else return 0; }
    virtual unsigned int getRegimeNbits() const; // get the number of regime bits (from which you can actually calculate the value)
    virtual long getRegimeBits() const; // get the regime bits (this is just for debug purposes... its not the value... just the bits)
    virtual unsigned int getExponentNbits() const; // get the number of bits in the encoded exponent value
    virtual unsigned int getExponentBits() const;  // get the raw bits for the exponent value
    virtual unsigned int getExponentValue() const; // get the actual value of the exponent (e.g. e^b) 
    virtual unsigned int getFractionNbits() const; // get number of active bits in the fraction that are not shifted off
    virtual unsigned int getFractionBits() const;  // get the raw bits for the numerator of the fraction.  Gets f for (1+f/maxf) part of posit
    virtual unsigned int getMaxFractionBits() const; // get the raw bits for the denominator of the fraction. Gets maxf in 1+f/maxf
    virtual int getUseed() const; // gets the useed, which is 2^2^es
    virtual void printBinary(); // convenient print for the binary value of the posit (debug)
    void printInfoCompact(); //
	static void increment(int32_t); 
    static void initializeCounter();
    static void printTraffic();
    static void writeAdvantage(ostream &os);
    static map<int, double> distribution();
    static void clearCounter();
    static double distillAdvantage();
    double toDouble();
    // not yet filled out
    
    friend Posit32 operator+ (const Posit32 &a, const Posit32 &b); // should we do ref to a or copy to stack?
    friend Posit32 operator- (const Posit32 &a, const Posit32 &b);
    friend Posit32 operator- (Posit32 a);
    friend Posit32 operator* (const Posit32 &a, const Posit32 &b);
    friend Posit32 operator/ (const Posit32 &a, const Posit32 &b);

    friend bool operator< (Posit32 a, Posit32 b);
    friend bool operator> (Posit32 a, Posit32 b);
    friend bool operator<= (Posit32 a, Posit32 b);
    friend bool operator>= (Posit32 a, Posit32 b);
    friend bool operator== (Posit32 a, Posit32 b);
    friend bool operator!= (Posit32 a, Posit32 b);


    Posit32& operator= (const Posit32& a);
    friend ostream& operator<< (ostream& os, Posit32 a);

    friend Posit32 sqrt(Posit32 a);
    friend Posit32 sin(Posit32 a);
    friend Posit32 cos(Posit32 a);
    friend Posit32 abs(Posit32 a);

    // Encode/Decode
    virtual void set(char *s);
    void set(string s);
    void set(mpf_class a); // set using exact version
    void setGNU(mpf_class a);
    void setFast(mpf_class a);
    void get(mpf_class &a) const; // transcode Posit to an MPF objectcd 
    void get(float &a);
    // void get(double &a);
    // void get(int32_t &a); // get to nearest integer?
};

//Default Posit32gmp will round to nearest and tie to even.
class Posit32gmp : public Posit32 {
public:
    // public:
    Posit32gmp(size_t _es,size_t maxsize=32):Posit32(_es,maxsize){}
    Posit32gmp(const Posit32gmp &a):Posit32(a){}
    Posit32gmp(const Posit32 &a):Posit32(a){}
    Posit32gmp(const int n):Posit32(n){}
    Posit32gmp(const double d):Posit32(d){}
    Posit32gmp(string s):Posit32(s){}
    Posit32gmp():Posit32(){}
    
    friend Posit32gmp operator+ (const Posit32gmp &a, const Posit32gmp &b);
    friend Posit32gmp operator- (const Posit32gmp &a, const Posit32gmp &b);
    friend Posit32gmp operator* (const Posit32gmp &a, const Posit32gmp &b);
    friend Posit32gmp operator/ (const Posit32gmp &a, const Posit32gmp &b);
    Posit32gmp& operator+= (const Posit32gmp& b);
    Posit32gmp& operator-= (const Posit32gmp& a);
    Posit32gmp& operator*= (const Posit32gmp& a);
    Posit32gmp& operator/= (const Posit32gmp& a);
    // upconvert and downconvert operations
};

#endif

