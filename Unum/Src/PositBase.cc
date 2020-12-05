#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "PositBase.hh"
using namespace std;

// perhaps use general bit-manip
// extract bits range
// mask bits

// TakeBits(int32_t &result, int32_t source,startbit, endbit)
// PutBits(int32_t &target, int32_t source,startbit, endbit); // put bits into target

/* f = Drop[RealDigits[y, 2, Max[2, nbits + 2 - Length[pbits]]],1]
 // (* Fraction bits *);
 pbits = Join[pbits, f];
 pbits = Take[pbits, nbits + 1] //  (* Keep the next bit *);
 */

//uint32_t getMaskMSB(){
//    return 1>>(this->size());
//}

void printBin(int32_t v){
    uint32_t mask=0x080000000; // is this a 32bit mask?
    // int32_t mask = 1<<31;
    int i;
    for(i=32;i>0;i--,mask>>=1){
        if(mask&v) putchar('1'); else putchar('0');
        if(!((i-1)%4)) putchar(' ');
    }
    putchar('\n');
    printf("\tIn Hex=[%0X]\n",v);
}


int64_t twosComp(int64_t v, size_t nbits){
    if (nbits == 64) return -v; //Special case, 1<<64 will wrap around otherwise. 

    int64_t twoscomp, mask;
    // nbits--;  // pull off the sign bit
    // do twos comp on remaining bits
    mask= 1lu << nbits; 
    mask--;
    return (-v) & mask;
}

int32_t twosComp(int32_t v,size_t nbits){
    if (nbits == 32) return -v; //Special case, 1<<32 will wrap around otherwise. 

    int32_t twoscomp, mask;
    // nbits--;  // pull off the sign bit
    // do twos comp on remaining bits
    mask=1<<nbits; 
    mask--;
    return (-v) & mask;
    
    //mask= 1<<nbits;
    //printf("\tmask "); printBin(mask);
    ////mask = 2**(num_bits - 1)
    //printf("\t -v"); printBin(-v)
    //printf("\t -v&mask="); printBin(-(v&mask));
    //printf("\t v & ~mask="); printBin(v&~mask);
    //return -(v & mask) + (v & ~mask);
}


PositBase::PositBase(size_t _es,size_t _size):es(_es),size(_size){
    // size_t useed1 = ipow(2,ipow(2,es)); // 2^2^es
    if (es==0) this->useed = 2;
    else       this->useed = 2L<<((2L<<(es-1L))-1L);
    //size_t useed2 = 2L<<((2L<<(es-1L))-1L); // equivalent to above expression, but faster
    // if(useed1!=useed2) { std::cout << "error: ipower() mismatch in PositBase for useed:  ipow() metho="<<useed1 <<" and shift method=" << useed2; exit(0);}
    //  else std::cout<<"useed is OK";
    //this->useed=useed2;
    
}

void PositBase::get(double &d){
    d=getValueDouble();
}

void PositBase::printInfo(){
    std::cout << "Stats:\n\tES=" << es << "\n\tsize=" << size << "\n\tuseed=" << useed << "\n";
    // std::cout << "Number Contained:";
    std::cout << "\t"; printBinary();
    std::cout << "\tsign=" << this->getSignBit();
    std::cout << "\n\tnregime=" << this->getRegimeNbits();
    std::cout << "\n\tregime=" << this->getRegimeBits();
    std::cout << "\n\texp_nbits=" << this->getExponentNbits();
    std::cout << "\n\texponent=" << this->getExponentBits();
    std::cout << "\n\tfraction nbits=" << this->getFractionNbits();
    std::cout << "\n\tfraction=" << this->getFractionBits();
    std::cout << "\n\tMax:Min ["  << this->getMax() << ":" << this->getMin() << "]\n";
}

double PositBase::getMax(){
    double duseed=(double)this->useed;
    double maxregime = (double)(this->size - 1);
    return pow(duseed,maxregime);
}

double PositBase::getMin(){
    double duseed=(double)this->useed;
    double maxregime = (double)(this->size - 1);
    return pow(duseed,-maxregime);
}

double PositBase::getValueDouble(){ // for sanity checking
    double v;
    double duseed = (double)this->useed;
    double des = (double)this->getExponentValue();
    double dregime = (double)this->getRegimeValue();
    double dfrac = (double)this->getFractionBits();
    int maxfracbits = this->getFractionNbits();
    double dmaxfrac; // = (double)this->getMaxFractionBits();
    if(isZero()) return 0.0;
    if(isInfinity()) { return 1.0/0.0;} // NaN
    
    //if(this->getRegimeBits() >= (this->size-2)){
    //    if(this->getSignBit()) return 1.0/0.0; // infinity
    //    else return 0.0;
    //}
   // printf("Useed[%lf]^dregime[%lf]->%lf  * 2^des[%lf]->%lf dfrac[%lf]/dmaxfrac[%lf]-->%lf\n",duseed,dregime,pow(duseed,dregime),des,pow(2.0,des),dfrac,dmaxfrac,(dmaxfrac>0)?(1.0+dfrac/dmaxfrac):0.0);
    v=pow(duseed,dregime) * pow(2.0,des) * ((this->getSignBit())?(-1.0):(1.0)); // fraction bits might not be duseed, but divide by max value of the dfrac
    if(maxfracbits>0) {
        // compute the represnetational range given the number of bits
        int maxfrac = 1<<maxfracbits; // max number representable with fraction
        dmaxfrac = (double)maxfrac;
       // printf("\tmaxfracbits=%u maxfrac=%u dmaxfrac=%lf dfrac=%lf frac = %lf\n", maxfracbits,maxfrac,dmaxfrac,dfrac, 1.0+dfrac/dmaxfrac);
        v*=(1.0 + dfrac/dmaxfrac);
    }
    // Fix Fraction Bits
    // Next do zero check and infinity check before returning a result
    return v;
}
