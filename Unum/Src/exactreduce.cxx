#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <stdint.h>
#define DEBUG

/**** when we get fancier, can do full inheritance from pure virtual base class

 class Accumulator{ // pure virtual base class
  
 };

 class FloatAccumulator: public Accumulator{ // inherit for single precision accumulator
  
 };

 class DoubleAccumulator: public Accumulator{ // inherit for double precision accumulator
 };

 class BigIntAccumulator: public Accumulator{ // inherit for 
 };
*/

void printBinary(uint64_t v){ // debug print
  int i;
  printf("v=%09lld \tb=",v);
  for(i=0;i<64;i++){
    if(!(i%4)) printf(" ");
    if((v<<i)&0x8000000000000000 ) printf("1");
    else printf("0");
  }
}

void printBinaryVerbose(uint64_t v){ // debug print
  int i;
  printf("v=%09lld \n",v);
  for(i=0;i<64;i++){
    if((v<<i)&0x8000000000000000 ) printf("%03u: 1\n",63-i);
    else printf("%03u: 0\n",63-i);
  }
}

void printBinary(uint32_t v){ // debug print
  int i;
  printf("v=%09u \tb=",v);
  for(i=0;i<32;i++){
    if(!(i%4)) printf(" ");
    if((v<<i)&0x8000000000000000 ) printf("1");
    else printf("0");
  }
}

void printBinaryVerbose(uint32_t v){ // debug print
  int i;
  printf("v=%09u \n",v);
  for(i=0;i<32;i++){
    if((v<<i)&0x8000000000000000 ) printf("%03u: 1\n",31-i);
    else printf("%03u: 0\n",31-i);
  }
}

/*
FP is (-1)^sign x 2^(exponent-exponentbias) x 1.mantissa
exponent bias for IEEE 754 is 1023
*/

union IntDouble {
  double d;
  uint64_t i;
};

class BigInt {
  //1000=8
  //0111=7 63-60
  //1111=f 59-56
  //1111=f  55-52
  //zmask=0x0000000000000000
  // building up masks for double precision number
  static const uint64_t smask=0x8000000000000000,emask=0x7FF0000000000000,mmask=0x0000FFFFFFFFFFFF,exp_bias=1023,normbit=0x0001000000000000; 
protected:
  // DP format
  // Sign: 63
  // Exponent: 62-52  (11 bits)
  // Mantissa: 0-51 (52 bits)
  uint64_t bigint[32]; // 2048 bits integer
  union IntAlias {
    uint64_t *u64;
    uint32_t *u32;
    unsigned char *u8;
  };
  IntAlias subword; // for accessing sub-words of big-int
  bool bigsign;
public: 
  BigInt(BigInt &b){memcpy(this->bigint,b.bigint,32*sizeof(uint64_t)); subword.u64=bigint;}
  BigInt(double d){this->set(d); subword.u64=bigint;}
  BigInt(){bzero(bigint,sizeof(uint64_t)*32); subword.u64=bigint;}
  void set(double d){
    // mask off exponent & mantissa of FP number
    IntDouble intdoub;
    uint64_t di;
    uint64_t mantissa[2], exponent; // write into lower word of mantissa 
    bool carry=0, sign;
    intdoub.d=d; // assign double to union
    di=intdoub.i; // safely extract long int from union (evade type checking safely)
    // first shift the mantissa by the LSB's of the exponent
    sign=(smask&di)?1:0;
    exponent = (di & emask)>>52; // shift it back down to zero (do not apply exponent bias)
    // need to deal with 2s complement "negative" exponents for values < 1.0 (currently not handled)
    mantissa[1]=mantissa[0]=(di & mmask)||(exponent?normbit:0); // don't need to shift (adding normbit)
    puts("Mantissa: "); printBinary(mantissa[0]); puts("");
    // and then shift up by round-off bits of exponent
    
    {
      uint64_t msb,lsb; // bigshift is placement in bigint and residual is how much to shift within mantissa
      msb = ((0x0FFFC0 & exponent)>>6);  // 0xFFFF-0x003F --> 0xFFC0
      lsb = (0x0000003F & exponent);  // LSB
      // find out what will shift off the top and mask, shift down, and put into next exp
      // then shift up the bottom (sign extend shift?)

      // MSB tells us which word to select in the large integer
#ifdef DEBUG
      puts("\n EXP: "); printBinary(exponent);
      puts("\n MSB: "); printBinary(msb);
      puts("\nMSB2: "); msb>>=6; printBinary(msb);
#endif


      // LSB tells us how much to shift up mantissa within its word (might need to shift with non-int offset)
#ifdef DEBUG
      puts("\n LSB: "); printBinary(lsb); puts("\n\n");
#endif
      mantissa[1]=mantissa[0]>>(64-lsb); // insert the items into mantissa that would have shifted off left side of word
      mantissa[0]<<=lsb;// and then left-shift the rest of the mantissa (and let the bits shift off end of word)

      bigint[msb]=mantissa[0]; // assign bottom-half of mantissa to the bigint
      bigint[msb+1]=mantissa[1]; // and also what got shifted off of the top of the mantissa (overflow) 
#ifdef DEBUG
      printf("\n original mantissa: "); printBinaryVerbose(di&mmask); puts("");
      printf("\n shifted mantissa0: "); printBinaryVerbose(mantissa[0]); puts("");
      printf("\n shifted mantissa1: "); printBinaryVerbose(mantissa[1]); puts("");
#endif
    }
    
    // Then determine which word to shift the mantissa into for the bigint[] (put into a temporary array)
    // perform the integer addition of the lower parts of the mantissa and continue if carry is set
  }
  double get(){
    int top;
    // figure top bit of bigint and then project it back down to denormalized double (need to figure out how to renormalize)
    // first, scan to find top occuped word
    for(top=31;top<=0;top--){ // from MSB to LSB
      if(bigint[top]) break; // exit if we have data values in mantissa
    }
    // repack as a denormalized float
    //  <<<<<<   this is empty... still on ToDo list >>>>
    // then need to renormalize
    return 0.0;
  }
  void add(double d){ // convert to bigint first then add
    BigInt tmp(d);
    this->add(tmp);
  }
  void add(BigInt &d){
    // go from right to left and carry if overflow
    bool carry=0;
    for(register int i=0;i<32;i++){
      register uint64_t a=bigint[i];
      register uint64_t b=d.bigint[i];
      register uint64_t c;
      c=a+b; // add the two 64bit ints together
      if(carry) c++; // increment additionally if there is an overflow from previous addition
      if(c<a || c<b) carry=1; // if result is less than either input, then it overflowed (so set overflow bit)
      else carry=0; // otherwise clear the carry bit
    }
  // CPP overflow wraps around
  // but if the result is less than the original, then it has overflowed (so carry)
  // perform the operation on all elements, and then do carry op for next pass
  }
  // debug
  uint64_t getMantissa(uint64_t i){ return i&mmask; }
  uint64_t getExponent(uint64_t i){ return (i&emask)>>52; }
  uint64_t getSign(uint64_t i){ return (i&smask)?1:0;}
};


int main(){
  BigInt b;
  IntDouble u;
  u.d=120000000000.31235416;

  printBinaryVerbose(u.i);
  printBinary(u.i);
  printf("\n\n mantissa: "); printBinary(b.getMantissa(u.i)); puts("");
  printf("\n\n exponent: "); printBinary(b.getExponent(u.i)); puts("");
  printf("\n\n sign: "); printBinary(b.getSign(u.i)); puts("");
  puts("\n done \n");

  b.set(u.d);
  b.add(10000.0);
  puts("final Answer");  printBinaryVerbose(u.i); puts("\n end \n");

  return 1;
}
