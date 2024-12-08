#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include "Posit32nc.hh"

Posit32nc::Posit32nc(size_t _es,size_t maxsize):PositBase(_es,maxsize),d(0){
        if(size>32) {std::cerr << "Error: allocated > max type size.  Abort!\n"; abort();}
}
    // Posit32(size_t _es):PositBase(sizeof(int32_t)*8,_es),d(0){}

//bool Posit32nc::getSignBit(){
//        register uint32_t msb=1L<<(size-1);
//        return (this->d & msb)?1:0;
//}


unsigned int Posit32nc::getRegimeNbits(){
        // positive vs negative regimes
        int n, nbits;
        uint32_t mask; // initialize the mask to point to regime start
        register bool start_bit = regimePolarity();
        register uint32_t msb=1L<<(size-1);
        // printf("GetRegimeBits startbit polarity=%u, msb[%0X]\n",start_bit,msb);
        // perhaps check for infinity and 0 cases first
        if(isInfinity() || isZero()) return 0; // zero
        for(n=size-1,nbits=0, mask=(msb)>>1; polarity(mask & this->d)==start_bit && n>0; mask>>=1,nbits++,n--){
           // printf("\t First time through n=%u nbits[%u] mask[%0X] polarity[%d]\n",n,nbits,mask,polarity(mask & this->d));
        } // if polarity of next bit
        // printf("NumRegimeBits=%u\n",nbits);
        return nbits;
}

long Posit32nc::getRegimeBits(){
        if(this->getRegimeNbits()<=0) return 0;
        if(regimePolarity()) return (this->getRegimeNbits()-1);
        else return -(this->getRegimeNbits());
}

unsigned int Posit32nc::getExponentNbits(){
        int nbits;
        if(isInfinity() || isZero()) return 0;
        nbits = size - this->getRegimeNbits() - 2; // why was it getSign() before?
        if(nbits<0) return 0;
        else { // return min of nbits or es
            if(nbits>es) return es;
            else return nbits;
        }
}

unsigned int Posit32nc::getExponentValue(){
    int32_t v;
    int ebits=this->getExponentNbits();
    //std::cerr << "num ebits=" <<ebits <<'\n';
    v=this->d; // else don't invert (positive number is not twoscomp
    //printf("\n twosCompCorrected v="); printBin(v); puts("\t ***********\n");
    
    //printf("ebits=%u and es=%lu\n",ebits, this->es);
    if(ebits<=0) {
        // puts("no ebits");
        return 0;
    }
    if(ebits < this->es){  // rep bits removes some bits but must correspondingly shift left the result
        uint32_t mask=(1<<ebits); // mask off the ebits
        int b;
        mask-=1;
        b = v & mask; // so mask off the correct bits
        // and then shift up by es-ebits
        b<<=(this->es - ebits);
        //printf("\n\t ebits<this->es mask[%0X] D&Mask[%0X]\n",mask,(v & mask));
        return b;
    }
    else if(ebits >= this->es){ // if number of exp bits is greater
        // now just need to worry about how much to right-shift the exponent down for mask
        uint32_t mask=1<<(this->es); // set mask for es
        int32_t nshift;
        mask -= 1;
        // debug("\tval before shift",v);
        nshift=(this->size - ebits - 2 - this->getRegimeNbits()); // check to see if this changes for neg regime?  (no, it doesn't change)
        //printf("RegimeNbits=%u ebits=%u total shift=%u\n",this->getRegimeNbits(),ebits,nshift);
        //if(this->getRegimeNbits()>0) nshift--; // shift off an extra bit
        //printf("\t *** Actually shift=%u\n",nshift);
        v>>=nshift;
        // debug("\tval post shift=",v);
        //printf("shift compute is this->size[%lu]-1-ebits[%u]-regimeNbits[%u]",this->size, ebits, this->getRegimeNbits());
        //printf("ebits>=this->es mask[%0X] val[%0X] D&Mask[%0X]\n",mask,v,(v & mask));
        //debug("\tv&mask=",v&mask);
        return (v & mask);
    }
    std::cerr<<"Error:  getExponent() should not be here!  [" << ebits << "]\n";
    return 0; // catch all
}

unsigned int Posit32nc::getExponentBits(){
    int32_t v;
    int ebits=this->getExponentNbits();
    v=this->d; // else don't invert (positive number is not twoscomp
    //printf("\n twosCompCorrected v="); printBin(v); puts("\t ***********\n");
    
    //printf("ebits=%u and es=%lu\n",ebits, this->es);
    if(ebits<=0) {
        // puts("no ebits");
        return 0;
    }
    if(ebits < this->es){  // rep bits removes some bits but must correspondingly shift left the result
        uint32_t mask=(1<<ebits); // mask off the ebits
        int b;
        mask-=1;
        b = v & mask; // so mask off the correct bits
        // and then shift up by es-ebits
        // b<<=(this->es - ebits); // correction
        //printf("\n\t ebits<this->es mask[%0X] D&Mask[%0X]\n",mask,(v & mask));
        return b;
    }
    else if(ebits >= this->es){ // if number of exp bits is greater
        // now just need to worry about how much to right-shift the exponent down for mask
        uint32_t mask=1<<(this->es); // set mask for es
        int32_t nshift;
        mask -= 1;
        // debug("\tval before shift",v);
        nshift=(this->size - ebits - 2 - this->getRegimeNbits()); // check to see if this changes for neg regime?  (no, it doesn't change)
        //printf("RegimeNbits=%u ebits=%u total shift=%u\n",this->getRegimeNbits(),ebits,nshift);
        //if(this->getRegimeNbits()>0) nshift--; // shift off an extra bit
        //printf("\t *** Actually shift=%u\n",nshift);
        v>>=nshift;
        // debug("\tval post shift=",v);
        //printf("shift compute is this->size[%lu]-1-ebits[%u]-regimeNbits[%u]",this->size, ebits, this->getRegimeNbits());
        //printf("ebits>=this->es mask[%0X] val[%0X] D&Mask[%0X]\n",mask,v,(v & mask));
        //debug("\tv&mask=",v&mask);
        return (v & mask);
    }
    std::cerr<<"Error:  getExponent() should not be here!  [" << ebits << "]\n";
    return 0; // catch all
}


unsigned int Posit32nc::getFractionNbits(){
        int nbits = size - this->getRegimeNbits() - this->getExponentNbits() - 2;
        if(nbits<0) return 0;
        else return nbits;
}

unsigned int Posit32nc::getFractionBits(){
        int fbits=this->getFractionNbits();
        if(fbits<=0) return 0;
        else {
            uint32_t mask=(1<<fbits)-1;
            return (this->d & mask);
        }
        std::cerr<<"Error:  getFraction() should not be here!  " << fbits << "\n";
        return 0; // catch all
}

unsigned int Posit32nc::getMaxFractionBits(){ // what is the max fraction?
        register int nbits=getFractionNbits();
        register unsigned int val=(1<<(nbits+1))-1;
        // debug("getMaxFractrionBits",val);
        if(nbits==0) return 0;
        else return val;
}

int Posit32nc::getUseed(){ // flaw... int may be two small of a type
        return this->useed;
}

void Posit32nc::printBinary(){
        uint32_t mask=1L<<(size-1);
        int i;
        std::cout << "Size is " << this->size << " [ ";
        for(i=size;i>0;i--,mask>>=1){
            if(mask&d) putchar('1'); else putchar('0');
            if(!((i-1)%4)) putchar(' ');
        }
        printf("]\tIn Hex=[%0X]",this->d);
}

void Posit32nc::printInfoCompact(){
    this->printBinary();
    printf("\t val[nbits]=bits r[%d]=%ld e[%u]=%u f[%u]=%u value=%lf\n",
                            this->getRegimeNbits(),this->getRegimeBits(),
                            this->getExponentNbits(),this->getExponentBits(),
                            this->getFractionNbits(),this->getFractionBits(),
                            this->getValueDouble());
}

    // not yet filled out
Posit32nc operator+ (Posit32nc a, Posit32nc b){return a;}
Posit32nc operator- (Posit32nc a, Posit32nc b){return a;}
Posit32nc operator* (Posit32nc a, Posit32nc b){return a;}
Posit32nc operator/ (Posit32nc a, Posit32nc b){return a;}

Posit32nc &Posit32nc::operator= (const Posit32nc& a){
        if(this!=&a){
            // copy stuff over from a to this
            this->set(a.d);
        }
        return *this;
}

// Incomplete... need to finish with MPF class set(mpf_class) implementation.
void Posit32nc::set(char *s){
    mpf_class v(s); // decimal MPZ class
    this->set(v);
    // what are the operators for MPZ?
    // need sgn(), sqrt,  mpz_sqrt (mpz t rop, const mpz t op)
    
    // power, +, -, /, *, and
    // int cmp(op1, op2) <=, >= (Compare op1 and op2. Return a positive value if op1 > op2,
    // zero if op1 = op2,
    // or a negative value if op1 < op2.  (can use case -1, 0, 1) with default of "fail"
    // has mpz_popcount() and mpz_scan0() and scan1
    //Scan op, starting from bit starting bit, towards more significant bits,
    // until the first 0 or 1 bit (respectively) is found. Return the index of the found bit.
    
    // int mpz_tstbit (const mpz t op, mp bitcnt t bit_index) [Function]
    // Test bit bit index in op and return 0 or 1 accordingly
}


// now here comes the fun part (must twoscomp this result)
void Posit32nc::set(mpf_class v){
    // OK, for this we have to estimate best regime, then best exponent to encode
    // at each stage, also need to determine how many bits we have to work with to encode the next part of the equation
    mpf_class e,r,p;
    mpf_class two=2.0;
    int k=0,b=0; // regime
    int es_nbits=this->es;
    int32_t posit=0; // build up this posit
    int npositbits=0; // have accumulated zero bits
    int ebase;
    
    mpf_class abs_v,minposit,maxposit,m_useed;
    // Sign Bit
    //puts("Append Sign Bit --------------");
    if(m_sgn(v)<=0) { AppendBitTo(posit,1); } // its a positive number
    else { AppendBitTo(posit,0); } // its negative or infinity
    npositbits++; // either way we appended a bit
    
    // Regime:  Need to break this down into 5 cases (really is positable() + 3 classes)
    //      abs(v) is absolute value of the incoming infinite-precision GMP number
    //      abs(v) > maxPosit (max representable number)  --> set Posit to infinity
    //      abs(v) < minPosit -->  set Posit to zero
    //      abs(v) >= 1  --> Regime is 0-maxregime
    //      abs(v) < 1 --> Regime is negative (need to iterate to recip of abs(v))
    // printf("after sign nbits[%u]\n",npositbits); printBin(posit);
    m_useed=this->useed;
    m_abs(abs_v,v); // get the absolute value of 'a'
    // need to create a POSITABLE function or equivalent to determine if its in range of minPosit to maxPosit
    // minposit is useed^(-nbits+2)
    // maxposit is useed*(nbits-2)
    m_pow(maxposit,m_useed,(this->size-2)); // useed ^ (nbits-2) -> maxposit
    m_recip(minposit,maxposit); // recip of the max Posit is minposit
    // now do a range check to see if it is zero or an infinity (need to figure out infinity with GMP
    if(abs_v < minposit){ // if absv is < minimum representable posit (is it < or <= ??)
        // set to 0 and exit this function... its min posit so clamp value to zero...
        this->setToZero(); // we are done since min representable sets to zero (can exit)
        npositbits=this->size; // clamp to size
        return; // nothing more to do
    }
    if(abs_v > maxposit){ // is value great than the maximum representable posit?  (is it > or >= ??)
        // this is above the representable range of this posit (set posit value to infinity)
        this->setToInfinity();
        npositbits=this->size;
        return;
    }
    
    if(abs_v >= 1){ // if value is greater than or equal to 1, then we are in positive regime
        // so we will iterate to find best positive power for the useed
        k=0; // starting regime is 0
        r=m_useed; // start with 1.0
        AppendBitTo(posit,1); npositbits++; // will always have at least one bit
        while(r <= abs_v && npositbits<(this->size)){ // while value is less than current useed^k, then increase k (bounded by size of posit)
            r*=m_useed;
            k++;
            AppendBitTo(posit,1); npositbits++;
        }
        if(npositbits<this->size) { AppendBitTo(posit,0); npositbits++; } // cap off the regime
    }
    else { // otherwise, we are in the negative regime space (but need to iterate using
        k=1;
        mpf_class iabs_v;
        r=m_useed;
        m_recip(iabs_v,abs_v);
        AppendBitTo(posit,0); npositbits++;
        while(r < iabs_v && npositbits < this->size){ // r < iabs_v and bounded by max size of posit
            r *= m_useed; // increase regime  (* useed)
            k++;
            AppendBitTo(posit,0); npositbits++;
        }
        
        if( npositbits<this->size) { AppendBitTo(posit,1); npositbits++; } // cap off the regime (except if already max)
        k=-k; // invert the result (because we were iterating on the recip)
    }
    //printf("\tDebug:  I have k=%d and npositbits=%u regime=%u out of %lu\n",k,npositbits,npositbits-1,this->size);
    //  ********** Done with computing the Regime
    //  ********** Now compute the exponent for the ES
    // Fortunately the exponent for the ES is always positive
    //puts("Now Computing the Exponent=============="); // something has gone wrong with computing b
    
    { // push the stack
        b=0; // virtual exponent count (e_b is exponent value)
        m_pow(e,two,b);
        m_pow(r,m_useed,k);
        p = r * e; // build up our current posit (r is precomputed)
        while(p <= abs_v){ // keep increasing exponent until it is greater than target value
            b++;
            m_pow(e,two,b);
            p = r * e; // compute new interim posit
        }
        if(b>0) b--; // b always reflects one iteration past target (unless it exited before first iteration)
        //printf("****** End Iterations b=%u\n",b);
    }
    // ************ done with computing exponent
    // need to see how many bits of the exponent to push onto the end starting with MSB of computed exponent
    //
    //printf("Debug: Computed Exponent max es size=%lu\t b=%u  npositbits=%u\n",this->es,b,npositbits);
    // perhaps shift b down by nbits available
    if(this->es>0) { // make sure we have at least one ES bit designated
        int i,emask=1;
        emask <<= (this->es-1);
        // b = the exponent value... we start with the MSB and push bits onto the posit until we run out of ES or run out of space in posit to cram bits
        for(i=0; (i<this->es) && (npositbits<this->size); i++, emask>>=1){
            if(emask & b) { AppendBitTo(posit,1); npositbits++; } // append a 1
            else {          AppendBitTo(posit,0); npositbits++; } // append a zero if ES bit is 0
        }
        
    }
    
    //puts("Now Computing Fraction================");
    if(npositbits < this->size){
        unsigned long numerator=0,denominator=1,numerator_accumulate=0;
        mpf_class fraction; // keep re-initializing with frac
        mpf_class one=1;
        // just keep shifting ones into the denominator and then compare
        // put these into fractional form.
        while(npositbits < this->size){
            AppendBitTo(denominator,0); // really its just a left shift (but this fits with our bit manip paradigm)
            numerator=numerator_accumulate; // set to most recent numerator value
            AppendBitTo(numerator,1);// first try with a 1.. if it is too big, then append a zero
            m_div(fraction,numerator,denominator);
            fraction += one; // add one to the fraction;
            m_pow(r,m_useed,k);// regime = unum^k
            m_pow(e,two,b); // exponent = base&b (base accounts for shifting of exp off of posit)
            p=r*e*fraction; // build up our unum
            //std::cout << "\tFracIteration State p=" << p << "\t r:e:frac="<< r << ':' << e << ':' << fraction << "\tFrac num/denom=" << numerator << "/" << denominator << "\n";
            if(p>abs_v){ // then we went too far (so make this fraction bit a zero instead of a 1
                AppendBitTo(numerator_accumulate,0);
                AppendBitTo(posit,0); npositbits++;
            }
            else { // it worked, so we can append a 1
                AppendBitTo(numerator_accumulate,1);
                AppendBitTo(posit,1); npositbits++;
            }
            // keep doing this until we run out of bits in the posit
        }
        // Fraction is now incorporated into the posit encoding and we should be DONE!
    }
    // printf("debug: check for npositbits overflow %u vs size=%lu\n",npositbits,this->size);
    //printBin(posit);
    //****** OK, we are done.
    // write result back to the data item
    // now we need to twos-comp the posit if it is negative
    this->d = posit;
    // now we need to correct the polarity
    if(this->isInfinity() || this->isZero()) return;
    
}



void Posit32nc::get(mpf_class &a){
    mpf_class regime,exponent,frac; // accumulate values into this infinite precision class
    mpf_class one,two;
    mpf_class m_useed,m_es,m_regime,m_frac,m_maxfrac; // intermediates
    two=2; one=1;
    m_useed = this->useed; // hopefully that translates correctly
    m_es=this->getExponentBits();
    m_regime=this->getRegimeBits();
    m_frac=this->getFractionBits();
    m_maxfrac=this->getMaxFractionBits();
    if(isZero()) { a=0; return;}
    if(isInfinity()) { a=1000000000.0; return; } // how do I set an infinity for GMP???? (set to absurdly high number to prevent a floating point exception
    //  void mpf_pow_ui (mpf_t rop, const mpf_t op1, unsigned long int op2): Set rop to op1 raised to the power op2.
    //op1=m_useed.get_mpf_t();
    // can't do negative powers, so need to do the recip if raised to negative power
    // mpf_ui_div() to do the reciprocal
    if(regimePolarity()){
        std::cout << "Doing MPF Power on Regime\n";
        mpf_pow_ui(regime.get_mpf_t(),m_useed.get_mpf_t(),(unsigned long)this->getRegimeBits());
    }
    else{ // negative regime
        // must treat this differently since mpf_pow_ui doesn't understand negative exponents
        int32_t tmp = -this->getRegimeBits();
        mpf_class regime_tmp;
        mpf_pow_ui(regime_tmp.get_mpf_t(),m_useed.get_mpf_t(),(unsigned long)tmp);
        mpf_ui_div(regime.get_mpf_t(),(unsigned long)1,regime_tmp.get_mpf_t());
    }
    // should do a check to see if there are any exponent bits... if 0 then don't bother.
    std::cout << "now the exponent bits es=" << this->getExponentBits() << "\n";
    mpf_pow_ui(exponent.get_mpf_t(),two.get_mpf_t(),(unsigned long)this->getExponentBits());
   // frac=one+m_frac/m_maxfrac;
    // a=regime*exponent*frac; // put it all together
    a=regime*exponent;
    if(m_frac>0){ // bring in the fraction bits
        a=a*(1.0+m_frac/m_maxfrac);
    }
    if(this->getSignBit()) {
        mpf_class neg=a;
        mpf_neg(a.get_mpf_t(),neg.get_mpf_t());
    }
}

void Posit32nc::get(float &a){
    double v; // compute it all in double and then downconvert to float (don't round prematurely)
    double duseed = (double)this->useed;
    double des = (double)this->getExponentBits();
    double dregime = (double)this->getRegimeBits();
    double dfrac = (double)this->getFractionBits();
    double dmaxfrac = (double)this->getMaxFractionBits();
    if(isZero()) {a=0.0; return;}
    if(isInfinity()) { a=1.0/0.0; return;} // Cheap way to encode a NaN
    v=pow(duseed,dregime) * pow(2.0,des) * ((this->getSignBit())?(-1.0):(1.0));
    if(dmaxfrac>0.0) v*=(1.0 + dfrac/dmaxfrac); // dfrac == 0 if there are no fraction bits, so it works)
    a=(float)v; // downconvert the double to float
}

//Posit32nc::Posit32nc(const Posit32nc &other) : PositBase(other) {
    // Copy all relevant members from `other` to `this`.
//    this->d = other.d;  // Example: Adjust as per your class members.
    // Add similar assignments for other data members, if any.
//}

