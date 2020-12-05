
#include <math.h>
#include <iostream>
#include <bitset>
#include <vector>
#include "Posit32.hh"
#include "gmpxx.h"
using namespace std;

/*
experiment(){
    e=(useed >> 1);
    while(e >= 2) {
        if(y >= e) {
            AppendTo(pbits, 1);
            y /= e;
        }
        else {
           AppendTo(pbits, 0);
        }
    }
    e >>= 1;
}
*/

map<pair<int32_t, int32_t>, long> Posit32::uses = map<pair<int32_t, int32_t>, long>();
bool Posit32::counterInitialized = false;

//Initialize/re-initialize the counter.
void Posit32::initializeCounter() {
    
    if (counterInitialized) { clearCounter(); return; } //re-initialize.
    
    Posit32gmp p;
    p.d = 1 << ES;
    Posit32gmp pH = 2*p;

    while (p != pH) {
        uses.insert(pair<pair<int32_t, int32_t>, long>(pair<int32_t, int32_t>(p.d, pH.d), 0));
        p = pH;
        pH *= 2;
    }
    
    counterInitialized = true;
}

//Adds value to traffic counter.
void Posit32::increment(int32_t p) {
    if (!counterInitialized || p==0) return;
	for (pair<pair<int32_t, int32_t>, long> pair : uses) {
        if (p >= pair.first.first && p < pair.first.second) { uses[pair.first]++; break; }
    }
}

//Displays traffic within set of intervals.
void Posit32::printTraffic() {
    if (!counterInitialized) return;
    Posit32 lower, upper;
    puts("\nDisplaying traffic.\n");
    for (pair<pair<int32_t, int32_t>, long> pair : uses) {
        lower.d = pair.first.first;
        upper.d = pair.first.second;
        cout << pair.second << " in range[" << lower << ", " << upper << "]\n";
    }
}


//Arranges traffic based on fraction bit advantage for posit. Pass 'cout' as argument
//to print data to the console. 
void Posit32::writeAdvantage(ostream& os) {
    if (!counterInitialized) return;

    int fESbits = 8; //how many exp bits floats have. Necessary for determining posit advantage.
    map<int, double> m;
    Posit32 p;
    int advantage, regime;
    long total=0;
    for (pair<pair<int32_t, int32_t>, long> pair : uses) total += pair.second; //get total uses.
    
    for (pair<pair<int32_t, int32_t>, long> pair : uses) {
        p.d = pair.first.first;
        advantage = fESbits - (1 + ES + p.getRegimeNbits()); 

        if (m.count(advantage)==0) m[advantage] = (((double) pair.second)/((double) total))*100;
        else m[advantage] += (((double) pair.second)/((double) total))*100;
    }

    string separator = "";
    for (pair<int, double> p : m) 
        { os << separator << p.first;  separator=","; }
    
    os << endl;
    separator = "";
    for (pair<int, double> p : m) 
        { os << separator << p.second; separator=","; }
    os << endl;
} 

map<int, double> Posit32::distribution() {
    if (!counterInitialized) exit(1);

    int fESbits = 8; //how many exp bits floats have. Necessary for determining posit advantage.
    map<int, double> m;
    Posit32 p;
    int advantage, regime;
    long total=0;
    for (pair<pair<int32_t, int32_t>, long> pair : uses) total += pair.second; //get total uses.
    
    for (pair<pair<int32_t, int32_t>, long> pair : uses) {
        p.d = pair.first.first;
        advantage = fESbits - (1 + ES + p.getRegimeNbits()); 

        if (m.count(advantage)==0) m[advantage] = (((double) pair.second)/((double) total))*100;
        else m[advantage] += (((double) pair.second)/((double) total))*100;
    }

    return m;
} 

double Posit32::distillAdvantage() {
    if (!counterInitialized) return 0;

    int fESbits = 8; //how many exp bits floats have. Necessary for determining posit advantage.
    map<int, double> m;
    Posit32 p;
    int advantage, regime;
    long total=0;
    for (pair<pair<int32_t, int32_t>, long> pair : uses) total += pair.second; //get total uses.
    
    for (pair<pair<int32_t, int32_t>, long> pair : uses) {
        p.d = pair.first.first;
        advantage = fESbits - (1 + ES + p.getRegimeNbits()); 

        if (m.count(advantage)==0) m[advantage] = (((double) pair.second)/((double) total));
        else m[advantage] += (((double) pair.second)/((double) total));
    }

    double avg=0;
    for (pair<int, double> p : m) avg += p.first * p.second;
    return avg;
}

void Posit32::clearCounter() {
    for (auto it=uses.begin(); it != uses.end(); it++) (*it).second = 0;
}

Posit32::Posit32():PositBase(ES, SIZE),d(0) {}

Posit32::Posit32(size_t _es,size_t maxsize):PositBase(_es,maxsize),d(0) {
        if(size>32) {std::cerr << "Error: allocated > max type size.  Abort!\n"; abort();}
}
    // Posit32(size_t _es):PositBase(sizeof(int32_t)*8,_es),d(0){}

//For working with mixture of standard numeric values in code.
Posit32::Posit32(const int n):PositBase(ES, SIZE) {
    mpf_class m(n);
    this->set(m);
}

Posit32::Posit32(const double d):PositBase(ES, SIZE) {
    mpf_class m(d);
    this->set(m);
}

Posit32::Posit32(string s):PositBase(ES, SIZE) {
    this->set(s);
}

// Retrieve the MSB that determines the sign fo the POSIT
bool Posit32::getSignBit() const {
    register uint32_t msb=1L<<(size-1);
    return (this->d & msb)?1:0;
}

unsigned int Posit32::getRegimeNbits() const  {
    // positive vs negative regimes
    int n, nbits;
    uint32_t mask; // initialize the mask to point to regime start
    uint32_t v;
    register bool start_bit = regimePolarity();
    register uint32_t msb=1L<<(this->size-1);
    if(this->getSignBit()) v=twosComp(this->d,this->size-1);  // invert two's comp if negative
    else v=this->d; // else don't invert (positive number is not twoscomp
    // printf("GetRegimeBits startbit polarity=%u, msb[%0X]\n",start_bit,msb);
    // perhaps check for infinity and 0 cases first
    if(isInfinity() || isZero()) return 0; // zero
    for(n=size-1,nbits=0, mask=(msb)>>1; polarity(mask & v)==start_bit && n>0; mask>>=1,nbits++,n--){
        // printf("\t First time through n=%u nbits[%u] mask[%0X] polarity[%d]\n",n,nbits,mask,polarity(mask & this->d));
    } // if polarity of next bit
    // printf("NumRegimeBits=%u\n",nbits);
    return nbits;
}

long Posit32::getRegimeBits() const {
    int32_t n=this->getRegimeNbits();
    if(n<=0) return 0;
    if(regimePolarity()) return (n-1);
    else return -n;
}

unsigned int Posit32::getExponentNbits() const { // works the same wheter 2s-comp or not
    int nbits;
    if(isInfinity() || isZero()) return 0;
    nbits = size - this->getRegimeNbits() - 2; // why was it getSign() before?
    if(nbits<0) return 0;
    else { // return min of nbits or es
        if(nbits>es) return es;
        else return nbits;
    }
}

unsigned int Posit32::getExponentValue() const {
    int32_t v;
    int ebits=this->getExponentNbits();
    //std::cerr << "num ebits=" <<ebits <<'\n';
    if(this->getSignBit()) v=twosComp(this->d,this->size-1);  // invert two's comp if negative
    else v=this->d; // else don't invert (positive number is not twoscomp
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
        return b;  // should this be a power of ebits???  < possible bug >
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
        
        v>>=nshift; // raising the exponent value to the power nshift
        
        // debug("\tval post shift=",v);
        //printf("shift compute is this->size[%lu]-1-ebits[%u]-regimeNbits[%u]",this->size, ebits, this->getRegimeNbits());
        //printf("ebits>=this->es mask[%0X] val[%0X] D&Mask[%0X]\n",mask,v,(v & mask));
        //debug("\tv&mask=",v&mask);
        return (v & mask);
    }
    std::cerr<<"Error:  getExponent() should not be here!  [" << ebits << "]\n";
    return 0; // catch all
}

unsigned int Posit32::getExponentBits() const {
    int32_t v;
        int ebits=this->getExponentNbits();
        //std::cerr << "num ebits=" <<ebits <<'\n';
        if(this->getSignBit()) v=twosComp(this->d,this->size-1);  // invert two's comp if negative
        else v=this->d; // else don't invert (positive number is not twoscomp
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

unsigned int Posit32::getFractionNbits() const {
    int nbits = size - this->getRegimeNbits() - this->getExponentNbits() - 2;
    if(nbits<0) return 0;
    else return nbits;
}

unsigned int Posit32::getFractionBits() const {
    int32_t v;
        int fbits=this->getFractionNbits();
    if(this->getSignBit()) v=twosComp(this->d,this->size-1);  // invert two's comp if negative
    else v=this->d; // else don't invert (positive number is not twoscomp
    
        if(fbits<=0) return 0;
        else {
            uint32_t mask=(1<<fbits)-1;
            return (v & mask);
        }
        std::cerr<<"Error:  getFraction() should not be here!  " << fbits << "\n";
        return 0; // catch all
}

unsigned int Posit32::getMaxFractionBits() const { // what is the max fraction?
        register int nbits=getFractionNbits();
        register unsigned int val=(1<<(nbits+1))-1;
        // debug("getMaxFractrionBits",val);
        if(nbits==0) return 0;
        else return val;
}

int Posit32::getUseed() const { // flaw... int may be two small of a type
        return this->useed;
}

void Posit32::printBinary(){
        uint32_t mask=1L<<(size-1);
        int i;
        std::cout << "Size is " << this->size << " [ ";
        for(i=size;i>0;i--,mask>>=1){
            if(mask&d) putchar('1'); else putchar('0');
            if(!((i-1)%4)) putchar(' ');
        }
        printf("]\tIn Hex=[%0X]",this->d);
}

void Posit32::printInfoCompact(){
    this->printBinary();
    printf("\t val[nbits]=bits r[%d]=%ld e[%u]=%u f[%u]=%u value=%lf 1/value = 1/%lf\n",
                            this->getRegimeNbits(),this->getRegimeBits(),
                            this->getExponentNbits(),this->getExponentBits(),
                            this->getFractionNbits(),this->getFractionBits(),
                            this->getValueDouble(),1.0/this->getValueDouble());
}

    // not yet filled out
Posit32 operator+ (Posit32 a, Posit32 b){return a;}
Posit32 operator- (Posit32 a, Posit32 b){return a;}
Posit32 operator* (Posit32 a, Posit32 b){return a;}
Posit32 operator/ (Posit32 a, Posit32 b){return a;}

Posit32& Posit32::operator= (const Posit32& a){
        if(this!=&a){
            // copy stuff over from a to this
            this->set(a.d);
        }
        return *this;
}


// Incomplete... need to finish with MPF class set(mpf_class) implementation.
void Posit32::set(char *s){
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

void Posit32::set(string s){
    mpf_class v(s); 
    this->set(v); 
}

void Posit32::set(mpf_class a){
    if (FAST) setFast(a);
    else      setGNU (a);
}

void Posit32::setFast(mpf_class v) { 
    if (v == 0) { this->d=0; return; }
        
    bool sign = v < 0;
    m_abs(v, v);
    
    mpf_class two=2.0;
    int es_nbits=this->es;
    
    mpf_class minposit,maxposit,m_useed;
    

    m_useed=this->useed; 
    m_pow(maxposit,m_useed,(this->size-2)); // useed ^ (nbits-2) -> maxposit
    m_recip(minposit,maxposit); // recip of the max Posit is minposit
    // now do a range check to see if it is zero or an infinity (need to figure out infinity with GMP
    if(v < minposit){ // if absv is < minimum representable posit (is it < or <= ??)
        // set to 0 and exit this function... its min posit so clamp value to zero...
        //this->setToZero(); // we are done since min representable sets to zero (can exit)
        if (sign) this->set(-minposit);
        else      this->set(minposit);
        //this->setToZero();
        return;
    }
    if(v > maxposit){ // is value great than the maximum representable posit?  (is it > or >= ??)
        // this is above the representable range of this posit (set posit value to infinity)
        //cerr << "Posit set to infinity" << endl;
        //this->setToInfinity();
        if (sign) this->set(-maxposit);
        else      this->set(maxposit);
        return;
    }
    
    bool regime, lsb, guard, sticky; 
    int64_t posit; //Leaving extra room. Will reduce to <= 32 bits later. 
    double d = v.get_d(); 
    
    regime = v >= 1; //positive regime space. 
    
    //if we lost some bits when converting to double, then record for rounding purposes. 
    if (mpf_class(d) < v) sticky=1; 
    
    //IEEE double uses biased exponent (bias of 1023), so need to subtract to get twos complement integer.
    uint64_t bias = 1023l << 52; 

    uint64_t bits; 
    memcpy(&bits, &d, sizeof(bits));

    //converting from biased exponent to twos complemented exponent.
    bits -= bias;  

    //remove sign bit if necessary. (i.e exp was less than 0 and subtracting the bias flipped the sign.)
    bits &= 0x7FFFFFFFFFFFFFFF; 
    
    if (regime) { //positive regime case. 
        
        //get posit started. lower es bits of double exponent are the same as for posit.
        //after this block, posit will be {1 (temporary bit),1,0 (regime terminator),esbits,0,0,...}
        //temporary bit is so we can arithmetic shift our regime into the posit.
        uint64_t exp = (bits << (12-this->es)) >> (64-this->es); 
        posit = 0xC000000000000000+(exp<<(64-(this->es+3))); 
        
        //get regime value from double and remove temporary bit.
        //after this posit will be {0, 1...1 (regime), 0(regime terminator), esbits, 0, 0, ...}
        int scale = bits >> (64-(12-this->es)); //get regime.  
        posit >>= scale; //apply regime
        posit &= 0x7FFFFFFFFFFFFFFF; //remove temporary bit.  
        
        //get as many fraction bits from double as possible.
        //and append them to posit. posit is now {0, 1...1 (regime), 0(regime terminator), esbits, fraction...}
        uint64_t fraction = (bits << 12) >> (this->es+scale+3); 
        posit += fraction; 

        //at this point posit has more fraction bits than is necessary. 
        //will take advantage of this for easy rounding.  
        
    } else { //negative regime case. similar but subtle differences. 
        int64_t scaling = bits >> 52; //double exp component (11 bits)

        //twos complement it, so that it gives positive value. 
        //will be useful later. 
        scaling = twosComp(scaling, 11); 

        //remove exp for posit the same as before, but in this case we will twos complement the exp. 
        //this subtlety is due to the fact that es is always positive
        int64_t exp = (scaling << (64-this->es)) >> (64-this->es); 
        exp = twosComp(exp, this->es); 

        //posit is now {0, 1 (regime terminator), esbits, 0...}
        posit = 0x4000000000000000+(exp<<(64-(this->es+2)));

        //if our scaling component is a power of useed then the case is simple. otherwise we must 
        //add one extra regime bit.
        //posit is now {0, 0...0 (regime), 1 (regime terminator), esbits, 0...}
        int32_t scale = scaling >> this->es;
        if (exp!=0) scale++;
        posit >>= scale;

        //add fraction same as before. 
        //posit is now {0, 0...0 (regime), 1 (regime terminator), esbits, fraction...}
        uint64_t fraction = (bits << 12) >> (this->es+scale+2);
        posit += fraction;
    }
    
    //now we must round.
    lsb    = (posit << this->size-1) >> 63; //least significant bit. 
    guard  = (posit << this->size) >> 63; //one bit off the edge.
    sticky = sticky || ((posit << this->size+1) != 0); //is there something after guard?
    //if guard = 1, then its at least at the half way point, either sticky will push it over the edge. or
    //if lsb is 1 then that means we should tie to even by rounding the odd lsb. 
    if (guard && (sticky || lsb)) posit += 1l << (64-this->size);

    //put posit into the right spot and convert to 32 bit integer. Also twos complement if necessary. 
    posit >>= (64-this->size);
    int64_t mask = (1l << this->size)-1;
    posit &= mask;
    if (USES) increment(posit);
    
    if (sign) posit = twosComp(posit, this->size);
    
    this->d = posit;
}

// now here comes the fun part (must twoscomp this result)
void Posit32::setGNU(mpf_class v){
    
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
    if(m_sgn(v)<=0) { AppendBitTo(posit,0); } // changed from 1!! its a negative number. We will two's complement at end. 
    else { AppendBitTo(posit,0); } // its postive or infinity
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
    if (abs_v==0) { this->d = 0; return; }
    
    if(abs_v < minposit){ // if absv is < minimum representable posit (is it < or <= ??)
        // set to 0 and exit this function... its min posit so clamp value to zero...
        //this->setToZero(); // we are done since min representable sets to zero (can exit)
        if (m_sgn(v)<=0) this->set(-minposit);
        else this->set(minposit);
        npositbits=this->size; // clamp to size
        return; // nothing more to do
    }
    if(abs_v > maxposit){ // is value great than the maximum representable posit?  (is it > or >= ??)
        // this is above the representable range of this posit (set posit value to infinity)
        //cerr << "Posit set to infinity" << endl;
        //this->setToInfinity();
        if (m_sgn(v)<=0) this->set(-maxposit);
        else this->set(maxposit);
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

    //Now we need to ensure proper rounding. (Round to nearest tie to even)
    
    Posit32 lower(this->es, this->size), upper(this->es, this->size);
    lower.d = posit; 
    upper.d = posit+1; //Next representable
    mpf_class lowerMPF, upperMPF, mean;
    
    lower.get(lowerMPF);
    upper.get(upperMPF); 

    if (upper.isInfinity()) this->d = posit;
    else {
        if (lower.getExponentNbits() < ES) { //then the guard bit is an exp bit, so we use geometric mean.
            mean = lowerMPF*upperMPF;
            m_sqrt(mean, mean);
        } else {
            mean = (lowerMPF+upperMPF)/2;
        }
        if      (abs_v < mean) this->d = posit;
        else if (abs_v > mean) this->d = posit+1;
        else                   this->d = ((posit << 31) >> 31) ? posit+1 : posit;
    }

    // now we need to correct the polarity
    //if(this->isInfinity() || this->isZero()) return;
    if      (this->isInfinity()) this->d -= 1;
    else if (this->isZero())     this->d += 1;
    
    if (USES)               increment(this->d);
    if (m_sgn(v)<=0) this->d = twosComp(this->d, this->size); //two's complement if negative. 
}
    

/*
void PositBase::CountBitsInPosit(int &r, int &e, int &f){
     // sanity check
        // how many bits for regime ?
        int regimebits = 2+k; // regime + sign bit + end-regime-bit
        int esbits = this->size - regimebits; // see how many bits we have left for the exponent
        int fracbits;
        if(esbits > this->es) esbits = this->es;
        // so how many bits do we have left for the fraction?
        fracbits = this->size - (esbits + regimebits); // this is how many bits are left for the fraction?
        if(fracbits>0){
            // compute the value of the maximum fraction
            int maxfrac = 1<<fracbits; // shift a 1 over by fracbits places
            maxfrac--; // and then subtract one from result to get actual max value (check)
        }
        // so now we know how many bits are left for the fraction... need to encode the fraction
        // use the shift on (choose closest if MSB is on vs MSB off and then right shift both top and bottom of the fraction
    f=fracbits;
    e=esbits;
    r=regimebits;
}
*/

void Posit32::get(mpf_class &a) const {
    mpf_class regime,exponent,frac; // accumulate values into this infinite precision class
    mpf_class one,two;
    mpf_class m_useed,m_es,m_regime,m_frac,m_maxfrac; // intermediates
    int maxfracbits;
    two=2; one=1;
    m_useed = this->useed; // hopefully that translates correctly
    
    // these have been corrected to properly twos-comp the binary rep if the posit is negative
    // so it should extract the correct values
    m_es=this->getExponentValue();
    m_regime=this->getRegimeValue();
    m_frac=this->getFractionBits();
    maxfracbits=this->getFractionNbits();
    if(isZero()) { a=0; return;}
    if(isInfinity()) { a=1000000000.0; return; } // how do I set an infinity for GMP???? (set to absurdly high number to prevent a floating point exception
    //  void mpf_pow_ui (mpf_t rop, const mpf_t op1, unsigned long int op2): Set rop to op1 raised to the power op2.
    //op1=m_useed.get_mpf_t();
    // can't do negative powers, so need to do the recip if raised to negative power
    // mpf_ui_div() to do the reciprocal
    if(regimePolarity()){
        // std::cout << "Doing MPF Power on Regime\n";
        mpf_pow_ui(regime.get_mpf_t(),m_useed.get_mpf_t(),(unsigned long)this->getRegimeValue());
    }
    else{ // negative regime
        // must treat this differently since mpf_pow_ui doesn't understand negative exponents
        int32_t tmp = -this->getRegimeBits();
        mpf_class regime_tmp;
        mpf_pow_ui(regime_tmp.get_mpf_t(),m_useed.get_mpf_t(),(unsigned long)tmp);
        mpf_ui_div(regime.get_mpf_t(),(unsigned long)1,regime_tmp.get_mpf_t());
    }
    // should do a check to see if there are any exponent bits... if 0 then don't bother.
    //std::cout << "now the exponent bits es=" << this->getExponentBits() << "\n";
    mpf_pow_ui(exponent.get_mpf_t(),two.get_mpf_t(),(unsigned long)this->getExponentValue()); // accounts for shift in exponent value
   // frac=one+m_frac/m_maxfrac;
    // a=regime*exponent*frac; // put it all together
    a=regime*exponent; // computing a new exp
    if(m_frac>0){ // bring in the fraction bits
        int maxfrac = 1<<maxfracbits; // max number representable with fraction
        m_maxfrac=maxfrac; // convert to infinite precision
        // OK, this is a bogus way of calculating this, but it seems to work
        a=a*(1.0+m_frac/m_maxfrac);
    }
    if(this->getSignBit()) { // tack on the sign bit
        mpf_class neg=a;
        mpf_neg(a.get_mpf_t(),neg.get_mpf_t());
    }
}

void Posit32::get(float &a){
    double v; // compute it all in double and then downconvert to float (don't round prematurely)
    double duseed = (double)this->useed;
    // these get/set routines correctly twos-comp the numbers if the sign is negative
    // so it should calculate the correct result.
    double des = (double)this->getExponentValue();
    double dregime = (double)this->getRegimeValue();
    double dfrac = (double)this->getFractionBits();
    double dmaxfrac = (double)this->getMaxFractionBits();
    if(isZero()) {a=0.0; return;}
    if(isInfinity()) { a=1.0/0.0; return;} // Cheap way to encode a NaN
    v=pow(duseed,dregime) * pow(2.0,des) * ((this->getSignBit())?(-1.0):(1.0));
    if(dmaxfrac>0.0) v*=(1.0 + dfrac/dmaxfrac); // dfrac == 0 if there are no fraction bits, so it works)
    a=(float)v; // downconvert the double to float
}

Posit32gmp operator+(const Posit32gmp &a, const Posit32gmp &b){
    if      (a.isInfinity()) return a;
    else if (b.isInfinity()) return b;
    else if (a.d == 0) return b;
    else if (b.d == 0) return a;
    mpf_class aa,bb,rr;
    Posit32gmp r(a); // use a as a prototype to set es and size for r
    a.get(aa); // decode posit to an GMP number
    b.get(bb); // decode posit to an GMP number
    rr = aa+bb; // perform the arithmetic exactly using GMP
    r.set(rr); // and then re-encode as a Posit
    return r;
}

Posit32gmp operator-(const Posit32gmp &a, const Posit32gmp &b){
    if      (a.isInfinity()) return a;
    else if (b.isInfinity()) return b;
    else if (a.d == 0) return -b;
    else if (b.d == 0) return a;
    mpf_class aa,bb,rr;
    Posit32gmp r(a); // use a as a prototype to set es and size for r
    a.get(aa); // decode posit to an GMP number
    b.get(bb); // decode posit to an GMP number
    rr = aa-bb; // perform the arithmetic exactly using GMP
    r.set(rr); // and then re-encode as a Posit
    return r;
}

Posit32gmp operator*(const Posit32gmp &a, const Posit32gmp &b){
    if      (a.isInfinity()) return a;
    else if (b.isInfinity()) return b;
    else if (a.d == 0) return a;
    else if (b.d == 0) return b;
    mpf_class aa,bb,rr; 
    Posit32gmp r(a); // use a as a prototype to set es and size for r
    a.get(aa); // decode posit to an GMP number
    b.get(bb); // decode posit to an GMP number
    rr = aa*bb; // perform the arithmetic exactly using GMP
    r.set(rr); // and then re-encode as a Posit
    return r;
}

Posit32gmp operator/(const Posit32gmp &a, const Posit32gmp &b){
    if      (a.isInfinity()) return a;
    else if (b.isInfinity()) return b;
    else if (b.isZero())     { Posit32gmp c; c.setToInfinity(); return c; }
    else if (a.d == 0) return a;
    mpf_class aa,bb,rr;
    Posit32gmp r(a); // use a as a prototype to set es and size for r
    a.get(aa); // decode posit to an GMP number
    b.get(bb); // decode posit to an GMP number
    rr = aa/bb; // perform the arithmetic exactly using GMP
    r.set(rr); // and then re-encode as a Posit
    return r;
}

bool operator< (Posit32 a, Posit32 b) {
    //return a.d < b.d; //Because Posit32.d is a signed integer, this works. Nevermind only works when size=32.
    mpf_class aa,bb,rr;
    a.get(aa); // decode posit to an GMP number
    b.get(bb); // decode posit to an GMP number
    return aa < bb;
}

bool operator<= (Posit32 a, Posit32 b) {
    return (a < b) || (a==b);
}

bool operator> (Posit32 a, Posit32 b) {
    return !(a <= b);
}

bool operator>= (Posit32 a, Posit32 b) {
    return !(a < b);
}

bool operator== (Posit32 a, Posit32 b) {
    return a.d == b.d;
}

bool operator!= (Posit32 a, Posit32 b) {
    return a.d != b.d;
}

Posit32gmp& Posit32gmp::operator+= (const Posit32gmp& b) {
    this->d = ((*this) + b).d;
    return *this;
}

Posit32gmp& Posit32gmp::operator-= (const Posit32gmp& b) {
    this->d = ((*this) - b).d;
    return *this;
}

Posit32gmp& Posit32gmp::operator*= (const Posit32gmp& b) {
    this->d = ((*this) * b).d;
    return *this;
}

Posit32gmp& Posit32gmp::operator/= (const Posit32gmp& b) {
    this->d = ((*this) / b).d;
    return *this;
}

Posit32 operator- (Posit32 a) {
    a.d = twosComp(a.d, a.size);
    return a;
}

ostream& operator<<(ostream& os, Posit32 a) {
    return os << a.getValueDouble(); 
} 

Posit32 sqrt(Posit32 a) {
    mpf_class p;
    a.get(p);
    m_sqrt(p, p);
    a.set(p);
    return a;
}

Posit32 sin(Posit32 a) {
    mpf_class m;
    a.get(m);
    m_sin(m, m, 20);
    a.set(m);
    return a;
}

Posit32 cos(Posit32 a) {
    mpf_class m;
    a.get(m);
    m_cos(m, m, 20);
    a.set(m);
    return a;
}

Posit32 abs(Posit32 a) {
    if (a < 0) {
        a.d = twosComp(a.d, a.size);
        return a;
    }
    else return a;
}

double Posit32::toDouble() {
    mpf_class a;
    this->get(a);
    return a.get_d();
}

