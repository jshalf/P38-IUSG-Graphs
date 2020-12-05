#ifndef __BFLOAT16_HH_
#define __BFLOAT16_HH_

#include <string>
using namespace std;
#define max_bfloat16 0x7F7F0000
#define min_bfloat16 0x00010000
#define prevent_overflow 0

class bfloat16 {
public:
    float f; //bfloat16 object takes on same value as f, which is limited to 16 msb. 

    bfloat16(int n      );
    bfloat16(float f    );
    bfloat16(double d   );
    bfloat16(mpf_class m);
    bfloat16(string s   );
    bfloat16(           );

    //Round the float 
    float round(float f, bool sticky=0);

    friend bfloat16 operator+ (bfloat16 a,bfloat16 b); 
    friend bfloat16 operator- (bfloat16 a,bfloat16 b);
    friend bfloat16 operator- (bfloat16 a           );
    friend bfloat16 operator* (bfloat16 a,bfloat16 b);
    friend bfloat16 operator/ (bfloat16 a,bfloat16 b);
    
    friend bool operator<  (bfloat16 a, bfloat16 b);
    friend bool operator>  (bfloat16 a, bfloat16 b);
    friend bool operator<= (bfloat16 a, bfloat16 b);
    friend bool operator>= (bfloat16 a, bfloat16 b);
    friend bool operator== (bfloat16 a, bfloat16 b);
    friend bool operator!= (bfloat16 a, bfloat16 b);
    
    bfloat16& operator+= (const bfloat16& b);
    bfloat16& operator-= (const bfloat16& a);
    bfloat16& operator*= (const bfloat16& a);
    bfloat16& operator/= (const bfloat16& a);
    bfloat16& operator=  (const bfloat16& a);

    float getMin();
    float getMax();
    
    friend ostream& operator<< (ostream& os, bfloat16 a);

    friend bfloat16 sqrt(bfloat16 a);
    friend bfloat16 sin (bfloat16 a);
    friend bfloat16 cos (bfloat16 a);
    friend bfloat16 abs (bfloat16 a);
};

#endif