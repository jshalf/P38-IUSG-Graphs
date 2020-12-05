
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <gmpxx.h>
#include <iostream>
#include "PositBase.hh"
#include "Posit32.hh"
#include "Posit32nc.hh"
/* 
int main(int argc, char*argv[]){
    Posit32 p(2,6); // time to start doing experiments with smaller numbers
    mpf_class a,b,c;
    
    a=1.00;
    b=-1.00;
    
    p.set(a); // assign positiveMPF to the posit
    if(p.getSignBit()==1) printf("Postive has negative sign bit incorrect\n");
    else printf("correct\n");
    
    exit(0);
    
    {
        unsigned int i,tc;
        size_t maxbits = 6;
        for(i=0;i<maxbits;i++){
            printf("Number:");
            printBin(i);
            tc=twosComp(i,maxbits);
            printf("TwosComp");
            printBin(tc);
        }
    }
    // return 1;
    
    a = 1234.5;
    b = "-5678.5888";
    c = a+b;
    std::cout << "a=" << a << " b=" << b << " c=" << c << "\n";
     std::cout << "sum is " << c << "\n";
     std::cout << "absolute value is " << abs(c) << "\n";
    p.set(0b00101);
    p.get(a); // overwrite a
    std::cout << "GMP value out is " << a << "\n";
    // if(p.isZero()) puts("Zero");
    //if(p.isInfinity()) puts("Infinity");
    //p.set(0xFFFFFFFF);
    p.printInfo(); putchar('\n');
    printf("Double Precision Value Conversion=%5.5lf\n",p.getValueDouble());

    { // run through the full range of posit bit representations
        unsigned int i,maxbits=(1<<(p.size));
        // maxbits>>=1; // half of maxbits (just do the positive numbers first
        printf("MaxSize=%u\n",maxbits);
        for(i=0;i<maxbits;i++){
            printf("TestLoop i=%u\n",i);
            puts("=========================================");
            printf("Encode a Raw Binary into the POSIT and then interpret it as a decimal FP number\n");
            p.set(i);
            p.printInfoCompact();
            // printf("Undo Twos Comp="); printBin(p.undoTwosComp()); puts("\n");
            printf("\tNow I get an infinite precision value from the posit");
            p.get(a);
            std::cout << "\tMPEvalue=" << a << "\n";
            
            // need to fix SET again to handle twos comp
            printf("Now I will SET the posit from that infinite precision value and see if it is the same \n");
            p.set(a); // now set it and see if we get the same result out
            p.printInfoCompact();
            //p.get(a);
            //std::cout << "\t Re-Evaluate=" << a << "\n\n";
            puts("=========================================\n\n\n");
        }
        
        // Now we do arithmetic in all possible combinations
        {
            int i,j;
            Posit32gmp a(p),b(p),r(p);
            mpf_class ma,mb,mr;
            float fa,fb,fr;
            double da,db,dr;
            unsigned int maxbits=(1<<(p.size));
            
            for(i=0;i<maxbits;i++){
                a.set(i);
                a.get(ma);
                da=a.getValueDouble();
                fa=(float)da;
                if(a.isInfinity()) break;
                for(j=0;j<maxbits;j++){
                    b.set(j);
                    b.get(mb); // pull out the exact rep (will compare exact rep to posit rep
                        // will do the same for double and single precision
                    if(b.isInfinity()) break;
                    db=b.getValueDouble();
                    fb=(float)db;
                    r=a+b; // do it in exact posit
                    mr= ma + mb; // do for GMP precision
                    dr = da + db; // double
                    fr = fa + fb; // single
                    // printf("da+db==>DR  %le+%le==>%le\n",da,db,dr);
                    {
                        // now we promote posits to infinite precision
                        mpf_class a1,b1,r1,mdr,mfr,diff,ddiff,fdiff;
                        a.get(a1); b.get(b1); r.get(r1);
                        mdr = dr; // convert double to infinite precision to compute error relative to exact
                        mfr=fr; // convert single to infinite precision to compute error relative to exact
                        m_diff(diff,r1,mr); // error for posit precision
                        m_diff(ddiff,mdr,mr); // error for double precision
                        m_diff(fdiff,mfr,mr); // error for single precision
                        // now need to add interval
                        // Output is Posit32 a,b,r < r=a+b >  and GMP ma,mb,mr  < mr=ma+b > and difference
                        std::cout << a1 << ',' << b1 << ',' << r1 << '\t' << ma << ',' << mb << ',' << mr << '\t' << diff // << '\n';
                           << '\t' << mdr << ',' << ddiff << '\t' << mfr << ',' << fdiff << '\n';
                    }
                }
            }
        }
    }
    /*
        have encode and decode functions as members
        encodePostFromFloat32, encodePositFromMPF()
        and complementary decode functions
        decodePositToFloat32 decodePositToFloat64(), decodePositToMPF().
     
        don't forget to do descriptions of PathForward work
     */
    
//    return 1;
//}

