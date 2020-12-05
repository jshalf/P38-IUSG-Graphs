#ifndef __Multiplication_hh__
#define __Multiplication_hh__

/**
 * Multiplication.cpp
 * 
 * Used by Matrix.cc to convert two decimal values into two posit values
 * 	and multiply them together then return the result of the multiplication.
 *
 * Sanjana Shah
 * August 13, 2018
 */


#include <iostream>
#include <cmath>
#include <string>
#include <bitset>
#include <sstream>

using namespace std;

//In temporary registers, 24th and 25th bit hold the carry-out, and the implicit 1
//Extra two bits are used to check for implicit 2 case or implicit 0 case
#define FRACBITS 25

//Fixed 6 bits for regime
#define REGBITS 6

//Fixed 2 bits for exponent value
#define ESBITS 2

#define USEED (pow(2, pow(2, ESBITS)))

/*
 * This class is used to convert a decimal (base 10 representation) value into
 * a 32-bit Posit number.
 */
class MultDecToPosit {

	private:
		double closestTwoPwr;

	public:

		//Absolute value of input
		double decRep;

		//Raw input
		double decRepUsr;

		//Actual length needed to represent regime
		int regLength;

		/*
		 * Constructor
		 */
		MultDecToPosit() {
			decRep = 0.0;
			decRepUsr = 0.0;
			regLength = 0;
			closestTwoPwr = 0.0;
		}

		inline double getClosest() {
			return closestTwoPwr;
		}

		/*
		 * Calculate the sign
		 */
    int calcSign(int decRepUsr);

		/*
		 * Calculate k value
		 */
    int calcK();
		/*
		 * Calculate exponent value
		 */
    int calcExp(int k);

		/*
		 * Calculate decimal (base 10 representation) of fraction
		 * by dividing input by closest two power
		 */
    double calcFrac();

		/*
		 * Converts regime, exponent, mantissa parts into Posit format.
		 * Joins the four parts into one 32-bit Posit number
		 */
    bitset<32> positFormat(short s, int k, int expVal, double frac);

		/*
		 * Driver method for converting decimal (base 10 representation)
		 * number into Posit 32-bit bitset
		 */
    bitset<32> control(double userInput);

		//Getters
		inline int getRegLength() {return regLength;}
		inline int getSign() {return calcSign(decRepUsr);}
		inline int getK() {return calcK();}
		inline int getExp() {return calcExp(calcK());}
		inline double getFrac() {return calcFrac();}
};



class MultExtraction {

	bitset<1> s1;
	bitset<REGBITS> reg1;
	bitset<ESBITS> expVal1;
	bitset<FRACBITS> frac1;
	bitset<1> s2;
	bitset<REGBITS> reg2;
	bitset<ESBITS> expVal2;
	bitset<FRACBITS> frac2;

	private:

	public:
    void extract(bitset<32> x, int regLengthX, bitset<32> y, int regLengthY);

		//Getters
		inline bitset<1> 			getSign1() {return s1;}
		inline bitset<REGBITS> 		getReg1()  {return reg1;}
		inline bitset<ESBITS> 		getExp1()  {return expVal1;}
		inline bitset<FRACBITS> 	getFrac1() {return frac1;}
		inline bitset<1> 			getSign2() {return s2;}
		inline bitset<REGBITS> 		getReg2()  {return reg2;}
		inline bitset<ESBITS> 		getExp2()  {return expVal2;}
		inline bitset<FRACBITS>		getFrac2() {return frac2;}
};

/*
 * Multiplies the separate parts of the given two IEEE decimal values (sign, regime, exponent, mantissa)
 */
 
class Multiplication {

	bitset<1> s;
	bitset<REGBITS> reg;
	bitset<ESBITS> expVal;
	bitset<FRACBITS> frac;

	bitset<1> s1;
	bitset<REGBITS> reg1;
	bitset<ESBITS> expVal1;
	bitset<FRACBITS> frac1;
	bitset<1> s2;
	bitset<REGBITS> reg2;
	bitset<ESBITS> expVal2;
	bitset<FRACBITS> frac2;

	public:
		inline void setValues(bitset<1> n1,
				bitset<REGBITS> r1,
				bitset<ESBITS> e1,
				bitset<FRACBITS> f1,
				bitset<1> n2,
				bitset<REGBITS> r2,
				bitset<ESBITS> e2,
				bitset<FRACBITS> f2) {
			s1 = n1;
			reg1 = r1;
			expVal1 = e1;
			frac1 = f1;
			s2 = n2;
			reg2 = r2;
			expVal2 = e2;
			frac2 = f2;
		}

    bitset<REGBITS> createReg(int k);

		/*
		 * Converts fraction to binary
		 */
    bitset<FRACBITS> mantToBin(double num);

};

/*
 * Organized class used to call methods and run through entire multiplication given two
 * decimal values in IEEE notation.
 */
 
class MultiplyControl{
	
private:

public:

    bitset<32> mult(double decRepUsr, double decRepUsr2);
};

class MultPositToDec {
	
	bitset<1> posSign;
	bitset<REGBITS> posRegime;
	bitset<ESBITS> posExp;
	bitset<FRACBITS> posFrac;
	
	private:
	
	public:
	
	MultPositToDec(){}
	
    int calcKPosit(bitset<REGBITS> reg);
	
    void posExtract(bitset<32> positVal);
	
    int control(bitset<32> positVal);
	
};

/*
int main() {
	MultiplyControl m;
	m.mult(3, -3);
}	*/

#endif

