#ifndef __Addition_hh__
#define __Addition_hh__
/**
 * Addition.hh
 * 
 * Used by Matrix.cc to add either two decimal values (by converting them
 *	into two posit values and adding parts together), or two posit values
 * 	(by extracting the parts and adding them together).
 * 
 * Since multiplication outputs values in posit format, which are then added
 * together in matrix multiplication, the latter will be more useful for matrix
 * algebra.
 * 	
 * The result of the addition is returned.
 *
 * Sanjana Shah
 * August 13, 2018
 */

#include <iostream>
#include <cmath>
#include <string>
#include <bitset>

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
class AddDecToPosit {

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
    AddDecToPosit();

    double getClosest();

		/*
		 * Calculate the sign
		 */
    int calcSign();

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
		inline int getSign() {return calcSign();}
		inline int getK() {return calcK();}
		inline int getExp() {return calcExp(calcK());}
		inline double getFrac() {return calcFrac();}
};

/*
 * This class extracts sign, regime, exponent, and mantissa from two posit values
 */
class AddExtraction {

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
    void extract(bitset<32> x, bitset<32> y);

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


class Operations {
	short c;

	public:
		inline short getC() {return c;}

		//sums two fraction values, could replace with Adder circuit in hardware
    bitset<32> sumTwoBin(bitset<32> num1, bitset<32> num2, bitset<32> ctemp);
};


class Addition {

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

		/*
		 * Compares regimes of two Posit values and divides one of the two fractions respectively
		 */
    bitset<FRACBITS> compareReg(AddDecToPosit dtp, AddDecToPosit dtp2);

		/*
		 * Compares exponents of two Posit values and divides one of the two fractions respectively
		 */
    bitset<FRACBITS> compareExpVal(bitset<FRACBITS> frac1, bitset<FRACBITS> frac2);
};

/*
 * Organized class used to call methods and run through entire addition given two
 * decimal values in IEEE notation.
 */

class AddControl{
	
private:

public:

    void add();
};

/*
 * Organized class used to call methods and run through entire addition given two
 * decimal values in IEEE notation.
 */
class GivenTwoPosit {
	private:
	
		bitset<1> s1;
		bitset<REGBITS> reg1;
		bitset<ESBITS> expVal1;
		bitset<FRACBITS> frac1;
		bitset<1> s2;
		bitset<REGBITS> reg2;
		bitset<ESBITS> expVal2;
		bitset<FRACBITS> frac2;
	
	public:
	
	GivenTwoPosit(){}

    int calcKPosit(bitset<REGBITS> reg);
    bitset<FRACBITS> compareRegPosit(int k1, int k2, bitset<32> x, bitset<32> y);
    bitset<FRACBITS> compareExpValPosit();
    bitset<32> addPosits(bitset<32> x, bitset<32> y);
};

class AddPositToDec {
    
	bitset<1> posSign;
	bitset<REGBITS> posRegime;
	bitset<ESBITS> posExp;
	bitset<FRACBITS> posFrac;
	
	private:
	
	public:
	
	AddPositToDec(){}
	
    void posExtract(bitset<32> positVal);
	
    int control(bitset<32> positVal);
};
/*

int main() {
	double decRepUsr, decRepUsr2;
	cout << "\n-----------Addition-----------\nEnter 1st Decimal Value: ";
	cin >> decRepUsr;
	cout << "Enter 2nd Decimal Value: ";
	cin >> decRepUsr2;
		
	AddDecToPosit dtp, dtp2;
	bitset<32> positOne = dtp.control(decRepUsr);
	bitset<32> positTwo = dtp2.control(decRepUsr2);
	
	cout << "One: " << positOne << endl;
	cout << "Two: " << positTwo << endl;
	
	GivenTwoPosit gtp;
	AddPositToDec ptd;
	cout << ptd.control(gtp.addPosits(positOne, positTwo));
}
*/
#endif

