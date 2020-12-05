/**
 * Multiplication.cpp
 * Four Arithmetic Operations using Posit Numbers
 *
 * TBD: Multiplication (two's complement) for negative*positive numbers.
 * Also, in some cases, fraction*wholeNumber does not seem to be working
 * 	ex. 0.05 * 1234 doesn't work, but 0.0032*35 works
 *
 * Sanjana Shah
 * August 13, 2017
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
class DecToPosit {

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
		DecToPosit() {
			decRep = 0.0;
			decRepUsr = 0.0;
			regLength = 0;
			closestTwoPwr = 0.0;
		}

		double getClosest() {
			return closestTwoPwr;
		}

		/*
		 * Calculate the sign
		 */
		int calcSign() {
			//If positive, set value and return 0
			if(decRepUsr > 0) {
				decRep = decRepUsr;
				return 0;
			} else {
				//If negative, set value to absolute value and return 1
				decRep = abs((double)decRepUsr);
				return 1;
			}
		}

		/*
		 * Calculate k value
		 */
		int calcK() {
			int k = 0;
			if(decRep >= 1) {	//greater than or equal to 1
				int temp = 1;
				//start at 1, keep multiplying until closest useed power is reached without surpassing
				while(decRep >= temp*USEED) {
					temp *= USEED;
					k++;
				}
			} else {				//greater than 0, less than 1
				double temp = 1.0;
				k = -1;
				cout << decRep << " " << temp/USEED << endl;
				//start at 1, keep dividing until closest useed power is reached without surpassing
				while(decRep <= temp/USEED) {
					temp /= USEED;
					k--;
					cout << temp/USEED << endl;
				}
			}
			return k;
		}
		/*
		 * Calculate exponent value
		 */
		int calcExp(int k) {
			int expVal = 0;
			int multiple = 2;
			double temp = pow(USEED, k);
			//using useed^k, keep multiplying by 2 until closest 2 power is reached without surpassing
			while(decRep >= temp*multiple) {
				temp *= multiple;
				expVal++;
			}
			return expVal;
		}

		/*
		 * Calculate decimal (base 10 representation) of fraction
		 * by dividing input by closest two power
		 */
		double calcFrac() {
			double frac = (double)(decRep)/closestTwoPwr;
			return frac;
		}

		/*
		 * Converts regime, exponent, mantissa parts into Posit format.
		 * Joins the four parts into one 32-bit Posit number
		 */
		bitset<32> positFormat(short s, int k, int expVal, double frac) {
			//Keep adding to a string, then go character by character and turn it into a 32-bit bitset

			string posString = "";

			//SIGN
			if(s == 0) posString += "0";
			else if(s == 1) posString += "1";

			//REGIME
			int regLength = 0;
			int shift = 5;
			bitset<REGBITS> temp2;

			if(decRep >= 0 && decRep < 1) {		//If the value is less than 1
				for(int i = k; i < 0; i++) {
					temp2[shift] = 0;			//Add k number of zeros
					shift--;
					regLength++;
				}
				temp2[shift] = 1;				//Append terminating 1
				shift--;
				regLength++;
			} else {
				for(int i = 0; i < k+1; i++) {
					temp2[shift] = 1;			//Add k+1 number of ones
					regLength++;
					shift--;
				}
				temp2[shift] = 0;				//Append terminating 0
				shift--;
				regLength++;
			}
			for(int i = (REGBITS - 1); i >= 0; i--) {	//Append regime to the string
				if(temp2[i] == 0)
					posString += "0";
				else
					posString += "1";
			}

			//EXPONENT
			string str;
			bitset<ESBITS> exp (expVal);
			for(int i = (ESBITS-1); i >= 0; i--) {		//Append exponent to the string
				if(exp[i] == 0)
					posString += "0";
				else
					posString += "1";
			}

			//MANTISSA
			frac -= 1.0;
			string mant = "";
			shift = 23;
			bitset<23> temp4, temp;
			double divResult = frac;
			int rem;

			if(divResult != 0) {							//Convert decimal representation into binary bitset
				while(divResult != 1 && shift > 0) {
					divResult *= 2;
					if(divResult > 1) {
						divResult -= 1;
						rem = 1;
					} else {
						rem = 0;
					}
					temp4[shift] = rem;					//Append remainder to fraction bitset
					shift--;
					cout << temp4 << " " << divResult << " " << shift << endl;
				}
				temp = 1 << (shift+1);
				for(int i = 0; i < 32; i++)
					temp4[i] = temp4[i] || temp[i];

				for(int i = 23; i > 0; i--){
					mant = mant + to_string(temp4[i]);
				}
			} else {
				for(int i = 23; i > 0; i--){
					mant = mant + "0";
				}
			}

			posString += mant;							//Append fraction to the string

			cout << endl << posString << endl;
			bitset<32> y;

			cout << endl << endl << endl;

			//Convert string into a 32-bit bitset
			for(int i = 31; i >= 0; i--) {
				if(posString.at(i) == '0') {
					y[31-i] = 0;
				} else {
					y[31-i] = 1;
				}
			}

			return y;
		}

		/*
		 * Driver method for converting decimal (base 10 representation)
		 * number into Posit 32-bit bitset
		 */
		bitset<32> control(double userInput) {
			// Get sign, k, expVal, frac
			decRepUsr = userInput;
			int s = calcSign();
			int k = calcK();
			int expVal = calcExp(k);
			cout << "k " << k << "e " << expVal << endl;
			cout << "closest 2 pwr " << pow(USEED, k)*pow(2, expVal) << endl;
			closestTwoPwr = pow(USEED, k)*pow(2, expVal);
			double frac = calcFrac();
			cout << frac << endl;

			//Format in Posit notation
			return positFormat(s, k, expVal, frac);
		}

		//Getters
		int getRegLength() {return regLength;}
		int getSign() {return calcSign();}
		int getK() {return calcK();}
		int getExp() {return calcExp(calcK());}
		double getFrac() {return calcFrac();}
};



class Extraction {

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
		void extract(bitset<32> x, int regLengthX, bitset<32> y, int regLengthY) {
			//SIGN
			s1[0] = x[32];
			s2[0] = y[32];

			//REGIME
			// preparing regime mask
			bitset<32> regMask = ~(~0 << REGBITS);
			regMask = regMask << (31-REGBITS);

			//extracting regime from posit
			for(int i = 0; i < REGBITS; i++) {
				reg1[REGBITS - i] = x[31 - i] && regMask[31 - i];		//obtain regime bits for first number
			}
			for(int i = 0; i < REGBITS; i++) {
				reg2[REGBITS - i] = y[31 - i] && regMask[31 - i];		//obtain regime bits for second number
			}

			//EXPONENT
			// preparing exponent mask
			bitset<32> expMask = ~(~0 << ESBITS);
			expMask = expMask << (25-ESBITS);

			//and them
			for(int i = 0; i < ESBITS; i++) {
				expVal1[ESBITS - 1 - i] = x[24-i] && expMask[24-i];		//obtain exponent bits for first number
			}
			for(int i = 0; i < ESBITS; i++) {
				expVal2[ESBITS - 1 - i] = y[24-i] && expMask[24-i];		//obtain exponent bits for second number
			}

			//FRACTION
			// preparing fraction mask
			bitset<32> fracMask = ~(~0 << FRACBITS);

			//and them
			for(int i = 0; i < FRACBITS; i++) {
				frac1[FRACBITS - 2 - i] = x[23-i] && fracMask[23-i];		//obtain fraction bits for first number
			}
			for(int i = 0; i < FRACBITS; i++) {
				frac2[FRACBITS - 2 - i] = y[23-i] && fracMask[23-i];		//obtain fraction bits for second number
			}

			frac1[23] = 1;
			frac2[23] = 1;

			cout << s1 << " " << reg1 << " " << expVal1 << " " << frac1 << endl;
			cout << s2 << " " << reg2 << " " << expVal2 << " " << frac2 << endl;
			cout << endl;
		}

		//Getters
		bitset<1> 			getSign1() {return s1;}
		bitset<REGBITS> 		getReg1()  {return reg1;}
		bitset<ESBITS> 		getExp1()  {return expVal1;}
		bitset<FRACBITS> 	getFrac1() {return frac1;}
		bitset<1> 			getSign2() {return s2;}
		bitset<REGBITS> 		getReg2()  {return reg2;}
		bitset<ESBITS> 		getExp2()  {return expVal2;}
		bitset<FRACBITS>		getFrac2() {return frac2;}
};

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
		void setValues(bitset<1> n1,
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

		bitset<REGBITS> createReg(int k) {
			int regLength = 0;
			int shift = REGBITS-1;
			bitset<REGBITS> temp2;

			if(k < 0) {		//If the value is less than 1
				k = abs(k);
				for(int i = 0; i < k; i++) {
					temp2[shift] = 0;	//Add k number of zeros
					shift--;
					regLength++;
				}
				temp2[shift] = 1;		//Append terminating 1
				shift--;
				regLength++;
			} else {
				for(int i = 0; i < k+1; i++) {
					temp2[shift] = 1;		//Add k+1 number of ones
					regLength++;
					shift--;
				}
				temp2[shift] = 0;			//Append terminating 0
				shift--;
				regLength++;
			}
			return temp2;
		}

		/*
		 * Converts fraction to binary
		 */
		bitset<FRACBITS> mantToBin(double num) {
			string mant = "";
			int shift = FRACBITS-1;
			bitset<FRACBITS> temp4;
			int rem;

			if(num != 0) {
				while(num != 1.0 && shift > 0) {
					if(num > 1.0) {
						num -= 1.0;
						rem = 1;
					} else {
						rem = 0;
					}
					temp4[shift] = rem;
					shift--;
					num *= 2.0;
				}
				temp4[shift] = 1;
			}
			return temp4;
		}

};

int main() {

	bitset<1> s1;
	bitset<REGBITS> reg1;
	bitset<ESBITS> expVal1;
	bitset<FRACBITS> frac1;
	bitset<1> s2;
	bitset<REGBITS> reg2;
	bitset<ESBITS> expVal2;
	bitset<FRACBITS> frac2;

	// Get decimal value from user
	double decRepUsr, decRepUsr2;
	cout << "Enter 1st Decimal Value: ";
	cin >> decRepUsr;
	cout << endl;
	cout << "Enter 2nd Decimal Value: ";
	cin >> decRepUsr2;

	//Two objects for each of the posit values
	DecToPosit dtp;
	DecToPosit dtp2;

	int regLengthX = 0, regLengthY = 0;
	bitset<32> x, y;

	//Calls DecToPosit for each value, creating a 32-bit posit value
	x = dtp.control(decRepUsr);
	regLengthX = dtp.getRegLength();
	y = dtp2.control(decRepUsr2);
	regLengthY = dtp2.getRegLength();

	//extract values and store locally in main
	Extraction ext;
	ext.extract(x, regLengthX, y, regLengthY);

	s1 = ext.getSign1();
	reg1 = ext.getReg1();
	expVal1 = ext.getExp1();
	frac1 = ext.getFrac1();
	s2 = ext.getSign2();
	reg2 = ext.getReg2();
	expVal2 = ext.getExp2();
	frac2 = ext.getFrac2();

	Multiplication mult;
	mult.setValues(s1, reg1, expVal1, frac1, s2, reg2, expVal2, frac2);

	bitset<32> temp1;
	bitset<32> temp2;
	bitset<FRACBITS> frac;

	for(int i = 32; i > 0; i--) {
		temp1[i] = frac1[i];
		temp2[i] = frac2[i];
	}

	cout << endl << endl;

	//SIGN
		bitset<1> sign = (s1 | s2);
		cout << "Or for sign " << sign << endl;

	//REGIME
		bitset<REGBITS> regime = mult.createReg(dtp.getK() + dtp2.getK());
		cout << "Add k's for regime " << regime << endl;

	//EXPONENT
		bitset<ESBITS+1> exponent = expVal1.to_ulong() + expVal2.to_ulong();
		cout << "Add expVals for exponent " << exponent << endl;

		int exp = exponent.to_ulong();
		if(exp > (pow(2, ESBITS)-1)) {
			cout << "1. Exp is greater than representable range" << endl;
			regime = mult.createReg(dtp.getK() + dtp2.getK() + 1);
			exponent = exponent.to_ulong() - USEED;
		}

		cout << endl;

	//FRACTION

		//create decimal representations for the two fractions
		//multiply them - can do later with adder, full adder, etc. circuits
		double sum = 0, sum2 = 0;
		for(int i = 0; i < 23; i++) {
			sum += frac1[i] * (pow(2, -23+i));
		}
		sum += 1;

		for(int i = 0; i < 23; i++) {
			sum2 += frac2[i] * (pow(2, -23+i));
		}
		sum2 += 1;
		sum = sum * sum2;

		//Implicit number is 2. Divide fraction by 2 and shift to exponent,
		//if exponent carries out, shift to regime.
		if(sum >= 2) {
			cout << "Frac Greater than 2, Shift to Exp" << endl;
			sum /= 2;
			exp += 1;
			exponent = exp;
		}
		bitset<FRACBITS> fracTemp = mult.mantToBin(sum);
		bitset<FRACBITS-2> fracFinal;
		for(int i = 0; i < 23; i++) {
			fracFinal[i] = fracTemp[i+1];
		}

		//if exponent is not representable
		if(exp > (pow(2, ESBITS)-1)) {
			cout << "2. Exp is greater than representable range" << endl;
			regime = mult.createReg(dtp.getK() + dtp2.getK() + 1);
			exponent = exp - pow(2, ESBITS);
		}


		//PRINT

		cout << endl << endl;

		bitset<ESBITS> expFinal;

		for(int i = 0; i < ESBITS; i++) {
			expFinal[i] = exponent[i];
		}

		cout << sign << " ";

		if(regime[REGBITS-1] == 1) {
			int regCount = 5;
			while(regCount >= 0 && regime[regCount] != 0) {
				cout << regime[regCount];
				regCount--;
			}
			cout << "0";
		} else if(regime[REGBITS-1] == 0) {
			int regCount = 5;
			while(regCount >= 0 && regime[regCount] != 1) {
				cout << regime[regCount];
				regCount--;
			}
			cout << "1";
		}

		cout << " " << expFinal << " " << fracFinal << endl;
}


