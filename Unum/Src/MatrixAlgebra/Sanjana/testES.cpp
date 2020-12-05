/**
 * testES.cpp
 * 
 * Used to test ideal es bits for 32 bit and 64 bit posit values.
 * ES is independent variable - along a numberline, goal is to find trend
 * 	between ES bits and the total number of bits required to represent
 * 	the decimal number
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
#define FRACBITS 24

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
				//start at 1, keep dividing until closest useed power is reached without surpassing
				while(decRep <= temp/USEED) {
					temp /= USEED;
					k--;
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
			int size = 23;
			shift = size;
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
				}
				temp = 1 << (shift+1);
				for(int i = 0; i < 32; i++)
					temp4[i] = temp4[i] || temp[i];

				for(int i = size; i > 0; i--)
					mant = mant + to_string(temp4[i]);
										
			} else {
				for(int i = size; i > 0; i--){
					mant = mant + "0";
				}
			}

			posString += mant;							//Append fraction to the string
			bitset<32> y;

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
			bitset<32> ifZero = 0;
			if (userInput == 0) { return ifZero;}
			// Get sign, k, expVal, frac
			decRepUsr = userInput;
			int s = calcSign();
			int k = calcK();
			int expVal = calcExp(k);
			// cout << "k " << k << "e " << expVal << endl;
			// cout << "closest 2 pwr " << pow(USEED, k)*pow(2, expVal) << endl;
			closestTwoPwr = pow(USEED, k)*pow(2, expVal);
			double frac = calcFrac();

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

	bitset<1> s;
	bitset<REGBITS> reg;
	bitset<ESBITS> expVal;
	bitset<FRACBITS> frac;

	private:

	public:
		void extract(bitset<32> x, int regLengthX) {
			//SIGN
			s[0] = x[32];

			//REGIME
			// preparing regime mask
			bitset<32> regMask = ~(~0 << REGBITS);
			regMask = regMask << (31-REGBITS);

			//extracting regime from posit
			for(int i = 0; i < REGBITS; i++) {
				reg[REGBITS-i] = x[31 - i] && regMask[31 - i];		//obtain regime bits for first number
			}
			
			int regCount = 5;
			int rs = 0;
			string regime = "";
			if(reg[REGBITS-1] == 1) {
				int regCount = 5;
				while(regCount >= 0 && reg[regCount] != 0) {
					regCount--;
					rs++;
					regime += "1";
				}
				regime += "0";
			} else if(reg[REGBITS-1] == 0) {
				while(regCount >= 0 && reg[regCount] != 1) {
					regCount--;
					rs++;
					regime += "0";
				}
				regime += "1";
			}
			rs++;
			cout << regime << endl;
			
			//EXPONENT
			// preparing exponent mask
			bitset<32> expMask = ~(~0 << ESBITS);
			expMask = expMask << (25-ESBITS);
			
			//and them
			for(int i = 0; i < ESBITS; i++) {
				expVal[ESBITS - 1 - i] = x[24-i] && expMask[24-i];		//obtain exponent bits for first number
			}

			//FRACTION
			// preparing fraction mask
			bitset<32> fracMask = ~(~1 << abs(FRACBITS+(REGBITS-rs)-3));

			//and them
			for(int i = 0; i < FRACBITS+(REGBITS-rs); i++) {
				frac[FRACBITS - 2 - i] = x[23-i] && fracMask[23-i];		//obtain fraction bits for first number
			}

			//cout << s << " " << regime << " " << expVal << " " << frac << endl;
		}

		//Getters
		bitset<1> 			getSign() {return s;}
		bitset<REGBITS> 	getReg()  {return reg;}
		bitset<ESBITS> 		getExp()  {return expVal;}
		bitset<FRACBITS> 	getFrac() {return frac;}
};

class testES{
	private:
	
	public:
	int control(double decRepUsr) {
		DecToPosit dtp;
		bitset<32> combine = dtp.control(decRepUsr);
		//cout << combine << endl;
		int j = 0;
		for(int i = 0; i < 32; i++){
			if(combine[32-i] == 1){
				j=i;
			}
		}
		cout << j << endl;
		
		return j;
	}
};

int main() {
	testES tes;
	int j;
	for(double i = 0.0000001; i < 1 ; i += 0.00001){
		cout << "Number: " << i << endl;
		j = tes.control(i);
		if(j == 31)
			cout << "----------------------------------------------------------31" << endl;
	}
}
