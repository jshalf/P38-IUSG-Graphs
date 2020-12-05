/**
 * Addition.cpp
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
		AddDecToPosit() {
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
			
			if(userInput == 0) {
				return positFormat(0,0,0,0);
			}
			
			// Get sign, k, expVal, frac
			decRepUsr = userInput;
			int s = calcSign();
			int k = calcK();
			int expVal = calcExp(k);
			//cout << "k " << k << "e " << expVal << endl;
			//cout << "closest 2 pwr " << pow(USEED, k)*pow(2, expVal) << endl;
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
		void extract(bitset<32> x, bitset<32> y) {
			//SIGN
			s1[0] = x[31];
			s2[0] = y[31];

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
			bitset<32> fracMask = ~(~1 << abs(FRACBITS-3));

			//and them
			for(int i = 0; i < FRACBITS; i++) {
				frac1[FRACBITS - 2 - i] = x[23-i] && fracMask[23-i];		//obtain fraction bits for first number
			}
			for(int i = 0; i < FRACBITS; i++) {
				frac2[FRACBITS - 2 - i] = y[23-i] && fracMask[23-i];		//obtain fraction bits for second number
			}

			frac1[23] = 1;
			frac2[23] = 1;

			//cout << s1 << " " << reg1 << " " << expVal1 << " " << frac1 << endl;
			//cout << s2 << " " << reg2 << " " << expVal2 << " " << frac2 << endl;
			//cout << endl;
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


class Operations {
	short c;

	public:
		short getC() {return c;}

		//sums two fraction values, could replace with Adder circuit in hardware
		bitset<32> sumTwoBin(bitset<32> num1, bitset<32> num2, bitset<32> ctemp) {
			bitset<32> sum1;

			for(int i = 0; i < 32; i++){
			   sum1[i] = ((num1[i] ^ num2[i]) ^ c); // c is carry
			   c = ((num1[i] & num2[i]) | (num1[i] & c)) | (num2[i] & c);
			}

			if((num1[0] == 0 && num2[0] == 0) || (num1[0] == 1 && num2[0] == 1))
				sum1[0] = 0;
			c = sum1[FRACBITS];

			return sum1;
		}
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

			if(s == 1) {			//If the value is negative
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
					temp4[shift] = rem;			//Append remainder to bitset
					shift--;
					num *= 2.0;
				}
				temp4[shift] = 1;
			}
			return temp4;
		}

		/*
		 * Compares regimes of two Posit values and divides one of the two fractions respectively
		 */
		bitset<FRACBITS> compareReg(AddDecToPosit dtp, AddDecToPosit dtp2) {
			cout << "Comparing Regimes" << endl;
			long r1 = dtp.getK();
			long r2 = dtp2.getK();

			//decimal (base 10 representations) of the two fractions
			//TBD: Future work is to convert to bitset (binary format) then shift
			//	   Here I am just converting to int, dividing, then converting to bitset
			long f1 = frac1.to_ulong();
			long f2 = frac2.to_ulong();
			int count = 0;

			//Always shift to higher
				//If first regime is greater than the second, make second regime equal to first by
				//multiplying second regime by useed (adding to k). Then divide second fraction by useed
			if(r1 > r2) {
				while(r2 != r1) {
					cout << "reg1 > reg2" << endl;
					r2 += 1;
					cout << "Dividing frac2 by useed" << endl;
					count++;
					f2 = f2 >> int(log(USEED)/log(2));
				}
				frac2 = bitset<FRACBITS>(f2);
				cout << "frac2 " << frac2 << endl;

				double x = 1.0/(pow(USEED, count));
				bitset<FRACBITS> y = mantToBin(x);
				y = y >> 2;
				for(int i = 0; i <= FRACBITS; i++){
					frac1[i] = frac1[i] || y[i];
				}

				reg2 = createReg(r2);
				return frac2;
			} else if(r2 > r1) {
				//If second regime is greater than the first, make first regime equal to second by
				//multiplying first regime by useed (adding to k). Then divide first fraction by useed
				while(r1 != r2) {
					cout << "reg2 > reg1" << endl;
					r1 += 1;
					cout << "Dividing frac1 by useed" << endl;
					count++;
					f1 = f1 >> int(log(USEED)/log(2));
				}
				frac1 = bitset<FRACBITS>(f1);

				double x = 1.0/(pow(USEED, count));
				bitset<FRACBITS> y = mantToBin(x);
				y = y >> 2;
				for(int i = 0; i <= FRACBITS; i++){
					frac1[i] = frac1[i] || y[i];
				}

				reg1 = createReg(r1);
				return frac1;
			}
			return frac1;
		}

		/*
		 * Compares exponents of two Posit values and divides one of the two fractions respectively
		 */
		bitset<FRACBITS> compareExpVal(bitset<FRACBITS> frac1, bitset<FRACBITS> frac2) {
			cout << endl << "Comparing Exponents"<< endl;

			//TBD: Instead of converting to integer and comparing, natively compare withing bitsets
			long exp1 = expVal1.to_ulong();
			long exp2 = expVal2.to_ulong();

			int count = 0;

			if(exp1 > exp2) {
				cout << "exp1 > exp2" << endl;
				//TBD: Instead of iteration, future work is to subtract smaller exponent value from larger
				while(exp2 != exp1) {
					count++;
					exp2 += 1;
					//cout << "Dividing frac2 by 2" << endl;
					frac2 = frac2 >> 1;
					count++;
				}
				return frac2;
			} else if(exp2 > exp1) {
				cout << "exp2 > exp1" << endl;
				while(exp1 != exp2) {
					exp1 += 1;
					//cout << "Dividing frac1 by 2" << endl;
					frac1 = frac1 >> 1;
					count++;
				}
				return frac1;
			}
			//The two exponents are now the same, default return can be any of the two fractions
			return frac1;
		}
};

/*
 * Organized class used to call methods and run through entire addition given two
 * decimal values in IEEE notation.
 */

class AddControl{
	
private:

public:

	void add() {
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
		cout << "Addition\nEnter 1st Decimal Value: ";
		cin >> decRepUsr;
		cout << endl;
		cout << "Enter 2nd Decimal Value: ";
		cin >> decRepUsr2;

		//Two objects for each of the posit values
		AddDecToPosit dtp;
		AddDecToPosit dtp2;

		if(decRepUsr2 == 0) {
			cout << dtp.control(decRepUsr) << endl;
		} else if(decRepUsr == 0) {
			cout << dtp.control(decRepUsr2) << endl;
		} else {

			int regLengthX = 0, regLengthY = 0;
			bitset<32> x, y;

			//Calls DecToPosit for each value, creating a 32-bit posit value
			x = dtp.control(decRepUsr);
			regLengthX = dtp.getRegLength();

			y = dtp2.control(decRepUsr2);
			regLengthY = dtp2.getRegLength();

			//extract values and store locally in main
			AddExtraction ext;
			ext.extract(x, y);

			s1 = ext.getSign1();
			reg1 = ext.getReg1();
			expVal1 = ext.getExp1();
			frac1 = ext.getFrac1();
			s2 = ext.getSign2();
			reg2 = ext.getReg2();
			expVal2 = ext.getExp2();
			frac2 = ext.getFrac2();

			Addition add;
			add.setValues(s1, reg1, expVal1, frac1, s2, reg2, expVal2, frac2);

			//Make reg same, divide frac
			if(dtp.getK() > dtp2.getK()) {
				reg2 = reg1;
				frac2 = add.compareReg(dtp, dtp2);
			} else if(dtp.getK() < dtp2.getK()){
				reg1 = reg2;
				frac1 = add.compareReg(dtp, dtp2);
			}

			//Make exp same, divide frac
			if(expVal1.to_ulong() > expVal2.to_ulong()) {
				expVal2 = expVal1;
				frac2 = add.compareExpVal(frac1, frac2);
			} else if(expVal1.to_ulong() < expVal2.to_ulong()){
				expVal1 = expVal2;
				frac1 = add.compareExpVal(frac1, frac2);
			}

			/*
			cout << endl << "Regimes and ExpVals are same" << endl;
			cout << s1 << " " << reg1 << " " << expVal1 << " " << frac1 << endl;
			cout << s2 << " " << reg2 << " " << expVal2 << " " << frac2 << endl;
			cout << endl;
			*/

			//exponent and regime values are now same for both Posit numbers

			bitset<32> temp1;
			bitset<32> temp2;
			bitset<32> ctemp = 0;
			bitset<FRACBITS> frac;

			for(int i = 32; i > 0; i--) {
				temp1[i] = frac1[i];
				temp2[i] = frac2[i];
			}

			Operations op;


			//SIGN
			bitset<1> sign;
			if(s1 == s2){sign = s1;}
			else if(decRepUsr > decRepUsr2){sign = s1;}
			else {sign = s2;}

			//REGIME
			bitset<REGBITS> regime = reg1;

			//EXPONENT
			bitset<ESBITS> exponent = expVal1;

			//FRACTION
			bitset<32> sumFrac = op.sumTwoBin(temp1, temp2, ctemp);

			// preparing fraction mask
			bitset<32> fracMask = ~(~0 << FRACBITS);

			// extracting fraction bits for sum
			for(int i = 0; i < FRACBITS; i++) {
				frac[FRACBITS - i] = sumFrac[FRACBITS-i] && fracMask[FRACBITS-i];		//keep &&-ing to get last ESBITS values
			}

			if(frac[FRACBITS-1] == 1 && frac[FRACBITS-2] == 0) {
				//Implicit number is 2. Divide fraction by 2 and shift to exponent,
				//if exponent carries out, shift to regime.
				int temp = expVal1.to_ulong() +1;
				if(temp > (pow(2, ESBITS)-1)) {
					regime = add.createReg(dtp.getK() + 1);
					exponent = bitset<ESBITS>(0);
				} else {
					exponent = bitset<ESBITS>(temp);
				}
				frac = frac >> 1;
			} else if(frac[FRACBITS-1] == 0 && frac[FRACBITS-2] == 0) {
				//Implicit number is 0. Multiply fraction until implicit number becomes 1.
				//If exponent has enough bits to pull from, do so. Else, pull from regime
				int count = 0;
				while(frac.to_ulong() != 0 && frac[FRACBITS-2] != 1) {
					frac = frac << 1;
					count++;
				}
				if(expVal1.to_ulong() >= count) {	//Since both expVals are same, we can compare either to count
					//If subtracting count from expVal does not result in a negative number
					exponent = bitset<ESBITS>(expVal1.to_ulong() - count);
				} else if (expVal1.to_ulong() < count) {
					//Else, subtract from regime
					regime = add.createReg(dtp.getK() - 1);
					exponent = bitset<ESBITS>(exponent.to_ulong() + int(log(USEED)/log(2)));
				}
			}

			bitset<FRACBITS-2> fracFinal;

			for(int i = 0; i < 23; i++)
				fracFinal[i] = frac[i];

			fracFinal[FRACBITS-2] = 0;

			cout << "Addition: " << sign << " " << regime << " " << exponent << " " << fracFinal << endl;
		}
	}
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

	int calcKPosit(bitset<REGBITS> reg) {
		int negative = 0;
		int positive = 0;
		int k = 0;
		int regCount = 5;
		
		if(reg[REGBITS-1] == 1) {
			while(regCount >= 0 && reg[regCount] != 0) {
				regCount--;
				positive++;
			}
			k = positive-1;
		} else if(reg[REGBITS-1] == 0) {
			while(regCount >= 0 && reg[regCount] != 1) {
				regCount--;
				negative++;
			}
			k = 0-negative;
		}
		
		return k;
	}
	
	bitset<FRACBITS> compareRegPosit(int k1, int k2, bitset<32> x, bitset<32> y){
				
		//cout << s1 << " " << reg1 << " " << expVal1 << " " << frac1 << endl;
		//cout << s2 << " " << reg2 << " " << expVal2 << " " << frac2 << endl;
		
		long r1 = calcKPosit(reg1);
		long r2 = calcKPosit(reg2);

		//decimal (base 10 representations) of the two fractions
		//TBD: Future work is to convert to bitset (binary format) then shift
		//	   Here I am just converting to int, dividing, then converting to bitset
		//long f1 = frac1.to_ulong();
		//long f2 = frac2.to_ulong();
		int count = 0;

		//Always shift to higher
			//If first regime is greater than the second, make second regime equal to first by
			//multiplying second regime by useed (adding to k). Then divide second fraction by useed
			
		if(r1 > r2) {
			//cout << "r1 > r2" << endl;
			while(r2 != r1) {
				r2 += 1;
				count++;
				frac2 = frac2 >> int(log(USEED)/log(2));
			}
			//frac2 = bitset<FRACBITS>(f2);

			double x = 1.0/(pow(USEED, count));
			Addition a;
			bitset<FRACBITS> y = a.mantToBin(x);
			y = y >> 2;
			for(int i = 0; i <= FRACBITS; i++){
				frac1[i] = frac1[i] || y[i];
			}

			return frac2;
		} else if(r2 > r1) {
			//cout << "r2 > r1" << endl;
			//If second regime is greater than the first, make first regime equal to second by
			//multiplying first regime by useed (adding to k). Then divide first fraction by useed
			while(r1 != r2) {
				r1 += 1;
				count++;
				frac1 = frac1 >> int(log(USEED)/log(2));
			}
			//frac1 = bitset<FRACBITS>(f1);
			double x = 1.0/(pow(USEED, count));
			Addition a;
			bitset<FRACBITS> y = a.mantToBin(x);
			y = y >> 2;
			for(int i = 0; i <= FRACBITS; i++){
				frac1[i] = frac1[i] || y[i];
			}
			return frac1;
		}
		return frac1;
	}
	
	bitset<FRACBITS> compareExpValPosit() {
		//TBD: Instead of converting to integer and comparing, natively compare within bitsets
		
		//cout << s1 << " " << reg1 << " " << expVal1 << " " << frac1 << endl;
		//cout << s2 << " " << reg2 << " " << expVal2 << " " << frac2 << endl;
						
		long exp1 = expVal1.to_ulong();
		long exp2 = expVal2.to_ulong();

		int count = 0;

		if(exp1 > exp2) {
			//cout << "exp1 > exp2" << endl;
			//TBD: Instead of iteration, future work is to subtract smaller exponent value from larger
			while(exp2 != exp1) {
				count++;
				exp2 += 1;
				frac2 = frac2 >> 1;
				count++;
			}
			return frac2;
		} else if(exp2 > exp1) {
			//cout << "exp2 > exp1" << endl;
			while(exp1 != exp2) {
				exp1 += 1;
				frac1 = frac1 >> 1;
				count++;
			}
			return frac1;
		}
		//The two exponents are now the same, default return can be any of the two fractions
		return frac1;
	}
	
	
	bitset<32> addPosits(bitset<32> x, bitset<32> y) {
		AddExtraction ext;
		ext.extract(x, y);

		s1 = ext.getSign1();
		reg1 = ext.getReg1();
		expVal1 = ext.getExp1();
		frac1 = ext.getFrac1();
		s2 = ext.getSign2();
		reg2 = ext.getReg2();
		expVal2 = ext.getExp2();
		frac2 = ext.getFrac2();
		
		Addition add;
		add.setValues(s1, reg1, expVal1, frac1, s2, reg2, expVal2, frac2);
		
		int k1 = calcKPosit(reg1);
		int k2 = calcKPosit(reg2);

		//Make reg same, divide frac
		if(k1 > k2) {
			frac2 = compareRegPosit(k1, k2, x, y);
			reg2 = reg1;
		} else if(k1 < k2){
			frac1 = compareRegPosit(k1, k2, x, y);
			reg1 = reg2;
		}

		//Make exp same, divide frac
		if(expVal1.to_ulong() > expVal2.to_ulong()) {
			frac2 = compareExpValPosit();
			expVal2 = expVal1;
		} else if(expVal1.to_ulong() < expVal2.to_ulong()){
			frac1 = compareExpValPosit();
			expVal1 = expVal2;
		}

		//exponent and regime values are now same for both Posit numbers

		bitset<32> temp1;
		bitset<32> temp2;
		bitset<32> ctemp = 0;
		bitset<FRACBITS> frac;

		for(int i = 32; i > 0; i--) {
			temp1[i] = frac1[i];
			temp2[i] = frac2[i];
		}

		Operations op;
		
		string addRsltStr = "";

		//SIGN
		bitset<1> sign;
		if(s1 == s2){sign = s1;}
		else if((frac1 == frac2) && (s1 != s2)) {sign = 0;}
		else{
			int counter1 = 300, counter2 = 300;
			for(int i = 0; i < FRACBITS; i++){
				if(frac1[i] == 1)
					counter1 = i;
				if(frac2[i] == 1)
					counter2 = i;
				if(counter1 < counter2){
					sign = s1;
					break;
				} else if(counter2 < counter1){
					sign = s2;
					break;
				}
			}
		}
		
		addRsltStr += sign.to_string();

		//REGIME
		bitset<REGBITS> regime = reg1;
		
		addRsltStr += regime.to_string();

		//EXPONENT
		bitset<ESBITS> exponent = expVal1;
		
		addRsltStr += exponent.to_string();

		//FRACTION
		bitset<32> sumFrac = op.sumTwoBin(temp1, temp2, ctemp);

		// preparing fraction mask
		bitset<32> fracMask = ~(~0 << FRACBITS);
		
		// extracting fraction bits for sum
		for(int i = 0; i < FRACBITS; i++) {
			frac[FRACBITS - i] = sumFrac[FRACBITS-i] && fracMask[FRACBITS-i];		//keep &&-ing to get last ESBITS values
		}

		if(frac[FRACBITS-1] == 1 && frac[FRACBITS-2] == 0) {
			//Implicit number is 2. Divide fraction by 2 and shift to exponent,
			//if exponent carries out, shift to regime.
			int temp = expVal1.to_ulong() +1;
			if(temp > (pow(2, ESBITS)-1)) {
				regime = add.createReg(k1 + 1);
				exponent = bitset<ESBITS>(0);
			} else {
				exponent = bitset<ESBITS>(temp);
			}
			frac = frac >> 1;
		} else if(frac[FRACBITS-1] == 0 && frac[FRACBITS-2] == 0) {
			//Implicit number is 0. Multiply fraction until implicit number becomes 1.
			//If exponent has enough bits to pull from, do so. Else, pull from regime
			int count = 0;
			while(frac.to_ulong() != 0 && frac[FRACBITS-2] != 1) {
				frac = frac << 1;
				count++;
			}
			if(expVal1.to_ulong() >= count) {	//Since both expVals are same, we can compare either to count
				//If subtracting count from expVal does not result in a negative number
				exponent = bitset<ESBITS>(expVal1.to_ulong() - count);
			} else if (expVal1.to_ulong() < count) {
				//Else, subtract from regime
				regime = add.createReg(k1 - 1);
				exponent = bitset<ESBITS>(exponent.to_ulong() + int(log(USEED)/log(2)));
			}
		}

		bitset<FRACBITS-2> fracFinal;

		for(int i = 0; i < 23; i++)
			fracFinal[i] = frac[i];

		fracFinal[FRACBITS-2] = 0;

		//cout << sign << " " << regime << " " << exponent << " " << fracFinal << endl;
		
		addRsltStr = sign.to_string() + regime.to_string() + exponent.to_string() + fracFinal.to_string();
		//cout << addRsltStr << endl;
		
		//convert string to bitset
			
		bitset<32> addRslt(addRsltStr);
		//cout << addRslt <<  " j " << endl;
		return addRslt;
	}
};

class AddPositToDec {
	
	bitset<1> posSign;
	bitset<REGBITS> posRegime;
	bitset<ESBITS> posExp;
	bitset<FRACBITS> posFrac;
	
	private:
	
	public:
	
	AddPositToDec(){}
	
	void posExtract(bitset<32> positVal) {
		//extract
		
		//SIGN
		posSign[0] = positVal[31];
		
		//REGIME
		// preparing regime mask
		bitset<32> regMask = ~(~0 << REGBITS);
		regMask = regMask << (31-REGBITS);
		
		//extracting regime from posit
		for(int i = 0; i < REGBITS; i++) {
			posRegime[REGBITS - i] = positVal[31 - i] && regMask[31 - i];		//obtain regime bits
		}
		
		//EXPONENT
		// preparing exponent mask
		bitset<32> expMask = ~(~0 << ESBITS);
		expMask = expMask << (25-ESBITS);
		
		//and them
		for(int i = 0; i < ESBITS; i++) {
			posExp[ESBITS - 1 - i] = positVal[24-i] && expMask[24-i];		//obtain exponent bits
		}

		//FRACTION
		// preparing fraction mask
		bitset<32> fracMask = ~(~1 << abs(FRACBITS-3));

		//and them
		for(int i = 0; i < FRACBITS; i++) {
			posFrac[FRACBITS - 2 - i] = positVal[23-i] && fracMask[23-i];	//obtain fraction bits
		}

		posFrac = posFrac << 2;

		//cout << "PTD val: " << posSign << " " << posRegime << " " << posExp << " " << posFrac << endl;
		//cout << endl;
	}
	
	int control(bitset<32> positVal) {
		AddDecToPosit dtp;
		bitset<32> zeroPosit = dtp.control(0);
		if(positVal == zeroPosit) {return 0;}
		
		posExtract(positVal);
		
		//do sign last to check if less than 0 and assign negative through subtraction of value from 0
		
		//Calculate k
		GivenTwoPosit gtp;
		int k = gtp.calcKPosit(posRegime);
		
		//Calculate expVal
		double expVal = posExp.to_ulong();
		
		//Calculate fraction sum
		double divide = 0.5;
		double sum = 0;
		for(int i = 0; i < FRACBITS; i++) {
			if(posFrac[(FRACBITS-1)-i] == 1) {
				sum += divide;
			}
			divide /= 2;
		}
		sum += 1;			//Add implicit one
		
		double useed = pow(2, (pow(2, ESBITS)));
		double exponent = pow(2, expVal);
				
		double decimal = pow(useed, k) * exponent * sum;
		
		return decimal;
	}
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
