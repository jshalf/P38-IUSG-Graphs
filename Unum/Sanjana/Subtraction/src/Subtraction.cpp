/**
 * Subtraction.cpp
 * Four Arithmetic Operations using Posit Numbers
 *
 *	TBD: Was currently working on implementing subtraction of
 *	values between 0 and 1 (ex. 0.007 - 67)
 *
 * Sanjana Shah
 * August 13, 2017
 */

#include <iostream>
#include <cmath>
#include <string>
#include <bitset>
#include <vector>

using namespace std;

//In temporary registers, 24th and 25th bit hold the carry-out, and the implicit 1
//However, these two bits are not represented in the final form
//TBD: Maybe don't need 25 bits for fraction, only 24, as the only case
//is to check if implicit bit (bit 24) is 0 or 1 (Subtraction will never result in implicit 2)
//Need to check ^ if 01-10 -> needs impicit one for carry
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
			double temp = pow(USEED, k);
			//using useed^k, keep multiplying by 2 until closest 2 power is reached without surpassing
			while(decRep >= temp*2) {
				temp *= 2;
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
			bitset<32> posBitset;               //final bitset with running place
            int runCount = 31;

			//SIGN
			if(s == 0) posBitset[runCount] = 0;
			else if(s == 1) posBitset[runCount] = 1;
            runCount--;

			//REGIME
			int shift = REGBITS-1;
			//vector<bool> temp (regLength, 0);
            bitset<REGBITS> temp;
            
			if(decRep >= 0 && decRep < 1) {		//If the value is less than 1
				for(int i = k; i < 0; i++) {
					temp[shift] = 0;			//Add k number of zeros
					shift--;
					regLength++;
				}
				temp[shift] = 1;				//Append terminating 1
				shift--;
				regLength++;
			} else {
				for(int i = 0; i < k+1; i++) {
					temp[shift] = 1;			//Add k+1 number of ones
					regLength++;
					shift--;
				}
				temp[shift] = 0;				//Append terminating 0
				shift--;
				regLength++;
			}
			for(int i = REGBITS-1; i >= REGBITS-regLength; i--) {	//Append regime to the bitset
				posBitset[runCount] = temp[i];
                runCount--;
            }
            
			//EXPONENT
			string str;
			bitset<ESBITS> exp (expVal);
			for(int i = (ESBITS-1); i >= 0; i--) {		//Append exponent to the bitset
                posBitset[runCount] = exp[i];
                runCount--;
			}

			//MANTISSA
            
            int extra = REGBITS - regLength;
            
			frac -= 1.0;
			shift = 22;
            vector<bool> temp2 (23+extra, 0);
            vector<bool> temp3 (23+extra, 0);
			double divResult = frac;
			int rem;

			if(divResult != 0) {						//Convert decimal representation into binary bitset
				while(divResult != 1 && shift > 0) {
					divResult *= 2;
					if(divResult > 1) {
						divResult -= 1;
						rem = 1;
					} else {
						rem = 0;
					}
					temp2[shift] = rem;					//Append remainder to fraction bitset
					shift--;
				}
                temp3[shift+1] = 1;
				for(int i = 0; i < 23; i++)                                                                             //??why32 not 23
					temp2[i] = temp2[i] || temp3[i];

                for(int i = 22; i >= 0; i--){
					posBitset[runCount] = temp2[i];
                    runCount--;
				}
			} else {
                for(int i = 22; i >= 0; i--){
                    posBitset[runCount] = 0;            //mant = bitset of 23 0's - simplify
                    runCount--;
                }
			}

			return posBitset;
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
			//cout << "k " << k << "e " << expVal << endl;
			//cout << "closest 2 pwr " << pow(USEED, k)*pow(2, expVal) << endl;
			closestTwoPwr = pow(USEED, k)*pow(2, expVal);
			double frac = calcFrac();
			//cout << frac << endl;

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
 * Given two Posit 32-bit bitsets, extracts the four parts - sign, regime, exponent, mantissa
 */
class Extraction {

	bitset<1> s1;
	bitset<ESBITS> expVal1;
    bitset<1> s2;
	bitset<ESBITS> expVal2;

    vector<bool> reg1;
    vector<bool> reg2;
    vector<bool> frac1;
    vector<bool> frac2;
    
	private:

	public:
		void extract(bitset<32> x, int regLengthX, bitset<32> y, int regLengthY) {
			//SIGN
			s1[0] = x[31];
			s2[0] = y[31];

			//REGIME
            reg1 = vector<bool> (regLengthX, 0);                                //NEED TO MAKE CONSTANT - can use extra in frac (6-3 for 79, 38)
            reg2 = vector<bool> (regLengthY, 0);                                //NEED TO MAKE CONSTANT
            
			// preparing regime mask for first number
			bitset<32> regMaskX = ~(~0 << regLengthX);
			regMaskX = regMaskX << (31-regLengthX);
            
			//extracting regime from posit
			for(int i = 1; i <= regLengthX; i++) {
				reg1[regLengthX - i] = x[31 - i] && regMaskX[31 - i];		//obtain regime bits for first number
			}

            //preparing regime mask for second number
            bitset<32> regMaskY = ~(~0 << regLengthY);
            regMaskY = regMaskY << (31-regLengthY);
            
            //extracting regime from posit
			for(int i = 1; i <= regLengthY; i++) {
				reg2[regLengthY - i] = y[31 - i] && regMaskY[31 - i];		//obtain regime bits for second number
			}
            
			//EXPONENT
			// preparing exponent mask for first number
			bitset<32> expMaskX = ~(~0 << ESBITS);                           //es bits fixed, cannot be variable
			expMaskX = expMaskX << ((FRACBITS+(REGBITS-regLengthX))-ESBITS);
            
            // preparing exponent mask for first number
            bitset<32> expMaskY = ~(~0 << ESBITS);                           //es bits fixed, cannot be variable
            expMaskY = expMaskY << ((FRACBITS+(REGBITS-regLengthY))-ESBITS);

			//and them
			for(int i = ESBITS; i > 0; i--) {
				expVal1[ESBITS-i] = x[(FRACBITS+(REGBITS-regLengthX))-i] && expMaskX[(FRACBITS+(REGBITS-regLengthX))-i];		    //obtain exponent bits for first number
                expVal2[ESBITS-i] = y[(FRACBITS+(REGBITS-regLengthY))-i] && expMaskY[(FRACBITS+(REGBITS-regLengthY))-i];            //obtain exponent bits for second number  25+6-3-(1)  25+6-3-(0)
			}
            
			//FRACTION
			// preparing fraction mask
            int extra1 = REGBITS - regLengthX;
            int extra2 = REGBITS - regLengthY;

            vector<bool> fracMaskX (32, 0);
            vector<bool> fracMaskY (32, 0);
            
            for(int i = 0; i < FRACBITS+extra1-2; i++){
                fracMaskX[i] = 1;
            }
            for(int i = 0; i < FRACBITS+extra2-2; i++){
                fracMaskY[i] = 1;
            }
            
            frac1 = vector<bool> (FRACBITS+extra1-2, 0);
            frac2 = vector<bool> (FRACBITS+extra2-2, 0);
            
			//and them
			for(int i = 0; i < FRACBITS+extra1-2; i++) {
				frac1[i] = x[i] && fracMaskX[i];		//obtain fraction bits for first number
			}
			for(int i = 0; i < FRACBITS+extra2-2; i++) {
				frac2[i] = y[i] && fracMaskY[i];		//obtain fraction bits for second number
			}
            
			//add implicit 1 to fraction in 24th bit
			frac1[FRACBITS+extra1-1] = 1;
			frac2[FRACBITS+extra2-1] = 1;

            cout << s1 << " ";
            for(int i = reg1.size()-1; i >= 0; i--) cout << reg1[i];
            cout << " " << expVal1 << " ";
            for(int i = frac1.size()-1; i >= 0; i--) cout << frac1[i];
            cout << endl;
            
            cout << s2 << " ";
            for(int i = reg2.size()-1; i >= 0; i--) cout << reg2[i];
            cout << " " << expVal2 << " ";
            for(int i = frac2.size()-1; i >= 0; i--) cout << frac2[i];
            cout << endl;
            
        }

		//Getters
		bitset<1> 			getSign1() {return s1;}
		vector<bool> 	    getReg1()  {return reg1;}
		bitset<ESBITS> 		getExp1()  {return expVal1;}
		vector<bool>     	getFrac1() {return frac1;}
		bitset<1> 			getSign2() {return s2;}
		vector<bool>     	getReg2()  {return reg2;}
		bitset<ESBITS> 		getExp2()  {return expVal2;}
		vector<bool>    	getFrac2() {return frac2;}
};

class Operations {
    
public:
    vector<bool> subTwoBin(vector<bool> num1, vector<bool> num2, int extra1, int extra2)
    {
        vector<bool> diff (num1.size(), 0);
        int temp = 0;
        for(int i = 0; i < num1.size(); i++){
            if(num1[i] == 0 && num2[i] == 1){
                int back = i;
                while (num1[back] != 1) {
                    temp++;
                    back++;
                }
                num1[back] = 0;
                for(int j = 1; j <= temp; j++) {
                    num1[back-j] = 1;
                }
                diff[i] = 1;
            } else {
                diff[i] = num1[i]-num2[i];
            }
        }
        return diff;
    }
};


class Subtraction {

	bitset<1> s;
	bitset<REGBITS> reg;
	bitset<ESBITS> expVal;
	bitset<FRACBITS> frac;

	bitset<1> s1;
    bitset<ESBITS> expVal1;
	bitset<1> s2;
	bitset<ESBITS> expVal2;
    
    vector<bool> reg1;
    vector<bool> reg2;
    
    vector<bool> frac1;
    vector<bool> frac2;

	public:

		void setSign(bitset<1> sign) {
			s = sign;
		}

		void setValues(bitset<1> n1,
				vector<bool> r1,
				bitset<ESBITS> e1,
				vector<bool> f1,
				bitset<1> n2,
				vector<bool> r2,
				bitset<ESBITS> e2,
				vector<bool> f2) {
			s1 = n1;
			reg1 = r1;
			expVal1 = e1;
			frac1 = f1;
			s2 = n2;
			reg2 = r2;
			expVal2 = e2;
			frac2 = f2;
		}

		vector<bool> createReg(int k) {
			int regLength = 0;
			int shift = REGBITS-1;
			vector<bool> temp2 (6, 0);

			if(k < 0 || s == 1) {	//If the value is less than 1 or is negative
				k = abs(k);
				for(int i = 0; i < k; i++) {
					temp2[shift] = 0;		//Add k number of zeros
					shift--;
					regLength++;
				}
				temp2[shift] = 1;			//Append terminating 1
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
		vector<bool> mantToBin(double num, int size) {
			int shift = size - 1;
			vector<bool> temp4 (size, 0);
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
		vector<bool> compareReg(DecToPosit dtp, DecToPosit dtp2) {
			cout << "Comparing Regimes" << endl;
			long r1 = dtp.getK();
			long r2 = dtp2.getK();

			//decimal (base 10 representations) of the two fractions
			//TBD: Future work is to convert to bitset (binary format) then shift
			//	   Here I am just converting to int, dividing, then converting to bitset
            
            int num = 0;
            for(int i = 0; i < frac1.size(); i++) if(frac1[i] == 1) { num += 2^i; }
            cout << "f1 got " << num << endl;
            long f1 = num;
            
            num = 0;
            for(int i = 0; i < frac2.size(); i++) if(frac2[i] == 1) { num += 2^i; }
            cout << "f2 got " << num << endl;
			long f2 = num;
            
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
				//frac2 = f2;

				double x = 1.0/(pow(USEED, count));                 //replace with division

				vector<bool> y = mantToBin(x, frac1.size());
                
				for(int i = 0; i <= frac1.size(); i++){
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
                cout << f1 << endl;
				//frac1 = f1;

				double x = 1.0/(pow(USEED, count));
				vector<bool> y = mantToBin(x, frac1.size());
               
                for (int i = 0; i < y.size(); i++) {
                    cout << y[i];
                }
                cout << endl;
                for (int i = 0; i < y.size(); i++) {
                    cout << y[i];
                }
                
				for(int i = 0; i <= frac1.size(); i++){
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
		vector<bool> compareExpVal(vector<bool> frac1, vector<bool> frac2) {
			cout << endl << "Comparing Exponents" << endl;

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
					cout << "Dividing frac2 by 2" << endl;
                    frac2.insert(frac2.end(), 0);                 //INSERTING 0
                    frac2.erase(frac2.begin());                   //DELETING 0
					count++;
				}
				return frac2;
			} else if(exp2 > exp1) {
				cout << "exp2 > exp1" << endl;
				while(exp1 != exp2) {
					exp1 += 1;
					cout << "Dividing frac1 by 2" << endl;
                    frac1.insert(frac1.end(), 0);                 //INSERTING 0
                    frac1.erase(frac1.begin());                   //DELETING 0
					count++;
				}
				return frac1;
			}
			//The two exponents are now the same, default return can be any of the two fractions
			return frac1;
		}
};

int main() {

	bitset<1> s1;
	vector<bool> reg1;
	bitset<ESBITS> expVal1;
	vector<bool> frac1;
	bitset<1> s2;
	vector<bool> reg2;
	bitset<ESBITS> expVal2;
	vector<bool> frac2;

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

	bitset<32> x, y;

	//Calls DecToPosit for each value, creating a 32-bit posit value
	x = dtp.control(decRepUsr);
	y = dtp2.control(decRepUsr2);

	//extract values and store locally in main
	Extraction ext;
	ext.extract(x, dtp.getRegLength(), y, dtp2.getRegLength());

	s1 = ext.getSign1();
	reg1 = ext.getReg1();
	expVal1 = ext.getExp1();
	frac1 = ext.getFrac1();
	s2 = ext.getSign2();
	reg2 = ext.getReg2();
	expVal2 = ext.getExp2();
	frac2 = ext.getFrac2();

	Subtraction subt;
	subt.setValues(s1, reg1, expVal1, frac1, s2, reg2, expVal2, frac2);

	bitset<32> temp1;
	bitset<32> temp2;
	bitset<FRACBITS> frac;

	for(int i = 32; i > 0; i--) {
		temp1[i] = frac1[i];
		temp2[i] = frac2[i];
	}

	//SIGN
	bitset<1> sign;
	if(decRepUsr2 > decRepUsr)
		sign = 1;
	else
		sign = 0;

	subt.setSign(sign);

	//set regimes equal to each other, then call compareReg, updating fraction
	if(dtp.getK() > dtp2.getK()) {
		reg2 = reg1;
		frac2 = subt.compareReg(dtp, dtp2);

	} else if(dtp.getK() < dtp2.getK()){
		reg1 = reg2;
		frac1 = subt.compareReg(dtp, dtp2);
	}

	//set exponents equal to each other, then call compareReg, updating fraction
	if(expVal1.to_ulong() > expVal2.to_ulong()) {
		expVal2 = expVal1;
		frac2 = subt.compareExpVal(frac1, frac2);
	} else if(expVal1.to_ulong() < expVal2.to_ulong()){
		expVal1 = expVal2;
		frac1 = subt.compareExpVal(frac1, frac2);
	}

   
    cout << s1 << " ";
    for(int i = reg1.size()-1; i >= 0; i--) cout << reg1[i];
    cout << " " << expVal1 << " ";
    for(int i = frac1.size()-1; i >= 0; i--) cout << frac1[i];
    cout << endl;
    
    
    cout << s2 << " ";
    for(int i = reg2.size()-1; i >= 0; i--) cout << reg2[i];
    cout << " " << expVal2 << " ";
    for(int i = frac2.size()-1; i >= 0; i--) cout << frac2[i];
    cout << endl;
    
	//exponent and regime values are now same for both Posit numbers

	//SIGN -> done above (due to createReg method needing updated sign value)

	//REGIME
	vector<bool> regime = reg1;

	//EXPONENT
	bitset<ESBITS> exponent = expVal1;

	//FRACTION
	//bitset<FRACBITS> one = 1;

	Operations op;

	//cout << frac1 << endl;
	//cout << frac2 << endl;

	//frac2 = op.subTwoBin(~(frac2), one);
    //cout << frac2 << endl;

    //check if second is smaller or bigger than first - if bigger then switch fractions
    //and sign of final number
    bool change = false;
    
    for(int i = frac1.size(); i >= 0; i--){
        if(frac1[i] == 1)
            break;
        if(frac2[i] == 1)
            change = true;
    }
    
    vector<bool> temp (frac1.size(), 0);
    
    if(change){
        //switch frac1 and frac2 because the second number is greater than the first
        for(int i = 0; i < frac1.size(); i++){
            temp[i] = frac1[i];
            frac1[i] = frac2[i];
            frac2[i] = temp[i];
        }
        //switch the final sign
        if(sign == 0) sign = 1; else sign = 0;
    }
    
	vector<bool> subFrac = op.subTwoBin(frac1, frac2, REGBITS-dtp.regLength, REGBITS-dtp2.regLength);
    subFrac[0] = 0;                         //arithmetic mean round
    
	//TBD: Maybe implicit 2 case is not needed, only implicit 0 case is.
	//	   Could get rid of top half of if condition.
	if(subFrac[FRACBITS-2] == 0 && subFrac[FRACBITS-1] == 1) {
		//Implicit number is 2. Divide fraction by 2 and shift to exponent,
		//if exponent carries out, shift to regime.
		int temp = expVal1.to_ulong() +1;
		if(temp > (pow(2, ESBITS)-1)) {
			cout << "Carry to reg2" << endl;
			regime = subt.createReg(dtp.getK() + 1);
			exponent = bitset<ESBITS>(0);
		} else {
			cout << "Carry to expVal2" << endl;
			exponent = bitset<ESBITS>(temp);
		}
		frac = frac >> 1;
	} else if (subFrac[FRACBITS-2] == 0 && subFrac[FRACBITS-1] == 0) {
		//Implicit number is 0. Multiply fraction until implicit number becomes 1.
		//If exponent has enough bits to pull from, do so. Else, pull from regime
		cout << "Left Shifting until 1" << endl;
		int count = 0;
		while(frac.to_ulong() != 0 && frac[FRACBITS-2] != 1) {
            subFrac.erase(subFrac.end());
            subFrac.insert(subFrac.begin(), 0);
			count++;
            for(int i = subFrac.size()-1; i >= 0; i--) cout << subFrac[i];
            cout << endl;
		}
		//Both exponent values are same. Can compare either to count
		if(expVal1.to_ulong() >= count) {
			//If subtracting count from the exponent does not result in a negative number,
			//then subtract count from exponent
			exponent = bitset<ESBITS>(expVal1.to_ulong() - count);
		} else if (expVal1.to_ulong() < count) {
			//Else, subtract from regime
			regime = subt.createReg(dtp.getK() - 1);
			exponent = bitset<ESBITS>(exponent.to_ulong() + int(log(USEED)/log(2)));
		}
	}
    

    bool exist = false; //check if there is a 1 that we can bring to the beginning
    
    for(int i = 0; i < subFrac.size(); i++) {
        if(subFrac[i] == 1) exist = true;
    }

    int count = 0;
    if(exist == true && subFrac.back() == 0){
        while(subFrac.back() != 1) {
            subFrac.erase(subFrac.end());
            subFrac.insert(subFrac.begin(), 0);             //make frac 23*2 bits - losing precision
            count++;
        }
    }
    
    cout << "count " << count << endl;
    //fraction pushed forward, count stores
    
    cout << "this is k: " << dtp.getK() << endl;
    
    for(int i = 0; i < subFrac.size(); i++) cout << subFrac[i];
    cout << endl;
    
    //Both exponent values are same. Can compare either to count
    if(expVal1.to_ulong() >= count) {
        //If subtracting count from the exponent does not result in a negative number,
        //then subtract count from exponent
        exponent = bitset<ESBITS>(expVal1.to_ulong() - count);
        cout << "exponent " << exponent << endl;
    } else if (expVal1.to_ulong() < count) {
        //Else, subtract from regime
        regime = subt.createReg(dtp.getK() - 1);
        exponent = bitset<ESBITS>(exponent.to_ulong() + int(log(USEED)/log(2)));
        cout << "borrow from regime ";
        for(int i = regime.size(); i > 0; i--) cout << regime[i];
        cout << ". exponent: " << exponent << endl;
    }
    /*
    if(subFrac[FRACBITS-2] == 0) {
        cout << "Left Shifting until 1." << endl;
        //subFrac = subFrac << 1;
        int count = 0;
        while(subFrac[FRACBITS-1] != 1) {
            subFrac = subFrac << 1;
            count++;
            cout << subFrac << endl;
            cout << "subFrac " << subFrac << " " << subFrac[FRACBITS-1] << subFrac[FRACBITS-2] << endl;
        }
        //Both exponent values are same. Can compare either to count
        if(expVal1.to_ulong() >= count) {
            //If subtracting count from the exponent does not result in a negative number,
            //then subtract count from exponent
            exponent = bitset<ESBITS>(expVal1.to_ulong() - count);
        } else if (expVal1.to_ulong() < count) {
            //Else, subtract from regime
            cout << "there " << endl;
            regime = subt.createReg(dtp.getK() - 1);
            exponent = bitset<ESBITS>(exponent.to_ulong() + int(log(USEED)/log(2)));
        }
    }*/


	//PRINT

	cout << endl;

	bitset<ESBITS> expFinal;

	for(int i = 0; i < ESBITS; i++) {
		expFinal[i] = exponent[i];
	}

	cout << sign << " ";
    for(int i = regime.size()-1; i >= 0; i--) cout << regime[i];

    //1110 or 0001
	/*if(regime[REGBITS-1] == 1) {
		int regCount = regime.size();
		while(regCount >= 0 && regime[regCount] != 0) {
			cout << regime[regCount];
			regCount--;
		}
		cout << "0";
	} else if(regime[REGBITS-1] == 0) {
		int regCount = regime.size();
		while(regCount >= 0 && regime[regCount] != 1) {
			cout << regime[regCount];
			regCount--;
		}
		cout << "1";
	}*/
    
    cout << " " << expFinal << " ";
    for(int i = subFrac.size()-1; i >= 0; i--) cout << subFrac[i];
    cout << endl;
}