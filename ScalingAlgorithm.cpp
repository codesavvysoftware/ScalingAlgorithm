// ScalingAlgorithm.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <math.h>

double TwoPiR = 2 * 3.14159;

unsigned long long GetScaledVal64(double);

unsigned int CalcBeta(unsigned int ScaleFactor, unsigned int Freq, unsigned int SampleRate)
{
	static bool KCalculated;

	static unsigned int K;

	if (!KCalculated) {

		double ApproxForBeta1HZ1US = 0.000003117109216600;

		unsigned long long temp = GetScaledVal64(ApproxForBeta1HZ1US);

		K = temp;

		KCalculated = true;
	}

	unsigned int shift = 26 - ScaleFactor;

	unsigned int ApproxScaled = K;

	ApproxScaled >>= 6;

	unsigned int FreqSampleScaled = Freq * SampleRate;

	unsigned int BetaScaled = FreqSampleScaled * ApproxScaled;

	BetaScaled >>= shift;

	return BetaScaled;
}

unsigned int GetScaledVal(double FracPart, unsigned int ScaleFactor = 16)
{
	double accum = FracPart;

	while (accum > 1.0)
	{
		accum -= 1.0;
	}

	unsigned int scaledval = 0;

	if (accum == 0.0)
	{
		return 0;
	}

	int i = 0;

	double PowerOf2;


	for (i = 0; i < ScaleFactor; i++)
	{
		PowerOf2 = 1.0;

		PowerOf2 /= (2 << i);

		double interVal = accum - PowerOf2;

		scaledval <<= 1;

		if (interVal >= 0.0)
		{
			scaledval |= 1;

			accum = interVal;
		}

		if (interVal == 0.0) break;

	}

	return scaledval;

}
unsigned long long GetScaledVal64(double FracPart)
{
	double accum = FracPart;

	while (accum > 1.0)
	{
		accum -= 1.0;
	}

	unsigned long long scaledval = 0;

	if (accum == 0.0)
	{
		return 0;
	}

	int i = 0;

	double PowerOf2;


	for (i = 0; i < 32; i++)
	{
		PowerOf2 = 1.0;

		PowerOf2 /= (2 << i);

		double interVal = accum - PowerOf2;

		scaledval <<= 1;

		if (interVal >= 0.0)
		{
			scaledval |= 1;

			accum = interVal;
		}

		if (interVal == 0.0) break;

	}

	return scaledval;

}

unsigned int Freq(unsigned int CornerFreqHZ, unsigned int SamplePeriodUS, unsigned int RC)
{
	double Period = SamplePeriodUS;

	Period /= 1000000.0;

	double Freq = CornerFreqHZ;

	Freq *= 6.28318;

	double Numerator = Freq * Period;

	double Denominator = RC;

	Denominator += Numerator;

	double Beta = Numerator / Denominator;

	return GetScaledVal(Beta);
}

unsigned int TakeRecip(unsigned int val)
{
	unsigned int OnePointZero = static_cast<unsigned int>(-(static_cast<int>(val)));

	unsigned int accum = 1;

	while (OnePointZero >= val)
	{
		accum += 1;

		OnePointZero -= val;
	}

	return accum;
}

unsigned long long TakeRecip64(unsigned long long val)
{
	unsigned long long OnePointZero = static_cast<unsigned long long>(-(static_cast<long long>(val)));

	unsigned long long accum = 1;

	while (OnePointZero >= val)
	{
		accum += 1;

		OnePointZero -= val;
	}

	return accum;
}
unsigned int ScaledFreq(unsigned int CornerFreqHZ, unsigned int SamplePeriodUS, unsigned int RC)
{
	unsigned int TwoPiScaled = 0x6487e;

	unsigned int RadiansPerSec = TwoPiScaled * CornerFreqHZ;

	unsigned int RadiansRecip = TakeRecip(RadiansPerSec);

	unsigned int CornerFreqRecip = TakeRecip(CornerFreqHZ);

	unsigned int Numerator = RC << 16;

	Numerator *= (1000000 << 16);

	Numerator *= RadiansRecip;

	Numerator >>= 16;

	Numerator *= CornerFreqRecip;

	Numerator >>= 16;
	
	unsigned int Scaled = 0;

	return Scaled;
}

double ConvertScaledNumToReal(unsigned int ScaledNum, unsigned int ScaleFactor = 16)
{
	unsigned int Mask = 1 << ScaleFactor;

	Mask -= 1;

	unsigned int Mask0 = Mask ^ 0xffffffff;


	unsigned int WholeNum = (ScaledNum & Mask0) >> ScaleFactor;

	double WholeNumPart = WholeNum;
	
	unsigned int Frac = ScaledNum;
	
	Frac &= Mask;

	double FracPart = Frac;

	double denominator = pow(2.0, (double)ScaleFactor);

	FracPart *= (1.0 / denominator);

	double RealVal = WholeNumPart + FracPart;

	return RealVal;
}

double ConvertScaledNumToReal64(unsigned long long ScaledNum)
{
	double WholeNumPart = (ScaledNum & 0xffffffff00000000) >> 32;

	double FracPart = (ScaledNum & 0xffffffff);

	FracPart *= (1.0 / 0x100000000);

	double RealVal = WholeNumPart + FracPart;

	return RealVal;
}
#define U16REC1(A,M,S) (unsigned int)( ( ( (unsigned int)(A) * (unsigned int)(M) )  >> 16u) >> (S))
#define U16REC2(A,M,S) (unsigned int)( ( ( ( ( (unsigned int)(A) * (unsigned int)(M) )  >> 16u) + (A)) >> 1u ) >> (S))


int main()
{
	unsigned int Freq1 = 1;

	unsigned int Freq10 = 10;

	unsigned int Freq20 = 20;

	unsigned int Freq50 = 50;

	unsigned int Freq100 = 100;

	unsigned int SampleRate = 50;

	unsigned int ScaleFactor = 26;

	unsigned int Beta50Scaled26 = CalcBeta(ScaleFactor, Freq50, SampleRate);

	double Beta50RealScaled26 = ConvertScaledNumToReal(Beta50Scaled26, ScaleFactor);

	ScaleFactor = 16;
	
	unsigned int Beta50Scaled16 = CalcBeta(ScaleFactor, Freq50, SampleRate);

	double Beta50RealScaled16 = ConvertScaledNumToReal(Beta50Scaled16, ScaleFactor);


	unsigned int SampleRateTest = 50;

	/*while (SampleRate > 0) {
		Beta1 = SampleRate * ApproxScaled;

		Beta1 >>= 10;

		if (Beta1 == 0) {
			break;
		}

		SampleRate--;

	}*/
	
	double FreqNew = 1;

	double USConv = 1000000;

	double MSConv = 1000;

	double CyclesPerUS = FreqNew / USConv;

	double PeriodNewUS = 1;

	double PiNew = 3.14159;

	double PiRecipNew = 1.0 / PiNew;

	double BetaNew = CyclesPerUS * PeriodNewUS;

	BetaNew /= (PiRecipNew + BetaNew);
	
	unsigned long long BaseFactor = GetScaledVal64(0.0001570791392);

	double DenomNew = CyclesPerUS * PeriodNewUS * PiNew;

	DenomNew = 1.0 / DenomNew;

	DenomNew += 1.0;

	DenomNew = 1.0 / DenomNew;


	double Freq50Rate50 = DenomNew * 2500.0;


	unsigned long BaseFactorX = BaseFactor;

	unsigned int F50 = 2500;

	F50 *= (BaseFactorX >> 4);

	unsigned int F50Shift = F50 >> 16;

	double F50Real = ConvertScaledNumToReal(F50Shift);

	unsigned int n = 10 << 16;
	
	unsigned int q, r;
	
	q = (n >> 1) + (n >> 2);

	q += (q >> 4);

	q += (q >> 8);

	q += (q >> 16);

	q >>= 3;

	r = n - (((q << 2) + q) << 1);

	q += (r > 9);
	
	unsigned int OneHalfTest = 3183;



	unsigned int OneTenthScaled = U16REC1(OneHalfTest, 0xcccdu, 2u);

	double OneTenth = ConvertScaledNumToReal(OneTenthScaled);

	
	double OneOverPiFracPart = 0.098861838;

	unsigned int OneOverPix10000ScaledFrac = GetScaledVal(OneOverPiFracPart);

	unsigned int OneOverPix10000ScaledWhole = (3183 << 16);

	unsigned int OneOverPix10000Scaled = OneOverPix10000ScaledWhole | OneOverPix10000ScaledFrac;

	double OneOverPix10000Check = ConvertScaledNumToReal(OneOverPix10000Scaled);

	unsigned int CheckVal = 3183 << 16;

	unsigned int OneOverPix10000IntPart = (OneOverPix10000Scaled & 0xffff0000);// >> 16;

	unsigned int TwoOverTauScaledWhole = U16REC2(OneOverPix10000IntPart, 0x47afu, 5u);

//	TwoOverTauScaledWhole <<= 16;

	double TwoOverTauCheck = ConvertScaledNumToReal(TwoOverTauScaledWhole);

	//unsigned int TwoOverTauScaledFrac = U16REC1(O)

	TwoOverTauScaledWhole = U16REC1(OneOverPix10000Scaled, 0xcccdu, 3u);

	double Period = 0.5;

	unsigned int PeriodScaled = GetScaledVal(Period);

	double TwoOverTau = ConvertScaledNumToReal(TwoOverTauScaledWhole);

	double Omega = Period / (TwoOverTau + Period);


	unsigned long long TwoPi64 = GetScaledVal64(2.0*.14159);

	TwoPi64 |= 0x600000000;

	double TwoPiReal = 3.14159 * 2;

	double TwoPiTest = ConvertScaledNumToReal64(TwoPi64);

	unsigned long long FreqHz = 1;

	unsigned long long SamplePeriodUS = 50;

	unsigned long long Microseconds = 1000000;

	Microseconds <<= 32;

	unsigned long long SamplePeriodSeconds = TakeRecip64( Microseconds );

	SamplePeriodSeconds *= SamplePeriodUS;

	unsigned long long TwoPiFreqHz = TwoPi64 * FreqHz;

	unsigned long long IntComponent = (TwoPiFreqHz >> 32) * SamplePeriodSeconds;

	unsigned long long FracComponent = ((TwoPiFreqHz & 0xffffffff) * SamplePeriodSeconds) >> 32;

	unsigned long long RadiansPerSec = IntComponent + FracComponent;

	double RadiansPerSecReal = (2 * 3.14159 * 50 * 1) / 1000000;

	double RadiansPerSecTest = ConvertScaledNumToReal64(RadiansPerSec);

	unsigned long long RC64 = 0x200000000;

	unsigned long long Denom = TakeRecip64(RC64 + RadiansPerSec);

	unsigned long long Beta = (RadiansPerSec * Denom ) >> 32;

	double DenomReal = ( RC64 >> 32 ) + RadiansPerSecReal;
	double BetaReal = RadiansPerSecReal / DenomReal;
	double BetaTest = ConvertScaledNumToReal64(Beta);

	unsigned long long BetaInter = 0x8000;

	BetaInter += Beta;

	BetaInter >>= 16;

	unsigned int Beta32 = BetaInter;

	double Beta32Real = ConvertScaledNumToReal(Beta32);





	unsigned int TwoPointZero = 0x20000;
	
	unsigned int RecipOfTwoPointZero = TakeRecip(TwoPointZero);

	double FracPart = TwoPiR - 6.0;



	unsigned int scaledVal = GetScaledVal(FracPart);

	unsigned int TwoPiScaled = 0x60000;

	TwoPiScaled |= scaledVal;


	unsigned int OneOverTwoPi = TakeRecip(TwoPiScaled);

	unsigned int FreqScaled = 50;

	unsigned int RadiansHzScaled = TwoPiScaled * FreqScaled;

	double FreqReal = 50;

	double RadiansHzReal = (2.0 * 3.14159 * 50);

	double RadiansHzTest = ConvertScaledNumToReal(RadiansHzScaled);

	unsigned int DenomRightTermDenom = TakeRecip(RadiansHzScaled);

	double DenomRightTermDenomReal = 1.0 / RadiansHzReal;

	double DenomRightTermDenomTest = ConvertScaledNumToReal(DenomRightTermDenom);

	unsigned int RCScaled = 2;

	double RCReal = 2.0;

	unsigned int DenomRightTermScaled = RCScaled * DenomRightTermDenom;

	double DenomRightTermTest = ConvertScaledNumToReal(DenomRightTermScaled);

	double DenomRightTermReal = RCReal * DenomRightTermDenomReal;

	unsigned int BetaDenomScaled = 1000000 << 16;

	BetaDenomScaled += DenomRightTermScaled;

	double BetaDenomReal = 1000000.0 + DenomRightTermReal;

	double BetaDenomTest = ConvertScaledNumToReal(BetaDenomScaled);

	unsigned int BetaDenomRecipScaled = TakeRecip(BetaDenomScaled);

	unsigned int BetaDenomRecipReal = 1.0 / BetaDenomReal;

	double BetaDenomRecipTest = ConvertScaledNumToReal(BetaDenomRecipScaled);

	unsigned int BetaScaled = 1000000 * BetaDenomRecipScaled;

	BetaReal = 1000000.0 * BetaDenomRecipReal;

	BetaTest = ConvertScaledNumToReal(BetaScaled);






	double x = (TwoPiR * 50.0 * 50.0)/1000000.0;

	unsigned int rightAnswwer = GetScaledVal(x);

	double y = (TwoPiR * 50.0);

	unsigned int Another = GetScaledVal(y);

	unsigned int freq = Freq(50, 50, 2);

	unsigned int scaledFreq = ScaledFreq(50, 50, 2);


}

