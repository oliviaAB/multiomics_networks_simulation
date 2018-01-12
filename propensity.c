// Begin propensity.c
#include <stdio.h>      
#include <math.h>
#include <R.h>

#define RNA 0
#define P_NA 1
#define P_A 2

//#define k_TC 5
#define k_TL 1
#define p0_DR 0.02
#define p0_DP 0.0083
#define n 2
#define theta 1000.

double Transcription(double time, double *y) {
	return k_TC*(1-pow(y[P_A],n)/(pow(y[P_A],n)+pow(theta,n)));
//	return k_TC*(1-y[P]/(y[P]+(QTL_TC*theta)));
}

double Translation(double time, double *y){
	return k_TL*y[RNA];
}

double RNAdecay(double time, double *y){
	return p0_DR*y[RNA];
}

double ProteinNAdecay(double time, double *y){
	return p0_DP*y[P_NA];
}

double ProteinAdecay(double time, double *y){
	return p0_DP*y[P_A];
}

double ProteinActivation(double time, double *y){
	return y[P_NA];
}