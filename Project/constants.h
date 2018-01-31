/*
	Author:Aneirin John Baker 
	Date : 12/12/17
	Description: Script to keep all of the constants for this program in. Makes it alot eaiser
		to have them all in one spot. 

	To use oscillator flip it needs to be set at compile time not at run time.

	the extern keyword used to define f makes it available to all other scripts which
	this class is a header in. 
*/
#include <fftw3.h>

#define REAL 0
#define IMAG 1
#define I complex<double>(0,1)
#define ONE complex<double>(1,0)
#define ZERO complex<double>(0,0)
#define PI 3.14159265
#define PI2 6.28318530718
#define Oscillator_flip 1 //1 - harmonic 0-Anharmonic

//extern double f=0;
