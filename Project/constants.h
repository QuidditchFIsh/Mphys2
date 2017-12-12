/*
	Author:Aneirin John Baker 
	Date : 12/12/17
	Description: Script to keep all of the constants for this program in. Makes it alot eaiser
		to have them all in one spot. 

	To use oscillator flip it needs to be set at compile time not at run time.
*/
#define REAL 0
#define IMAG 1
#define I complex<double>(0,1)
#define ONE complex<double>(1,0)
#define ZERO complex<double>(0,0)
#define PI 3.14159265
#define Oscillator_flip 0; //1 - harmonic 0-Anharmonic

double f =0;

void set_f(double input)
{
	f = input;
}