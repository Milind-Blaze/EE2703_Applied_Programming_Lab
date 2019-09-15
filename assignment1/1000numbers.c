// Program to print 1000 numbers with four digits using the given algorithm

#include <stdio.h>
#include <math.h>

#define pi 3.141596

int main(){

	int  i= 0;
	double array[1000]={0},j;
	
	array[0]= 0.2000;

	for(i=1; i<1000; i++){
		array[i]= modf((array[i-1] + M_PI)*100, &j);
		printf("%0.4f \n", array[i]);
	}

}