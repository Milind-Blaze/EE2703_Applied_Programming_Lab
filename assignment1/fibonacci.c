// program to print the first 10 fibonacci numbers

#include <stdio.h>

int main(){

	int i=0, n=1, nold=1, new= 0;
	printf("1  %d \n", nold);
	printf("2  %d \n", n);

	for(i=3; i<=10; i++){

		new= nold+ n;
		nold= n;
		n= new;
		printf("%d  %d\n",i,n);
	}

return 0;

}