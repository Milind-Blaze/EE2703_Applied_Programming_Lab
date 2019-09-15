# program to make a list of thousand numbers using the algorithm given in the problem
from math import pi as alpha
 #alpha= round(alpha,4)

def frac(number):
	return number- int(number)

def main(no):
	n = [0]*1000
	n[0] = 0.2000
	for i in range(1,no):
		n[i]= (frac((n[i-1]+alpha)*100))
	
	for i in range(0, no):
		n[i] = round(n[i],4)
		print n[i]
		


main(1000)