from __future__ import division

import numpy as np
from scipy.integrate import quad


def function1(a):
	a= np.array(a)
	return 1.0/(1.0+a**2)

def taninv(a):
	a= np.array(a)
	array= []
	for i in a:
		array.append(quad(function1, 0.0, i))
	array= np.array(array)
#	array[:,0] = [round(i,5) for i in array[:,0]]	
	return (array)

def trapezoid_integrate(a,b,h,f):
	avbl_points= np.linspace(a,b,((b-a)/h)+1)					# making the assumption that b-a/h is integral					
	values= f(avbl_points)
	integral= []
	integral= h*(np.cumsum(values)- 0.5*( [values[0]]*len(values) + values ))
	return np.array(integral)


def trapezoid_integrate_for_loop(a,b,h,f):
	avbl_points= np.linspace(a,b,((b-a)/h)+1)					# making the assumption that b-a/h is integral					
	values= f(avbl_points)
	integral= []
	for i in range(0,len(values)):
		integral.append(h*(np.sum(values[:i+1])- 0.5*(values[0]+values[len(values[:i+1])-1])))
	return np.array(integral)




		
