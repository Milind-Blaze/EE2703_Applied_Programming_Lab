from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from functions import function1, taninv, trapezoid_integrate, trapezoid_integrate_for_loop


# Step 2


start= 0.0
stop= 5.0
num= 51												# to include the last element and avoid obo error
vector= np.linspace(start, stop, num)				# calling x vector


values= function1(vector)

##########
# Step 3

plt.plot(vector, values, "bo")
plt.xlabel("x")
plt.ylabel("$f(x)$")
plt.title("""Plot of $f(x) = 1/(1+x^2) $""")
plt.show()
plt.close()


###########
# Step 4

taninv_calculated= np.array(taninv(vector)[:,0])							# remember that the quad function returns (integral, error)
taninv_actual= np.arctan(vector)
#with open("table.txt","w") as f:
#	for i in range(0,51):
#		onevalue= str(vector[i])+ ' & ' + str(round(taninv_calculated[i],5))+ ' & ' + str(round(taninv_actual[i],5))+ '\\\ \n '
#		f.write(onevalue) 
error= taninv_calculated- taninv_actual




plt.figure(1)

plt.subplot(211)
plt.plot( vector, taninv_calculated, "ro", label= "quad")
plt.plot( vector, taninv_actual, "k-", label="$tan^{-1}(x)$", linewidth=2)
plt.xlabel("x")
plt.ylabel("$\int^x_0 du/(1+u^2) $")
plt.title("A comparision of quad and np.arctan()")
plt.legend(loc= "lower right")

plt.subplot(212)
plt.semilogy(vector, error, "ro")
plt.xlabel("x")
plt.ylabel("Error")
plt.title("Error in $\int^x_0 dx/(1+x^2) $")
plt.tight_layout()
plt.show()
plt.close()



plt.semilogy(vector, (taninv(vector))[:,1], "ro")
plt.xlabel("x")
plt.ylabel("Error from quad")
plt.title("Error returned by the quad function")

plt.tight_layout()												# ensures reasonable gaps between the first subplot and the second				
plt.show()

plt.close()




###########
# Step 5

import time as t

# vectorised code
t1= t.time()
print trapezoid_integrate(0,5,0.1,function1)
t2= t.time()
print t2-t1
# code with the for trapezoid_integrate_for_loop
t1= t.time()
print trapezoid_integrate_for_loop(0,5,0.1,function1)
t2= t.time()
print t2-t1

# vectorised code
t1= t.time()
print trapezoid_integrate(0,5,0.0001,function1)
t2= t.time()
print t2-t1
# code with the for trapezoid_integrate_for_loop
t1= t.time()
print trapezoid_integrate_for_loop(0,5,0.0001,function1)
t2= t.time()
print t2-t1

#to determine h
interval_start= 0
interval_end= 1
h= 0.1
hlist=[]

tolerance= 1
integrals=[]
i=1
error=[]

integrals.append(trapezoid_integrate( interval_start, interval_end, h, function1))

while tolerance > 10**(-8):
	h= h/2
	hlist.append(h)
	integrals.append(trapezoid_integrate( interval_start, interval_end, h, function1))
	feed= np.linspace(interval_start, interval_end,((interval_end- interval_start)/h)+1)
	actual_error= max(abs(integrals[i]- np.arctan(feed)))
	error.append( [max( abs( integrals[i][::2]- integrals[i-1])), actual_error])
	i+=1
	tolerance= error[len(error)-1][0]

error= np.array(error)
plt.loglog(hlist, error[:,0],"ro",label="estimated error")
plt.loglog(hlist, error[:,1],"b+",label="actual error")
plt.legend(loc="lower right")
plt.show()
plt.close()









