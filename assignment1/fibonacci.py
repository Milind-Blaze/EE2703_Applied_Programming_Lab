# program to print the first ten fibonacci numbers

n = 1
nold = 1
print 1,nold
print 2, n
for i in range(3,11):
	new = n+nold
	nold = n
	n = new
	print i, n

	
