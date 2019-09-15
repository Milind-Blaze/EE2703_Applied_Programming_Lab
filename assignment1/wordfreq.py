# program to read in the words of a given file and then determine the word frequency

import sys
import string

inputfile= sys.argv[1]



with open(inputfile, 'r') as f:
	file= f.read()
	exclude = set(string.punctuation)
	file = ''.join(ch for ch in file if ch not in exclude)
	words= file.strip().split()
	wordcount= {}
	for word in words:
		if word in wordcount:
			wordcount[word] +=1
		else:
			wordcount[word] = 1

##for word in wordcount:
##	print word, wordcount[word]

#for word in sorted(wordcount, key= wordcount.get):			# key argument tells sorted the basis on which its 	
#	print word, wordcount[word]								# to sort wordcount

for word in sorted(wordcount, key= wordcount.get, reverse =True):			# key argument tells sorted the basis on which its 	
	print word, wordcount[word]												# to sort wordcount

y=[]
x= sorted(wordcount, key= wordcount.get, reverse= True)[:20]
for word in x:
	y.append(int(wordcount[word]))

#print x,y

import matplotlib.pyplot as plt 

#plt.hist(y,range(20))
mybar = plt.bar(range(len(x)), y, color='green', alpha=0.4)
plt.xticks(range(len(x)),x)
plt.title("word frequency")
plt.xlabel("word")
plt.ylabel("frequency")
plt.show()

