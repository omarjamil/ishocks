#!/disks/jets/cosmic/o.jamil/software/python/bin/ env python

from scipy import *
from pylab import *

print "Logarithmic/Geometric re-binning program:", '\n'

# var = int(raw_input("Enter 1 for logarithmic and 2 for geometric: "))
# #if var == 2:
# #    geo = float(raw_input("Enter the geometric ratio: "))
# infile = raw_input("Enter the input file name: ")
# oufile = raw_input("Enter the output file name: ")
# binNums = float(raw_input("Enter the number of bins: "))
# xLoBound = float(raw_input("Enter the x-lower bound: "))
# xUpBound = float(raw_input("Enter the x-upper bound: "))


parArray=[]
inpar = open ("rebPar.par","r")
for line in inpar:
    columns = line.split("=")
    parArray.append( columns[1] )

var = int(parArray[0])
infile = parArray[1]
infile = infile.replace("\n","")
oufile = parArray[2]
oufile = oufile.replace("\n","")
binNums = float(parArray[3])
xLoBound = float(parArray[4])
xUpBound = float(parArray[5])
    
Y=[]
X=[]
inp = open (infile,"r")
for line in inp:
    columns = line.split("\t")
    X.append( float(columns[0]) )
    Y.append( float(columns[1]) )

#temp for adding up the values that fall in the bin selected
temp=[]
#For logarithmic spacing
if var == 1:
    X1=logspace(log10(xLoBound), log10(xUpBound), binNums)
#For geometric series
if var == 2:
    geo=(xUpBound/xLoBound)**(1./(binNums-1.))
    X1=[]
    for g in range(1, int(binNums)+1, 1):
        X1.append(xLoBound*(geo**(g-1.)))

Y1=[]
count = 1.0
sum = 0.0
prevSum = 0.0

#For logarithmic progression
if var == 1:
    outfile=open(oufile,'wr')
    for i in range(size(X1)):
        for k in range(size(Y)):
            if i==size(X1)-1:
                temp.append( Y[k] )
                sum += Y[k]
                count += 1.0
            elif i < size(X1)-1:
                if X[k] >= X1[i] and X[k] < X1[i+1]:
                    temp.append( Y[k] )
                    sum += Y[k]
                    count += 1.0
                else:
                    temp.append( 0.0 )
        Y1.append( (sum/count) )
        #prevSum = sum
        sum = 0.0
        count = 1.0
        del temp[:]

    for m in range(size(X1)):
        if Y1[m] != 0.0:
            s1=str(X1[m])
            s2=str(Y1[m])
            outfile.write(s1+"\t"+s2+"\n")


#For geometric progression
if var == 2:
    outfile=open(oufile,'wr')
    for i in range(size(X1)):
        for k in range(size(Y)):
            if i==size(X1)-1:
                temp.append( Y[k] )
                sum += Y[k]
                count += 1.0
            elif i < size(X1)-1:    
                if X[k] >= X1[i] and X[k] < X1[i+1]:
                    temp.append( Y[k] )
                    sum += Y[k]
                    count += 1.0
                else:
                    temp.append( 0.0 )
        Y1.append( (sum/count) )
        #prevSum = sum
        sum = 0.0
        count = 1.0
        del temp[:]

    #output the x and y values
    for m in range(size(X1)):
        if Y1[m] != 0.0: #this does not output the zero values
            s1=str(X1[m])
            s2=str(Y1[m])
            outfile.write(s1+"\t"+s2+"\n")

inp.close()
outfile.close()
               
