#!/usr/bin/env python
import sys

def readfile(filename):
	f=open(filename,'r')
	content=[]
	for line in f:
		if len(line.strip())!=0:
			content.append(line.strip())
	return content

#reshapes a list of form [[a,b,c,d],[e,f,g,h]] into form [[a,e],[b,f],[c,g],[d,h]]
def reshape(twoDim):
	out=[]
	for i in range(len(twoDim[0])):
		out.append([twoDim[0][i],twoDim[1][i]])
	return out

fileNames=sys.argv[1:3]
print "comparing", fileNames[0], "and", fileNames[1] 

files=[]
for name in fileNames:
	files.append(readfile(name))

if len(files[0])!=len(files[1]):
	print "files differ in length:", fileNames[0], ":", len(files[0]), "lines,", fileNames[1], ":", len(files[1]), "lines"
	sys.exit(1)

files=reshape(files)

maxRelativeDeviation, maxAbsoluteDeviation = 1e-2, 1e-5
highestRelDev, sumRelDev, numComparisons = .0, .0, 0

for line in files:
	a, b = line[0].replace(":","").split(), line[1].replace(":","").split()
	threeDim=(line[0].find("z")!=-1)
	if int(a[0].split("=")[1]) != int(b[0].split("=")[1]) or int(a[1].split("=")[1]) != int(b[1].split("=")[1]) or threeDim and int(a[2].split("=")[1]) != int(b[2].split("=")[1]):
		print r"lines don't equal:", a, b
		continue
	if threeDim:
		a,b=float(a[3]), float(b[3])
	else:
		a,b=float(a[2]), float(b[2])
	if a==b:
		numComparisons += 1
		continue
	absDev = abs(a - b)
	if absDev > 1e-6:
		if a<1e-30 or b<1e-30:
			print "one value was too small", line
		relDev = max(.0, absDev / min(abs(a), abs(b)))
		highestRelDev = max(highestRelDev, relDev)
		sumRelDev += relDev
	numComparisons+=1
	if absDev > maxAbsoluteDeviation and relDev > maxRelativeDeviation:
		print "line", line, "differs by", relDev * 100, "%"
	
print "comparison of", numComparisons, "entries finished\nmaximum relative deviation was", highestRelDev * 100, "%, average relative deviation was", sumRelDev/numComparisons * 100, "%"
