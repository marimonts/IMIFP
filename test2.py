#!usr/bin/python2.7

import sys 
import re
import numpy as np
import pylab as pl
#import for libraries, sys several functions, incl arguments usage

threeletter = ["ALA", "CYS", "CYX", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
oneletter = ['A', 'C', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

""" Stores a matrix from your input file, stores each line containing ifps, and stores header and ifps in strings of 7 bits """
def storeData(text):
	# Split text in lines
	mtrx = text.split('\n')
	matrix = []
	ifp = []
	header = "" #string instead of list!
	for line in mtrx[1:]:
		if "|" in line:
			ifp = []
			temp = line.split('|')
			#header.append(temp[0])
			# Convert and append header
			header = headerConv(temp[0])
			ifps = temp[1]
			ifps = ifps.strip()
			if len(ifps) % 7 == 0:
				while ifps:
					ifp.append(ifps[:7])
					ifps = ifps[7:]				
			else:
				print "IFPs no multiple of 7!"
			matrix.append([header, ifp])
		else:	
			break

# 	print matrix	
	return matrix	

"""Convert every line header to list containing a single letter aa and it's index"""
def headerConv(header):	
	head = header
	
	# Get three letter code
	aa3 = head.strip(' ')
	aa3 = aa3[-3:]
	
	number = head.strip(' ')
	number = number.strip('-')
	number = number[1:-3]

	# Replace three letter code for single letter
	if aa3 in threeletter:
		aa3 = oneletter[threeletter.index(aa3)]
	else:
		aa3 = aa3
	return aa3 + '-' + str(number)
	
""" Convert 7-bit IFP into decimal number """
def binTodec(ifp):	
	decimal = 0
	int(ifp)
	for i in range(6):  
		decimal = decimal + (int(ifp[i]) * (2**(6-i)))
	return decimal

""" Get index from user """
def getNumber(max = 2048, min = 0):
	index = "NaN"
	while(index.isdigit() == False):
		if index != "NaN":
			print "Please provide an integer value within " + str(min) + " and " + str(max)
		index = raw_input("Give an index: ").strip()
	return int(index)

IFP = 1
HEADER = 0

				
"""Begin"""

filename = sys.argv[1]
input = open(filename, "r")
text_input = input.read()

""" Split text in blocks, each containing one matrix, skip first """
text_split = text_input.split("All atom comparison\n")[1:]

""" For every block of text process and store in matrix """
mtrices = []
for mtrx in text_split:
	matrix = []
	matrix = storeData(mtrx)
	mtrices.append(matrix)

#mtrices has one element/matrix in the input. Each matrix is also a list, first element full headers, second element 7-bit ifps

""" Interaction bit masks"""
hydrophobic = 6
aromatic2 = 5
aromatic1 = 4
hydrogen2 = 3
hydrogen1 = 2
ionic2 = 1
ionic1 = 0

""" Set bitmask for specific interaction """
#mask = aromatic2

for mtrx in mtrices:
	decilines = []
	headers = []
	decimals = []
	for line in mtrx:
		for ifps in line[IFP]:
			decifps = binTodec(ifps)
			decimals.append(decifps)
			#decimals.append((decifps & (1<<hydrogen1)) >> hydrogen1)
		headers.append(line[HEADER])
		decilines.append([line[HEADER], decimals])

#i1 = getNumber(len(mtrices[0]),1)
#i2 = getNumber(len(mtrices[0]),1)

#print mtrices[0][i1-1][1][i2-1]


"""Converting to numpy arrays for plotting + heatmap plotting"""
size = len(mtrices[3])
arr = np.array(decimals).reshape(size, size)

arr_common = np.delete(arr, [range(50,100)], 0)
arr_common = np.delete(arr_common, [range(100,150)], 1)

pl.pcolor(arr_common)
pl.title('IMIFPs')
pl.colorbar()
pl.axis([0, size, 0, size])
pl.show()


input.close()
