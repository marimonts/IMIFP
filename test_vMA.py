#!usr/bin/python2.7

import sys 
import re
import argparse
import numpy as np
import pylab as pl
import math
from sets import Set
import html_generator

#import for libraries, sys several functions, incl arguments usage

threeletter = ["ALA", "CYS", "CYX", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
oneletter = ['A', 'C', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

""" Stores a matrix from your input file, stores each line containing ifps, and stores header and ifps in strings of 7 bits """
def storeData(text):
	# Split text in lines
	mtrx = text.split('\n')
	matrix = []
	header = "" #string instead of list!
	for line in mtrx[1:]:
		ifp = []
		if "|" in line:
			temp = line.split('|')
			# Convert and append header

#			content = headerConv(temp[0])
#			aa1 = content[0]
#			number = content[1]
#			header = aa1 + '-' + str(number)

			aa, number = headerConv(temp[0])

			ifps = temp[1]
			ifps = ifps.strip()
			if len(ifps) % 7 == 0:
				while ifps:
					ifp.append(ifps[:7])
					ifps = ifps[7:]				
			else:
				print "IFPs no multiple of 7!"
			matrix.append([aa, number, ifp])
		else:	
			break
	
	return matrix	

"""Convert every line header to list containing a single letter aa and it's index"""
def headerConv(header):	
	
	# Get aminoacid letter code
	aa = header.strip(' ')
	aa = aa[-3:]

	# Get residue index from the header	
	number = header.strip(' ')
	number = number.strip('-')
	number = number[1:-3]

	# Replace three letter code for single letter
	if aa in threeletter:
		aa = oneletter[threeletter.index(aa)]

	return aa, number

def headerConv_short(header):	
	aa1 = header[0:1]
	number = int(header[-3:])

	return aa1, number
	
""" Convert 7-bit IFP into decimal number """
def binTodec(ifp):	
	decimal = 0
	int(ifp)
	for i in range(7):  
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


def deleteNonCommon(positions, aa, only_data):
	# Find positions which are not in all matrices
	total = set(positions[0])
	only_curated_data = []
	for i in range(1, len(positions)):
		total = total | set(positions[i])
	uncommon = set()
	for s in positions:
		uncommon.update(total ^ set(s))
	for i in range(len(positions)):
		for element in uncommon:
			if element in positions[i]:
				idx = positions[i].index(element)
				del positions[i][idx]
				del aa[i][idx]
				deleteX = np.delete(only_data[i], [idx] , 0)
				only_data[i] = np.delete(deleteX, [idx] , 1)
		only_curated_data.append(only_data[i])
	return list(positions), list(aa), only_curated_data

# Returns how many bits are set in the binary representation of a number of any length
def getSetBits(val):
	n = 0
	while(val > 0):
		n+=val&1
		val = val>>1

	return n

# Returns the tanimoto score of two ifps (default:0x7f takes all bits into acount)
def getTanimoto(ifp1, ifp2, mask=0x7f):
	a_sum = 0
	b_sum = 0
	c_sum = 0

	# For every ifp in the provided lists:
	for i in range(1, len(ifp1)):
		for j in range(0, i):
			# Perform bitmask operation
			a = ifp1[i][j] & mask
			b = ifp2[i][j] & mask

			# Add new values
			a_sum += getSetBits( a )
			b_sum += getSetBits( b )
			c_sum += getSetBits( a & b )	# Determine number of bits set (1) in both ifps
			
	# Catch division by zero errors if both ifps are 0!
	try:
		tanimoto = float(c_sum) / abs( a_sum + b_sum - c_sum )
	except ZeroDivisionError:
		# Make tanimoto 1.0 as ifps are both 0 and thus equal
		tanimoto = 0.0

	return tanimoto

			
				
"""Begin"""

filename = sys.argv[1]
input = open(filename, "r")
output_ifps = "IMIFP_output_ifps.txt"
text_input = input.read()
output1 = open(output_ifps, "w")


""" Split text in blocks, each containing one matrix, skip first """
text_split = text_input.split("All atom comparison\n")[1:]

""" For every block of text process and store in matrix """
mtrices = []
for mtrx in text_split:
	mtrices.append( storeData(mtrx) )

#mtrices has one element/matrix in the input. Each matrix is also a list, first element full headers, second element 7-bit ifps

only_data = []
aa_list = []
number_list = []
for mtrx in mtrices:
	decimals = []
	aa = []
	number = []
	for line in mtrx:
		aa.append(line[0])
		number.append(line[1])
		print >> output1, '\t' + line[0] + line[1],
	aa_list.append(aa)
	number_list.append(number)
	print >> output1, '\n',
	for line in mtrx:
		print >> output1, line[0]+ line[1],
		for ifps in line[2]:
			decifps = binTodec(ifps)
			decimals.append(decifps)
			print >> output1, '\t' + ifps,	
		print >> output1, '\n',
	size = len(mtrx)
	only_data.append(np.array(decimals).reshape(size,size))

positions_curated, aa_curated, only_curated_data = deleteNonCommon(number_list, aa_list, only_data)

print positions_curated	
print getTanimoto(only_curated_data[0], only_curated_data[1])
# For every ifp in the provided lists:
#print only_curated_data[0]
#print only_curated_data[2]

#print getTanimoto(only_curated_data[0], only_curated_data[2])

#print getTanimoto(only_curated_data[0], only_curated_data[1])

"""Converting to numpy arrays for plotting + heatmap plotting"""

# arr = np.array(decimals).reshape(size, size)

# pl.pcolor(arr)
# pl.title('IMIFPs')
# pl.colorbar()
# pl.axis([0, size, 0, size])
# pl.show()

#Creating html output file for data visualization
# html_file = open("ifp_output.html", "w")

# html_file.write("<html>\n\r<heading><title>IFP</title></heading>\n\r<font face='verdana'><body><center>\n\r")
# html_file.write("<h1>Interaction Fingerprint Viewer</h1>\n\r")
# html_file.write("<table border='0' cellpadding='1'>\n\r")

# # n = number of matrices	
# n = len(only_curated_data)

# hydrophobic = 0b1000000
# aromatic = 0b0110000
# hbond = 0b0001100
# ionic = 0b0000011
# mask = hydrophobic + aromatic + hbond + ionic

# # IFP matrices are square so rows = columns
# for row in range(len(only_curated_data[0])):
# 	for col in range(len(only_curated_data[0][row])):
	
# 		cell = []
# 		header = []
# 		for i in range(n):
# 			cell.append(only_curated_data[i][row][col] & mask)
			
# 		if sum(cell) > 0:
# 			print str(positions_curated[0][row]) + "\t" + str(positions_curated[0][col]) + "\t",
# 			for i in range(len(cell)):
# 				print aa_curated[i][row] + "-" + aa_curated[i][col] + "\t" + str(cell[i]) + "\t",
# 				header.append(aa_curated[i][row] + "-" + aa_curated[i][col])
# 			print " "

# 			html = html_generator.ifpTable(cell, header)
# 			html = html_generator.makeTableRow( [str(positions_curated[0][row]) + "-" + str(positions_curated[0][col]), html] )
# 			html_file.write(html)


# html_file.write("</table>")
# html_file.write("</body></font>\n</html>")

# html_file.close()

output1.close()	


#print mtrices[0][i1-1][1][i2-1]



input.close()
