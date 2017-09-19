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

def preprocessMatrix(text):
	# Split text in lines
	mtrx = text.split('\n')
	
	positions = []
	aminoacids = []
	ifps = []

	header = ""
	for line in mtrx[1:]:		# Skip the first line (header)
		if "|" in line:
			temp = line.split('|')

			# Append the aminoacid and position to their respective list
			aa, position = headerConv(temp[0])
			positions.append(position)
			aminoacids.append(aa)

			# Convert an entire line ifp bitstrings to numbers and append to list
			ifp_row = []
			ifp_line = temp[1].strip()
			if len(ifp_line) % 7 == 0:
				while ifp_line:
					ifp_row.append( int(ifp_line[:7],2) )	# Convert ifp bitstring to decimal number
					ifp_line = ifp_line[7:]					# Delete processed from the line
				# Save the this row of ifps as list of lists
				ifps.append(ifp_row)	
			else:
				print "IFPs no multiple of 7!"	
		else:
			break

	return list([aminoacids, positions, ifps])

"""Convert every line header to list containing a single letter aa and it's index"""
def headerConv(header):	
	
	# Get aminoacid letter code
	aa = header.strip(' ')
	aa = aa[-3:]

	# Get residue index from the header	
	number = header.strip(' ')
	number = number.strip('-')
	number = int(number[1:-3])

	# Replace three letter code for single letter
	if aa in threeletter:
		aa = oneletter[threeletter.index(aa)]

	return aa, number

# Creates a sorted list with the positions that are common in all matrices
def getCommonPositions(matrix_container):
	# Index of the positions element in the matrices
	pos = 1

	common = set(matrix_container[0][pos])			# Init common with first matrix positions
	for matrix in matrix_container[1:]:	# Go over every matrix, starting from the second
		common = common & set(matrix[pos])

	return sorted(list(common))

# Returns how many bits are set in the binary representation of a number of any length
def getSetBits(val):
	n = 0
	while(val > 0):
		n+=val&1
		val = val>>1

	return n

def similarityMatrix(matrices, positions, mask=0x7f):
	# Determine the number of matrices
	n_matrices = len(matrices)

	# Create list lists of dimension n_matrices initialized to 1.0
	similarity_matrix = []
	for row in range(0, n_matrices):
		similarity_matrix.append( [1.0 for i in range(n_matrices)] )

	# For every unique matrix combination in the provided list of matrices calculate Tc 
	for row in range(1, n_matrices):
		for col in range(0, row):
			score = getTanimoto(matrices[row], matrices[col], positions, mask)
			similarity_matrix[row][col] = score
			similarity_matrix[col][row] = score

	return similarity_matrix

# Returns the tanimoto score between two matrices (default:0x7f takes all bits into acount)
def getTanimoto(m1, m2, positions, mask=0x7f):
	a_sum = 0
	b_sum = 0
	c_sum = 0

	# Perform once for every unique interaction between the matrices:
	for i in range(1, len(positions)):
		for j in range(0, i):
			# Select the interaction between positions at indeces i and j
			p1 = positions[i]
			p2 = positions[j]

			# Obtain the indeces of both matrices that refer to the positions p1 and p2 in each
			i1_m1 = m1[1].index(p1)
			i2_m1 = m1[1].index(p2)
			i1_m2 = m2[1].index(p1)
			i2_m2 = m2[1].index(p2)

			# Get the ifp values of both matrices at positions p1 and p2 and perform bitmask
			a = m1[2][i1_m1][i2_m1] & mask
			b = m2[2][i1_m2][i2_m2] & mask

			# Obtain the set number of bits in each bitstring and sum
			a_sum += getSetBits( a )
			b_sum += getSetBits( b )
			c_sum += getSetBits( a & b )	# Determine number of bits set (1) in both ifps
	

	# Catch division by zero errors if both ifp matrices are 0!
	try:
		tanimoto = float(c_sum) / abs( a_sum + b_sum - c_sum )
	except ZeroDivisionError:
		# Make tanimoto 0.0 as ifp matrices are both 0, no similarity
		tanimoto = 0.0
	return tanimoto

def similarityMatrixPlot(tanimoto_similarity, filename, filepath=""):
	size = len(tanimoto_similarity[0])
	arr = np.array(tanimoto_similarity).reshape(size, size)
	pl.pcolor(arr)
	pl.title('Tc Similarity Matrix')
	pl.colorbar()
	pl.axis([1, size, size, 1])
	# pl.set_aspect(1)
	# width, height = arr.shape
	# for x in xrange(width):
	#     for y in xrange(height):
	#         pl.annotate(round(arr[x][y], 2), xy=(y, x))
	pl.show()	
	pl.savefig(filepath+filename)
	pl.clf()
	fig.savefig(filepath+filename)

def cutMatrix(matrix, positions):
	indeces_uncommon = []
	indeces_common = []
	uncommon = set(matrix[1]) - set(positions) # get a list of uncommon positions
	for p in uncommon:
		indeces_uncommon.append(matrix[1].index(p)) # Convert uncommon positions to corresponding indeces
	for q in positions:
		indeces_common.append(matrix[1].index(q)) # Convert common positions to corresponding indeces
	residues_curated = []
	for i in indeces_common:
		residues_curated.append(matrix[0][i])

	# Convert to a numpy array
	size = len(matrix[1])
	matrix = np.array(matrix[2]).reshape(size,size)
	matrix = np.delete(matrix, indeces, 0) # delete rows of uncommon indeces
	matrix = np.delete(matrix, indeces, 1) # delete columns of uncommon indeces

	return matrix, residues_curated

def frequencyMatrix(matrices, positions, mask=0x7f, exactFilter="TRUE"):
	# Determine the number of positions
	size = len(matrices[0][2])
	sumarrays = 0
	total = len(matrices)
	for matrix in matrices:
		arr = cutMatrix(matrix, positions)[0]
		if(exactFilter == "TRUE"):
			arr[(arr & mask) != mask] = 0
		else:
			arr[(arr & mask) == 0] = 0
		arr[arr > 0] = 1
		sumarrays = sumarrays + arr
	sumarrays = (sumarrays / (float(total))*100)
	return sumarrays.round(2)

def frequencyMatrixPlot(matrix, positions, filename, filepath=""):
	fig, ax = pl.subplots()

	step = len(positions) / 5
	indeces = range(0, len(positions), step)
	for i in range(0, len(indeces)):
		indeces[i] = positions[indeces[i]]

	ax.set_xticklabels(indeces)
	ax.set_yticklabels(indeces)
	pl.pcolor(matrix)
	pl.title('Frequency Matrix')
	pl.colorbar()
#	pl.axis([positions[0], positions[-1], positions[0], positions[-1]])
# pl.set_aspect(1)
# width, height = arr.shape
# for x in xrange(width):
#     for y in xrange(height):
#         pl.annotate(round(arr[x][y], 2), xy=(y, x))
	pl.show()
	pl.savefig(filepath+filename)
	pl.clf()

	fig.savefig(filepath+filename)

				
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
	mtrices.append( preprocessMatrix(mtrx) )

# mtrices has one element/matrix in the input. Each matrix is also a list, first element AA, second element number, third element 7-bit ifps

common = getCommonPositions(mtrices)
# tanimoto_similarity = similarityMatrix(mtrices, common, 0x7f)
# similarityMatrixPlot(tanimoto_similarity, "Tanimoto_similarity_matrix.png")
# FM = frequencyMatrix(mtrices, common, mask=0x7f, exactFilter="FALSE")
# frequencyMatrixPlot(FM, common, "frequency_matrix.png")
for matrix in mtrices:
	temp = cutMatrix(matrix, common)
	"""UGLYYYYYYYYYYY"""

#Creating html output file for data visualization
html_file = open("ifp_output.html", "w")

html_file.write("<html>\n\r<heading><title>IFP</title></heading>\n\r<font face='verdana'><body><center>\n\r")
html_file.write("<h1>Interaction Fingerprint Viewer</h1>\n\r")
html_file.write("<table border='0' cellpadding='1'>\n\r")

# n = number of matrices	
n = len(mtrices)

hydrophobic = 0b1000000
aromatic = 0b0110000
hbond = 0b0001100
ionic = 0b0000011
mask = hydrophobic + aromatic + hbond + ionic

# IFP matrices are square so rows = columns
for row in range(len(only_curated_data[0])):
	for col in range(len(only_curated_data[0][row])):
	
		cell = []
		header = []
		for i in range(n):
			cell.append(only_curated_data[i][row][col] & mask)
			
		if sum(cell) > 0:
			for i in range(len(cell)):
				header.append(aa_curated[i][row] + "-" + aa_curated[i][col])

			html = html_generator.ifpTable(cell, header)
			html = html_generator.makeTableRow( [str(positions_curated[0][row]) + "-" + str(positions_curated[0][col]), html] )
			html_file.write(html)


html_file.write("</table>")
html_file.write("</body></font>\n</html>")

html_file.close()

output1.close()	

input.close()
