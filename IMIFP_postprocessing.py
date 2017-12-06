#!usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys 
import re
import argparse
import numpy as np
import pylab as pl
import math
from sets import Set
import html_generator
from sklearn.cluster import KMeans
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectPercentile
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
	for line in mtrx[2:]:		# Skip the first two lines (title and header)
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
					rev_ifp = str(ifp_line[:7])[::-1]
					unified_ifp = (rev_ifp[0] or rev_ifp[1]) + (rev_ifp[2] or rev_ifp[3]) + (rev_ifp[4] or rev_ifp[5]) + rev_ifp[6]
					print int(unified_ifp,2)
					ifp_row.append( int(unified_ifp,2) )	# Convert ifp bitstring to decimal number
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

# def similarityMatrix(matrices, positions, mask=0x7f):
# 	# Determine the number of matrices
# 	n_matrices = len(matrices)

# 	# Create list lists of dimension n_matrices initialized to 1.0
# 	similarity_matrix = []
# 	for row in range(0, n_matrices):
# 		similarity_matrix.append( [1.0 for i in range(n_matrices)] )

# 	# For every unique matrix combination in the provided list of matrices calculate Tc 
# 	for row in range(1, n_matrices):
# 		for col in range(0, row):
# 			score = getTanimoto(matrices[row], matrices[col], positions, mask)
# 			print matrices[row]
# 			similarity_matrix[row][col] = score
# 			similarity_matrix[col][row] = score

# 	return similarity_matrix

# # Returns the tanimoto score between two matrices (default:0x7f takes all bits into acount)
# def getTanimoto(m1, m2, positions, mask=0x7f):
# 	a_sum = 0
# 	b_sum = 0
# 	c_sum = 0

# 	# Perform once for every unique interaction between the matrices:
# 	for i in range(1, len(positions)):
# 		for j in range(0, i):
# 			# Select the interaction between positions at indeces i and j
# 			p1 = positions[i]
# 			p2 = positions[j]

# 			# Obtain the indeces of both matrices that refer to the positions p1 and p2 in each
# 			i1_m1 = m1[1].index(p1)
# 			i2_m1 = m1[1].index(p2)
# 			i1_m2 = m2[1].index(p1)
# 			i2_m2 = m2[1].index(p2)

# 			# Get the ifp values of both matrices at positions p1 and p2 and perform bitmask
# 			a = m1[2][i1_m1][i2_m1] & mask
# 			b = m2[2][i1_m2][i2_m2] & mask

# 			# Obtain the set number of bits in each bitstring and sum
# 			a_sum += getSetBits( a )
# 			b_sum += getSetBits( b )
# 			c_sum += getSetBits( a & b )	# Determine number of bits set (1) in both ifps
	

# 	# Catch division by zero errors if both ifp matrices are 0!
# 	try:
# 		tanimoto = float(c_sum) / abs( a_sum + b_sum - c_sum )
# 	except ZeroDivisionError:
# 		# Make tanimoto 0.0 as ifp matrices are both 0, no similarity
# 		tanimoto = 0.0
# 	return tanimoto

def similarityMatrixPlot(tanimoto_similarity, filename, filepath=""):
	size = len(tanimoto_similarity[0])
	arr = np.array(tanimoto_similarity).reshape(size, size)
	print arr
	pl.pcolor(arr)
	pl.title('Tc Similarity Matrix')
	pl.colorbar()
	pl.axis([1, size, size, 1])
	# pl.set_aspect(1)
	# width, height = arr.shape
	# for x in xrange(width):
	#     for y in xrange(height):
	#         pl.annotate(round(arr[x][y], 2), xy=(y, x))
	pl.savefig(filepath+filename)
	#pl.clf()

def cutMatrixAA(matrix, positions):
	indeces_common = []
	
	for q in positions:
		indeces_common.append(matrix[1].index(q)) # Convert common positions to corresponding indeces
	residues_curated = []
	for i in indeces_common:
		residues_curated.append(matrix[0][i])

	return residues_curated

def cutMatrixIFP(matrix, positions):
	indeces = []
	uncommon = set(matrix[1]) - set(positions) # get a list of uncommon positions
	for p in uncommon:
		indeces.append(matrix[1].index(p)) # Convert uncommon positions to corresponding indeces
	# Convert to a numpy array
	size = len(matrix[1]) 
	matrix = np.array(matrix[2]).reshape(size,size)
	matrix = np.delete(matrix, indeces, 0) # delete rows of uncommon indeces
	matrix = np.delete(matrix, indeces, 1) # delete columns of uncommon indeces

	return matrix

def frequencyMatrix(matrices, positions, mask=0x7f, exactFilter="TRUE"):
	sumarrays = 0
	total = len(matrices)
	for matrix in matrices:
		arr = cutMatrixIFP(matrix, positions)
		if(exactFilter == "TRUE"):
			arr[(arr & mask) != mask] = 0
		else:
			arr[(arr & mask) == 0] = 0
		arr[arr > 0] = 1
		sumarrays = sumarrays + arr
		break
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
	pl.savefig(filepath+filename)
	#pl.clf()
def getTanimotoSpecificInteractions(m1, m2, positions, list_interactions, mask=0x7f):
	a_sum = 0
	b_sum = 0
	c_sum = 0

	interactions = []
	#from the list given by user, make a list with the pair split 
	for item in list_interactions:
		item = item.split('-')
		interactions.append(item)


	# Perform once for every unique interaction between the matrices:
	i = 0
	for interaction in list_interactions:
		# Select the interaction between positions in list interactions
		p1 = int(interactions[i][0])
		p2 = int(interactions[i][1])
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
		i = i+1	

	# Catch division by zero errors if both ifp matrices are 0!
	try:
		tanimoto = float(c_sum) / abs( a_sum + b_sum - c_sum )
	except ZeroDivisionError:
		# Make tanimoto 1.0 as ifp matrices are both 0
		tanimoto = 1.0
	return tanimoto

def clusterXinteraction(matrices, positions, list_interactions, mask=0x7f, exactFilter="TRUE"):
	# Determine the number of matrices
	n_matrices = len(matrices)

	# Create list lists of dimension n_matrices initialized to 1.0
	similarity_matrix = []
	for row in range(0, n_matrices):
		similarity_matrix.append( [1.0 for i in range(n_matrices)] )

	# For every unique matrix combination in the provided list of matrices calculate Tc 
	for row in range(1, n_matrices):
		for col in range(0, row):
			score = getTanimotoSpecificInteractions(matrices[row], matrices[col], positions, list_interactions, mask) #do it only in the chosen interaction pairs
			similarity_matrix[row][col] = score
			similarity_matrix[col][row] = score

	return similarity_matrix


def conservedOut(matrices, positions, mask=0x7f):
	bin_matrices = []
	bin_tmp = []
	pos_x = []
	pos_y =  []
	bin_arrays = []
	# formatting each matrix and storing it in bin with only IFPs
	for matrix in matrices:
		bin_matrices.append(cutMatrixIFP(matrix, positions))
	n = len(bin_matrices)
	# change data ordering from receptor-based grouped to interaction pair grouped
	for row in range(len(bin_matrices[0])):
		for col in range(len(bin_matrices[0][row])):
			cell = []
			
			for i in range(n):
				cell.append(bin_matrices[i][row][col] & mask)
			bin_tmp.append(cell)
			pos_x.append(positions[row])
			pos_y.append(positions[col])

	# merge everything together by making a list of lists, in which every sublist contains the same element position of the original lists, e.g. ([64, 64, 0, 64, 64...], 332, 336)	

	# zip positions x and y together (zip is NOT append, each element is zipped together with the element of the same index in the other list(s))
	interactions = zip(pos_x, pos_y)
	# create two separate arrays and merge them
	bin_arrays_positions = np.array(interactions)
	bin_arrays_ifps = np.array(bin_tmp)
	bin_nparray = np.concatenate((bin_arrays_ifps, bin_arrays_positions), axis=1)
	# "delete" fully conserved interactions by summing up the values of the ifps / amount of receptors, it checks if the result equals the first value (which happens when the interaction is fully conserved)
	# np.delete mixes up in iterative mode, therefore I copy the non-conserved into a list, and back to nparray > alt: store index of rows to delete and give the list to np.delete
	tmp_list = []
	for row in bin_nparray:
		if sum(row[:-2])/n != row[0]:
			tmp_list.append(row)
	bin_nparray = np.array(tmp_list)
	np.savetxt(output2, bin_nparray, fmt='%d',delimiter='	', newline='\n')
	return bin_nparray 

def similarityMatrix(bin_nparray, mask=0x7f):
	# Determine the number of matrices, which equals number of columns in the array -2 coordinates
	n_matrices = (np.shape(bin_nparray)[1] ) - 2

	# Create list lists of dimension n_matrices initialized to 1.0
	similarity_matrix = []
	for row in range(0, n_matrices):
		similarity_matrix.append( [1.0 for i in range(n_matrices)] )
	
	# For every unique matrix combination in the provided list of matrices calculate Tc 
	for row in range(1, n_matrices):
		for col in range(0, row):
			score = getTanimoto(bin_nparray[:,row], bin_nparray[:,col], mask)
			# to make it square, mirror images:
			similarity_matrix[row][col] = score
			similarity_matrix[col][row] = score

	return similarity_matrix

# Returns the tanimoto score between two matrices (default:0x7f takes all bits into acount)
def getTanimoto(m1, m2, mask=0x7f):
	a_sum = 0
	b_sum = 0
	c_sum = 0

	# Perform once for every unique interaction between the matrices:
	for i in range(0, len(m1)):

			# Get the ifp values of both matrices at position i and perform bitmask
			a = m1[i] & mask
			b = m2[i] & mask

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

def npTanimotoXgroup(npmatrix, mask=0x7f):
	nr_col_matrix = np.shape(npmatrix)[1]
	column_Tcs = []
 	for row in npmatrix:
 		row = row[1:]
 		Tcs = []
 		len_row = len(row)
 		for i in range(1, len_row):
 			for j in range (0, i):

			# Get the ifp values of both matrices at position i and perform bitmask
				a = row[i] & mask
				b = row[j] & mask
 				a_sum = 0
				b_sum = 0
				c_sum = 0
				a_sum += getSetBits( a )
				b_sum += getSetBits( b )
				c_sum += getSetBits( a & b )
				try:
					tanimoto = float(c_sum) / abs( a_sum + b_sum - c_sum )
				except ZeroDivisionError:
					# Make tanimoto 0.0 as ifp matrices are both 0, no similarity
					tanimoto = 0.0
				Tcs.append(tanimoto)
		
		column_Tcs.append(sum(Tcs))
	npmatrix = np.insert(npmatrix, nr_col_matrix, column_Tcs, axis=1)	

	return npmatrix

def featureSelection(bin_nparray, mask=0x7f):
 	nr_columns = np.shape(bin_nparray)[1] #get the number of columns in array to automatically split the array in subarrays
 	list_columns = np.split(bin_nparray, range(1, nr_columns), axis =1) #creates one list in which each element is one column of the original array
 	
 	# create a unique column for the positions with a int format (residues contained in two last columns)
 	positionX = []
 	for x in list_columns[-2]:
 		positionX.append(x)
 	positionY = []
  	for y in list_columns[-1]:
 		positionY.append(y)	
 	interacting_positions = []
 	for i in range(len(positionX)): 
 			interacting_positions.append((str(positionX[i][0])+ str(positionY[i][0])))
 	
 	#merging the columns that cluster together... HARD CODING! 
 	inactives = np.concatenate((list_columns[0:5]), axis=1)
 	ni = np.shape(inactives)[1] # how many receptors in the group for future parameters
 	actives = np.concatenate((list_columns[5:10]), axis=1)
 	na = np.shape(actives)[1] # how many receptors in the group for future parameters
 	test_set = np.concatenate((list_columns[0:10]), axis=1)
	nt = np.shape(test_set)[1] # how many receptors in the group for future parameters

 	inactives = np.insert(inactives, 0, interacting_positions, axis =1)
 	actives = np.insert(actives, 0, interacting_positions, axis =1)
 	test_set = np.insert(test_set, 0, interacting_positions, axis=1)

 	#Adding a column to the end with the sum of all Tc possible combination
 	inactives = npTanimotoXgroup(inactives)
 	actives = npTanimotoXgroup(actives)
 	
 	#sort each array based on the Tc sum and store it back
 	tmp = np.argsort(inactives[:,(ni+1)])[::-1] # ni + 1 represents the last column after the IFPs of this group
 	tmp2 = np.argsort(actives[:,(na+1)])[::-1] # na + 1 represents the last column after the IFPs of this group
 	inactives = inactives[tmp]
 	actives = actives[tmp2]

 	# sum IFP values (contained in [0:5] and insert a new column in the end with the result) AND add more value if the IFP has a polar bit (0bXXX1111)
	sumi = ni+2 # will represent the index of the column of Tc sum of this group
	sum_inactives = []
	polarmask = 0b0001111 # select here mask to assign extra value to the IFP, 0b0001111 = extra weight to any polar interaction
	for row in inactives:
		counter = 0
		#for item in row[1:-1]:
		#	if (item & polarmask) > 0:
		#		counter = counter + 40 # assign here how much extra weight to give to polar interactions
		sum_R = sum(row[1:-1])
		sum_inactives.append((sum_R + counter)/ni)
	sum_inactives = np.array(sum_inactives)
	inactives = np.insert(inactives, sumi, sum_inactives, axis=1)

	suma = na+2 # will represent the index of the column of Tc sum of this group
	sum_actives = []
	for row in actives:
		counter = 0
		#for item in row[1:-1]:
		#	if (item & polarmask) > 0:
		#		counter =+ 40 # assign here how much extra weight to give to polar interactions
		sum_R = sum(row[1:-1])
		sum_actives.append((sum_R + counter)/na)
	sum_actives = np.array(sum_actives)
	actives = np.insert(actives, suma, sum_actives, axis=1)
	
	# creating files with all inactive and inactive ordered based on the similarity per group, calculate difference of the Tcs sum
	commoninactives = []
	commonactives = []
	sum_difference = []
	for rowinactives in inactives:
		for rowactives in actives:
			if rowinactives[0] == rowactives[0]: #use index of interaction positions (always first column)
				sum_difference.append(abs(rowinactives[sumi] - rowactives[suma])) #use indeces of IFP sums
				commoninactives.append(rowinactives)
				commonactives.append(rowactives)
	sum_difference = np.array(sum_difference)
	commoninactives = np.array(commoninactives)
	commonactives = np.array(commonactives)
	np.savetxt(output3, commoninactives, fmt='%d',delimiter='	', newline='\n')
	np.savetxt(output4, commonactives, fmt='%d',delimiter='	', newline='\n')
	
	lastcol = np.shape(commoninactives)[1]
	commoninactives = np.insert(commoninactives, lastcol, sum_difference, axis=1) # add difference of IFP sums at the end

	# creating an array+txt file with IFPS of a group, their interaction positions, Tcsum, IFPs of the other group, their Tcsum, and IFP sum difference
	whole_set = np.array(commoninactives)
	whole_set = np.delete(whole_set, sumi, axis=1) #delete sum IFPs first group

	commonactives = np.delete(commonactives, 0, axis=1)
	for column in range((np.shape(commonactives)[1])-1): # range corresponding to number of receptors = columns -1 (last columns are not IFPs)
		whole_set = np.insert(whole_set, (sumi + column), commonactives[:,column], axis=1) # add receptors in position (sumi + column), so number of columns plus iterator, and in the end add sum Tcs group2

	#delete redundancy (delete rows that represent the same interaction)(e.g. 349350 = 350349)
	pos_todelete = []
	for pos in whole_set[:,0]: #use index of column interaction positions
		pos = str(pos)
		tmp1 = pos[0:3] 
		tmp2 = pos[-3:]
		if tmp1 > tmp2:
			pos_todelete.append(int(pos))
	
	clean_whole_set = []
	for row in whole_set:
		if row[0] in pos_todelete: #use index of interaction positions
			pass
		else:
			clean_whole_set.append(row)
	clean_whole_set = np.array(clean_whole_set)	

	#delete redundancy (delete rows that represent the same interaction)(e.g. 349350 = 350349)
	test_set_clean = []
	for row in test_set:
		if row[0] in pos_todelete:
			pass
		else:
			test_set_clean.append(row)
	test_set_clean = np.array(test_set_clean)
	
	np.savetxt(output5, clean_whole_set, fmt='%d',delimiter='	', newline='\n') # create tab separated text file 

	return clean_whole_set, test_set_clean, ni
#for filtering out Tc sums lower than a certain cutoff and save in txt files
	#rowfilter_inactives, rowfilter_actives = TcSumFilter(inactives, actives)

def TcSumFilter(arr1, arr2):
	#returns arrays (one x group) after filtering out rows under a certain cutoff of Tcsum
	rowfilter_arr1 = []
	for row in arr1:
		if row[6] >= 6:## hard coding, change if used 
			rowfilter_arr1.append(row)
	rowfilter_arr1 = np.array(rowfilter_arr1) 
	
	rowfilter_arr2 = []
	for row in arr2:
		if row[6] >= 6:## hard coding, change if used 
			rowfilter_arr2.append(row)
	rowfilter_arr2 = np.array(rowfilter_arr2)

	return rowfilter_arr1, rowfilter_arr2 

def clustering(clean_whole_set, test_set_clean, ni, nr_features):


	sorted_set = clean_whole_set[clean_whole_set[:,-1].argsort()[::-1]] #sort entries based on sum_differences

	#store the label of the features to select then in the test set
	selected_features = []
	for row in sorted_set[0:nr_features]: # how many features to consider with best sum_differences score
		selected_features.append(row[0]) # row[0] --> interaction positions label

	top50 = sorted_set[0:nr_features] # (*) get ranking features with best sum_differences score 
	shape_col =  np.shape(top50)[1]
	top50 = np.delete(top50, (shape_col-1, shape_col-2), axis=1) # make table suitable for Kmeans, delete non IFP data columns
	top50 = np.delete(top50, 0, axis=1) # make table suitable for Kmeans, delete non IFP data columns
	top50 = np.delete(top50, ni, axis=1) # make table suitable for Kmeans, delete non IFP data columns

	training_set = np.array(top50)

		# select only data of the selected features in the test set for Kmeans comparison 
	test_set_selectedFeatures = []
	for item in selected_features:
		for row in test_set_clean:
			if row[0] == item:
				test_set_selectedFeatures.append(row)
			else:
				pass

	test_set_selectedFeatures = np.array(test_set_selectedFeatures)
	
	test_set = np.delete(test_set_selectedFeatures, 0, axis=1) # make table suitable for Kmeans, delete non IFP data columns 

	top50 = np.transpose(top50) # make table suitable for Kmeans, axis = samples (x, rows), features (y, columns)
	training_set = np.transpose(training_set)
	test_set = np.transpose(test_set)
	kmeans_clustering =  KMeans(n_clusters=2).fit(training_set)
	print kmeans_clustering.labels_
	print kmeans_clustering.predict(test_set)


					
"""Begin"""

filename = sys.argv[1]
input = open(filename, "r")
output_PCA = "matrix_for_PCA.txt"
text_input = input.read()
output2 = open(output_PCA, "w")
output3 = open("inactives_conserved_interactions.txt", "w")
output4 = open("actives_conserved_interactions.txt", "w")
output5 = open("whole_set.txt", "w")


print "-------------------- IMIFP ---------------------\n\n"
sys.stdout.write("Reading file " + filename + "...")
sys.stdout.flush()
""" Split text in blocks, each containing one matrix, skip first """
text_split = text_input.split("File:	")[1:]


""" For every block of text process and store in matrix """
mtrices = []
for mtrx in text_split:
	mtrices.append( preprocessMatrix(mtrx) )


# mtrices has one element/matrix in the input. Each matrix is also a list, first element AA, second element number, third element 7-bit ifps

common_positions = getCommonPositions(mtrices)
# sys.stdout.write("DONE\n\n")
# sys.stdout.flush()


# # Create a similarity matrix between all matrix combinations based on the tanimoto score
# sys.stdout.write("Calculating tanimoto based similarity matrix...")
# sys.stdout.flush()


# #tanimoto_similarity = similarityMatrix(mtrices, common_positions, 0x7f)
# #similarityMatrixPlot(tanimoto_similarity, "tanimoto_similarity_matrix.png")

# sys.stdout.write("DONE\nSaved as: tanimoto_similarity_matrix.png\n\n")
# sys.stdout.flush()

# calculates similarity only based on interaction pairs in the list
list_interactions = ['550-644','246-640', '132-137']
#similarityXinteraction = clusterXinteraction(mtrices, common_positions, list_interactions, mask=0x7f, exactFilter="FALSE")
#similarityMatrixPlot(similarityXinteraction, "tanimoto_similarityXinteraction_matrix.png")

# # Create a frequency matrix of occuring interactions
# sys.stdout.write("Calculating frequency interaction matrix...")
# sys.stdout.flush()

# #FM = frequencyMatrix(mtrices, common_positions, mask=0x7f, exactFilter="FALSE")
# #frequencyMatrixPlot(FM, common_positions, "frequency_matrix.png")

# sys.stdout.write("DONE\nSaved as: frequency_matrix.png\n\n")
# sys.stdout.flush()

array_conservedout = conservedOut(mtrices, common_positions) 
clean_whole_set, test_set_clean, ni = featureSelection(array_conservedout)
clustering(clean_whole_set, test_set_clean, ni, 50)
#nonconserved_SM = similarityMatrix(array_conservedout)
#similarityMatrixPlot(nonconserved_SM, "tanimoto_similarity_nonconserved_matrix.png")

sys.stdout.write("Creating HTML output file...")
sys.stdout.flush()

aa_curated = []
only_curated_data = []
for matrix in mtrices:
	aa_curated.append(cutMatrixAA(matrix, common_positions))
	only_curated_data.append(cutMatrixIFP(matrix, common_positions))



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
			html = html_generator.makeTableRow( [str(common_positions[row]) + "-" + str(common_positions[col]), html] )
			html_file.write(html)


html_file.write("</table>")
html_file.write("</body></font>\n</html>")

html_file.close()

sys.stdout.write("DONE\nWritten to: ifp_output.html\n\n")
sys.stdout.flush()

output2.close()
output3.close()
output4.close()
output5.close()



input.close()
