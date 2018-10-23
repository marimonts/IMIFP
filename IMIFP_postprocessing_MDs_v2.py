#!usr/bin/python2.7
# -*- coding: utf-8 -*-

import sys 
import re
import argparse
import numpy as np
import pylab as pl
import math
import string
from sets import Set
import html_generator
from sklearn.cluster import KMeans
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectPercentile
from math import sqrt
#import for libraries, sys several functions, incl arguments usage

_threeletter = ["ALA", "CYS", "CYX", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"]
_oneletter = ['A', 'C', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
_watersymbol = ['WAT', 'HOH']


""" Extracts the name of the inout odb file of a matrix in a variable """
def Filename(text):

	# Split text in lines
	mtrx = text.split('\n') #split text in lines
	file = mtrx[0] #line 0 cntains file path
	filename = file.split('/')[-1] #last element of path is filename
	filename = filename[0:-4]
	return filename

""" Stores a matrix from your input file, stores each line containing ifps, and stores header and ifps in strings of 7 bits """
def preprocessMatrix(text):
	# Split text in lines
	mtrx = text.split('\n')
	
	positions = []
	aminoacids = []
	ifps = []
	file = mtrx[0]
	filename = file.split('/')[-1]
	#print filename
	header = ""
	for line in mtrx[2:]:		# Skip the first two lines (title and header)
		if "|" in line:
			temp = line.split('|')
			# Append the aminoacid and position to their respective list
			aa, position = headerConv(temp[0]) #Convert every line header to list containing a single letter aa and its index
			positions.append(position)
			aminoacids.append(aa)

			# Convert an entire line ifp bitstrings to numbers and append to list
			ifp_row = []
			ifp_line = temp[1].strip()


			if len(ifp_line) % 7 == 0:
				while ifp_line:
					ifp = str(ifp_line[:7])
					#turn on following to simplify 7-bit IFPs to 4-bit IFPs
					#unified_ifp = str(int(ifp[0])) + str(int(ifp[1]) or (int(ifp[2]))) + str(int(ifp[3]) or (int(ifp[4]))) + str(int(ifp[5]) or (ifp[6]))
					ifp_row.append(int(ifp,2))	# Convert ifp bitstring to decimal number
					ifp_line = ifp_line[7:]					# Delete processed from the line

				# Save the this row of ifps as list of lists
				ifps.append(ifp_row)



			else:
				print "IFPs no multiple of 7!"	
		else:
			'I am going to break'
			break

	return list([aminoacids, positions, ifps])

"""Convert every line header to list containing a single letter aa and it's index"""
def headerConv(header):	
	
	header = header.strip(' ')

	# Get aminoacid letter code

	aa = header[-3:]

	# Get residue index from the header	
	if header[0] in string.ascii_letters:
		number = header.strip('-')
		number = int(number[1:-3])
	else:
		number = header.strip('-')
		number = int(number[0:-3])
	
	#print number
	# Replace three letter code for single letter
	if aa in _threeletter:
		aa = _oneletter[_threeletter.index(aa)]

	return aa, number


# Creates a sorted list with all the positions that are present in any matrix (at least once)
def getAllPositions(matrix_container):


	common = set(matrix_container[0][1])			# Init common with first matrix [0] positions [0][1]

	for matrix in matrix_container[1:]:	# Go over every matrix, starting from the second

		common = common.union(set(matrix[1])) #index 1 contains positions
	return sorted(list(common))

#function that converts each matrix to numpy array of ifps 
def formatMatrix(matrix, mask=0x7f):

	bin_tmp = []
	pos_x = []
	pos_y =  []
	bin_arrays = []

	# formatting each matrix and storing it in bin with only IFPs
	matrixnp = np.array(matrix[2])

	# change data ordering from receptor-based grouped to interaction pair grouped
	for row in range(len(matrixnp)):

		for col in range(len(matrixnp[row])):
			cell = []
			cell.append(matrixnp[row][col] & mask)


			bin_tmp.append(cell)
			pos_x.append(matrix[1][row])
			pos_y.append(matrix[1][col])

	# merge everything together by making a list of lists, in which every sublist contains the same element position of the original lists, e.g. ([64, 64, 0, 64, 64...], 332, 336)	

	# zip positions x and y together (zip is NOT append, each element is zipped together with the element of the same index in the other list(s))
	interactions = zip(pos_x, pos_y)

	# create two separate arrays and merge them
	bin_arrays_positions = np.array(interactions)
	bin_arrays_ifps = np.array(bin_tmp)
	bin_nparray = np.concatenate((bin_arrays_positions, bin_arrays_ifps), axis=1)
	
	#deleting all interaction pairs with IFP = 0
	bin_array = bin_nparray.tolist()
	bin_array_nozeros = []
	for row in bin_array:
		if row[-1] != 0:
			bin_array_nozeros.append(row)


	#np.savetxt(output2, bin_nparray, fmt='%d',delimiter='	', newline='\n') #look later into printing into file if useul

	return bin_array_nozeros

#delete redundancy (delete rows that represent the same interaction)(e.g. 349350 = 350349)
def deleteDoubles(bin_nparray, mask=0x7f):

	pos_tokeep = []
	for row in bin_nparray: #use column index
		if row[0] < row[1]: #to delete doubles, only keep the ones where the first residue is a lower poistion than the second, and not the othwe way around (e.g. 349 350, but NO 350 349)
			pos_tokeep.append(row)

	fullSet= np.array(pos_tokeep) # convert to array

	return fullSet
"""Function that calculates basic water statistics and adds a 8th bit to water mediated interaction pairs. Then creates a new ifp list with only those interaction pairs with a water mediated HBond"""
def waterStatistics(aa, positions, ifps):
	aa_x = []
	aa_y = []

	nr_water_interactions = 0
	nr_water_hbonds = 0

	ifps = ifps.tolist()

	hbond_interactions = []
	for row in ifps: #iterate over each row of the ifps array, which current format is col0: position x, col1: position y, col2+: ifps
		row.insert(2, aa[positions.index(row[0])]) # inserts in position 2 the AA (in aa[0]) in the index of position with index of the element row[0]
		row.insert(3, aa[positions.index(row[1])]) # inserts in position 3 the AA (in aa[0]) in the index of position with index of the element row[1]

	#calculate basic statistics
		if row[-1] != 0:
			if row[-3] in _watersymbol or row[-2] in _watersymbol:
				nr_water_interactions = nr_water_interactions + 1 # calculates total number of (any)interactions that involve water
				if bin(row[-1])[5] or bin(row[-1])[6] == "1":
					#row[-1] = row[-1] + (1 << 7) # when water in interaction and H-bond bits on, add additional 8th bit to the least significant side
					nr_water_hbonds = nr_water_hbonds + 1 # calculates total number of hbonds that involve water
					hbond_interactions.append(row[0:4])	

	return ifps, hbond_interactions

"""Creates a dictionary with all waters every residue interacts with, to start analysis"""
def waterAdvanced(aa, positions, hbond_interactions):

	residueX = []	
	residueY = []
	
	for interaction in hbond_interactions: #iterate over the list that contains all info lists (ifp, posx, posy, aax, aay)

		residueX.append(str(interaction[0])+'-'+interaction[2]) # put together posx and aax
		residueY.append(str(interaction[1])+'-'+interaction[3]) #put together posy and aay


	hbond_interactions = zip(residueX,residueY) #put together residues interacting 
	interaction_dictionary = {} #create empty dictionary
	interactors_query = set(residueX).union(set(residueY)) #create reference list with all resisudes of interest (selected before, meaning involving at least 1 water AND a hydrogen bond)
	for residue in interactors_query: #iterate over the reference list with all residues
		interaction_dictionary[residue] = [] #create a key of the residue with no values 
		for row in hbond_interactions: #iterate over the data list (only labels now)
			if residue in row: #following will check both residues in the interaction to add as a value for each key only the interacting residue to not repeat the key
				for element in row:
					if element == residue:
						pass
					else:
						interaction_dictionary[residue].append(element) #add the non-repeat residue as value for the key

	return interaction_dictionary

def waterNetRWR(interaction_dictionary): # to retrieve all interacting trios with a water involved
	trios = []
	for key in interaction_dictionary.keys(): #iterate over every key in the dictionary (reference)

		if str(key) not in _watersymbol:
			for partner in interaction_dictionary[key]: # iterate over all the residues with which your initial search interacts with
				
				for keys, values in interaction_dictionary.iteritems(): #iterate over all items in dictionary, keys and values separately
					if partner in keys: #look for matches
						
						for p in interaction_dictionary[partner]: #look for values of each partner 
							trios_tmp = []
							if (p != key) and (str(p) not in _watersymbol): #do not repeat the original search + look for residue - water - residue combinations
								#print str(key) + ' - ' + str(partner) + ' - ' + str(p)
								trios_tmp.append(key)
								trios_tmp.append(partner)
								trios_tmp.append(p)
							else:
								pass
							trios.append(trios_tmp)
	trios = filter(None, trios) #deleting empty elements in list
	trios_tokeep = [] #to remove duplicates (each trio is twice in the list in opposite order!)
	for trio in trios:
		if int(trio[0].split('-')[0]) < int(trio[2].split('-')[0]): #compare the residue number of first and last residue of the trio
			trios_tokeep.append(trio) #keep only one combination

	for trio in trios_tokeep:
		print >>output6, trio


	return trios_tokeep

def waterNetRWWR(interaction_dictionary): # to retrieve all interacting quadruplets with the format: residue - water - water - residue
	quadruplets = []
	
	for key in interaction_dictionary.keys(): #iterate over every key in the dictionary (reference)

		if str(key) not in _watersymbol:
			for partner in interaction_dictionary[key]: # iterate over all the residues with which your initial search interacts with
				
				#for keys, values in interaction_dictionary.iteritems(): #iterate over all items in dictionary, keys and values separately
				if partner in interaction_dictionary.keys() and str(partner) in _watersymbol: #look for matches
					
					for p in interaction_dictionary[partner]: #look for values of each partner 
						
						if (p not in interaction_dictionary[key]) and (str(p) in _watersymbol): #do not repeat the original search
							for forth in interaction_dictionary[p]:
								quadruplets_tmp = []
								if (forth != key) and (forth != p) and (str(forth) not in _watersymbol) and (forth not in interaction_dictionary[partner]):
									quadruplets_tmp.append(key)
									quadruplets_tmp.append(partner)
									quadruplets_tmp.append(p)
									quadruplets_tmp.append(forth)
								else:
									pass
								quadruplets.append(quadruplets_tmp)
	
	
	quadruplets = filter(None, quadruplets) #deleting empty elements in list

	quadruplets_tokeep = [] #to remove duplicates (each trio is twice in the list in opposite order!)
	for quadru in quadruplets:
		if int(quadru[0].split('-')[0]) < int(quadru[3].split('-')[0]): #compare the residue number of first and last residue of the trio
			quadruplets_tokeep.append(quadru) #keep only one combination
	for t in quadruplets_tokeep:
		print >>output7, t

def waterClusters(interaction_dictionary): # to retrieve all interacting trios with a water involved
	clusters = []
	for key in interaction_dictionary.keys(): #iterate over every key in the dictionary (reference)

		if str(key) not in _watersymbol: #look for the interaction partners of every water


					#interaction_dictionary[key].remove(partner) #remove waters from clusters
			if len(interaction_dictionary[key]) >= 2: #only look into clusters that connect more than one residue
				clusters.append(interaction_dictionary[key])


	for cluster in clusters:
		print >>output8, cluster

def waterClustersExtended(interaction_dictionary): # to retrieve all interacting trios with a water involved
	
	clusters = []
	for key in interaction_dictionary.keys(): #iterate over every key in the dictionary (reference)
		extended = []
		if str(key) not in _watersymbol: #look for the interaction partners of every water
			extended.append(interaction_dictionary[key])
		
			for partner in interaction_dictionary[key]:

				if str(partner) in _watersymbol and len(interaction_dictionary[partner]) >= 2 and (partner != key):
		
					

					extended.append(interaction_dictionary[partner])
				

		clusters.append(extended)

	clusters = filter(None, clusters)

	for cluster in clusters:
	 	print >>output9, cluster

"""Transform the list of RWR interactions into the same format as the data tables. E.g. 'position x, position y, aa x, aa y' """
def transformRWR(RWR):

	RWRtransformed = []
	for RWRset in RWR:
		RWR_tmp = []
		RWRset.pop(1)

		RWR_tmp.append(RWRset[0].split('-')[0])
		RWR_tmp.append(RWRset[1].split('-')[0])
		RWR_tmp.append(RWRset[0].split('-')[1])
		RWR_tmp.append(RWRset[1].split('-')[1])
		RWR_tmp.append(128)
		
		RWRtransformed.append(RWR_tmp)


	return RWRtransformed



def ifpsWaterMediatedInteractions(aas, positions, ifps_with_waters, RWRtransformed):


	ifps_tokeep = []
	for row in ifps_with_waters:
		if row[2] not in _watersymbol and row[3] not in _watersymbol:
			ifps_tokeep.append(row)


	a = [str(row[0]) + str(row[1]) for row in ifps_tokeep]
	b = [str(row[0]) + str(row[1]) for row in RWRtransformed] 

	setintersect = set(a) & set(b)


	indecesA = [a.index(x) for x in setintersect]
	indecesB = [b.index(x) for x in setintersect]

	for i in indecesA:
		ifps_tokeep[i][-1] = ifps_tokeep[i][-1] + (1 << 7)

	return ifps_tokeep

def combineMatrices(matrices):
	
	aax_col = 2
	aay_col = 3
	ifp_col = 4 				# Position in the matrix where ifps are stored
	interactions = list() 		# A list of lists with the interactions found in each matrix in 'posx_posy'
	all_interactions = set() 	# A unique referenceset containing all interactions found in all matrices

	for matrix in matrices:
		interactions.append([str(row[0]) + "_" + str(row[1]) for row in matrix]) # Make a string combining posx and posy in the format 'posx_posy' and append to a list of lists (a list for every matrix)

		all_interactions = all_interactions.union( set(interactions[-1]) ) 		 # For every matrix processed, add de posx_posy that did not exist in the set, to it (for reference)

	output = list()																 # Output will be a new list (we do not modify initial matrices!)
	for interaction in list(all_interactions):									 # !!!Iterate over every possible 'posx_posy' present in the reference set!!!
		row = [interaction] + [int(x) for x in interaction.split("_")] 			 # Start row with 'posx_posy' 

		for m in range(len(matrices)):											 # Iterate over every matrix
			if(interaction in set(interactions[m])):							 # If 'interaction' (initial 'for' loop) is present in the specific list of interactions of this matrix...
				i = interactions[m].index(interaction)							 # ... add the ifp value for that matrix to the row (looking for its index in the list of interactions of its matrix)
				matrix_row = matrices[m][i]
				row = row + [matrix_row[aax_col] + "-" + matrix_row[aay_col], matrix_row[ifp_col]]							 	
			else:																 # Otherwise, add a value 0 for the given matrix
				row = row + ["X-X",0]

		output.append(row)

	# Sort the list of lists by the second element (posx) in each entry list using a lambda
	output.sort(key=lambda output:output[1])

	return output



"""Begin"""

filename = sys.argv[1]
if sys.argv[2] == 'waters':
	watersTrue = sys.argv[2]

input = open(filename, "r")
text_input = input.read()


print "-------------------- IMIFP ---------------------\n\n"
sys.stdout.write("Reading file " + filename + "...")
sys.stdout.flush()

""" Split text in blocks, each containing one matrix, skip first """
text_split = text_input.split("File:	")[1:]


""" For every block of text process and store in matrix """
mtrices = []
filenames = []
for mtrx in text_split:
	filenames.append(Filename(mtrx)) #extract the name of the file where the matrix comes from, to append to output name files
	mtrices.append( preprocessMatrix(mtrx) ) # initial formatting




#all_positions = getAllPositions(mtrices) # Creates a sorted list with all the positions that are present in any matrix (at least once)


matrix_reformat = []
aa_compilation = []
position_compilation = []

for matrix in mtrices:
	aa_compilation.append(matrix[0])
	position_compilation.append(matrix[1])
	matrix_reformat.append(formatMatrix(matrix)) 


matrix_noDoubles = []
for matrix in matrix_reformat:
	matrix_noDoubles.append(deleteDoubles(matrix)) #delete redundancy (delete rows that represent the same interaction)(e.g. 349350 = 350349)

ifps_with_waters = []
hbond_interactions = []
water_interaction_dictionary = []
RWR = []

for i in range(len(matrix_noDoubles)):
	ifps_with_waters_tmp, hbond_interactions_tmp = waterStatistics(aa_compilation[i], position_compilation[i], matrix_noDoubles[i]) #calculate basic stats and format
	ifps_with_waters.append(ifps_with_waters_tmp)
	hbond_interactions.append(hbond_interactions_tmp)

	water_interaction_dictionary.append(waterAdvanced(aa_compilation[i], position_compilation[i], hbond_interactions[i])) #creates the water mediated interaction dictionary needed for later analysis
	
	output6 = open(filenames[i]+"_waterNetRWR.txt", "w")
	output7 = open(filenames[i]+"_waterNetRWWR.txt", "w")
	output8 = open(filenames[i]+"_waterClusters.txt", "w")
	output9 = open(filenames[i]+"_waterClustersExtended.txt", "w")
	RWR.append(waterNetRWR(water_interaction_dictionary[i])) # trios of interactions of type residue - water - residue, saved into lists
	waterNetRWWR(water_interaction_dictionary[i]) # quadruplets of interactions of type residue - water - water - residue 
	waterClusters(water_interaction_dictionary[i]) # print all interacting residues per water molecule
	waterClustersExtended(water_interaction_dictionary[i]) # print all interacting residues per water molecule and if one of them is another water, print residues interactiong with the 2nd water
	
RWRtransformed = []
for RWRset in RWR:
	RWRtransformed.append(transformRWR(RWRset))
	

ifps_with_water_mediated = []
for i in range(len(ifps_with_waters)):
	ifps_with_water_mediated.append(ifpsWaterMediatedInteractions(aa_compilation[i], position_compilation[i], ifps_with_waters[i], RWRtransformed[i]))

combined_matrices = combineMatrices(ifps_with_water_mediated)


#Creating html output file for data visualization
html_file = open("ifp_output_test.html", "w")

html_file.write("<html>\n\r<head><title>IFP</title> <link rel='stylesheet' href='ifps.css'> </head>\n\r<font face='verdana'><body><center>\n\r")
html_file.write("<h1>Interaction Fingerprint Viewer</h1>\n\r")
html_file.write("<table border='0' cellpadding='1'>\n\r")


html = "<th class='rotate'>"
for f in filenames:
	html = html + "<div><span>" + f + "</span></div>"
html = html + "</th>"
html_file.write(html)


mask = 0b11111111

for interaction in combined_matrices:

	header = []
	cell = []

	for i in range(3, len(interaction), 2):
		header = header + [interaction[i]]
		cell = cell + [interaction[i+1] & mask]

	html = html_generator.ifpTable(cell, 8, header)
	html = html_generator.makeTableRow( [interaction[0], html] )
	html_file.write(html)


html_file.write("</table>")
html_file.write("</body></font>\n</html>")

html_file.close()

output6.close()
output7.close()
output8.close()
output9.close()