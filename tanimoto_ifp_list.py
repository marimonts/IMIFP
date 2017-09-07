
# invert (~ifp & 0x7f)

# Returns how many bits are set in the binary representation of a number of any length
def getSetBits(val):
	n = 0
	while(val > 0):
		n+=val&1
		val = val>>1

	return n

def getMyDistanceMetric(ifp1, ifp2, mask=0x7f):
	# Perform bitmask operation
	ifp1 &= mask
	ifp2 &= mask

	# Calculate the total number of bits
	#total_bits = getSetBits(mask) * len(ifp1)	# Once you start working with lists of ifps you need to do this
	total_bits = getSetBits(mask)

	# Get all equal bits (0-0 and 1-1)
	equal_bits = ~(ifp1 ^ ifp2) & mask

	return float( getSetBits(equal_bits)) / total_bits 

# Returns the tanimoto score of two ifps (default:0x7f takes all bits into acount)
def getTanimoto(ifp1, ifp2, mask=0x7f):
	a_sum = 0
	b_sum = 0
	c_sum = 0
	# For every ifp in the provided lists:
	for i in range(0,len(ifp1)):
		# Perform bitmask operation
		a = ifp1[i] & mask
		b = ifp2[i] & mask

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

# Get two bit strings from the user
ifp_list1 = []
ifp_string = raw_input("ifp list 1: ")
for ifp in ifp_string.split():
	ifp_list1.append( int(ifp, 2) )

ifp_list2 = []
ifp_string = raw_input("ifp list 2: ")
for ifp in ifp_string.split():
	ifp_list2.append( int(ifp, 2) )

mask = int(raw_input("Provide a mask: ") ,2)

# Print the tanimoto score
print "Tanimoto score: " + str(getTanimoto(ifp_list1, ifp_list2, mask))

#print "myDistanceMetric: " + str(getMyDistanceMetric(ifp1, ifp2))



