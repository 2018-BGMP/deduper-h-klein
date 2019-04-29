# Helena Klein
# Deduper, single-end
# October 31, 2018
# Usage: python3 deduper_draft.py -f <input sam>

import re 
import argparse
import random


def get_arguments():
	parser = argparse.ArgumentParser(description="Remove PCR duplicates from single-end reads with UMIs. At this time, we can only support single-end duplicate removal. Also, input sam files need to be sorted by RNAME before running this program. For updates, please contact hklein3@uoregon.edu.")
	parser.add_argument('-f', '--file', required=True, help='Indicate the absolute path to a SAM file from which to remove PCR duplicates')
	parser.add_argument('-u', '--umi', required=False, help='If UMIs are not randomers, indicate the absolute path to a file with a list of all known UMIs')
	parser.add_argument('-out', '--output_filename', required=False, help='Indicate the name of your output file')
	parser.add_argument('-p', '--paired', required=False, action='store_true', help='Indicates input SAM contains paired-end sequences.')
	parser.add_argument('-r', '--random', required=False, action='store_true', help='Output at random the pcr duplicate encountered')

	return parser.parse_args()
	
args = get_arguments()
in_sam = args.file
umi_file = args.umi
out_file = args.output_filename
random_out = args.random

out_format = '{}_deduped'
if out_file ==None:
	out_file = out_format.format(in_sam)
	
if umi_file == None:
	randomer=True
else:
	randomer=False
	
if args.paired == True:
	print('The current version of klein_deduper.py does not support this option')
	quit()

# Take input UMI file and convert it to a searchable set
def separate(filename):
	"""Converts the raw unique molecule identifier file into an set for later use."""
	with open(filename, 'rt') as fh:
		set_of_UMIs = []
		for line in fh:
			line = line.strip()
			set_of_UMIs.append(line)
		set_of_UMIs = set(set_of_UMIs)

	return set_of_UMIs


# This function should separate the important bits out of the raw read line
def read_IDs(line):
	"""Reads in the sequence line of a SAM file and returns the separated read identifiers for removing PCR duplicates later"""
	
	line = line.strip().split('\t')
	# separate the UMI from the QNAME
	UMI = line[0].split(":")[-1]
    
	# separate the Flag from the column position, convert to integer
	flag = int(line[1])
    
	#Select the chromosome from RNAME (column 3)
	chromosome = line[2]
    
	# Select position and convert to int
	position = int(line[3])
    
	CIGAR = line[5]
  
	return UMI, flag, chromosome, position, CIGAR


# Get the strand out of the flag
def strand(flag):
	'''Determine the strand of the input line by using the flag'''
	if (flag & 16):
		result = "reverse"
	else:
		result= "forward"
	return result


# These functions adjust the start position of each of these reads, 
#but they do not factor in any unexpected letters, like H or P


def adjust_forward(POS, CIGAR):
	"""Adjust the start position of the forward strand due to soft clipping"""
	if "S" in CIGAR:
		CIGAR = CIGAR.split('S')[0]
		if re.search('[A-Z]', CIGAR)== None:
			CIGAR = int(CIGAR)
			position = POS - CIGAR
		else:  
			position = POS
	else:
		position = POS
	return position


def adjust_reverse(POS, CIGAR):
	adjustment = 0
	# Add together all the integers in the cigar string
	for item in re.split('[A-Z]', CIGAR):
		if len(item)>0:
			adjustment += int(item)
	# Subtract the things that don't affect the actual start position of the reverse strand
	if "S" in CIGAR:
		# We don't want to adjust our start position if the soft clipping is to the left of the POS
		first_S = CIGAR.split('S')[0]
		if re.search('[A-Z]', first_S)==None:
			adjustment -= int(first_S)
	# Insertions do not affect the rightmost starting position on the reverse strand
	if "I" in CIGAR:
		first_I = CIGAR.split("I")[0]
		first_I = re.split('[A-Z]',first_I)[-1]
		adjustment -=int(first_I)
        
	# Do final adjustment (subtract 1 because it's zero based)
	POS += adjustment-1

	return POS


# write head of sam file to output file, and the set of not-pcr-duplicates

def write_file(in_file, out_file, set_of_lines_to_output):
	with open(in_file, 'rt') as fh, open(out_file, 'a') as out:
		output = '{}\n'
		LN=0
		for line in fh:
			line = line.strip()
			if '@' in line:
				out.write(output.format(line))
			elif '@' not in line:
				LN += 1
				if LN in set_of_lines_to_output:
					out.write(output.format(line))



# empty the umi dictionary either when a new RNAME is encountered, or the end of the file is reached

def empty_dict(dictionary ,in_set):
	duplicates = 0
	for entry in dictionary:
		if len(dictionary[entry]) == 1:
			in_set.add(dictionary[entry][0])
		else:
			duplicates += 1
			if random_out:
				choice = random.randint(0,len(dictionary[entry])-1)
				in_set.add(dictionary[entry][choice])
			else:
				in_set.add(dictionary[entry][0])
			
	return in_set, duplicates

# Counter initialization section for output stats

forward = 0
reverse = 0
PCR_duplicates = 0
low_quality_UMI = 0
chromosomes = 0

# This chunk executes the previously defined functions to remove PCR duplicates from single-end data
# Then writes them to a file

# Make set of known UMIs for later reference
if randomer == False:
	UMIs = separate(umi_file)

# create a dict for UMI pairs based on chromosome
umi_dict = {}
key_format= '{}:{}:{}:{}'

# These will be used in the following to clear the working memory of the previous chromosome when a new 
# one is encountered
chromosome = ""
old_chromosome = "none"
unique = set()

# Open sorted sam file
with open(in_sam, 'rt') as fh:
	LN = 0
	unique = set()
    
	for line in fh:
		line = line.strip()
		if "@" not in line:
            # Line number counter for sequence lines in sam file
			LN += 1
			if LN%1000000 == 0:
				working = 'Working on line {}'
				print(working.format(LN+1))
			old_chromosome = chromosome
			UMI, flag, chromosome, position, CIGAR = read_IDs(line)
            
            # empty and reset umi_dict when new RNAME is encountered
			if old_chromosome != chromosome:
				unique, duplicates = empty_dict(umi_dict, unique)
				PCR_duplicates += duplicates
				umi_dict = {}
				chromosomes += 1
                
            # If we have a list of UMIs, then low quality UMIs will 
            #be filtered out during this step, as will undetermined bases without further work. 
			if randomer == False:
				if UMI in UMIs:
                
               	 # Adjust start position depending on strand
					result = strand(flag)
					if result == 'forward':
						forward += 1
						position = adjust_forward(position, CIGAR)
					elif result == 'reverse':
						reverse += 1
						position = adjust_reverse(position, CIGAR)
                    
                
					key = key_format.format(UMI, flag, chromosome, position)
                
					if key in umi_dict:
						umi_dict[key].append(LN)
					else:
						umi_dict[key] = [LN]
				else:
					low_quality_UMI += 1
			else:
			# This means there are random UMIs. Since we have no information about the quality of the UMI, 
			# then filter out any UMI with an "N"
				if "N" not in UMI:
					result = strand(flag)
					if result == 'forward':
						forward += 1
						position = adjust_forward(position, CIGAR)
					elif result == 'reverse':
						reverse += 1
						position = adjust_reverse(position, CIGAR)
					key = key_format.format(UMI, flag, chromosome, position)
					if key in umi_dict:
						umi_dict[key].append(LN)
					else:
						umi_dict[key] = [LN]
				else:
					low_quality_UMI += 1
  # add most recent dictionary to output set
unique, duplicates = empty_dict(umi_dict, unique)
PCR_duplicates += duplicates 
# write to output
write_file(in_sam, out_file, unique)

print()
print("Results")
print("Total reads processed:", LN)
print("Total number of chromosomes:", chromosomes)
print("PCR duplicates removed:", PCR_duplicates)
print("Low-quality UMIs removed:", low_quality_UMI)
print('Reads on the forward strand:', forward)
print("reads on the reverse strand:", reverse)
print()


