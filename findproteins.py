import re
import sys
import os
import shutil

if not os.path.exists("./proteins"):
	os.mkdir("./proteins")

gencode = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'}

#transcribes dna
def transcribe(dna):
	base_replace = {'A':'u', 'T':'a' , 'C':'g', 'G':'c'}

	for i, o in base_replace.items():
		dna = dna.replace(i,o)
	rna = dna.upper()
	return rna


#returns a list of the indexes at which a codon is found in a DNA strand.
def list_indices(codon, rna):
	indices = []
	index = -1
	#begin at -1 so the first position to search from is 0 after one iteration
	while True:
		#find next index of substring, by starting after its last known position
		if rna == 0:
			break
		index = rna.find(codon, index + 1)
		if index == -1:  
			break  
			#all occurrences have been found
		indices.append(index)
	return indices

#returns a list of the lists of start and stop positions
def start_stop_fun(rna):
	starts = list_indices(r'AUG',rna)
	if starts == []:
		error = "sequence has no start codons"
		rna = 0

	amber = list_indices(r'UAG',rna)
	ochre = list_indices(r'UAA',rna)
	opal = list_indices(r'UGA',rna)
	stops = amber + ochre + opal

	if stops == []:
		error = "sequence is an open reading frame"
		rna = 0


	return [starts, stops]


#Translates input rna, returning a dictionary of every reading frame's coding sequence and encoded protein
def translate(rna):
	d = 0

	start_stop = start_stop_fun(rna)
	starts = start_stop[0]
	stops = start_stop[1]
	protein_dict = {}
	#these are LISTS

	for start in starts:
		pairable = []
		for stop in stops:
			d += 1
			protein = ""
			#initialize blank protein
			coding_rna = rna[start:stop]
			length = len(coding_rna)
			codon_num = float(length) / 3.0
			#finds the "number" of codons

			if rna != 0 and stop > start and codon_num.is_integer():
				for num in range(int(codon_num)):
				#this is numbers starting with 0 and ending with the number of codons in the rna
					error = ""
					#error is a blank string here since there is no error
					codon = coding_rna[(num*3):((num*3)+3)]
					#takes a single codon based on what iteration of the range it is on.
					residue = gencode.get(codon,"X")
					#gets a residue matching the codon, defaults to x with no matching codon.
					protein += residue
					#extends the protein by one residue
					data = [coding_rna, error, protein]
					protein_dict[d] = data
				#appends a protein string to the returned full dictionary of proteins, indexed by reading frame.
			else:
				error = "frame " + str(d) + " is not a viable reading frame."

			#does this mess up? We'll have to check.


	return protein_dict

#gets the file from whatever is typed after the program in cmd
infile = sys.argv[1]
file = open(infile)

#iterator for file name (sequence number)
i = 1
out_path = "./proteins/seq1.txt"


for line in file:
	#resets to fail case on each line
	rna = 0
	
	#fail a line if it's not DNA/RNA
	if r'[^ACTNGU]' in line.upper():
		error = "sequence is not nucleic acid!"
		data = {0:["",error,""]}
		rna = 0 #error marker

	#fails a line that has bad encoding.
	elif re.match(r"[ACG]*U+[ACG]*T+[ACG]*", line.upper())\
	 or re.match(r"[ACG]*T+[ACG]*U+[ACG]*", line.upper()):
		error = "sequence contains both Thymine and Uracil!"
		data = {0:["",error,""]}
		rna = 0 #error marker

	#transcribes the line if it's DNA
	elif 'T' in line.upper():
		rna = transcribe(line)
	else:
		rna = line.upper()

	if rna != 0:
		data = translate(rna)
		#this is a library by frame index d, which returns list data of the coding frame

	#finds a filename that isn't taken, and makes it the new output path
	while os.path.exists(out_path):
		out_path = "./proteins/seq" + str(i) + ".txt"
		i += 1

	#creates an output file based on set output path 
	output = open(out_path,"w+")
	output.write("RAW:     " + rna)
	output.close()
	#need to reopen here so that we append instead of overwriting during the loop.
	output = open(out_path,"a")

	#output path for non-errors
	if rna != 0:
		for d, datum in data.items():
			#writes frame number
			output.write("FRAME " + str(d) + ": \n")
			#writes coding sequence, error, and protein, if any of each.
			output.write("\t" + datum[0] + "\n\t" + datum[1] + datum[2] + "\n")

	#output path for errors
	else:
		output.write(error)
	
	output.close()


#do garbage collection
