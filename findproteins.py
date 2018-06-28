import re
import sys
import os
import shutil
#import these in one line

#move next to iteration over file 

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

#transcription
def transcribe(dna):
	base_replace = {'A':'u', 'T':'a' , 'C':'g', 'G':'c'}

	for i, o in base_replace.items():
		dna = dna.replace(i,o)
	dna_trans = dna.upper()
	return dna_trans

def list_indices(codon, dna):
	indices = []
	index = -1
	#begin at -1 so the first position to search from is 0 after one iteration
	while True:
		#find next index of substring, by starting after its last known position
		index = dna.find(codon, index + 1)
		if index == -1:  
			break  
			#all occurrences have been found
		indices.append(index)
	return indices

#gets start and stop codon positions as lists of indices.
def start_stop_fun(dna):
	if 'AUG' in dna:
		starts = list_indices(r'AUG',dna)
	else:
		print("line " + str(c) + " has no start codon")
		dna = 0
	if 'UAG' in dna or 'UAA' in dna or 'UGA' in dna:
		stops = list_indices(r'U(AG|GA|AA)',dna)
	else:
		print("line " + str(c) + " is an open reading frame")
		dna = 0
	return [start, stop]

#translates
def translate(dna):
	d = 0

	start_stop = start_stop_fun(dna)
	start = start_stop[0]
	stop = start_stop[1]
	protein_dic = {}
	#these are LISTS

	for start in starts:
		for stop in stops:
			d += 1
			protein = ""
			#initialize blank protein
			coding_dna = dna[start:stop]
			length = len(dna_coding)
			codon_num = float(length) / 3.0
			#finds the "number" of codons

			if dna != 0 and stop > start and codon_num.is_integer():
				for num in range(int(codon_num)):
				#this is numbers starting with 0 and ending with the number of codons in the dna
					error = ""
					#error is a blank string here since there is no error
					codon = coding_dna[(num*3):((num*3)+3)]
					#takes a single codon based on what iteration of the range it is on.
					residue = gencode.get(codon,"X")
					#gets a residue matching the codon, defaults to x with no matching codon.
					protein += residue
					#extends the protein by one residue
				#appends a protein string to the returned full dictionary of proteins, indexed by reading frame.
			else:
				error = "frame" + str(d) + "is not a viable reading frame."

			data = [coding_dna, error, protein]
			protein_dict.[d] = dat
			#does this mess up? We'll have to check.


	return protein_dict


infile = sys.argv[1]
file = open(infile)
i = 1
c = 1
out_path = "./proteins/protein1.txt"


for line in file:
	dna = 0
	#resets to non-fail case on each line
	if r'[^ACTNGU]' in line.upper():
		error = "sequence is not nucleic acid!"
		#fail a line if it's not DNA/RNA
	elif r'[A+T+]' or r'[T+A+]' in line.upper()
		error = "sequence contains both Thymine and Uracil!"
		#fails a line that has bad coding.
	elif 'T' in line.upper():
		dna = transcribe(line)
		#transcribes the line if it's DNA
	else:
		dna = line.upper()

	if dna != 0:
		data = translate(dna)

		if len(protein) != 0:

			#creates a non-existent protein file
			while os.path.exists(out_path):
				out_path = "./proteins/seq" + str(i) + ".txt"
				i += 1

			output = open(out_path,"w+")
			output.write(protein)
			output.close()

	c += 1