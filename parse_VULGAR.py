#!/usr/bin/env python

# usage: parse_VULGAR.py [genomic_sequence] [cigar_line]

import sys, itertools

def reverse_complement(sequence):

	# turn the sequencing into a reverse complement string

	# complement nuc dictionary
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

	# return the reverse complement sequence
	return ''.join([comp[nuc] for nuc in sequence[::-1]])


def convert_sequence(sequence, VULGAR):

	# splice the sequence according to the VULGAR code
	
	# set variables
	start, stop, new_sequence = 0, 0, ''

	# grap the mutation and length from the CIGAR list
	for item in itertools.izip(* [iter(VULGAR)] * 3):
		stop += int(item[2])
		if item[0] in 'MS': new_sequence += sequence[start:stop]
		start = stop

	return new_sequence


def get_sequence():

	# return the sequence
	return [line.strip() for line in open(sys.argv[1])]


def parse_VULGAR():

	# get sequence and raw CIGAR code
	sequence = get_sequence()
	VULGAR = sys.argv[2].strip().split('\n')
	temp_VULGAR = []

	# walk throught the alignments and get the orientation
	# and VULGAR code for each one.
	for code in VULGAR:
		
		# get the direction of the sequence and extract the
		# VULGAR code.
		if ' - ' in code:
			pos, code = code.split(' - ')
		elif ' + ' in code:
			pos, code = code.split(' + ')
		else: return
		pos, code = pos.split(' '), code.split(' ')
		code[0] = int(code[0])
		temp_VULGAR.append([code, pos])

	# sort the list if needed so that the highest scoring
	# alignment will be used
	if len(temp_VULGAR) > 1: temp_VULGAR.sort(reverse=True)

	# check the direction and convert the sequence if needed
	#print 'start {0} stop {1}'.format(temp_VULGAR[0][1][-2], temp_VULGAR[0][1][-1])
	if int(temp_VULGAR[0][1][-1]) < int(temp_VULGAR[0][1][-2]):
		sequence[1] = sequence[1][int(temp_VULGAR[0][1][-1]):int(temp_VULGAR[0][1][-2])]
		sequence[1] = reverse_complement(sequence[1])
	else:
		sequence[1] = sequence[1][int(temp_VULGAR[0][1][-2]):int(temp_VULGAR[0][1][-1])]

	# obtain the sequence introns
	sequence[1] = convert_sequence(sequence[1], temp_VULGAR[0][0][1:])

	# write sequence info to stdout
	for line in sequence: print line

parse_VULGAR()
