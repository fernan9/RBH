# RBH
All scripts belonging to this algorithm implementation
RBH for Ortholog Search 1.0 17/04/2015

CONTACT
------
Fernan Rodrigo Pérez Gálvez  +52 (442) 272 4272  fernan954@gmail.com

GENERAL DESCRIPTION
-------------------
This code implements Reciprocal Best Hit algorithm for ortholog search
by comparing one species, denominated 'node', against other 'tip' species.
Assumes a set of BLAST 'all-against-all' files with the exact names expresed
in a list with the 'node' species at the beginning.

Specifications
--------------

IN:	[genome_ls]

	genome_ls (list)		Contains the code for each
					species whose blast 'all-against-all'
					file exists in both directions.
					e.g. "A_B.fasta" and "B_A.fasta"

ARG:	[op, th, r]

	op (int) 	val [1,2]	Indicates the topology of the ortholog
					search analysis.

	th (float)	stdval (0.001)	Indicates the maximum E_VALUE accepted
					during the hit	list discrimination.

	r (float) 	stdval (0.5)	The minimum permisible fraction 
					representing the alignment length over 
					the query length.
		 
OUT: 	[orth_ls]

	orth_ls (txt)			A list containing names ID of all 
					requences that are ortholog candidates.

Internal Functions
------------------

file_verif(list)	Creates all necesary names of the blast output files
			depending in the model of the analysis. Will throw a
			warning if a file is missing, although it will continue
			with the available data.

			IN: list of code names for organisms
			OUT: list of the files to be filled into a dictionary

fill_dict(list)		Parses information from each BLAST output file in the
			given list. An E_VALUE threshold and a length ratio 'r'
			is evaluated to append a sequence name into the key's
			list



LICENCE
-------
Reciprocal Best Hit for Ortholog Search (c) by Fernan Rodrigo Pérez Gálvez

This work is licensed under the Creative Commons Attribution 4.0
International License. To view a copy of this license, visit
http://creativecommons.org/licenses/by/4.0/.
