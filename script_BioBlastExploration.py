#! /usr/bin/env python

# need to get this library
from Bio.Blast import NCBIXML
# this is the direction to the blast results
result_handle = open("meth_staph.xml")

# initiate the parser
blast_records = NCBIXML.parse(result_handle)
E_VALUE_THRESH = 0.04

# this instruction is for a looping task
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        print('****Alignment****')
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print(hsp.align_length)
                print(blast_record.query_length)
                #print('****Alignment****')
                #print('query id:', ((blast_record.query).split(' ')[0]).strip())
                #print('length:', blast_record.query_length)
                #print('align title:', ((alignment.title).split(' ')[0]).strip())
                #print('hit id:', alignment.hit_id)
                #print('hit def:', alignment.hit_def)
                #print('accession:', alignment.accession)
                #print('accession:', hsp.expect)
