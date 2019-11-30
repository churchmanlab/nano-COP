"""

Date : October 31, 2017

Author : Heather Landry Drexler

Count the number of reads spanning the 5'SS and 3'SS of all introns in the genome and record splicing status.
This script will produce a file that can be used to calculate and plot splicing index values.

use : python junctionCounts_for_SplicingIndex.py intronCoordinates.txt input.bam output.txt
            
            iIntrons (input file with splice site coordinates for introns) [1]
            iBAM (input BAM file with only unique alignments) [2]
            oFile (output file with gene name, intron number, chr, intron start, intron end, strand, intron length, 
                    5' unspliced read number, 3' unspliced read number, exon-exon junction read numbers) [3]
                                            
"""
import sys, pysam, os, numpy, re, time, datetime
from collections import Counter

# print the start time in output file
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print("start: "+str(st))

# import files 
iIntrons = open(sys.argv[1], 'r')
iBAM = pysam.Samfile(sys.argv[2], 'rb')
oFile = open(sys.argv[3], 'w')

# write header to houtput file
oFile.write('gene'+'\t'+'intron'+'\t'+'chr'+'\t'+'start'+'\t'+'end'+'\t'+'strand'+'\t'+'5SS_count'+'\t'+'3SS_count'+'\t'+'splice_count'+'\n')

# classifier for cigar string
CigarNumtoOp = {0 : 'M',
                1 : 'I',
                2 : 'D',
                3 : 'N',
                4 : 'S',             
                5 : 'H',
                6 : 'P', 
                7 : '=',
                8 : 'X'}


# loop through a file with intron coordinates
# record features about the introns
for line in iIntrons:
    chrom = line.split('\t')[0]                 # chromosome
    start = line.split('\t')[1]                 # start coordinate of intron (last base of exon)
    end = line.split('\t')[2]                   # end coordinate of intron (last base of intron)
    gene = line.split('\t')[3].split('_')[1]    # gene name
    intron = line.split('\t')[3].split('_')[3]  # intron number within the gene
    strand = line.split('\t')[5][0]             # strand of gene with intron
    length = int(end)-int(start)                # length of intron 
    unspliced5SS = 0                            # variable for measuring unspliced reads that span 5'SS
    unspliced3SS = 0                            # variable for measuring unspliced reads that span 3'SS
    splice5SS = 0                               # variable for measuring spliced reads that span 5'SS
    splice3SS = 0                               # variable for measuring spliced reads that span 3'SS  

    # build a set for read names at each junction
    found = set()

    # get 5' SS and 3' SS positions for introns
    if strand=='+':
        pos5SS = int(start)+1
        pos3SS = int(end)
    if strand=='-':
        pos5SS = int(end)
        pos3SS = int(start)+1    
        
    
    # if the intron is within 100nt from end of chromosome, skip to next one      
    if int(pos5SS)<100:
        continue

    
    # fetch reads that span the 5'SS
    for read in iBAM.fetch(chrom, int(pos5SS)-1, int(pos5SS)):

        # set up a dictionary to count all cigar operators from each read 
        CigarOp_Counts = {'M': 0,
              'I': 0,
              'D': 0,
              'N': 0,
              'S': 0,             
              'H': 0,
              'P': 0, 
              '=': 0,
              'X': 0}
        
        # take reads with one end <=100 nt from the junction
        # this gets rid of reads that are spliced far away from the junction
        if pos5SS-read.pos<=100 or read.aend-pos5SS<=100: 
            
            # take reads that have ends >=3 nt from the junction
            # this avoids reads have ends close to the junction and might be misaligned
            if pos5SS-read.pos>=3 and read.aend-pos5SS>=3:
                  
                # if read is in the sense direction for the intron    
                if (strand=='+' and read.is_reverse and read.is_read2
                    or (strand=='+' and not read.is_reverse and read.is_read1)
                    or (strand=='-' and read.is_reverse and read.is_read1)
                    or (strand=='-' and not read.is_reverse and read.is_read2)):
                    
                    # get read name and check if it is already in the set
                    key = str(read.qname)+"_5SS_"+str(pos5SS)
                    if key not in found:
                        found.add(key)

                        for cigar_Entry in read.cigar:
                            cigarOp = CigarNumtoOp[cigar_Entry[0]]
                            op_Length = cigar_Entry[1]
                            CigarOp_Counts[cigarOp] += op_Length
                            if cigarOp=='N' and length-1<=op_Length<=length+1:
                                splice5SS+=1
                        if CigarOp_Counts['N'] == 0:
                            unspliced5SS += 1
    

    # fetch reads that span the 3'SS
    for read in iBAM.fetch(chrom, int(pos3SS)-1, int(pos3SS)):

        # set up a dictionary to count all cigar operators from each read 
        CigarOp_Counts = {'M': 0,
              'I': 0,
              'D': 0,
              'N': 0,
              'S': 0,             
              'H': 0,
              'P': 0, 
              '=': 0,
              'X': 0}
        
        # take reads with one end <=100 nt from the junction
        # this gets rid of reads that are spliced far away from the junction
        if pos3SS-read.pos<=100 or read.aend-pos3SS<=100: 
            
            # take reads that have ends >=3 nt from the junction
            # this avoids reads have ends close to the junction and might be misaligned
            if pos3SS-read.pos>=3 and read.aend-pos3SS>=3:
                  
                # if read is in the sense direction for the intron    
                if (strand=='+' and read.is_reverse and read.is_read2
                    or (strand=='+' and not read.is_reverse and read.is_read1)
                    or (strand=='-' and read.is_reverse and read.is_read1)
                    or (strand=='-' and not read.is_reverse and read.is_read2)):
                    
                    # get read name and check if it is already in the set
                    key = str(read.qname)+"_3SS_"+str(pos3SS)
                    if key not in found:
                        found.add(key)

                        for cigar_Entry in read.cigar:
                            cigarOp = CigarNumtoOp[cigar_Entry[0]]
                            op_Length = cigar_Entry[1]
                            CigarOp_Counts[cigarOp] += op_Length
                            if cigarOp=='N' and length-1<=op_Length<=length+1:
                                splice3SS+=1
                        if CigarOp_Counts['N'] == 0:
                            unspliced3SS += 1


    # write a file to output with gene/intron information 
    # and counts for reads that span the splice junctions
    oFile.write(str('NM_'+gene)+'\t'+str(int(intron)+1)+'\t'+str(chrom)+'\t'+str(start)+'\t'+str(end)+'\t'+str(strand)+'\t'+str(unspliced5SS)+'\t'+str(unspliced3SS)+'\t'+str(splice5SS)+'\n')


# print the end time in output file
ts = time.time()
st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
print("end: "+str(st))

# close all files
iIntrons.close()
iBAM.close()
oFile.close()





