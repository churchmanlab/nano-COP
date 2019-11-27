HELP_STRING = """

Date : April 21, 2018

Author : Heather Landry Drexler

This script will output all constitutively spliced introns in a sample based on mature RNA sequencing data. 
Using pysam, the script measures the percentage of correctly spliced reads spanning an intron's splice sites. 
This script also determines which polyA sites (or gene ends in a gtf file) are associated with constitutively 
spliced introns. The initial pass at making this script takes in unstranded RNA-seq reads, but can be easily
modified for stranded data in the future.

use : python intronCoords_from_RNAseq.py transcriptomeCoordinates.txt input.bam output.txt --min_overlap --min_cov --min_splice
            
            -i 				iTrx (input file with exon and intron coordinates e.g. NCBI_RefSeq_hg38_merge_parsed_sortByNameCoord.bed)
            -b 				iBAM (input BAM file with only unique alignments for unstranded sequencing data)
            -o 				oFile (tab delimited file with: chromosome, intron start, intron end, polyAsites("_" delimited))
            --min_overlap 	minimumm overlap required over the splice site
            --min_cov  		minimum coverage required over the splice site
            --min_splice	minimum splicing percentage required over the splice site (will be as a percent out of 100  e.g. 85)

output file (tab delimited):

chrom 	start 	end 	geneName_polyAsite (split by "_" if multiple)	strand

"""

import sys, pysam, os, re
import numpy as np
from getopt import getopt

def main(argv=None):
	if argv is None:
		argv = sys.argv

	inputTrx = ""
	inputBAM = ""
	outputFile = ""
	min_overlap = ""
	min_cov = ""
	min_splice = ""

	try:
		optlist, args = getopt(argv[1:], "hi:b:o:", ['min_overlap=','min_cov=','min_splice='])
	except:
		print ""
		print HELP_STRING
		sys.exit(1)
	   
	if len(optlist) == 0:
		print ""
		print HELP_STRING
		sys.exit(1)
	   
	for (opt, opt_arg) in optlist:

		if opt == '-h':
			print ""
			print HELP_STRING
			sys.exit(1)
		elif opt == '-i':
			inputTrx = opt_arg
		elif opt == '-b':
			inputBam = opt_arg
		elif opt == '-o':
			outputFile = opt_arg
		elif opt == '--min_overlap':
			min_overlap = int(opt_arg)
		elif opt == '--min_cov':
			min_cov = int(opt_arg)	
		elif opt == '--min_splice':
			min_splice = float(opt_arg)
			
	if (inputTrx == "" or inputBam == "" or outputFile == "" or min_overlap == "" or min_cov == "" or min_splice == ""):
		print HELP_STRING
		print "here"
		sys.exit(1)

##################################################################################################
	# open files and define features for script

	iTrx = open(inputTrx, 'r') 					# open the transcript file
	iBAM = pysam.Samfile(inputBam, 'rb') 		# read in the bam file with pysam

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

	# list of cigar strings that match to query (read) and reference (genome)
	Consumes_Query = ["M", "I", "S", "=", "X"]
	Consumes_Reference = ["M", "D", "N", "=", "X"]

	# save the gene name in a set if it's not already in it
	gene_polyAcoords = {}

	# make a dictionary with 3'SS coordinates and potential polyA site coordinates
	junction_polyAcoords = {}

##################################################################################################
	# loop through a file with feature coordinates sorted by gene name
	# and record information about constitutively spliced introns and their associated
	# polyA/gene end coordinates to a new file
	
	# record information about each feature
	for line in iTrx:
	    chrom = line.split('\t')[0]                     # chromosome
	    start = int(line.split('\t')[1])                # start coordinate of feature
	    end = int(line.split('\t')[2])                  # end coordinate of feature
	    gene = line.split('\t')[3]                      # gene name
	    feature = line.split('\t')[4].split('_')[0]     # feature: gene, exon, or intron
	    strand = line.split('\t')[5][0]                 # strand of gene with intron
	    intron_size = end-start                         # get size of intron from coords
	        
	    # get coordinates for intron  
	    coord=chrom+"_"+str(start)+"_"+str(end)+"_"+strand        

	    # get polyA sites for genes from gtf file
	    if feature=="gene":
	        if strand=='+':
	            polyA = int(end)
	        if strand=='-':
	            polyA = int(start)+1 
	            
	        # check if gene name is in the set, if not save it 
	        # start variables for counting the number of reads in features
	        if gene not in gene_polyAcoords.keys():
	        
	            # add gene name to set
	            gene_polyAcoords[gene] = polyA
	            
	        elif polyA != gene_polyAcoords[gene]:
	                print("ERROR: multiple polyA sites for gene name "+str(gene))

	        
	    # use intron coordinates to measure the percentage of molecules spliced 
	    # at that intron in RNA-seq data
	    if feature=="intron":  
	        
	        # set up variables for splicing counts
	        spliced = 0.0 
	        unspliced = 0.0
	        alt5SS = 0.0
	        alt3SS = 0.0
	        skipped = 0.0
	        incomplete = 0.0
	        incorrect = 0.0
	        
	        if coord in junction_polyAcoords.keys():
	            # add polyA site to the dictionary
	            if (str(gene_polyAcoords[gene]) not in junction_polyAcoords[coord]): 
	                junction_polyAcoords[coord].append(str(gene))
	                junction_polyAcoords[coord].append(str(gene_polyAcoords[gene]))
	                
	        elif coord not in junction_polyAcoords.keys():
	        
	            # build a sets for read names at each junction
	            # will be used to prevent overcounting reads from the same transcript
	            found = set()

	            # loop through reads that span the 5'SS
	            for read in iBAM.fetch(chrom, start-min_overlap, end+min_overlap):

	                # set up variables for measuring the length of cigar string operators
	                # for each aligned read to intron
	                end_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
	                start_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}
	                intron_counts = {'M': 0, 'I': 0, 'D': 0, 'N': 0, 'S': 0, 'H': 0, 'P': 0, '=': 0, 'X': 0}

	                # count for position within read
	                currentloc = int(read.pos) 

	                # read through the cigar string
	                for cigar_Entry in read.cigar:

	                    # measure features from the cigar string
	                    cigarOp = CigarNumtoOp[cigar_Entry[0]] # get the cigar operator
	                    op_Length = int(cigar_Entry[1]) # get length of cigar operator
	                    cigarOp_start = currentloc # get start position of operator
	                    if (cigarOp in Consumes_Reference):
	                        currentloc=currentloc+op_Length # add the cigar operator length to the current location coordinate 
	                    cigarOp_end = currentloc # get end position of operator

	                    # check if read is spliced at intron start (can be either 5' or 3' SS)
	                    if (cigarOp_start<start-min_overlap and cigarOp_end>=start-min_overlap):
	                        if (cigarOp_end>=start+min_overlap):
	                            count=min_overlap*2
	                        else:
	                            count=cigarOp_end-(start-min_overlap)
	                        start_counts[cigarOp] += count
	                    elif (cigarOp_start>=start-min_overlap and cigarOp_end<start+min_overlap):
	                        count=op_Length
	                        start_counts[cigarOp] += count
	                    elif (cigarOp_start<start+min_overlap and cigarOp_end>=start+min_overlap):
	                        if (cigarOp_start<=start-min_overlap):
	                            count=min_overlap*2
	                        else:
	                            count=(start+min_overlap)-cigarOp_start
	                        start_counts[cigarOp] += count

	                    # check if read is spliced at intron end (can be either 5' or 3' SS)
	                    if (cigarOp_start<end-min_overlap and cigarOp_end>=end-min_overlap):
	                        if (cigarOp_end>=end+min_overlap):
	                            count=min_overlap*2
	                        else:
	                            count=cigarOp_end-(end-min_overlap)
	                        end_counts[cigarOp] += count
	                    elif (cigarOp_start>=end-min_overlap and cigarOp_end<end+min_overlap):
	                        count=op_Length
	                        end_counts[cigarOp] += count
	                    elif (cigarOp_start<end+min_overlap and cigarOp_end>=end+min_overlap):
	                        if (cigarOp_start<=end-min_overlap):
	                            count=min_overlap*2
	                        else:
	                            count=(end+min_overlap)-cigarOp_start
	                        end_counts[cigarOp] += count


	                currentloc = int(read.pos) 
	                for cigar_Entry in read.cigar:
	                    cigarOp = CigarNumtoOp[cigar_Entry[0]] # get the cigar operator
	                    op_Length = int(cigar_Entry[1]) # get length of cigar operator
	                    cigarOp_start = currentloc
	                    if (cigarOp in Consumes_Reference):
	                        currentloc=currentloc+op_Length # add the cigar operator length to the current location coordinate 
	                    cigarOp_end = currentloc

	                    # check if spliced intron is the correct size               
	                    if (start_counts['M']==min_overlap and start_counts['N']==min_overlap and
	                       end_counts['M']==min_overlap and end_counts['N']==min_overlap):

	                        if (cigarOp_start<start and cigarOp_end>=start):
	                            if (cigarOp_end>=end):
	                                count=end-start
	                            else:
	                                count=cigarOp_end-start
	                            intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

	                        elif (cigarOp_start>=start and cigarOp_end<end):
	                            count=op_Length
	                            intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary

	                        elif (cigarOp_start<end and cigarOp_end>=end):
	                            if (cigarOp_start<=start):
	                                count=end-start
	                            else:
	                                count=end-cigarOp_start
	                            intron_counts[cigarOp] += count # add the cigar operator length to the counts dictionary          


	                # count splicing features for each intron
	                total_count = start_counts['M'] + start_counts['N'] + end_counts['M'] + \
	                end_counts['N'] + intron_counts['M'] + intron_counts['N']

	                # count number of events for each type of splicing
	                if(total_count>=min_overlap*2 and read.qname not in found):

	                    # add the read name to set to prevent overcounting reads coming
	                    # from the same transcript
	                    found.add(read.qname)

	                    if (start_counts['M']==min_overlap and start_counts['N']==min_overlap and 
	                        end_counts['M']==min_overlap and end_counts['N']==min_overlap):

	                        # count number of reads with "correct splicing" (spliced)
	                        if (intron_counts['N']==intron_size):
	                            spliced+=1

	                        # count number of reads with "incomplete splicing" (incomplete) 
	                        if (intron_counts['N']!=intron_size):
	                            incomplete+=1

	                    # count number of reads with "no splicing / retained intron" (unspliced)
	                    elif (start_counts['N']==0 and end_counts['N']==0):
	                        if (start_counts['M']==min_overlap*2 or end_counts['M']==min_overlap*2):
	                            unspliced+=1

	                    # count number of reads with "alternative 5'SS" (alt_5SS)
	                    elif (end_counts['M']==min_overlap and end_counts['N']==min_overlap and start_counts['N']!=min_overlap):
	                        if(strand=="+"):
	                            alt5SS+=1
	                        if(strand=="-"):
	                            alt3SS+=1

	                    # count number of reads with "alternative 3'SS" (alt_3SS)
	                    elif (start_counts['M']==min_overlap and start_counts['N']==min_overlap and end_counts['N']!=min_overlap):
	                        if(strand=="+"):
	                            alt3SS+=1
	                        if(strand=="-"):
	                            alt5SS+=1

	                    # count number of reads with "intron skipping / misalignment" (skipped)
	                    elif (start_counts['N']==min_overlap*2 and end_counts['N']==min_overlap*2):
	                        skipped+=1

	                    # count number of reads with "sequencing / alignment mistakes" (incorrect)
	                    else:
	                        incorrect+=1
	            
	            # count the reads that span the splice junction
	            total_noskip = spliced + unspliced + alt5SS + alt3SS + incomplete + incorrect 

	            # 0.001 is added to avoid division by zero for cases with no coverage
	            if total_noskip == 0: 
	            	total_noskip = 0.001
	            
	            # select introns that have at least the minimum coverage reads that span the splice junction
	            # and more than minimum percent of spliced transcripts at the junction
	            # these paramters can be made to be variable
	            if (total_noskip >= min_cov):
	                if (spliced/total_noskip >= min_splice):
	                    
	                    # add coordinate and polyA site to the junction dictionary
	                    junction_polyAcoords[coord] = [str(gene)]
	                    junction_polyAcoords[coord].append(str(gene_polyAcoords[gene]))

	# close the input file
	iTrx.close()								
	
##################################################################################################
	# write constitutive intron information to an output file
	
	# open the output file for writing
	oFile = open(outputFile, 'w')				
	
	# write a header line
	#oFile.write("\t".join(["chrom","start","end","polyA_sites\n"]))

	# loop through junction dictionary and write ouput to a new file
	for coord in junction_polyAcoords: 
	    chrom = coord.split("_")[0]
	    start = coord.split("_")[1]
	    end = coord.split("_")[2]
	    strand = coord.split("_")[3]
	    polyA = '_'.join(str(i) for i in junction_polyAcoords[coord])

	    # write coordinate and polyA site information to a new file
	    oFile.write(chrom+"\t"+start+"\t"+end+"\t"+polyA+"\t"+strand+"\n")

	# close the output file
	oFile.close()
				 
##################################################################################################
if __name__ == "__main__":
	sys.exit(main())







