#!/usr/bin/env python3
"""Removes PCR introduced duplicate sequences from a Sam formated file. Will adjust alignement based on soft clipping amount

For SINGLE END/PAIRED END reads, duplicates are removed if they fulfill the following: 
a) start at the same genomic coordinate (ie. same left most maping postion adjusted for soft clipping) 
b) have the same molecular id tag (index, barcode, umi, randomer)
c) Have the same strand orientation (+,-). 
The read with the highest quality score is kept as the non-duplicate read.
All reads that are not uniquely mapped ( ie. mapped and no secondary alignments) will not be retained.


How to use this tool:
- (Standard): User supplies One input file,
 1) A correctly formatted SAM file. The sam file should have molecular id tag in the header with no spaces(see example sam headers) 
 2) User must supply the length of the index used (ie -l 6) in the options command.
 3) if sam file is unsorted use the --sort option * assumes samtools 1.5 or better is installed.
 
 Options for dual indexes:
 1) If using dual indexes the following separtors are allowed (^ or -) to separate the two indexes in the sam header (NO SPACES).
 2) supply the combined length of both indcies in the ---length option (ie 2 indces of 6 length= 12)
    
Runtime Optimized: User supplies only one input file,
 1) SAM file that a) contains unique alignments only b) is sorted c) has a fixed length sequence containing the
    molecular tag appended to each read name.

Example sam headers:
 -(indexes can only contain the folowing characters [ACTGN^-] preferably at the end of the header but can be anywhere. *for now MUST BE CAPITAL lETTERS)
 ex1: HWI-ST354R:351:C0UPMACXX:6:2301:2151:AAGCTC
 ex2: HWI-ST354R:351:C0UPMACXX:6:2301:2151:AAGCTC^GGCTAA (AAGCTC-GGCTAA is acebtable)
 ex3: AACTTG:HWI-ST354R:351:C0UPMACXX:6:2301:2151
    
"""
__author__ = 'Devin Dinwiddie'
__contact__ = "devin-d.github.io"
__version__ = '1.0'
__copyright__ = """Copyright (C) 2017 Devin Dinwiddie.
   
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License 
    along with this program. If not, see <http://www.gnu.org/licenses/>."""
    
import subprocess
import logging
import tempfile
import sys
import os




# set up logger
logging.basicConfig(filename='dup_remove.log', level=logging.DEBUG, 
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger=logging.getLogger(__name__)

class sam_info(object):
    """ Filters a sam file line adjusts alignment start postion based on soft clipping amount, 
    converts quality into phred score, checks for strandeness and unique mapping, gets molecular ID from sam header
    
    Attributes:
        line -- single line from a sam file (str)
        start_pos -- lines left most mapping postion soft clip adjusted (int)
        qual_score --line quality score (int)
        strand -- lines strandeness (str)
        BC -- unique index found in sam header (str)
        
    """
    
    def __init__(self,line):
        """ init
        
        arguments
        line-- a single line from a sorted sam file
        """
        self.line=line
        self.start_pos=self.start_pos() 
        self.qual_score=self.convert_phred()
        self.strand=self.bit_check()
        self.BC=self.barcode_get()
        
    def soft_check(self,align):
        """ function to check for soft clipping *returns Soft clip amt"""
        soft_amt=0
        if align[1]=="S" or align[1]=="s":
            soft_amt=align[0]
        elif align[2]=="S" or align[2]== "s":
            soft_amt=align[:2]
    
        return(int(soft_amt))

    def convert_phred(self):
        """Converts a quality score line into a phred score *returns the sum of individual character scores* 
        assumes phred score ASCII_BASE 33"""
        qual=(self.line.split("\t")[10])
        sum_qual=0
        for letter in qual:
            n = ord(letter) - 33
            sum_qual+=n
            return sum_qual
        
    def start_pos(self):
        """ returns the left most mapping postion of an alignment
        adjusted for soft clipping if soft clipping is found"""
        start_pos=int(self.line.split("\t")[3])
        align=self.line.split("\t")[5]
        soft=self.soft_check(align)
        if soft  > 0:
            start_pos-=soft
        return start_pos
    
    def barcode_get(self):
        """ looks in the sam header for the (index,barcode,umi, randomer) assumes the barcode is at the end of the header with no spaces 
        and user input barcode length as it is in the header
        ex: NS500451:204:HH7GHBGXY:1:11101:10281:1048:AACCCG -barcode length 6
        NS500451:204:HH7GHBGXY:1:11101:10281:1048-AACCCG^ACCGCN -barcode length 13 (include ^)
        returns the barcode. if no barcode found read will be logged in the log file"""
#         lengthBC=6
#         BC=self.line.split("\t")[0][-6:]
        global COUNT
        BC=""
        try:
            BC = subprocess.check_output("cat {sam} |cut -f 1 | grep -v '^@'|sed -n '{NR}'p |grep -o '{seq}'".format(sam=sam,NR=COUNT,seq=seq), shell=True,universal_newlines=True)
            
        except:
            logger.error("Error@ line{co}: Could not find the indetifier in sam header,will not be included in output: {head}\n".format(co=COUNT,head=self.line.split("\t")[0]))
            pass
        COUNT +=1
        
        return BC.strip("\n")
        
    def bit_check(self):
        """checks that sequence is mapped and and no secondary allignment
        also gets the strandness (+ or -)
        *returns strandeness and empty string if uniq mapping
        *non uniq mapped will be returned but ignored for duplicates (ie not output)"""
        flag=int(self.line.split("\t")[1])
        strand="+"
        if ((flag & 4)!=4 and (flag & 256)!=256): #mapped and unique allignment
            umap= ""
        else:
            umap= None
            logger.error("Error@ line{co}: Read is unmapped or non unique alligned, will not be included in output: {head}".format(co=COUNT,head=self.line.split("\t")[0]))
        if ((flag & 16))==16: #strand is -
            strand="-"
        return (strand,umap)
    
class rmdup(object):
    """ interates a sam file collects Unique non pcr duplicate reads with best quality"""
    def __init__(self):
        self.sorted_sam=sorted_sam
        """ init
        
        arguments
        sorted sam-- a sam file that has been sorted by alignment postion (ie *samtools sort)
        """
        
    
    def inter_sam(sorted_sam):
        
        """ interates a sam file and removes pcr duplicates if removal conditions are met"""
         
        global NR
        bestScore=0 # storage for best quality score

        with open(sorted_sam)as sam:


            for line in sam:

                if not line.startswith('@'):
                    line=sam_info(line)
                    NR+=1
                    val=line.BC,line.start_pos,line.strand[0] #(barcode,start_pos, +or-)

                    if line.strand[1] is None or line.BC is "": # not mapped or secondary allignment or no index in header *ignore
                        continue
                    elif val in place.keys() and bestScore < line.qual_score : # if val is in dict but current line quality score is highest replace with current line

                        for key in place.keys():
                            if key == val:
                                place[val] = line.line
                            break
                    elif val not in place.keys(): #if values is not in dict put it in with line as value

                        place[val]=line.line
                        bestScore = line.qual_score
                    else: #for cases where duplicates with same quality score 1st read will be retained
                        continue
            return 
        
#class sam_sort(object):
    #""" sorting a sam file"""
    
    #def __init__(self,sam,out):
        
        #""" init
        
        #arguments
        #sam -- a sam file that has not been sorted (str)
        #out -- name of output file (str)
        #"""
        #self.sam=sam
        
        #self.out=out
        
def sort_sam(sam,out):
    """ sorts a sam file by left most mapping postion * defaults to 3M temp memory storage and 28 nodes (see samtools manual -m,-@) may need adjustments based on user system"""
    with tempfile.TemporaryFile() as f: #uses a temp file
        try:
            output=subprocess.check_output("samtools view -bS {sam} |samtools sort -m 3M -@ 28 -o {out} ".format(sam=sam,out=out),shell=True,stderr=subprocess.STDOUT) #runs samtools
        except subprocess.CalledProcessError as er :
            logger.error("error sorting sam file:\n{error}\n".format(error=er.output)) #logs any errors
            exit(1)
                
def file_check(parser, arg):
    """ checks if input files exist"""
    if not os.path.exists(arg):
        parser.error("The file {0} does not exist!".format(arg))
    else:
        return str(arg)
    
if __name__ == '__main__':
    import argparse, sys, logging
    default_tmp_dir = '/tmp'

    parser = argparse.ArgumentParser(description=__doc__.format(author=__author__, contact=__contact__), formatter_class=argparse.RawDescriptionHelpFormatter, add_help=False)

    pgroup = parser.add_argument_group("Input")
    pgroup.add_argument('sam', metavar='IN.sam', type=lambda x:file_check(parser, x), help='input sorted/unsorted SAM')

    ogroup = parser.add_argument_group("Options")
    #ogroup.add_argument('-2','--paired-end', dest='pe', action='count', default=0, help="use paired end deduping with template. SAM/BAM alignment must contain paired end reads. Degenerate read pairs (alignments for one read of pair) will be discarded.")
    #ogroup.add_argument('-f', dest='fq', metavar='INDEX.fq|READ.fq', type=lambda x:file_check(parser, x), default=None, help='FASTQ file containing the molecular tag sequence for each read name in the corresponding SAM/BAM file (required only for CASE 1 detailed above)')
    ogroup.add_argument('-o','--out', dest='out_prefix', help='prefix of output file for sorted sam',type=str)
    ogroup.add_argument('-s','--sort', dest='sort',action='store_true', default=False, help="input sam file needs to be sorted")
    ogroup.add_argument('-l','--length', dest='length', type=int, help="length of molecular tag sequence",required=True)
    #ogroup.add_argument('-T', dest='tmp_prefix', metavar='TEMP_DIR', type=lambda x:tmp_dir_check(parser, x), default=default_tmp_dir, help='directory for reading and writing to temporary files and named pipes (default: {0})'.format(default_tmp_dir))
    #ogroup.add_argument('--old-samtools', dest='old_samtools', action='store_true', default=False, help="required for compatibility with samtools sort style in samtools versions <=0.1.19")
    #ogroup.add_argument('--rmdup-only', dest='rmdup_only', action='store_true', default=False, help="required for only outputting duplicates removed file")
    #ogroup.add_argument('--debug', dest='debug', action='store_true', default=False, help=argparse.SUPPRESS)
    #ogroup.add_argument('-log', help="log file to write statistics to (optional)")
    ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
    ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')

    args = parser.parse_args()
    sam=args.sam #name of sam file
    ID=args.length #length of index
    COUNT=1 # sam line counter
    NR=0 # record counter
    seq="[ACGTN^-]\{%d,\}" %ID #target index
    place={}
    
    if args.sort:
        sort_sam(args.sam,args.out_prefix)
        args.sam=str(args.out_prefix)
        rmdup.inter_sam(str(args.sam))
    else:
        rmdup.inter_sam(str(args.sam))
        
    with open(args.sam.split(".")[0]+"_deduped.sam",'w') as dedup:
        for value in place.values():
            dedup.write('{}'.format(value))