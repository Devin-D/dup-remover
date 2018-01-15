#!/usr/bin/env python3
"""Removes PCR introduced duplicate sequences from a Sam formated file. Will adjust alignement based on soft clipping amount

Requires (python 3+)

For SINGLE END/PAIRED END reads, duplicates are removed if they fulfill the following: 
a) start at the same genomic coordinate (ie. same left most maping postion adjusted for soft clipping) 
b) have the same molecular id tag (index, barcode, umi, randomer)
c) Have the same strand orientation (+,-). 
d) mapped to the same contig.
The read with the highest quality score is kept as the non-duplicate read if -q flag is set.
Reads will be checked against 96 known umi's if -u flag is set otherwise assumes randomers
All reads that are not uniquely mapped ( ie. not mapped or have secondary alignments) will not be retained.




How to use this tool:
- (Standard): User supplies One input file,
 1) A correctly formatted SAM file. The sam file should have molecular id tag in the header with no spaces(see example sam headers) 
 2) User must supply the length of the index used (ie -l 6) in the options command.
 3) if sam file is unsorted use the --sort option with -o to specify the name of the output file * assumes samtools 1.5 or better is installed.
 
 Options for dual indexes:
 1) If using dual indexes the following separtors are allowed (^ or -) to separate the two indexes in the sam header (NO SPACES).
 2) supply the combined length of both indcies in the ---length option (ie 2 indces of 6 length= 12)
    
Runtime Optimized: User supplies only one input file,
 1) SAM file that a) contains unique alignments only b) is sorted c) has a fixed length sequence containing the
    molecular tag appended to each read name.

Example sam headers:
 -indexes must be located at the end of the header With NO SPACES (see examples below)
 ex1: HWI-ST354R:351:C0UPMACXX:6:2301:2151:AAGCTC
 ex2: HWI-ST354R:351:C0UPMACXX:6:2301:2151:AAGCTC^GGCTAA (AAGCTC-GGCTAA is ok too)
 
    
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
logging.basicConfig(filename='dup_remove.log', level=logging.INFO, 
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger=logging.getLogger(__name__)

class sam_info(object):
    """ Filters a sam file line adjusts alignment start postion based on soft clipping amount, 
    converts quality into phred score, checks for strandeness and unique mapping, gets molecular ID from sam header
    
    Attributes:
        line -- single line from a sam file (str)
        start_pos -- lines left most mapping postion soft clip adjusted (int)
        qual_score --line mapq score (int)
        strand -- lines strandeness (str)
        BC -- unique index found in sam header (str)
        
    """
    
    def __init__(self,line):
        """ init
        
        arguments
        line-- a single line from a sorted sam file
        """
        self.line=line
        self.contig=self.line.split('\t')[2]
        self.start_pos=self.start_pos()
        self.qual_score=self.line.split('\t')[4]
        self.strand=self.bit_check()
        self.BC=self.barcode_get()
        
    def soft_check(self,align):
        """ function to check for soft clipping *returns Soft clip amt"""
        soft_amt=0
        if align[1]=="S" :
            soft_amt=align[0]
        elif align[2]=="S" :
            soft_amt=align[:2]
    
        return(int(soft_amt))

    #def convert_phred(self):
        #"""Converts a quality score line into a phred score *returns the sum of individual character scores* 
        #assumes phred score ASCII_BASE 33"""
        #qual=(self.line.split("\t")[10])
        #sum_qual=0
        #for letter in qual:
            #n = ord(letter) - 33
            #sum_qual+=n
            #return sum_qual
        
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
        """ looks in the sam header for the (index,barcode,umi, randomer) assumes the barcode is at the end of  the header with no spaces 
        ex: NS500451:204:HH7GHBGXY:1:11101:10281:1048:AACCCG -barcode length 6
        NS500451:204:HH7GHBGXY:1:11101:10281:1048-AACCCG^ACCGCN -barcode length 13 (include ^)
        returns the barcode. if no barcode found read will be logged in the log file"""
        global NR
        global umi_list
        BC=self.line.split()[0].split(':')[-1]
        #try:
            #BC = subprocess.check_output("cat {sam} |cut -f 1 | grep -v '^@'|sed -n '{NR}'p |grep -o '{seq}'".format(sam=args.sam,NR=NR,seq=seq), shell=True,universal_newlines=True)
            #BC=BC.strip()#get rid of newline char
        #except:
            #logger.error("Error@ line{co}: Could not find the indetifier in sam header,will not be included in output: {head}\n".format(co=NR,head=self.line.split("\t")[0]))
            #pass
        
        if args.umi:
            if BC in umi_list:
                
                BC=BC
            
            else:
                
                BC=""
                logger.error("Error@ line{co}: Not a valid umi,will not be included in output: {head}\n".format(co=NR,head=self.line.split("\t")[0]))
                
        return BC
        
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
    
#class rmdup(object):
    #""" interates a sam file collects Unique non pcr duplicate reads with best quality"""
    #def __init__(self):
        #self.self=self
        #""" init
        
        
        #"""
        
    
def inter_sam(sorted_sam):
    
    """ interates a sam file and removes pcr duplicates if removal conditions are met. assumes sorted sam file
    arguments
    sorted sam-- a sam file that has been sorted by alignment postion (ie *samtools sort, -s flag)
    """
        
    global NR
    global UnMap
    global place
    global tot
    global badBc
    logger.info("Checking sorted sam file {sam} for pcr duplicates\n".format(sam=sorted_sam))
    with open(sorted_sam)as sam:


        for line in sam:

            if line.startswith('@'): #write @ line info to output file
                with open(args.sam+"_deduped",'a') as dedup:
                    dedup.write(line)
            else:
                
                NR+=1
                #if NR%10==0: #test to make sure its running
                    #print(NR)
                line=sam_info(line)
                key=line.BC,line.start_pos,line.strand[0],line.contig #(barcode,start_pos, +or-,contig)

                if line.strand[1] is None  : # not mapped or secondary allignment or no index in header *ignore
                    UnMap +=1
                    continue
                
                elif line.BC is "": # if no barcode or unknown umi
                    badBc+=1
                    continue
                
                elif key in place: #look in dictionary 
                    
                    if args.qual: #if keeping highest mapq score
                    
                        if place[key][1] < line.qual_score : # if val is in dict but current line quality score is highest replace with current line, qual score
                            
                            place[key] = (line.line,line.qual_score)
                        
                        else:#for cases where duplicates with same quality score 1st read will be retained 
                            
                            continue
                    else: # if args.qual not set 1st read of duplicates encountered will be retained
                        continue
    
                else:  #if values is not in dict put it in with line, qual score as value 
                    
                    place[key]=(line.line,line.qual_score)
                    
        return 

class dedup_writer(object):
    """object for writing output duplicates removed sam file"""
    def __init__(self):
        self.self=self
        """__init__"""
        
    def single_write(place):
        
        """ writer for single end data. outputs contents of place dict to "inputsam"_deduped.sam in directory where script is run"""
        logger.info("writing sam formatted file with PCR duplicates removed\n")
        #with open(args.sam.split(".")[0]+"_deduped.sam",'w') as dedup:
        with open(args.sam+"_deduped",'a') as dedup:
            for value in place.values():
                dedup.write('{}'.format(value[0]))
            
    def pair_write(place):
        """ writer for paired end data. outputs contents of place dict to "inputsam"_PE_deduped.sam in directory where script was run. assumes ordered dict, reades with no matching pair will be ignored.
        checks qname of 1st entry matches rname of 2nd entry""" 
        logger.info("writing Paired End sam formatted file with PCR duplicates removed\n")
        pecount=0
        global peRemo
        #with open(args.sam.split(".")[0]+"_PE_deduped.sam",'w') as dedup:
        with open(args.sam+"_deduped",'a') as dedup:
            firstvals = [v for v in sorted(list(place.values()))] # sort the dict values by header name (paired end will be inline with each other)
            while pecount < len(place)-1:
            
                 
                first2vals=firstvals[pecount:pecount+2] #grab the values from the sorted dict 2 at a time
            
                if first2vals[0][0].split("\t")[3] == first2vals[1][0].split("\t")[7]: #pos == pnext
                    dedup.write(first2vals[0][0]+first2vals[1][0]) #write the lines (dict value is line,qual score)
                    pecount+=2
                    peRemo+=2
                    continue
                else:
                    pecount+=1
                    logger.info("non Paired end@ line{co}: Could not find the matching pair ,will not be included in output: {head}\n".format(co=(pecount),head=first2vals[0]))
                    continue
        
        
        
        
def sort_sam(sam,out):
    """ sorts a sam file by left most mapping postion * defaults to 3M temp memory storage and 28 nodes (see samtools manual -m,-@) may need adjustments based on user system"""
    with tempfile.TemporaryFile() as f: #uses a temp file
        try:
            logger.info("sorting input sam file {sam}\n".format(sam=sam))
            output=subprocess.check_output("samtools view -bS {sam} |samtools sort -m 3M -@ 28 -o {out}.sam ".format(sam=sam,out=out),shell=True,stderr=subprocess.STDOUT) #runs samtools
            logger.info("sorted input sam file {sam} and output to {out}.sam\n".format(sam=sam,out=out))
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
    pgroup.add_argument('sam', metavar='IN.sam', type=lambda x:file_check(parser, x), help='input sorted/unsorted SAM. If unsorted must specify --sort')

    ogroup = parser.add_argument_group("Options")
    ogroup.add_argument('-p','--paired-end', dest='pe', action='store_true', default=False, help="use paired end deduping. SAM alignment must contain paired end reads. read pairs  with alignments for one read of pair  will be discarded.")
    ogroup.add_argument('-o','--out', dest='out_prefix', help='prefix of output file for sorted sam (optional if input is already sorted)',type=str)
    ogroup.add_argument('-s','--sort', dest='sort',action='store_true', default=False, help="input sam file needs to be sorted")
    ogroup.add_argument('-u','--umi', dest='umi',action='store_true', default=False, help="check index against 96 know umi's")
    ogroup.add_argument('-q','--qual', dest='qual',action='store_true', default=False, help="keep the read with the highest mapq score when duplicates found")
    #ogroup.add_argument('-l','--length', dest='length', type=int, help="length of molecular tag sequence",required=True)
    ogroup.add_argument('-v','--version', action='version', version='%(prog)s '+ __version__)
    ogroup.add_argument('-h','--help',action='help', help='show this help message and exit')

    args = parser.parse_args()
    #sam=args.sam #name of sam file
    #ID=args.length #length of index
    NR=0 # record counter
    UnMap=0 #unmapped counter
    peRemo=0 # non paired end counter
    tot=0 # total reads in dict
    badBc=0 # unidentified barcode
    #seq="[ACGTN^-]\{%d,\}" %ID #target index
    place={}
    umi_list=[]
    
    if args.umi:  # make a list of known umis
        
        with open("STL96.txt") as u:
            for line in u:
                line=line.strip()
                umi_list.append(line)
                
    if args.sort: #  sorting a sam file
        sort_sam(args.sam,args.out_prefix)
        args.sam=str(args.out_prefix +'.sam')
        inter_sam(str(args.sam))
    else:
        inter_sam(str(args.sam))
      
    if not args.pe: #writing output (dups removed)
        dedup_writer.single_write(place)
        tot+=len(place)
        logger.info("Dup Remover Summary Statistics\nTotal reads processed: {NR}\nTotal unmapped or secondary allignments: {unmap}\nTotal unidentified barcodes: {badbc}\nTotal duplicates removed: {DR}\nTotal reads retained: {D}\n".format(NR=NR,D=tot,unmap=UnMap,DR=NR-UnMap-tot,badbc=badBc))
    else:
        dedup_writer.pair_write(place)
        tot+=len(place)
        logger.info("Dup Remover Summary Statistics(Paired End)\nTotal reads processed: {NR}\nTotal unmapped or secondary allignments: {unmap}\nTotal unidentified barcodes: {badbc}\nTotal duplicates removed: {DR}\nTotal non paired reads: {rm}\nTotal P.E reads retained: {D}\n".format(NR=NR,D=peRemo,unmap=UnMap,DR=NR-UnMap-tot,rm=tot-peRemo,badbc=badBc))
        
