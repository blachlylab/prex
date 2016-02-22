#!/usr/bin/env python

# PrEx
# Promoter Extraction
# 
# James S Blachly, MD

from __future__ import print_function

import argparse
import json
import re
import os
import sys
import subprocess
import tempfile
from distutils.version import LooseVersion 
from pdb import set_trace as stop
import modGTF as gtf 

#import bio

# CONSTANTS
ENST = 1
ENSG = 2
UCSC = 3
REFSEQ = 4
SYMBOL = 5

id_descriptions = { ENST: 'ensembl! transcript',
                    ENSG: 'ensembl! gene',
                    UCSC: 'UCSC transcript id',
                    REFSEQ: 'NCBI Refseq id',
                    SYMBOL: 'Official gene symbol' }

# TO DO: figure out what column to pull in case of UCSC or REFSEQ identifiers
id_gff3_names    = { ENST: 'transcript_id',
                    ENSG: 'gene_id',
                    UCSC: 'TBD',
                    REFSEQ: 'TBD',
                    SYMBOL: 'gene_name' }

def load_config():
    """
    Read prex.json for default FASTA and GFF3
    """
    if os.path.isfile('prex.json'):
        with open('prex.json','r') as f:
            config = json.load(f)
        return config
    else:
        return dict()

def validate_file(filename):
    '''
    Check whether the input file exists and is readable
    '''
    normpath = os.path.abspath(os.path.expanduser(filename))
    if not os.path.isfile(normpath):
        abort("File not found: " + normpath)
    else:
        try:
            open(normpath)
            info(normpath)
        except Exception as e:
            abort(str(e))
    return normpath

def decode_id(identifier):
    """
    Take gene identifier and guess whether it is 
    gene symbol, or ensembl, refseq, or UCSC gene id
    """

    if re.search('^ENST[0-9]{11}', identifier): return ENST
    elif re.search('^ENSG[0-9]{11}', identifier): return ENSG
    elif re.search('^uc[0-9]{3}[a-z]{3}\.', identifier): return UCSC
    elif re.search('^[NX][GM]_', identifier): return REFSEQ
    elif re.search('[A-Z0-9][A-Za-z0-9]{1,}', identifier): return SYMBOL
    else:
        warn("I was unable to understand your gene id: " + identifier)
        return None

def using_new_bedtools():
    '''
    determine which version of bedtools is installed.
    return true if the newest version is being used, i.e. '-fo' flag removed from getfasta
    return false if using an older version that still uses '-fo'
    '''
    user_version = subprocess.check_output(['bedtools', '--version']).strip().split()[1]
    if LooseVersion(user_version) < LooseVersion('v2.25'):
        return False 
    else:
        return True


def bedtools_cmd(chrom, start, end, identifier, fasta_in, fasta_out):
    """
    TBD

    bedtools getfasta usage/brief summary
    
    bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>

    Option  Description
    -name   Use the name column in the BED file for the FASTA headers in the output FASTA file.
    -tab    Report extract sequences in a tab-delimited format instead of in FASTA format.
    -s      Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.
    -split  Given BED12 input, extract and concatenate the sequences from the BED blocks (e.g., exons)
    """
    bed_name = identifier + ";promoter;" + chrom + ":" + str(start) + "-" + str(end)
    bed_fields = [chrom, start, end, bed_name]
    bed_line   = '\t'.join(bed_fields) + '\n'       # fails without newline
    
    with tempfile.NamedTemporaryFile(mode='w') as bedfileptr:
        bedfileptr.write(bed_line)
        bedfileptr.flush()
        cmd = ['bedtools', 'getfasta', '-name', '-s', '-fi', fasta_in, '-bed', bedfileptr.name, '-fo', fasta_out]
        if using_new_bedtools():
            # new version of getfasta doesn't have -fo option. 
            # strip the option from the command and specify output file name
            outfile = open(fasta_out,'w')
            cmd = cmd[:-2]
        else:
            # old version of bedtools requires -fo option
            # set outfile to be none and pass output file option to bedtools call
            outfile = None

        info("Running " + ' '.join(cmd))
        # subprocess call
        subprocess.call(cmd, stdout=outfile)
        if outfile: 
            outfile.close()
    return True

def do_gff3_stuff(annot, id_column, identifier, up, down):
    # Do GFF3 stuff here
    
    pAnnot = annot.dropna(subset=['tag'])
    tmp = pAnnot[pAnnot[id_column]==identifier][['feature','seqname','start','strand']].drop_duplicates()
    if len(tmp) == 0 and len(annot[annot[id_column]==identifier]) == 0:
        abort("no {0} known by this identifier: {1}".format(id_column, identifier))
    
    
    # a check for whether any features have more than 1 primary start codon or no primary start codon. 
    # code will not work correctly in these cases
    if len(tmp[tmp["feature"]=="start_codon"]) < 1:
        warn("no principal isoform found for {0}".format(identifier))
        return Region()
    elif len(tmp[tmp["feature"]=="start_codon"]) > 1:
        warn("too many primary isoforms for {0}".format(identifier))
        return Region()
    
    
    CDS = tmp[tmp["feature"]=="start_codon"]["start"].values[0]
    chrom = tmp[tmp["feature"]=="start_codon"]["seqname"].values[0]
    strand = tmp[tmp["feature"]=="start_codon"]["strand"].values[0]
    if strand == "+":
        bed_start = CDS - up 
        bed_stop  = CDS + down 
    elif strand == "-":
        bed_start = CDS - down 
        bed_stop  = CDS + up
    return Region(chrom, bed_start, bed_stop)

class Region(object):
    def __init__(self, chrom=None, start=None, stop=None):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        # can also add name, score, and strand later on if necessary
    def __bool__(self):
        if any([self.chrom, self.start, self.stop]):
            return True 
        else:
            return False 
    # to override bool in python 2.7 and/or 3
    __nonzero__=__bool__

def main():
    parser = argparse.ArgumentParser(description='Return promoter sequence for given gene', 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('identifiers', nargs='+', help='Gene identifier: Gene symbol, ensembl! gene/transcript id, Refseq gene id, UCSC gene id')
    parser.add_argument('-f', '--fasta',metavar='filename', help='(multi)FASTA file')
    parser.add_argument('-g', '--gff3', metavar='filename', help='GFF3 formatted annotation')
    parser.add_argument('-u', '--up',  metavar='nt', type=int, default=1000, help='Bases upstream of TSS')
    parser.add_argument('-d', '--down',metavar='nt', type=int, default=500,  help='Bases downstream of TSS')

    args = parser.parse_args()

    # load default config from prex.json, if one exists
    config = load_config()  # no config file returns empty dict

    # overwrite defaults with cmdline parms
    if args.fasta: config['fasta'] = args.fasta
    if args.gff3: config['gff3'] = args.gff3
    
    if 'fasta' not in config or 'gff3' not in config:
        abort("Please specify a FASTA file and GFF3 annotation\n(or define defaults in your prex.json config file)")
        
    config['fasta'] = validate_file(config['fasta'])
    config['gff3']  = validate_file(config['gff3'])   
    info("loading gff3")
    annot = gtf.dataframe(config['gff3'])
    # cast start as int
    annot.loc[annot.index ,'start'] = annot['start'].astype(int)
    info("probing genes")     

    for identifier in args.identifiers:
        # Autodetect gene identifier
        id_type = decode_id(identifier)
        # get gff3 column while we're at it
        id_column = id_gff3_names[id_type]
        if id_type:
            info(identifier + " => " + id_descriptions[id_type])
        region = do_gff3_stuff(annot, id_column, identifier, args.up, args.down)
        if not bool(region):
            print()
            continue
        chrom = region.chrom
        bed_start = region.start
        bed_stop = region.stop 
        bedtools_cmd(chrom,str(bed_start),str(bed_stop), identifier, config['fasta'],str(identifier) + '_fastaout.fa')
        print()

#
# Print warning / exit messages according to template
#
def abort(msg, retval=1):
    print("[!!] " + msg)
    exit(retval)

def warn(msg):
    print("[* ] " + msg)

def info(msg):
    print("[ok] " + msg)

if __name__ == "__main__":
    main()
