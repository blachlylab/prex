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

def load_config():
    """Read prex.json for default FASTA and GFF3"""
    if os.path.isfile('prex.json'):
        with open('prex.json','r') as f:
            config = json.load(f)
        return config
    else:
        return dict()

def validate_file(filename):
    # I suppose we should also see if file is openable
    if not os.path.isfile(os.path.abspath(os.path.expanduser(filename))):
        abort("File not found: " + filename)
    else:
        info(filename)
    return os.path.abspath(os.path.expanduser(filename))

def decode_id(identifier):
    """Take gene identifier and guess whether it is 
    gene symbol, or ensembl, refseq, or UCSC gene id"""

    if re.search('^ENST[0-9]{11}', identifier): return ENST
    elif re.search('^ENSG[0-9]{11}', identifier): return ENSG
    elif re.search('^uc[0-9]{3}[a-z]{3}\.', identifier): return UCSC
    elif re.search('^[NX][GM]_', identifier): return REFSEQ
    elif re.search('[A-Z0-9][A-Za-z0-9]{1,}', identifier): return SYMBOL
    else:
        warn("I was unable to understand your gene id: " + identifier)
        return None

def bedtools_cmd(chr, start, end, fasta_in, fasta_out):
    """TBD

    bedtools getfasta usage/brief summary
    
    bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>

    Option  Description
    -name   Use the name column in the BED file for the FASTA headers in the output FASTA file.
    -tab    Report extract sequences in a tab-delimited format instead of in FASTA format.
    -s      Force strandedness. If the feature occupies the antisense strand, the sequence will be reverse complemented. Default: strand information is ignored.
    -split  Given BED12 input, extract and concatenate the sequences from the BED blocks (e.g., exons)
"""
    
    bed_fields = [chr, start, end, 'SOMEKINDOFNAME']
    bed_line   = '\t'.join(bed_fields) + '\n'       # fails without newline
    
    with tempfile.NamedTemporaryFile(mode='w') as bedfileptr:
        bedfileptr.write(bed_line)
        bedfileptr.flush()
        cmd = ['bedtools', 'getfasta', '-name', '-s', '-fi', fasta_in, '-bed', bedfileptr.name, '-fo', fasta_out]
        info("Running " + ' '.join(cmd))
        # subprocess call
        subprocess.check_output(cmd)
    return True

def main():
    parser = argparse.ArgumentParser(description='Return promoter sequence for given gene', 
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('identifier', help='Gene identifier: Gene symbol, ensembl! gene/transcript id, Refseq gene id, UCSC gene id')
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

    # Autodetect gene identifier
    id_type = decode_id(args.identifier)
    if id_type:
        info(args.identifier + " => " + id_descriptions[id_type])

    # Do GFF3 stuff here

    bedtools_cmd('chr1','1','100', config['fasta'],'fastaout.fa')
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
