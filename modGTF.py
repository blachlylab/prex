#!/usr/bin/env python
"""
GTF.py
Kamil Slowikowski
December 24, 2013

Read GFF/GTF files. Works with gzip compressed files and pandas.

    http://useast.ensembl.org/info/website/upload/gff.html
"""


from collections import defaultdict
import gzip
import pandas as pd
import re

from pdb import set_trace as stop


GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            # karl modification
            if line.split('\t')[2] not in ["gene", "transcript", "start_codon"]:
                continue
            # end karl modification 
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    # karl modification
    # how to identify the principal isoforms
    PRINCIPAL_TAG = "appris_principal_1"
    # end karl modification 
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            #karl mod
            # modify how tags are parsed. 
            #     enter None into dataframe when the principal tag isn't present
            #     this permits us to drop dataframe rows where tag column is None
            #     as a way to extract the principal isoforms
            if key == "tag" and str(PRINCIPAL_TAG) not in value:
                result[key] = None
            else:
                result[key] = _get_value(value)
            #end karl mod    
            #result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None
    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value
