#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import sys
import logging
import argparse

import numpy as np

LOG = logging.getLogger(__name__)

__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


def read_tsv(file, sep=None):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def read_abundance(file):

    header = []
    data = {}
    n = 0

    for line in read_tsv(file, "\t"):
        n += 1
        if n ==1:
            header = line
            continue

        data[line[0]] = np.array(line[1::]).astype(float)

    return header, data


def format_list2float(alist):

    r = []

    for i in alist:
        if i == 0:
            r.append("0")
        else:
            r.append("{:.5f}".format(i))

    return r


def sort2head(file, lines=20):

    header, data = read_abundance(file)

    print("\t".join(header))
    n = 0
    for otuid, abund in sorted(data.items(), key=lambda x: sum(x[1]), reverse=True):
        abund = format_list2float(abund)
        print("%s\t%s" % (otuid, "\t".join(abund)))
        n += 1
        if n >= lines:
            break

    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input metabolic abundance table, metabolites.tsv")
    parser.add_argument("-n", "--lines",  metavar="INT", type=int, default=20,
        help="print the first K lines instead of the first 20, default=20")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
name:
    sort2head.py: Sort the abundance file and print out the first few lines.
For exmple:
    sort2head.py differential_metabolites.tsv -n 20 > differential_metabolites_top20.tsv

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    sort2head(args.input, args.lines)


if __name__ == "__main__":

    main()
