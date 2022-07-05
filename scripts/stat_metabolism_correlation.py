#!/usr/bin/env python
#coding:utf-8

import os
import re
import sys
import logging
import argparse

import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from collections import OrderedDict

try:
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    pass
LOG = logging.getLogger(__name__)

__version__ = "1.2.1"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []


COLORS = ["#0e72cc", "#6ca30f", "#f59311", "#fa4343", "#16afcc", "#85c021",
          "#d12a6a", "#0e72cc", "#6ca30f", "#f59311", "#fa4343", "#16afcc"]


def read_tsv(file, sep=None):

    for line in open(file):
        if isinstance(line, bytes):
            line = line.decode("utf-8")
        line = line.strip()

        if not line:
            continue

        yield line.split(sep)


def read_abundance(file):

    samples = []
    data = OrderedDict()
    n = 0

    for line in read_tsv(file, "\t"):
        n += 1
        if n ==1:
            samples = line[1::]
            continue

        data[line[0]] = np.array(line[1::]).astype(float)

    return samples, data


def sort_list(data, index):

    r = []

    for i in index:
        if i not in data:
            i = i.replace("_", "-")
        r.append(data[i])

    return np.array(r)


def group_variable(x, y, index):

    dx = {}
    dy = {}

    for i in range(len(index)):
        sid = index[i].replace("-", "_")
        if "_" in sid:
            group = sid.split("_")[0]
        elif "_" in sid:
            group = sid.split("-")[0]
        else:
            group = re.match("\D+", sid).group()

        if group not in dx:
            dx[group] = []
            dy[group] = []
        dx[group].append(x[i])
        dy[group].append(y[i])

    return dx, dy
       

def read_metabolites(file, samples):

    data = OrderedDict()
    n = 0

    for line in read_tsv(file, "\t"):
        n += 1
        if n ==1:
            temp_s = line[1::]
            continue

        temp_d = {}
        for i in range(len(temp_s)):
            temp_d[temp_s[i]] = float(line[i+1])
        data[line[0]] = sort_list(temp_d, samples)

    return data


def linear_fit(x, y):

    z = np.polyfit(x, y, 1)
    
    return z[0], z[1]


def mkdir(d):
    """
    from FALCON_KIT
    :param d:
    :return:
    """
    d = os.path.abspath(d)
    if not os.path.isdir(d):
        LOG.debug('mkdir {!r}'.format(d))
        os.makedirs(d)
    else:
        LOG.debug('mkdir {!r}, {!r} exist'.format(d, d))

    return d


def plot_fit(dx, dy, w, b, xl, yl, prefix="out", title=""):

    prefix = mkdir(prefix)    
    plt.rcParams['font.family'] = ['Times New Roman']
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.rc('font' , family='Times New Roman')
    
    fig = plt.figure(figsize=[6, 5])
    ax = plt.axes([0.18, 0.15, 0.70, 0.75])

    n = 0
    tx = []
    ty = []   
    for i in dx:
        x = dx[i]
        y = dy[i]
        ax.scatter(x, y, c=COLORS[n], label=i, alpha=0.2)
        n += 1
        tx += x 
        ty += y
    ax.legend(loc="center left",
        frameon=False,
        labelspacing=0.5,
        handletextpad=0.3,
        handlelength=0.2,
        columnspacing=3.0,
        bbox_to_anchor=(1.02, 0.5))

    ax.plot([min(tx), max(tx)], [w*min(tx)+b, w*max(tx)+b], 'red',
             linewidth=1.5)
    
    plt.xlabel(xl, size=12)
    plt.ylabel(yl, size=12)
    ax.set_title(title, size=10)
    plt.savefig("%s/%s_%s.pdf" % (prefix, xl, yl))
    plt.savefig("%s/%s_%s.png" % (prefix, xl, yl), dpi=300)
    plt.close()

    return 0


def stat_metabolism_correlation(abundance, metabolites, prefix, plot=False):

    samples, data1 = read_abundance(abundance)
    data2 = read_metabolites(metabolites, samples)

    fo = open("%s.link.tsv" % prefix, "w")
    fo.write("Source\tTarget\tStrand\tvalue\n")
    print("Index\tTaxonomy\tCorrelation\tP value")
    n = 0
    for i in data2:
        meta = data2[i]
        temp = 0
        for j in data1:
            abun = data1[j]
            #result = np.corrcoef(meta, abun)
            #print(result[0][1])
            rvalue, pvalue = stats.pearsonr(meta, abun)
            if abs(rvalue) < 0.8 or pvalue >= 0.05:
                continue
            print("{0}\t{1}\t{2:.6f}\t{3:.6f}".format(i, j, rvalue, pvalue))
            temp += 1
            if plot and max(abun)-min(abun) >= 0.04:
                dx, dy = group_variable(abun, meta, samples)
                w, b = linear_fit(abun, meta,)
                title = "r={0:.4f},p={1:.4f}".format(rvalue, pvalue)
                plot_fit(dx, dy, w, b, xl=j, yl=i, prefix=prefix, title=title)
            if n >= 151:
               continue
            strand = "+"
            if rvalue <= 0:
                strand = "-"
            fo.write("{0}\t{1}\t{2}\t{3:.6f}\n".format(i, j, strand, abs(rvalue)))
        if temp >= 1:
            n += 1
    fo.close()
    return 0


def add_hlep_args(parser):

    parser.add_argument("input", metavar="FILE", type=str,
        help="Input Differential Species Abundance Table, difference_genus.tsv")
    parser.add_argument("-m", "--metabolites", metavar="FILE", type=str, required=True,
        help="Input the Differential Metabolic Distribution Table, differential_metabolites.tsv")
    parser.add_argument("-p", "--prefix", metavar="STR", type=str, default="out",
        help="Output directory for plot files, default=out")
    parser.add_argument("--plot", action="store_true", default=False,
        help="Plot the fitted curve, default=False.")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
stat_metabolism_correlation.py: Species and Metabolic Correlation Analysis
For exmple:
    stat_metabolism_correlation.py difference_genus.tsv -m differential_metabolites.tsv
    stat_metabolism_correlation.py difference_genus.tsv -m differential_metabolites.tsv --plot -p out

version: %s
contact:  %s <%s>\
    ''' % (__version__, " ".join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    stat_metabolism_correlation(args.input, args.metabolites, args.prefix, args.plot)


if __name__ == "__main__":

    main()
