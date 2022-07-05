#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import logging
import argparse

from dagflow import DAG, Task, ParallelTask, do_dag

LOG = logging.getLogger(__name__)
__version__ = "1.0.0"
__author__ = ("Xingguo Zhang",)
__email__ = "invicoun@foxmail.com"
__all__ = []
QUEUE = "-q all.q,s01"
R_BIN = "/Work/pipeline/software/Base/miniconda/v4.10.3/bin/"
PYTHON_BIN = "/Work/pipeline/software/Base/miniconda/v4.10.3/bin/"
SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")

def check_path(path):

    path = os.path.abspath(path)

    if not os.path.exists(path):
        msg = "File not found '{path}'".format(**locals())
        LOG.error(msg)
        raise Exception(msg)

    return path


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


def read_tsv(file, sep="\t"):

    for line in open(file):
        line = line.strip()

        if not line or line.startswith("#"):
            continue

        yield line.split(sep)


def create_combine_heatmap_tasks(sorts, abunds, positive, negative, job_type,
                                       work_dir="", out_dir=""):
    id = "combine_heatmap"
    tasks = ParallelTask(
        id=id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
python {script}/sort2head.py {positive} -n 50 > {{sorts}}.top50.positive.tsv
python {script}/sort2head.py {negative} -n 50 > {{sorts}}.top50.negative.tsv
export PATH={rbin}:$PATH
Rscript {script}/combine_heatmap.R {{abunds}} {{sorts}}.top50.positive.tsv {{sorts}}.positive_top50
Rscript {script}/combine_heatmap.R {{abunds}} {positive} {{sorts}}.positive
cp {{sorts}}.positive.combine_heatmap.p* {{sorts}}.positive_top50.combine_heatmap.p* {out_dir}
rm {{sorts}}.top50.positive.tsv
Rscript {script}/combine_heatmap.R {{abunds}} {{sorts}}.top50.negative.tsv {{sorts}}.negative_top50
Rscript {script}/combine_heatmap.R {{abunds}} {negative} {{sorts}}.negative
cp {{sorts}}.negative.combine_heatmap.p* {{sorts}}.negative_top50.combine_heatmap.p* {out_dir}
rm {{sorts}}.top50.negative.tsv
""".format(rbin=R_BIN,
            script=SCRIPTS,
            positive=positive,
            negative=negative,
            out_dir=out_dir),
        abunds=abunds,
        sorts=sorts,
    )

    return tasks


def create_correlation_heatmap_tasks(sorts, abunds, positive, negative, job_type,
                                       work_dir="", out_dir=""):

    id = "correlation_heatmap"
    tasks = ParallelTask(
        id=id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
python {script}/sort2head.py {positive} -n 50 > {{sorts}}.top50.positive.tsv
python {script}/sort2head.py {negative} -n 50 > {{sorts}}.top50.negative.tsv

export PATH={rbin}:$PATH
Rscript {script}/correlation_heatmap.R {{sorts}}.top50.positive.tsv {{abunds}} {{sorts}}.positive_top50
Rscript {script}/correlation_heatmap.R {positive} {{abunds}} {{sorts}}.positive
cp {{sorts}}.positive.correlation_heatmap.p* {{sorts}}.positive_top50.correlation_heatmap.p* {out_dir}
Rscript {script}/correlation_heatmap.R {{sorts}}.top50.negative.tsv {{abunds}} {{sorts}}.negative_top50
Rscript {script}/correlation_heatmap.R {negative} {{abunds}} {{sorts}}.negative
cp {{sorts}}.negative.correlation_heatmap.p* {{sorts}}.negative_top50.correlation_heatmap.p* {out_dir}
rm {{sorts}}.top50.positive.tsv {{sorts}}.top50.negative.tsv
""".format(rbin=R_BIN,
            script=SCRIPTS,
            positive=positive,
            negative=negative,
            out_dir=out_dir),
        abunds=abunds,
        sorts=sorts,
    )

    return tasks


def create_ellipse_heatmap_tasks(sorts, abunds, positive, negative, job_type,
                                       work_dir="", out_dir=""):

    id = "ellipse_heatmap"
    tasks = ParallelTask(
        id=id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
python {script}/sort2head.py {positive} -n 20 > top20.positive.tsv
python {script}/sort2head.py {negative} -n 20 > top20.negative.tsv
export PATH={rbin}:$PATH
Rscript {script}/ellipse_heatmap.R top20.positive.tsv {{abunds}} {{sorts}}.positive
cp {{sorts}}.positive.ellipse_heatmap.p* {out_dir}
Rscript {script}/ellipse_heatmap.R top20.negative.tsv {{abunds}} {{sorts}}.negative
cp {{sorts}}.negative.ellipse_heatmap.p* {out_dir}
""".format(rbin=R_BIN,
            script=SCRIPTS,
            positive=positive,
            negative=negative,
            out_dir=out_dir),
        abunds=abunds,
        sorts=sorts,
    )

    return tasks


def create_spearman_correlation_tasks(sorts, abunds, positive, negative, job_type,
                                       work_dir="", out_dir=""):

    id = "spearman_correlation"
    tasks = ParallelTask(
        id=id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={python}:$PATH
python {script}/stat_metabolism_correlation.py {{abunds}} -m {positive} -p {{sorts}}.positive --plot >{{sorts}}.positive.correlation.tsv
cp {{sorts}}.positive.correlation.tsv {out_dir}
cp -rf {{sorts}}.positive {out_dir}
python {script}/stat_metabolism_correlation.py {{abunds}} -m {negative} -p {{sorts}}.negative --plot >{{sorts}}.negative.correlation.tsv
cp {{sorts}}.negative.correlation.tsv {out_dir}
cp -rf {{sorts}}.negative {out_dir}
export PATH={rbin}:$PATH
Rscript {script}/correlation_network.R {{sorts}}.negative.link.tsv {{sorts}}.negative
Rscript {script}/correlation_circos.R {{sorts}}.negative.link.tsv {{sorts}}.negative
Rscript {script}/correlation_network.R {{sorts}}.positive.link.tsv {{sorts}}.positive
Rscript {script}/correlation_circos.R {{sorts}}.positive.link.tsv {{sorts}}.positive
cp {{sorts}}.negative.network.p* {{sorts}}.positive.network.p* {out_dir}
cp {{sorts}}.negative.circos.p* {{sorts}}.positive.circos.p* {out_dir}
""".format(python=PYTHON_BIN,
            rbin=R_BIN,
            script=SCRIPTS,
            positive=positive,
            negative=negative,
            out_dir=out_dir),
        abunds=abunds,
        sorts=sorts,
    )

    return tasks


def create_cca_tasks(sorts, abunds, positive, negative, job_type,
                     work_dir="", out_dir="", group=""):
    id = "cca"
    tasks = ParallelTask(
        id=id,
        work_dir=work_dir,
        type=job_type,
        option="-pe smp 2 %s" % QUEUE,
        script="""
export PATH={rbin}:$PATH 
Rscript {script}/cca_analysis.R {{abunds}} {positive} {{sorts}}.positive
Rscript {script}/cca_analysis.R {{abunds}} {negative} {{sorts}}.negative
cp {{sorts}}.positive.plot_cca.p* {out_dir}
cp {{sorts}}.negative.plot_cca.p* {out_dir}
if [ -n "{group}" ]; then
  Rscript {script}/metabolism_cca.R {{abunds}} {positive} {group} {{sorts}}.positive >{{sorts}}.positive.log
  Rscript {script}/metabolism_cca.R {{abunds}} {negative} {group} {{sorts}}.negative >{{sorts}}.negative.log
  cp {{sorts}}.positive.cca.p* {out_dir}
  cp {{sorts}}.negative.cca.p* {out_dir}
fi
""".format(rbin=R_BIN,
            script=SCRIPTS,
            positive=positive,
            negative=negative,
            group=group,
            out_dir=out_dir),
        abunds=abunds,
        sorts=sorts,
    )

    return tasks


def run_bioyimamas(abunds, positive, negative, group, job_type,
                    concurrent=30, refresh=10, work_dir="", out_dir=""):

    abunds = check_path(abunds)
    positive = check_path(positive)
    negative = check_path(negative)
    work_dir = mkdir(work_dir)
    out_dir = mkdir(out_dir)

    if group:
        group = check_path(group)

    sorts = []
    abundans = []
    for line in read_tsv(abunds, "\t"):
        sorts.append(line[0])
        abundans.append(check_path(line[1]))

    work_dict = {
        "combine_heatmap": "01_combine_heatmap",
        "correlation_heatmap": "02_correlation_heatmap",
        "ellipse_heatmap": "03_ellipse_heatmap",
        "spearman_correlation": "04_spearman_correlation",
        "cca": "05_CCA"
    }
    for k, v in work_dict.items():
        mkdir(os.path.join(work_dir, v))
        mkdir(os.path.join(out_dir, v))

    dag = DAG("run_bioyimamas")

    combine_heatmap_tasks = create_combine_heatmap_tasks(
        sorts=sorts, 
        abunds=abundans,
        positive=positive,
        negative=negative,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["combine_heatmap"]),
        out_dir=os.path.join(out_dir, work_dict["combine_heatmap"])
    )

    dag.add_task(*combine_heatmap_tasks)

    correlation_heatmap_tasks = create_correlation_heatmap_tasks(
        sorts=sorts,
        abunds=abundans,
        positive=positive,
        negative=negative,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["correlation_heatmap"]),
        out_dir=os.path.join(out_dir, work_dict["correlation_heatmap"])
    )
    dag.add_task(*correlation_heatmap_tasks)

    ellipse_heatmap_tasks = create_ellipse_heatmap_tasks(
        sorts=sorts,
        abunds=abundans,
        positive=positive,
        negative=negative,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["ellipse_heatmap"]),
        out_dir=os.path.join(out_dir, work_dict["ellipse_heatmap"])
    )
    dag.add_task(*ellipse_heatmap_tasks)

    spearman_correlation_tasks = create_spearman_correlation_tasks(
        sorts=sorts,
        abunds=abundans,
        positive=positive,
        negative=negative,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["spearman_correlation"]),
        out_dir=os.path.join(out_dir, work_dict["spearman_correlation"])
    )
    dag.add_task(*spearman_correlation_tasks)

    cca_tasks = create_cca_tasks(
        sorts=sorts,
        abunds=abundans,
        positive=positive,
        negative=negative,
        job_type=job_type,
        work_dir=os.path.join(work_dir, work_dict["cca"]),
        out_dir=os.path.join(out_dir, work_dict["cca"]),
        group=group
    )
    dag.add_task(*cca_tasks)

    do_dag(dag, concurrent, refresh)

    return 0


def bioyimamas(args):

    run_bioyimamas(
        abunds=args.abunds,
        positive=args.positive,
        negative=args.negative,
        group=args.group,
        job_type=args.job_type,
        concurrent=args.concurrent,
        refresh=args.refresh,
        work_dir=args.work_dir,
        out_dir=args.out_dir
    )

    return 0


def add_hlep_args(parser):

    parser.add_argument("abunds", metavar="STR", type=str,
        help="Input the microbial differential species abundance table.")
    parser.add_argument("-p", "--positive", metavar="FILE", type=str, required=True,
        help="Input positive ion differential metabolites.")
    parser.add_argument("-n", "--negative",  metavar="FILE", type=str, required=True,
        help="Input negative ion differential metabolites.")
    parser.add_argument("-g", "--group", metavar="FILE", type=str, default="",
        help="Input sample grouping table, group.list.")
    parser.add_argument("--concurrent", metavar="INT", type=int, default=10,
        help="Maximum number of jobs concurrent  (default: 10)")
    parser.add_argument("--refresh", metavar="INT", type=int, default=30,
        help="Refresh time of log in seconds  (default: 30)")
    parser.add_argument("--job_type", choices=["sge", "local"], default="local",
        help="Jobs run on [sge, local]  (default: local)")
    parser.add_argument("--work_dir", metavar="DIR", type=str, default="work",
        help="Work directory (default: current directory)")
    parser.add_argument("--out_dir", metavar="DIR", type=str, default="./",
        help="Output directory (default: current directory)")

    return parser


def main():

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format="[%(levelname)s] %(message)s"
    )
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''
URL: https://github.com/zxgsy520/bioyimamas
name:
    bioyimamas.py: Microbiology and Metabolism Association Analysis

attention:
    bioyimamas.py

version: %s
contact:  %s <%s>\
        ''' % (__version__, ' '.join(__author__), __email__))

    args = add_hlep_args(parser).parse_args()

    bioyimamas(args)


if __name__ == "__main__":

    main()
