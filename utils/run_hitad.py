#!/usr/bin/env python
#! -*- coding:utf-8 -*-

"""
%prog rice_20000_iced.matrix rice_20000_abs.bed [options]

    to generate a script to execute hitad pipeline.
"""

import re
import os
import os.path as op
import sys


from optparse import OptionParser
from TDGP.apps.grid import Cluster


p = OptionParser(__doc__)
p.add_option('-t', '--threads', type=int, default=12,
                help='the number of threads [default: %default]')
p.add_option('-o', '--outdir', default='./', 
                help='the output directory of results [default: %default]')
p.add_option('--exclude', 
                help='which chromosome you want to exclude, split by comma ')
p.add_option('--no_qsub', action='store_false', default=True,
                help='no qsub this jobs, only print out. [default: %default]')


opts, args = p.parse_args()
if len(args) != 2:
    sys.exit(p.print_help())

matrix, bed= args


if not op.exists(matrix):
    sys.stderr.write('ERROR: No such file of %s'%matrix)
    sys.exit()
if not op.exists(bed):
    sys.stderr.write('ERROR: No such file of %s'%bed)
outdir = opts.outdir
if not op.exists(outdir):
    os.makedirs(outdir)
ncpu = opts.threads
if opts.exclude:
    if op.exists(opts.exclude):
        exclude = "--exclude " + " ".join(i.strip() for i in open(opts.exclude) if i.strip())
    else:
        exclude = "--exclude " + " ".join(opts.exclude.split(','))
else:
    exclude = ""


header = Cluster().get_header(name='run_{}.sh'.format(op.basename(matrix).replace('.matrix', '')), 
                                threads=ncpu)
out_prefix = """
bed={bed}
matrix={matrix}
outdir={outdir}
ncpu={ncpu}
exclude='{exclude}'
""".format(bed=bed, matrix=matrix, outdir=outdir, ncpu=ncpu, 
         exclude=exclude)

out = """
    

matrix_path=$(dirname $matrix)
bed_path=$(dirname $bed)
matrix=$(basename $matrix)
bed=$(basename $bed)
resolution=`echo $matrix | perl -lne '/_(\d+)[_.]/ && print $1'`

echo "Convert matrix to cool..."
hicConvertFormat -m $matrix_path/$matrix --bedFileHicpro $bed_path/$bed --inputFormat hicpro --outputFormat cool -o ${outdir}/${matrix%%.matrix}.cool
echo "Convert Done"
rm ${outdir}/${matrix%%.matrix}.ini
echo "res:$resolution
  rep1:`pwd`/${outdir}/${matrix%%.matrix}.cool
" > ${outdir}/${matrix%%.matrix}.ini

"""


hitad_command = """
hitad -O ${outdir}/${matrix%%.matrix}.hitad_out.txt -d ${outdir}/${matrix%%.matrix}.ini --logFile ${outdir}/hitad.log -p $ncpu -W RAW  """ + exclude

hitad_command = hitad_command + "\ntad_merge.py ${outdir}/${matrix%%.matrix}.hitad_out.txt > ${outdir}/${matrix%%.matrix}.hitad_out.merged.txt"
hitad_command = hitad_command + "\ncooler dump -t bins ${outdir}/${matrix%%.matrix}.cool > ${outdir}/${matrix%%.matrix}_DI.bg"
hitad_command = hitad_command + "\ncooler dump -t chroms ${outdir}/${matrix%%.matrix}.cool > ${outdir}/${matrix%%.matrix}.chromsizes"
hitad_command = hitad_command + "\ncut -f 1-3 ${outdir}/${matrix%%.matrix}.hitad_out.merged.txt > ${outdir}/${matrix%%.matrix}.hitad.domain"
hitad_command = hitad_command + "\npython -m TDGP.analysis.tad quickPlotTAD ${outdir}/${matrix%%.matrix}.cool "
hitad_command = hitad_command + "${outdir}/${matrix%%.matrix}.hitad.domain ${outdir}/${matrix%%.matrix}.chromsizes -o ${outdir}/quickTADPlot_results/ | parallel -j $ncpu {}"
#hitad_command = hitad_command + "\nsort -k1,1 -k2,2n ${matrix%%.matrix}_DI.bg > ${matrix%%.matrix}_DI.sorted.bg"
#hitad_command = hitad_commdan + "\nbedGraphToBigWig ${matrix%%.matrix}_DI.sorted.bg "
hitad_command = hitad_command + "\npython -m TDGP.analysis.tad plotSizeDist ${outdir}/${matrix%%.matrix}.hitad_out.merged.txt -o ${outdir}/${matrix%%.matrix}.hitad_out.merged_dist.pdf"

with open('run_{}.sh'.format(op.basename(matrix).replace('.matrix', '')), 'w') as f_out:
    f_out.write(header + out_prefix + out + hitad_command)

if opts.no_qsub:
    os.system('sed -i "s/cd $PBS_O_WORKDIR//" {}'.format(matrix.replace('.matrix', '')))
    os.system('qsub {}'.format('run_{}.sh'.format(matrix.replace('.matrix', ''))))

 