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

p = OptionParser(__doc__)
p.add_option('-t', '--threads', type=int, default=12,
                help='the number of threads [default: %default]')
p.add_option('--exclude', 
                help='which chromosome you want to exclude, split by comma ')
p.add_option('--no_qsub', action='store_false', default=True,
                help='no qsub this jobs, only print out. [default: %default]')


opts, args = p.parse_args()
if len(args) != 2:
    sys.exit(p.print_help())

matrix, bed = args

if not op.exists(matrix):
    sys.stderr.write('ERROR: No such file of %s'%matrix)
    sys.exit()
if not op.exists(bed):
    sys.stderr.write('ERROR: No such file of %s'%bed)



ncpu = opts.threads
if opts.exclude:
    if op.exists(opts.exclude):
        exclude = "--exclude " + " ".join(i.strip() for i in open(opts.exclude) if i.strip())
    else:
        exclude = "--exclude " + " ".join(opts.exclude.split(','))
else:
    exclude = ""


out_prefix = """
bed={bed}
matrix={matrix}
ncpu={ncpu}
exclude='{exclude}'
""".format(bed=bed, matrix=matrix, ncpu=ncpu, exclude=exclude)

out = """
    

matrix_path=$(dirname $matrix)
bed_path=$(dirname $bed)
matrix=$(basename $matrix)
bed=$(basename $bed)
resolution=`echo $matrix | perl -lne '/_(\d+)[_.]/ && print $1'`

echo "Convert matrix to cool..."
/public1/home/stu_wangyibin/software/anaconda2/envs/hicexplorer/bin/hicConvertFormat -m $matrix_path/$matrix --bedFileHicpro $bed_path/$bed --inputFormat hicpro --outputFormat cool -o ${matrix%%.matrix}.cool
echo "Convert Done"
rm ${matrix%%.matrix}.ini
echo "res:$resolution
  rep1:`pwd`/${matrix%%.matrix}.cool
" > ${matrix%%.matrix}.ini

"""


hitad_command = """
hitad -O ${matrix%%.matrix}.hitad_out.txt -d ${matrix%%.matrix}.ini --logFile hitad.log -p $ncpu -W RAW  """ + exclude

hitad_command = hitad_command + "\ntad_merge.py ${matrix%%.matrix}.hitad_out.txt > ${matrix%%.matrix}.hitad_out.merged.txt"
hitad_command = hitad_command + "\ncooler dump -t bins ${matrix%%.matrix}.cool > ${matrix%%.matrix}_DI.bg"
#hitad_command = hitad_command + "\nsort -k1,1 -k2,2n ${matrix%%.matrix}_DI.bg > ${matrix%%.matrix}_DI.sorted.bg"
#hitad_command = hitad_commdan + "\nbedGraphToBigWig ${matrix%%.matrix}_DI.sorted.bg "

with open('run_{}.sh'.format(matrix.replace('.matrix', '')), 'w') as f_out:
    f_out.write(out_prefix + out + hitad_command)

if opts.no_qsub:
    os.system('qsub -pe mpi {} -j y -q all.q -cwd -S /bin/bash {}'.format(ncpu, 'run_{}.sh'.format(matrix.replace('.matrix', ''))))

