# Title     : TODO
# Objective : TODO
# Created by: YiBin
# Created on: 2019/6/17


library(optparse)
option_list = list(
    make_option(c('-m','--matrix'),type='character',default=NULL,
                help='matrix file from hicpro',metavar='character'),
    make_option(c('-w','--window'),type='numeric', default=5,
                help='window_size of TopDom',metavar='numeric'),
    make_option(c('-o','--outprefix'),type='character',default=NULL,
                help='out file prefix')
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (length(opts) != 4){
    print_help(opt_parser)
    stop("Must input -,m,w,o arguments")
}

library(TopDom)
TopDom(data=opts$matrix,window.size=opts$window,outFile=opts$outprefix)