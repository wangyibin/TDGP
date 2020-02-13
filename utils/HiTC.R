#!/usr/bin/env Rscript

library(optparse)
option_list = list(
    make_option(c('-b','--bed'),type='character',default=NULL,
                    help="bed file from hicpro", metavar='character'),
    make_option(c('-m','--matrix',type='character'),default=NULL,
                help='matrix file from hicpro',metavar='character'),
    make_option(c('-n','--binsize',type='numeric',default=NULL,
                help='bin size of matrix'))
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (length(opts) != 4){
    print_help(opt_parser)
    stop("Must input -b,m,n arguments")
}
library(HiTC)
require(Cairo)
# import data
hic <- importC(opts$matrix, opts$bed)

# Quality control
CairoPDF('hitc_cqc.pdf')
par(cfow=c(1, 1))
CQC(hic, winsize = opts$binsize, dev.new=FALSE, hist.dist=FALSE)
dev.off()

# whole genome chromosome heatmap
hic_pair <- HTClist(mclapply(hic, binningC,binsize=100000, bin.adjust=FALSE, method='sum', step=1))
CairoPDF('hitc_whole_heatmap.pdf')
mapC(forcePairwise(hic_pair), maxrange=150)
dev.off()

# draw per chromosome heatmap
for (chrn in strsplit(seqlevels(hic)," ")){
    hic_chr <- reduce(hic,chr=chrn)
    hic_chr <- HTClist(mclapply(hic_chr,binningC,binsize=100000,bin.adjust=FALSE,step=1))
    CairoPDF(paste(chrn,"_hitc_heatmap.pdf",sep=""))
    mapC(forcePairwise(hic_chr))
    dev.off()
}

# A/B compartments

