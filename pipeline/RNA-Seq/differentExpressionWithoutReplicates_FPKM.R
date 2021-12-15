#!/usr/bin/env Rscipt

library(edgeR)
library(optparse)

# bcv 
bcv <- 0.2
option_list = list(
    make_option(c('-i','--input'),type='character',default=NULL,
                help='input TMM normalized FPKM (two columns)', metavar='input')
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (length(opts) < 2){
    print_help(opt_parser)
    stop("Must input -,i arguments")
}



expr_file <- opts$input 
out_prefix <- sub(pattern='(.*)\\..*$', replacement="\\1", basename(expr_file))
expr_data = read.table(expr_file, sep='\t', header=T, row.names=1)
group <- 1:2
y <- DGEList(counts=expr_data, group=group)

# remove cpm less than 1
keep <- rowSums(cpm(y)>1) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]

y_bcv <- y
et <- exactTest(y_bcv, dispersion=bcv^2)

de_gene <- decideTestsDGE(et, p.value=1, lfc=2)
de_summary <- summary(de_gene)

write.table(de_summary, file=paste(out_prefix, 'summary', sep='.'), quote=FALSE, sep=' ')

up_gene_names <- names(de_gene@.Data[de_gene@.Data == 1, ])
down_gene_names <- names(de_gene@.Data[de_gene@.Data == -1, ])
up_gene_data <- expr_data[up_gene_names, ]
down_gene_data <- expr_data[down_gene_names,]
print(colnames(up_gene_data))
write.table(up_gene_data, file=paste(out_prefix, 'up.tsv', sep='.'), col.names=NA, quote=FALSE, sep='\t')
write.table(down_gene_data, file=paste(out_prefix, 'down.tsv', sep='.'),  quote=FALSE, sep='\t')

