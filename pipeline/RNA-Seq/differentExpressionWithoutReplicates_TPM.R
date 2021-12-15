#!/usr/bin/env Rscipt
# suppress warning
options(warn=-1)

library(edgeR)
library(optparse)

option_list = list(
    make_option(c('-i','--input'),type='character',default=NULL,
                help='input raw read counts without normalized! (two columns)', metavar='input')
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (length(opts) < 2){
    print_help(opt_parser)
    stop("Must input -,i arguments")
}


# bcv 
bcv <- 0.2

expr_file <- opts$input 
out_prefix <- sub(pattern='(.*)\\..*$', replacement="\\1", basename(expr_file))
expr_data = read.table(expr_file, sep='\t', header=T, row.names=1)
group <- 1:2
y <- DGEList(counts=expr_data, group=group)
tpm <- cpm(y)
# remove cpm less than 1
keep <- rowSums(cpm(y)>2) >= 1
y <- y[keep, , keep.lib.sizes=FALSE]
# normalized lib by TMM(trimmend mean of M-values)
y <- calcNormFactors(y)
# calculate the logFC and Pvalue
y_bcv <- y
et <- exactTest(y_bcv, dispersion=bcv^2)


tTags <- topTags(et, n=NULL)
out_table <- tTags$table
out_table <- out_table[ order(row.names(out_table)), ]
for (n in colnames(expr_data)) {
    out_table[n] <- tpm[row.names(out_table), n]
}
out_table <- out_table[ , -which(names(out_table) %in% c('logCPM'))]
write.table(data.frame('gene'=rownames(out_table), round(out_table, digits=2)), 
                file=paste(out_prefix, 'edgeR.tsv', sep='.'), 
                quote=F, row.names=F, sep='\t')


de_gene <- decideTestsDGE(et, p.value=1, lfc=2)
de_summary <- summary(de_gene)

write.table(de_summary, file=paste(out_prefix, 'summary', sep='.'), quote=FALSE, sep=' ')

up_gene_names <- names(de_gene@.Data[de_gene@.Data == 1, ])
down_gene_names <- names(de_gene@.Data[de_gene@.Data == -1, ])
up_gene_data <- tpm[up_gene_names, ]
down_gene_data <- tpm[down_gene_names,]

write.table(data.frame('gene'=rownames(up_gene_data), round(up_gene_data, digits=2)), 
        file=paste(out_prefix, 'up.tsv', sep='.'), row.names=F, quote=FALSE, sep='\t')
write.table(data.frame('gene'=rownames(down_gene_data), round(down_gene_data, digits=2)), 
        file=paste(out_prefix, 'down.tsv', sep='.'),  quote=FALSE, sep='\t')