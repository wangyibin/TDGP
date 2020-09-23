#!/usr/bin/env Rscript


library(ggplot2)
library(optparse)


option_list = list(
    make_option(c('-i','--input'),type='character',default=NULL,
                help='input go enrichment results from omicshare', metavar='input')
)
opt_parser = OptionParser(option_list=option_list)
opts = parse_args(opt_parser)
if (length(opts) < 2){
    print_help(opt_parser)
    stop("Must input -,i arguments")
}

inputFile <- opts$input

dat <- read.table(inputFile, sep='\t', header=T)

dat20 <- dat[order(dat$Pvalue), ][1:20, ]

xlabel <- 'Gene count'
p <- ggplot(dat20, aes(x=num, y=reorder(Descrption, num), fill=-log10(Pvalue))) + 
    geom_bar(stat = 'identity') +
    xlab(xlabel) +
    ylab("") + 
    theme(axis.text=element_text(size=16), 
            axis.title = element_text(size=20, face='bold'),
            legend.key.size = unit(0.5, 'cm'),
            legend.title = element_text(size=14, face='bold'),
            legend.text = element_text(size=12)) +
    scale_fill_gradient(low='blue', high='red') +
    expand_limits(x=0, y=0) + 
    scale_x_continuous(expand = c(0, 0))



ggsave(gsub("txt", "barplot.pdf", inputFile), p, width=13, height=10, dpi=300)
ggsave(gsub("txt", "barplot.png", inputFile), p, width=13, height=10, dpi=300)