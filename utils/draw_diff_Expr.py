#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path as op
import seaborn as sns
import sys
import time


def draw_diff_expr(in_data, highlight_ids, out_pre, in_xlabel=""):
    if op.exists(highlight_ids):
        highlight_list = [ i.strip() for i in open(highlight_ids) if i.strip() ]
    else:
        highlight_list = highlight_ids.split(',')
    x_label = []
    start_time = time.time()
    print("Reading data")
    plt.figure(figsize=(15.2, 10.8), dpi=100)
    
    if in_xlabel :
        selected_xlabel = [ i.strip() for i in open(in_xlabel) ]
    else: 
        selected_xlabel = ""

    all_data = {}
    highlight_data = {}
    with open(in_data, 'r') as f_in:
        linen = 1
        for line in f_in:
            data = line.strip().split('\t')
            if linen == 1:
                linen += 1
                for name in data:
                    if name != '':
                        x_label.append(name)
                if selected_xlabel:
                    value_keys = selected_xlabel
                else:
                    value_keys = x_label
                        
            else:
                gene_name = data[0]
                value_dict = dict(zip(x_label,data[1:]))
                if '.'.join(gene_name.split('.')[:2]) in highlight_list:
                    highlight_data[gene_name] = []
                    for key in value_keys:
                        highlight_data[gene_name].append(float(value_dict[key]))
                else:
                    all_data[gene_name] = []
                    for key in value_keys :
                        all_data[gene_name].append(float(value_dict[key]))
    
    if selected_xlabel:
        x_label = selected_xlabel
    #plt.xlabel("TEST")
    plt.style.use('ggplot')
    plt.ylabel("centered log2(fpkm2+1)", fontsize=20)
    plt.ylim(-10, 10)
    plt.xlim(-0.2, len(x_label)-0.8)
    plt.xticks(range(0, len(x_label)), x_label, rotation=-45, ha='left')
    plt.yticks(list(range(-10, 11)))
    X_data = list(range(0, len(x_label)))

    for gene in all_data:
        plt.plot(X_data, all_data[gene], color='darkgrey')
    
    import matplotlib.colors as mcolors
    #color_list = mcolors.CSS4_COLORS
    color_list = sns.hls_palette(20, l=.5, s=.7)
    #color_list = ['red', 'green', 'blue','cyan','pink','salmon']
    i = 0
    for gene in highlight_data:
        plt.plot(X_data, highlight_data[gene], color=color_list[i], marker='.', label=gene, lw=2)
        i += 1
    plt.legend()
    
    print("Save picture")
    
    plt.savefig(out_pre+".png", bbox_inches='tight')
    plt.savefig(out_pre+".pdf", filetype='pdf', bbox_inches='tight')
    plt.savefig(out_pre+".svg", filetype='svg', bbox_inches='tight')    
    end_time = time.time()
    print("cost time: %d"%(end_time-start_time))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python "+sys.argv[0]+" <in_data> <highlight_ids> <out_pre> [in_xlabel]")
    else:
        if len(sys.argv) == 5:
            proc, in_data, highlight_ids, out_pre, in_xlabel = sys.argv
        else:
            proc, in_data, highlight_ids, out_pre = sys.argv
            in_xlabel = ""
        draw_diff_expr(in_data, highlight_ids, out_pre, in_xlabel)
