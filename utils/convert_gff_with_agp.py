#!/usr/bin/env python
import sys


def convert_file(in_gff3, in_agp, out_gff):
	ctg_on_chr = {}
	with open(in_agp, 'r') as fin:
		for line in fin:
			if line.strip() == '' or line[0] == "#":
				continue
			data = line.strip().split()
			if data[4] == 'U':
				continue
			chrn = data[0]
			start_pos = int(data[1])
			end_pos = int(data[2])
			ctg = data[5]
			direct = data[-1]
			ctg_on_chr[ctg] = [chrn, start_pos, end_pos, direct]
	gff_db = {}	
	with open(in_gff3, 'r') as fin:
		for line in fin:
			if line[0] == '#' or len(line.strip()) == 0:
				continue
			data = line.strip().split()
			ctg = data[0]
			sp = int(data[3])
			ep = int(data[4])
			if ctg not in ctg_on_chr:
				continue
			chrn, csp, cep, cdir = ctg_on_chr[ctg]
			if cdir == '+':
				nsp = csp+sp-1
				nep = csp+ep-1
			else:
				nep = cep-sp+1
				nsp = cep-ep+1
				data[6] = '+' if data[6] == '-' else '-'
			data[0] = chrn
			data[3] = str(nsp)
			data[4] = str(nep)
			if data[2] == 'gene':
				gsp = nsp
			info = '\t'.join(data)
			if chrn not in gff_db:
				gff_db[chrn] = {}
			if gsp not in gff_db[chrn]:
				gff_db[chrn][gsp] = []
			gff_db[chrn][gsp].append(info)
		
	with open(out_gff, 'w') as fout:
		for chrn in sorted(gff_db):
			for gsp in sorted(gff_db[chrn]):
				fout.write('\n'.join(gff_db[chrn][gsp]))
				fout.write("\n###\n")


if __name__ == "__main__":
	if len(sys.argv) < 3:
		print("Usage: python "+sys.argv[0]+" <in_gff3> <in_agp> <out_gff>")
	else:
		in_gff3, in_agp, out_gff = sys.argv[1:]
		convert_file(in_gff3, in_agp, out_gff)

