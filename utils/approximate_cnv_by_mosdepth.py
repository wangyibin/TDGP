#!/usr/bin/env python
import sys, os
import gzip
import multiprocessing
import pysam 

def help_message():
	print("Usage: python "+sys.argv[0]+" -bam <bam_list_file> -g <genome_size> -l <read_length>  -o <out_file> [-t <thread_nums>]")


def get_opts(ARGV):
	opts = {}
	if len(ARGV) < 3:
		help_message()
		sys.exit(0)

	for i in range(1, len(ARGV), 2):
		key = ARGV[i][1:]
		value = ARGV[i+1]
		if key not in opts:
			opts[key] = value
	return opts


def calc_mapped_reads_count(in_bam):
	fn = in_bam + '.read_counts.txt'
	counts = 0
	#res = os.popen("samtools view -F 12"+in_bam)
	sam = pysam.AlignmentFile(in_bam, 'rb')
	for record in sam:
		if record.is_unmapped:
			continue
		if record.is_paired:
			counts += 1
	with open(fn, 'w') as f_out:
		f_out.write(str(counts))


def calc_read_depth(in_bam, in_bed):
	fn = in_bam.split(".")[0] 
	os.system('modepth -b {} {} {}'.format(in_bed, fn, in_bam))



def calc_pipeline(in_bams):
	calc_mapped_reads_count(in_bams)
	#calc_read_depth(in_bams, in_bed)
	pass


def quick_CNV(opts):
	bam_list = []
	name_list = []
	bed_rows = 0
	print("Calculating read depth and counts")
	with open(opts['bam'], 'r') as f_in:
		for line in f_in:
			print("\tDealing %s"%line.strip())
			bam_list.append(line.strip())
			name_list.append(line.strip().split('/')[-1].split('\\')[-1].split('.')[0])
	if 't' in opts:
		t_n = int(opts['t'])
	else:
		t_n = 1

	print("Creating processes pool")
	# bed_file = opts['bed']

	pool = multiprocessing.Pool(processes=t_n)
	for bam_file in bam_list:
		res = pool.apply_async(calc_pipeline, (bam_file,))
	pool.close()
	pool.join()

	genome_size = int(opts['g'])
	read_length = int(opts['l'])
	gene_length_db = {}

	# print("Reading bed")
	# with open(bed_file, 'r') as f_in:
	# 	for line in f_in:
	# 		data = line.strip().split()
	# 		gene_name = data[3]
	# 		s_p = int(data[1])
	# 		e_p = int(data[2])
	# 		length = e_p - s_p
	# 		if length < 0:
	# 			length = -length
	# 		gene_length_db[gene_name] = length
	
	print("Approximating CNV")
	# bed_rows = len(data)
	mapped_rc = {}
	coverage_db = {}
	for i in range(0, len(bam_list)):
		fn = bam_list[i] + '.read_counts.txt'
		with open(fn, 'r') as f_in:
			for line in f_in:
				if line.strip() == '':
					continue
				mapped_rc[name_list[i]] = int(line.strip())
		#os.remove(fn)

		fn = bam_list[i].split('.')[0] + '.regions.bed.gz'
		if name_list[i] not in coverage_db:
			coverage_db[name_list[i]] = {}
		with gzip.open(fn, 'rb') as f_in:
			for line in f_in:
				line = line.decode('UTF-8')
				if line.strip() == '':
					continue
				data = line.strip().split('\t')
				s_p = int(data[1])
				e_p = int(data[2])
				gene_name = data[3]
				length = e_p - s_p
				if length < 0:
					length = -length
				gene_length_db[gene_name] = length
				
				if gene_name not in coverage_db[name_list[i]]:
					coverage_db[name_list[i]][gene_name] = float(data[-1])
		#os.remove(fn)
	
	print("Writing result")
	out_file = opts['o']
	with open(out_file, 'w') as f_out:
		f_out.write("#Sample\\CopyNumber\t")
		f_out.write('\t'.join(list(sorted(gene_length_db.keys())))+'\n')
		for name in sorted(mapped_rc.keys()):
			f_out.write(name)
			seq_depth = mapped_rc[name]*1.0*read_length/genome_size
			for gene in sorted(gene_length_db.keys()):
				copy_number = 1.0*coverage_db[name][gene]*read_length/gene_length_db[gene]/seq_depth
				f_out.write('\t'+str(copy_number))
			f_out.write('\n')
	
	print("Finished")


if __name__ == "__main__":
	opts = get_opts(sys.argv)
	necessary_paras = ['bam', 'g', 'l', 'o']
	for key in necessary_paras:
		if key not in opts:
			help_message()
			sys.exit(0)
	quick_CNV(opts)
