  fasta=$1
  fastq_folder=$2
  fastq_1_suffix=$3
  fastq_2_suffix=$4
  sample=$5
  workdir=$6

  
  group=$sample
  mkdir -p $workdir/$sample
  cd $workdir/$sample

  # ******************************************
  # 1. Mapping reads with BWA-MEM, sorting
  # ******************************************
  ( $release_dir/bin/bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_folder/${sample}$fastq_1_suffix $fastq_folder/${sample}$fastq_2_suffix || echo -n 'error' ) | $release_dir/bin/sentieon util sort -r $fasta -o sorted.bam -t $nt --sam2bam -i -
  
  # ******************************************
  # 2. Metrics
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
  $release_dir/bin/sentieon plot metrics -o metrics-report.pdf gc=gc_metrics.txt qd=qd_metrics.txt mq=mq_metrics.txt isize=is_metrics.txt
  
  # ******************************************
  # 3. Remove Duplicate Reads
  # ******************************************
  $release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
  $release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 
###use Python to find non-unique mapping reads
  #$release_dir/bin/sentieon driver  -t $nt -i sorted.bam --algo Python -- /public1/user_program/sentieon-genomics-201711/bin/non_unique_mapping.py  Algo read_names.txt 
###use shard to find non-unique mapping reads
#  sh $release_dir/bin/shard.sh $nt $fasta deduped.bam read_names.txt
###use non_unique_mapping_uniform.sh, faster than shard.sh
  sh $release_dir/bin/non_unique_mapping_uniform.sh $nt $fasta deduped.bam read_names.txt
  $release_dir/bin/sentieon driver -t $nt -i deduped.bam --algo Dedup --rmdup --dup_read_name read_names.txt unique.bam
  # ******************************************
  # 4. Indel realigner
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta  -t $nt -i unique.bam --algo Realigner realigned.bam
  
  # ******************************************
  # 5. Base recalibration
  # ******************************************
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam --algo QualCal recal_data.table
  $release_dir/bin/sentieon driver -r $fasta -t $nt -i realigned.bam -q recal_data.table --algo QualCal recal_data.table.post
  $release_dir/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
  $release_dir/bin/sentieon plot bqsr -o recal_plots.pdf recal.csv