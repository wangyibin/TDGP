#!/usr/bin/perl -w

use Getopt::Std;
getopts "c:t:s:a:f:";


if (!defined $opt_c) {
    die "************************************************************************
    Usage: filterGetaAnno -c config.txt -t threads -s INT
      -h : help and usage.
      -c : config.txt
      -t : threads
      -s : step 1 or 2
      -a : augustus support cutoff (default 30), optional
      -f : FPKM cutoff (default 3), optional
************************************************************************\n";
}else{
  print "************************************************************************\n";
  print "Version 1.1\n";
  print "Copyright to Tanger\n";
  print "RUNNING...\n";
  print "************************************************************************\n";
	}

my $THREADS          = (defined $opt_t)?$opt_t:24;
my $STEP             = (defined $opt_s)?$opt_s:1;
my $AUSUPORTCUTOFF   = (defined $opt_a)?$opt_a:30;
my $FPKMCUTOFF       = (defined $opt_f)?$opt_f:3;
die "No such step defined: step $STEP\n" if($STEP!=1 and $STEP!=2);


open(IN, $opt_c) or die"";
while(<IN>){
	chomp;
	next if(/#/);
	if(/triUtil=(\S+)/){
		$triUtil = $1;
	}elsif(/PfamScan=(\S+)/){
		$PfamScan = $1;
	}elsif(/HMMDB=(\S+)/){
		$HMMDB = $1;
	}elsif(/R1=(\S+)/){
		$R1 = $1;
	}elsif(/R2=(\S+)/){
		$R2 = $1;
	}elsif(/CDS=(\S+)/){
		$CDS = $1;
	}elsif(/PEP=(\S+)/){
		$PEP = $1;
	}elsif(/GENE=(\S+)/){
		$GENE = $1;
	}elsif(/GFF=(\S+)/){
		$GFF = $1;
	}elsif(/HOMO=(\S+)/){
		$HOMO = $1;
	}elsif(/PWD=(\S+)/){
		$PWD = $1;
		}
	}
close IN;



### 0. check files
print "Checking files ...\n";
die "$triUtil not found!\n" if(!(-e $triUtil));
die "$PfamScan not found!\n" if(!(-e $PfamScan));
die "$HMMDB not found!\n" if(!(-e $HMMDB));
die "$R1 not found!\n" if(!(-e $R1));
die "$R2 not found!\n" if(!(-e $R2));
die "$CDS not found!\n" if(!(-e $CDS));
die "$PEP not found!\n" if(!(-e $PEP));
die "$GENE not found!\n" if(!(-e $GENE));
die "$GFF not found!\n" if(!(-e $GFF));
die "$HOMO not found!\n" if(!(-e $HOMO));
print "Done ...\n";

if($STEP==1){
### 1. running PfamScan
print "Running PfamScan\n";
system("date -R");
my $cmd = "perl ".$PfamScan." -fasta "." ".$PEP." -dir ".$HMMDB." -cpu ".$THREADS." -outfile pfam.out";
print "$cmd\n";
open(OUT, "> run_PfamScan.sh") or die"";
print OUT "#!/bin/bash\n";
print OUT "cd \$PBS_O_WORKDIR\n";
print OUT "$cmd\n";
close OUT;
print "Please submit run_PfamScan.sh ...\n\n\n";

### 2. running RSEM
print "Running RSEM ...\n";

system("grep '>' $CDS|sed 's/>//'|awk '{print \$1\"\t\"\$1}' > gene.map");
system("ln -s $CDS ./reference.fasta");
$cmd = $triUtil." --transcripts $PWD/reference.fasta --seqType fq --left ".$R1." --right ".$R2." --est_method RSEM --aln_method bowtie2 --gene_trans_map $PWD/gene.map --prep_reference --output_dir $PWD/rsem_outdir --thread_count ".$THREADS;
print "$cmd\n";
open(OUT, "> run_RSEM.sh") or die"";
print OUT "#!/bin/bash\n";
print OUT "cd \$PBS_O_WORKDIR\n";
print OUT "$cmd\n";
close OUT;
print "Please submit run_RSEM.sh ...\n\n\n";

print "3. BLASTP against homologs ...\n";
print "Making database for homologs ...\n";
#system("makeblastdb -in $HOMO -dbtype prot -out dbname");

my $n_parts = 20;
print "Spliting fasta into $n_parts files\n";
system("perl /share/home/stu_wangwenling/software/my_script-master/splitFa2parts.pl -i $PEP -n $n_parts");
system("rm cmd.list") if(-e "cmd.list");
open(OUT, "> cmd.list") or die"";
foreach my $i(1..$n_parts){
	my $fa = "part_".$i.".fasta";
	my $out = "blast".$i.".out.tmp";
	$cmd = "blastp -query ".$fa." -db ".$HOMO." -out ".$out." -evalue 0.001 -outfmt 6 -num_threads 4 -num_alignments 1 && rm $fa";
	print OUT "$cmd\n";
	}
close OUT;

my $n_core = int $THREADS/4;
open(OUT, "> run_BLASTP.sh") or die"";
print OUT "#!/bin/bash\n";
print OUT "cd \$PBS_O_WORKDIR\n";
print OUT "ParaFly -c cmd.list -CPU $n_core\n";
print OUT "cat blast*.out.tmp > BLASTP.OUT.TMP\n";
print OUT "rm blast*.out.tmp\n";
close OUT;

print "Please submit run_BLASTP.sh ...\n\n\n";

}elsif($STEP==2){

die "BLASTP not finished!\n" if(!(-e "BLASTP.OUT.TMP"));
die "PfamScan not finished!\n" if(!(-e "pfam.out"));
die "RSEM not finished!\n" if(!(-e "rsem_outdir/RSEM.genes.results"));

my %infordb;
open(IN, "BLASTP.OUT.TMP") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	next if($data[2]<30);
	next if($data[3]<100);
	$infordb{$data[0]}->{'BLASTP'}++;
	}
close IN;
	
open(IN, "grep -v '#' pfam.out |") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $n = @data;
	next if($n<2);
	$infordb{$data[0]}->{'PFAM'}++;
	}
close IN;
	
open(IN, "rsem_outdir/RSEM.genes.results") or die"";
<IN>;
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	$infordb{$data[0]}->{'FPKM'} = int $data[6];
	}
close IN;

open(IN, $GFF) or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my $n    = @data;
	next if($n<3);
	next if($data[2] ne "gene");
	my $g = $1 if(/ID=(\S+);/);
	   $g =~ s/;.*//g;
	if(/Augustus_transcriptSupport_percentage=(\S+);/){
		$auSuport=$1; $auSuport =~ s/;.*//g;
		$infordb{$g}->{'auSuport'} = $auSuport;
	}else{
		$infordb{$g}->{'auSuport'} = 0;
		}
	   
	}
close IN;

my $NUMKEEP = 0;
my $TOTAL   = 0;
my %keepdb = ();
open(OUT, "> gene.validate.txt") or die"";
print OUT "GeneID	BLAST	PFAM	FPKM	auSuport	Label\n";
foreach my $gene (sort keys %infordb){
	my $BLAST = (exists($infordb{$gene}->{'BLASTP'}))?$infordb{$gene}->{'BLASTP'}:0;
	my $PFAM  = (exists($infordb{$gene}->{'PFAM'}))?$infordb{$gene}->{'PFAM'}:0;
	my $FPKM  = (exists($infordb{$gene}->{'FPKM'}))?$infordb{$gene}->{'FPKM'}:0;
	my $auSuport = (exists($infordb{$gene}->{'auSuport'}))?$infordb{$gene}->{'auSuport'}:0;
	if($BLAST>1 or $PFAM>1 or $FPKM>$FPKMCUTOFF or $auSuport>$AUSUPORTCUTOFF){
		$flag = "keep";
	}else{
		$flag = "remove";
		}
	print OUT "$gene	$BLAST	$PFAM	$FPKM	$auSuport	$flag\n";
	$NUMKEEP++ if($flag eq "keep");
	$TOTAL++;
	$keepdb{$gene}++ if($flag eq "keep");
	}

close OUT;
print "Total number of genes: $TOTAL\n";
print "Number of genes retained: $NUMKEEP\n";

system("rm -rf FILTER/");
system("mkdir FILTER");

open(OUT, "> FILTER/filter.gff3") or die"";
open(IN, $GFF) or die"";
while(<IN>){
	chomp;
	if($_ eq ""){
		print OUT "$_\n";
		}elsif(/ID=(\S+)/){
		$gene = $1;
		$gene =~ s/;.*//g;
		$gene =~ s/\..*//g;
		print OUT "$_\n" if(exists($keepdb{$gene}));
	  }
	}
close IN;
close OUT;

open(OUT, "> FILTER/filter.cds") or die"";
open(IN, $CDS) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($n,$seq) = split(/\n/,$_,2);
	next if(!exists($keepdb{$n}));
	print OUT ">$n\n$seq\n";
	}
close IN;
close OUT;

open(OUT, "> FILTER/filter.pep") or die"";
open(IN, $PEP) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($n,$seq) = split(/\n/,$_,2);
	next if(!exists($keepdb{$n}));
	print OUT ">$n\n$seq\n";
	}
close IN;
close OUT;


open(OUT, "> FILTER/filter.gene") or die"";
open(IN, $GENE) or die"";
$/='>';
<IN>;
while(<IN>){
	chomp;
	my ($n,$seq) = split(/\n/,$_,2);
	next if(!exists($keepdb{$n}));
	print OUT ">$n\n$seq\n";
	}
close IN;
close OUT;
print "Done\n";
}

