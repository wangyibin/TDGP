#!/usr/bin/perl -w

use Getopt::Std;
getopts "a:b:p:c:d:t:";


if ((!defined $opt_a)|| (!defined $opt_b)  || (!defined $opt_p)) {
    die "************************************************************************
    Usage: perl $0 -a A.fasta -b B.fasta -p blast_program 
      -h : help and usage.
      -a : A.fasta, cds/protein in species A
      -b : B.fasta, cds/protein in species B
      -p : blast program, blastn/blastp
      -d : identity (optional, default is 0.6)
      -c : coverage (optional, defalut is 0.6)
      -t : threads (optional, default is 6)
      -o : output
************************************************************************\n";
}

my $coverage 	   = (defined $opt_c) ? $opt_c : 0.6;
my $identity     = (defined $opt_d) ? $opt_d : 0.6;
my $threads      = (defined $opt_t) ? $opt_t : 6;

my $opt_c = $coverage;
my $opt_d = $identity;

system("simple_blast.pl -p $opt_p -i $opt_a -d $opt_b -o AvsB.blast.out -c 6");
system("blastn_parse.pl -i AvsB.blast.out -o EAvsB.blast.out -q $opt_a -b 1 -c $opt_c -d $opt_d");

system("simple_blast.pl -p $opt_p -i $opt_b -d $opt_a -o BvsA.blast.out -c 6");
system("blastn_parse.pl -i BvsA.blast.out -o EBvsA.blast.out -q $opt_b -b 1 -c $opt_c -d $opt_d");

system("rm dbname*");

open(IN, "EAvsB.blast.out") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my ($na,$nb) = sort ($data[0],$data[1]);
	my $key  = $na.":".$nb;
	$infordb{$key}++;
	}
close IN;


open(OUT, "|sort -k 1 > table.v1.txt") or die"";
open(IN, "EBvsA.blast.out") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	my ($na,$nb) = sort ($data[0],$data[1]);
	my $key  = $na.":".$nb;
	$infordb{$key}++;
	}
close IN;

my %dataSet_one = ();
foreach my $key (keys %infordb){
	next if($infordb{$key}==1);
	my ($aa,$bb) = split(/:/,$key);
	print OUT "$aa	$bb\n";
	}
close OUT;


my %Atdb;
my %Btdb;
open(IN, "grep '>' $opt_a|sed 's/>//'|") or die"";
while(<IN>){
	chomp;
	$Atdb{$_}++;
	}
close IN;

open(IN, "grep '>' $opt_b|sed 's/>//'|") or die"";
while(<IN>){
	chomp;
	$Btdb{$_}++;
	}
close IN;


my %Ardb;
my %Brdb;
open(IN, "table.v1.txt") or die"";
while(<IN>){
	chomp;
	my ($na,$nb) = split(/\s+/,$_);
	$Ardb{$na}++;
	$Brdb{$nb}++;
	}
close IN;

open(OUT, "> table.v2.txt") or die"";
open(IN, "AvsB.blast.out") or die"";
while(<IN>){
	chomp;
	my ($a,$b) = (split/\s+/,$_)[0,1];
	next if(exists($Ardb{$a}));
	next if(exists($Brdb{$b}));
	$Ardb{$a}++; $Brdb{$b}++;
	print OUT "$a	$b\n";
	}
close IN;
close OUT;

system("cat table.v1.txt table.v2.txt > reciprocalBLAST_out.txt");

open(ON, "> onTable.list") or die"";
open(OUT, "> offTable.list") or die"";
foreach my $gene (keys %Atdb){
	print OUT "$gene\n" if(!exists($Ardb{$gene}));
	print ON "$gene\n" if(exists($Ardb{$gene}));
	}

foreach $gene (keys %Btdb){
	print OUT "$gene\n" if(!exists($Brdb{$gene}));
	print ON "$gene\n" if(exists($Brdb{$gene}));
	}
close OUT;
close ON;

system("cat $opt_a $opt_b > all.fasta");

system("perl /share/home/zhangxt/software/script/getSeqFromList.pl -l offTable.list -d all.fasta -o off.fasta");
system("perl /share/home/zhangxt/software/script/getSeqFromList.pl -l onTable.list -d all.fasta -o on.fasta");
system("simple_blast.pl -p $opt_p -i on.fasta -d off.fasta -n 1 -o onvsoff.blast.out -c $opt_t");
system("rm dbname*");
system("rm all.fasta on.fasta off.fasta");
system("rm onTable.list offTable.list");

my %matchdb;
open(IN, "onvsoff.blast.out") or die"";
while(<IN>){
	chomp;
	my ($on, $off, $iden) = (split/\s+/,$_)[0,1,2];
	$matchdb{$on} .= $off."|".$iden.",";
	}
close IN;

my $lNum = 0;
my %idb = (); ### identity database
open(OUT, "> tmp.txt") or die"";
open(IN, "reciprocalBLAST_out.txt") or die"";
while(<IN>){
	chomp;
	$lNum++;
	print OUT "$lNum	";
	my ($na,$nb) = split(/\s+/,$_);
	print OUT "$na	$nb	Paralogs:	";
	print OUT "\n" if (!exists($matchdb{$na}) and !exists($matchdb{$nb}));
	next if(!exists($matchdb{$na}) and !exists($matchdb{$nb}));
	my $line = "";
	$line .= $matchdb{$na}."," if(exists($matchdb{$na}));
	$line .= $matchdb{$nb}     if(exists($matchdb{$nb}));
	$line =~ s/,/ /g;
	my @genedb = split(/\s+/,$line);
	foreach my $para (@genedb){
		print OUT "$para,";
		my ($pgene,$pid) = split(/\|/,$para);
	  if(!exists($idb{$pgene})){
	  	$idb{$pgene}->{'id'}   = $pid;
	  	$idb{$pgene}->{'lNum'} = $lNum;
	  }else{
	  	next if($pid<$idb{$pgene}->{'id'});
	  	$idb{$pgene}->{'id'}   = $pid;
	  	$idb{$pgene}->{'lNum'} = $lNum;
	  	}		
		}
	print OUT "\n";
	}
close IN;
close OUT;

open(OUT, "> reciprocalBLAST_plus.txt") or die"";
my %orthdb = ();
open(IN, "tmp.txt") or die"";
while(<IN>){
	chomp;
	my @data = split(/\s+/,$_);
	print OUT "$data[0]	$data[1]	$data[2]	$data[3]	";
	$orthdb{$data[1]}++; $orthdb{$data[2]}++;
	print OUT "\n" if(!exists($matchdb{$data[1]}) and !exists($matchdb{$data[2]}));
	next if(!exists($matchdb{$data[1]}) and !exists($matchdb{$data[2]}));
	next if(!defined $data[4]);
	my @tmpdb = split(/,/,$data[4]);
	my %pgenedb = ();
	foreach my $i (0..$#tmpdb){
		my ($pgene,$pid) = split(/\|/,$tmpdb[$i]);
		$orthdb{$pgene}++;
		my $exp_lNum = $idb{$pgene}->{'lNum'};
		my $exp_id   = $idb{$pgene}->{'id'};
		next if($data[0] ne $exp_lNum);
		next if($pid != $exp_id);
#		my $key = $pgene."|".$pid;
		my $key = $pgene;
		$pgenedb{$key}++;
		}
	map {print OUT "$_,"} keys %pgenedb;
	print OUT "\n";
	}
close IN;

foreach my $g (sort keys %Atdb){
	next if(exists($orthdb{$g}));
	$lNum++;
	print OUT "$lNum	$g	NA\n";
	}
foreach $g (sort keys %Btdb){
	next if(exists($orthdb{$g}));
	$lNum++;
	print OUT "$lNum	NA	$g\n";
	}

close OUT;



