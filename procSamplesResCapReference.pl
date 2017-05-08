#!/usr/bin/perl


###########################################
#
#	proSamplesResCapReference.pl is a modified pipeline for the analysis of ResCap set of samples using a uniq reference database. 
#	First it performs a mapping process with Bowtie 2, then parse the SAM file to DRS file
#	and finally performs the ./procSam.pl script to calculate the count table.
#
#	Usage:
#	  Batch Mode:	
#		./proSamplesResCapReference.pl *R1.fastq.gz		(Supose that all the reads file are in the same folder
#									 		 and the forward reads have the subfix R1)
#
#	  Single Mode:
#		./proSamplesResCapReference.pl.pl ReadsFile.R1.fastq
#		
###########################################

foreach $ar (@ARGV)
{
	$r1 = $ar;
	$r2 = $ar;
	$r2 =~ s/R1/R2/;   ### Modify if it is neccesary to fix with the raw data file names
	
	@c = split(/\./,$ar);
	$name = $c[0];
	
	system("bowtie2 -1 $r1 -2 $r2 --no-discordant --no-mixed -p 24 --fast -S mapping.tmp.sam -a -x ./db/ResCapReference --no-head --no-sq");
	system("sed '1~2d' mapping.tmp.sam | cut -f1,3,4,5,8 | grep -v '*' > $name.all.drs");
	system("./proSam.pl $name.all.drs ./db/ResCapReference.length > $name.all.csv");
	
}
