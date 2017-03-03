#!/usr/bin/perl


###########################################
#
#	proSamples.pl is the standard pipeline for the analysis of ResCap set of samples. 
#	First it performs a mapping process with Bowtie 2, then parse the SAM file to DRS file
#	and finally performs the ./procSam.pl script to calculate the count table.
#
#	Usage:
#	  Batch Mode:	
#		./proSamples.pl *R1.fastq.gz		(Supose that all the reads file are in the same folder
#									 		 and the forward reads have the subfix R1)
#
#	  Single Mode:
#		./procSample.pl ReadsFile.R1.fastq
#		
###########################################

foreach $ar (@ARGV)
{
	$r1 = $ar;
	$r2 = $ar;
	$r2 =~ s/R1/R2/;   ### Modify if it is neccesary to fix with the raw data file names
	
	@c = split(/\./,$ar);
	$name = $c[0];
	
	system("bowtie2 -1 $r1 -2 $r2 --no-discordant --no-mixed -p 24 --fast -S mapping.tmp.sam -a -x ./db/argannot --no-head --no-sq");
	system("sed '1~2d' mapping.tmp.sam | cut -f1,3,4,5,8 | grep -v '*' > $name.abr.drs");
	system("./proSam.pl $name.abr.drs ./db/argannot.lengths > $name.abr.csv");
	
	system("bowtie2 -1 $r1 -2 $r2 --no-discordant --no-mixed -p 24 --fast -S mapping.tmp.sam -a -x ./db/relax_trim --no-head --no-sq");
	system("sed '1~2d' mapping.tmp.sam | cut -f1,3,4,5,8 | grep -v '*' > $name.rel.drs");
	system("./proSam.pl $name.rel.drs ./db/relax_trim.length > $name.rel.csv");
	
	system("bowtie2 -1 $r1 -2 $r2 --no-discordant --no-mixed -p 24 --fast -S mapping.tmp.sam -a -x ./db/bacmet --no-head --no-sq");
	system("sed '1~2d' mapping.tmp.sam | cut -f1,3,4,5,8 | grep -v '*' > $name.bac.drs");
	system("./proSam.pl $name.bac.drs ./db/bacmet.length > $name.bac.csv");
}
