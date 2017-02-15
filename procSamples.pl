#!/usr/bin/perl


foreach $ar (@ARGV)
{
	$r1 = $ar;
	$r2 = $ar;
	$r2 =~ s/R1/R2/;
	
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
