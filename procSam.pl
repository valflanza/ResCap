#!/usr/bin/perl

################################################
# 	
#	procSam.pl perform the counting of reads per gene, reads per kilobase per gene,
#		number of unequivocally mapped reads and alignment horizontal coverage. The output
#		is a table with the described fields as columns and the founded genes as rows.
#	procSam.pl is implemented by procSamples.pl is the standart pipeline.
#
#	Usage:
#		./procSam.pl fileWithGeneLengths.txt OutputProcSamplesfile.drs > sampleCounts.cvs
#
#
#
################################################
open(B,$ARGV[1]);
@len = <B>;
close B;

foreach $l (@len)
{
	chomp $l;
	@c = split(' ',$l);
	$longitud{$c[1]}=$c[0];
}


open(A,$ARGV[0]);
	
@prev = [0,0,0,0,0];
@actual = split('\t',<A>);
@next = split('\t',<A>);
	
if($actual[0] eq @prev[0] | $actual[0] eq @next[0])
{
	$counts{$actual[1]}++;
}else{
	$counts{$actual[1]}++;
	$countsSpecific{$actual[1]}++;
}

while($l =<A>)
{
	chomp $l;
	@prev = @actual;
	@actual = @next;
	@next = split('\t',$l);
	
	push(@{$position{$actual[1]}},$actual[2]);
	push(@{$position{$actual[1]}},$actual[4]);
	
	if($actual[0] eq @prev[0] | $actual[0] eq @next[0])
	{
		$counts{$actual[1]}++;
	}else{
		$counts{$actual[1]}++;
		$countsSpecific{$actual[1]}++;
	}
}

foreach $k (keys(%position))
{
	#print "$k\t".join(',',@{$position{$k}})."\n";
	$seq{$k} =~ s/\s//;
	$suma=0;
	undef @gen;
	#@c = split(":",$k);
	#$length = $c[-1];
	
	if (exists($longitud{$k}))
	{
		$length = $longitud{$k};
	
		for($i=0;$i<=$length;$i++)
		{
			$gen[$i]=0;
		}
	
		foreach $p (@{$position{$k}})
		{
		
			for($j=$p;($j<=$p+100 && $j<$length);$j++)
			{
				$gen[$j]=1;
			}
		}
		for($i=0;$i<scalar(@gen);$i++)
		{
			$suma += $gen[$i];
		}
		#print "$k\n";
		$coverage{$k} = $suma/$length;
	}else{
		$coverage{$k} =0;
	}
	
}


foreach $k (keys(%counts))
{
	#print "$k\n";
	#@c = split(":",$k);
	#$length = $c[-1];
	
	
	$length = $longitud{$k};
	#print "$k\t$length\n";
	$abunRel = ($counts{$k}/($length+1))*1000;
	
	if(exists($countsSpecific{$k}))
	{
		print "$k\t$counts{$k}\t$abunRel\t$countsSpecific{$k}\t$coverage{$k}\n";
	}else{
		print "$k\t$counts{$k}\t$abunRel\t0\t$coverage{$k}\n";
	}
}

