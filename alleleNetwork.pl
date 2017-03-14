#!/usr/bin/perl


foreach $ar (@ARGV)
{
	undef @lecturas;
	open(A,$ar);
	$anterior="";
	$read="";

	while ($l = <A>)
	{
		@c = split('\t',$l);
		$read = $c[0];
		$gen = $c[1];
		
		if ($read eq $anterior)
		{
			push(@lecturas,$gen);
		}else{
			if(@lecturas)
			{
				for($i=0;$i<=scalar(@lecturas);$i++)
				{
					for($j=i;$j<=scalar(@lecturas);$j++) 
					{
						$hash{"$lecturas[$i]\t$lecturas[$j]"}++;
						$hash{"$lecturas[$j]\t$lecturas[$i]"}++;
					}
				}
			}else{
				$hash{"$gen\t$gen"}++;			
			}
			
			undef @lecturas;
			$lecturas[0] = $gen;
		}
		$anterior = $read;
	}

	close A;
}

foreach $k (keys(%hash))
{
	print "$k\t$hash{$k}\n";
}
