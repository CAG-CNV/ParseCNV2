#!/usr/bin/perl

@ARGV>=2 or die "Usage: $0 AtoB.txt As.txt\n";

open (FAM, $ARGV[0]) or die "AtoB file $ARGV[0] not found!\n";

while (<FAM>) {
	s/[\r\n]+$//;
	my @ind = split (/\s+/, $_);
	if(exists($AtoB{$ind[0]}))
	{
		$AtoB{$ind[0]}.=",".$ind[1];
	}
	else
	{
		$AtoB{$ind[0]}=$ind[1];
	}
}

open (FAM2, $ARGV[1]) or die "A file $ARGV[1] not found!\n";
while (<FAM2>) {
	s/[\r\n]+$//;
	my @ind = split (/\s+/, $_);
	if(exists($AtoB{$ind[0]})){print $_."\t".$AtoB{$ind[0]}."\n";}
	else{print $_."\tNA\n";}
}
