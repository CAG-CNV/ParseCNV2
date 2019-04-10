#!/usr/bin/perl
# Joe Glessner

if (@ARGV < 4){print "USAGE: perl FilterCNV.pl IDsToRemove.txt CNVCalls.txt column keep/remove\n";}

$defFile = $ARGV[0];

$pedFile = $ARGV[1];

$column = $ARGV[2] - 1;  ###Let user use base 1 numbering

$keepRemove = $ARGV[3];

unless( open( DEFFILE, $defFile))
    {
        print "File not found!\n";
        exit;
    }
while (defined ($line = <DEFFILE>))
{
	chomp($line);
	$line =~ s/[\r\n]//g;
	@Vals=split(/\s+/,$line);
	$ExcludeIDs{$Vals[0]}=1;
	#print "#$line#";
}
unless( open( PEDFILE, $pedFile))
    {
        print "File not found!\n";
        exit;
    }
if($keepRemove eq "remove")
{
	open (RESULTFILE, ">$pedFile"."_".$keepRemove."_"."$defFile");
}
elsif($keepRemove eq "keep")
{
	open (RESULTFILE, ">$pedFile"."_".$keepRemove."_"."$defFile");
}
else
{
	print "keep/remove option must be either \"keep\" or \"remove\"\n";
	exit;
}
$i=0;
while (defined ($line = <PEDFILE>))
{
	chomp($line);
	$line =~ s/[\r\n]//g;
	@Vals=split(/\s+/,$line);
	#$Vals[4]=~s/\.int\.csv//;
	#print "#$Vals[4]#";
	if (exists($ExcludeIDs{$Vals[ $column ]}))
	{
		if($keepRemove eq "keep")
		{
			print RESULTFILE "$line\n";
		}
	}
	else
	{
		if($keepRemove eq "remove")
		{
			print RESULTFILE "$line\n";
		}
	}

}

close PEDFILE;
close RESULTFILE;
done

