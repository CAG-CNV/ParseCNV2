#!/usr/bin/perl
# Joe Glessner

if (@ARGV<2)
{
	print "USAGE: perl PercentSamples.pl ID_Category.txt SNP_SampleIDs.txt\n";
}
$defFile= $ARGV[0];

unless( open( DEFFILE, $defFile))
    {
        print "File not found!\n";
        exit;
    }
while (defined ($line = <DEFFILE>))
{
	$line =~ s/[\r\n]//g;
	@Vals=split(/\s+/,$line);
	$IndexAdmix{$Vals[0]} = $Vals[1];
	if(!(exists($Seen{$Vals[1]})))
	{
		push(@CohortIDs, $Vals[1]);
	}
	$Seen{$Vals[1]}=1;
}
push(@CohortIDs, "Missing");
####@CohortIDs=("ARIC","CARDIA","CFS","CHS","CSSCD","FHS","JHS","MESA","Missing");

$pedFile = $ARGV[1];
unless( open( PEDFILE, $pedFile))
    {
        print "File not found!\n";
        exit;
    }


open (RESULTFILE, ">$ARGV[1]"."_CohortCounts.txt");
$i=0;
print RESULTFILE "SNP";
while (defined ($CohortIDs[$i]))
{
	print RESULTFILE "\t$CohortIDs[$i]";
	$i++;
}
print RESULTFILE "\n";

$i=0;
while (defined ($line = <PEDFILE>))
{
	chomp($line);
	$line =~ s/[\r\n]//g;
	@CNVRSamples=split(/\t/,$line);
	$countCases = () = $CNVRSamples[1] =~ / /g;
	@Vals=split(/ /,$CNVRSamples[1]);
	$CountAGRE=0;
	$Cohorts="";
	for($i=1;$i<=$countCases;$i++)
	{
		if(exists($IndexAdmix{$Vals[$i]}))
		{	
			$Cohorts.=" ".$IndexAdmix{$Vals[$i]}." ";
		}
		else
		{
			$Cohorts.=" Missing ";
		}
	}
	print RESULTFILE "$CNVRSamples[0]";
	for($j=0;$j< $#CohortIDs+1;$j++)
	{
		$countCases = () = $Cohorts =~ / $CohortIDs[$j] /g;
		print RESULTFILE "\t$countCases";
	}
		
	print RESULTFILE "\n";
	
}

close PEDFILE;
close RESULTFILE;
done
