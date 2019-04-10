#!/usr/bin/perl
# Joe Glessner

#CallRate  PCGC_PDR_7-19-12.csv
#Call Rate PCGC_Omni1_SamplesTable.txt
#F_MISS    PCGC_Omni1-Quadv1-0_BED_Miss.imiss

$pedFile = $ARGV[0];

unless( open( PEDFILE, $pedFile))
    {
        print "File not found!\n";
        exit;
    }
open(RESFILE,">$ARGV[0]_CallRates.txt");
print RESFILE "ChipID\tCallRate\n";
###Illumina LIMs->Reports->LabMangement->ProjectDetailReport CallRate  PCGC_PDR_7-19-12.csv 
###BeadChip,SSR,SSC
$line = <PEDFILE>; #Header
@Vals=split(/,/,$line);
$i=0;
while( ($Vals[$i] ne 'Call Rate') && ($i <= $#Vals) )
{
	$i++;
}
if ($Vals[$i] ne 'Call Rate') 
{
	#print "No Call Rate column header. Not PDR\n";
}
else
{
	print "MATCH PDR, ";
	$myCallRateColNum=$i;
	$myDelimit=",";
	$file="PDR";
	$i=0;
	while( ($Vals[$i] ne 'BeadChip') && ($i <= $#Vals) )
	{
	$i++;
	}
	if ($Vals[$i] ne 'BeadChip') 
	{
	print "No BeadChip column header. Not PDR";
	}
	else
	{
		$myIDColNum=$i;
		print "IDColNum=".$i.", ";
	}
}

@Vals=();
###IlluminaGenomeStudio->SamplesTable->ExportDisplayedDataToAFile Call Rate PCGC_Omni1_SamplesTable.txt
###Array Info.Sentrix ID   Array Info.Sentrix Position
@Vals=split(/\t/,$line);
$i=0;
while( ($Vals[$i] ne 'Call Rate') && ($i <= $#Vals) )
{
	$i++;
}
if ($Vals[$i] ne 'Call Rate') 
{
	#print "No Call Rate column header. Not GenomeStudio Samples Table\n";
}
else
{
	print "MATCH GenomeStudioSamplesTable, ";
	$myCallRateColNum=$i;
	$myDelimit="\t";
	$file="GSSampTable";
	$i=0;
	while( ($Vals[$i] ne 'Array Info.Sentrix ID') && ($i <= $#Vals) )
	{
	$i++;
	}
	if ($Vals[$i] ne 'Array Info.Sentrix ID') 
	{
	print "No Array Info.Sentrix ID column header. Not GenomeStudio Samples Table";
	}
	else
	{
		$myIDColNum=$i;
		print "IDColNum=".$i.", ";
	}
}

@Vals=();
###Plink: Run --missing .imiss result file (be careful about intensity only CN probes which are always no call. To be safe, run plink with --geno 0.9 option and --make-bed then run --missing on this pre-filtered bed) F_MISS    PCGC_Omni1-Quadv1-0_BED_Miss.imiss
###FID IID
@Vals=split(/\s+/,$line);
$i=0;
while( ($Vals[$i] ne 'F_MISS') && ($i <= $#Vals) )
{
	$i++;
}
if ($Vals[$i] ne 'F_MISS') 
{
	#print "No F_MISS column header. Not Plink .imiss\n";
}
else
{
	print"MATCH Plink imiss, ";
	$myCallRateColNum=$i;
	$myDelimit='\s+';
	$file="PlinkImiss";
	$i=0;
	while( ($Vals[$i] ne 'IID') && ($i <= $#Vals) )
	{
	$i++;
	}
	if ($Vals[$i] ne 'IID') 
	{
	print "No IID column header. Not Plink .imiss";
	}
	else
	{
		$myIDColNum=$i;
		print "IDColNum=".$i.", ";
	}
}
print "CallRateColNum$myCallRateColNum";


@Vals=();
###General format
###SampleID	CallRate
@Vals=split(/\s+/,$line);
$i=0;
while( ($Vals[$i] ne 'CallRate') && ($i <= $#Vals) )
{
	$i++;
}
if ($Vals[$i] ne 'CallRate') 
{
	#print "No CallRate column header. Not Plink .imiss\n";
}
else
{
	print"MATCH General format, ";
	$myCallRateColNum=$i;
	$myDelimit='\s+';
	$file="General";
	$i=0;
	while( ($Vals[$i] ne 'SampleID') && ($i <= $#Vals) )
	{
	$i++;
	}
	if ($Vals[$i] ne 'SampleID') 
	{
	print "No SampleID column header. Not General format";
	}
	else
	{
		$myIDColNum=$i;
		print "IDColNum=".$i.", ";
	}
}
print "CallRateColNum$myCallRateColNum";


@Vals=();
if(!(defined($file)))
{
	print" No match PDR GSSampTable PlinkImiss or General";
}

while (defined ($line = <PEDFILE>))
{
	chomp($line);
	$line =~ s/[\r\n]+$//g;
	@Vals=split(/$myDelimit/,$line);
	print RESFILE "$Vals[$myIDColNum]";
	if($file eq "PDR")
	{
		print RESFILE "_".$Vals[$myIDColNum+1].$Vals[$myIDColNum+2];
	}
	elsif($file eq "GSSampTable")
	{
		print RESFILE "_".$Vals[$myIDColNum+1];
	}
	elsif($file eq "PlinkImiss")
	{
		$Vals[$myCallRateColNum]=1-$Vals[$myCallRateColNum];
	}
	print RESFILE "\t$Vals[$myCallRateColNum]\n"
}
close PEDFILE;
close RESFILE;