#!/usr/bin/perl
##Joseph Glessner PhD
##Center for Applied Genomics, Children's Hospital of Philadelphia
##University of Pennsylvania Perelman School of Medicine
##Developed under advising of Hakon Hakonarson
##Copyright (C) GNU General Public License 2020 Joseph Glessner parsecnv@gmail.com
##Linux bash environment is best, but Windows Cygwin can also work.
if (@ARGV < 1)
{
print "\nUSAGE MESSAGE:\n
perl ParseCNV2_Insert.pl -i CNV_Calls.rawcnv -c Cases.txt -b hg19 -o myParseCNV2OutputPrefix\n";
exit;
}
GetOptions('i=s'=>\$input,'o=s'=>\$out,'b=s'=>\$build,'p=f'=>\$mergePVar,'d=i'=>\$mergeDist,'t'=>\$tdt,'m=f'=>\$maxPInclusion,'bfile=s'=>\$bfile,'qc'=>\$qc,'log=s'=>\$log,'c=s'=>\$cases,'q=s'=>\$quantitativeTrait,'batch=s'=>\$batch);
$c="sort -k1,1 -k2,2n temp/$out$input.rawcnv2 | sed 's/\\t\$//' > temp/$out$input.rawcnv2Sorted";$o=`$c`;

if (not $mergePVar)
{       
        $mergePVar=1;
}
if (not $mergeDist)
{       
        $mergeDist=1000000;
}
if (not $maxPInclusion)
{       
        $maxPInclusion=0.05;
}

#$prefix = $pedFile;
#$out .= $prefix."_PlinkP_";
#$pedFile="temp/".$prefix."CNV_Verbose.stats";

#$dir = 'temp';
#unless(-d $dir)
#{
#    mkdir $dir or die;
#}

#open (RESFILESigCaseEnrichDel, ">temp/$out"."RESFILESigCaseEnrichDel.txt");
#open (RESFILESigControlEnrichDel, ">temp/$out"."RESFILESigControlEnrichDel.txt");
#open (RESFILESigCaseEnrichDup, ">temp/$out"."RESFILESigCaseEnrichDup.txt");
#open (RESFILESigControlEnrichDup, ">temp/$out"."RESFILESigControlEnrichDup.txt");

#open (LOG, ">$out".".log");

#print $myCommandLine;
#print LOG $myCommandLine;

if (not $lowHighRisk)
{
        $lowHighRisk="high";
        #print "NOTICE: --lowHighRisk set to 'high' by default\n";
        #print LOG "NOTICE: --lowHighRisk set to 'high' by default\n";
}
if($cases ne "")
{
	open(DEL,"temp/$out"."del.assoc.fisher.model");
	open(DUP,"temp/$out"."dup.assoc.fisher.model");
}
elsif($quantitativeTrait ne "")
{
	open(DEL,"temp/$out"."del.qassoc");
        open(DUP,"temp/$out"."dup.qassoc");
}
if($tdt ne "")
{
	open(DEL,"temp/$out"."del.tdt");
	open(DUP,"temp/$out"."dup.tdt");
}
$line=<DEL>;
@Vals=split(/\s+/,$line);
$i=0;
while( ($Vals[$i] ne 'P') && ($i <= $#Vals) )
{       
        $i++;
}
if ($Vals[$i] ne 'P')
{
        print "ERROR: No P column header.\n";
        exit;
}
$myPColNum=$i;
$i=0;
while( ($Vals[$i] ne 'SNP') && ($i <= $#Vals) )
{
        $i++;
}
if ($Vals[$i] ne 'SNP')
{
        print "ERROR: No SNP column header.\n";
        exit;
}
$mySNPColNum=$i;
$i=0;
$myBetaTColNum=0;
$myORColNum=0;
while($Vals[$i] ne 'BETA' && ($Vals[$i] ne 'T' || $tdt ne "") && $Vals[$i] ne 'OR' && ($i <= $#Vals) )
{
        $i++;
}
if ($Vals[$i] ne 'BETA' && $Vals[$i] ne 'T' && $Vals[$i] ne 'OR' )
{
        print "ERROR: No BETA, T, or OR column header (at least 1 needed).\n";
        exit;
}
if($Vals[$i] eq 'BETA' || $Vals[$i] eq 'T')
{
        $myBetaTColNum=$i;#print"BETA or T matches\n";
}
if($Vals[$i] eq 'OR')
{
        $myORColNum=$i;#print"OR matches\n";
}
$myCaseEnrichStat=$Vals[$i];
$i=0;
while( ($Vals[$i] ne 'CHR') && ($i <= $#Vals) )
{
        $i++;
}
if ($Vals[$i] ne 'CHR')
{
        print "ERROR: No CHR column header.\n";
        exit;
}
$myChrColNum=$i;
$i=0;
while( ($Vals[$i] ne 'BP') && ($i <= $#Vals) )
{
        $i++;
}
if ($Vals[$i] ne 'BP')
{
        print "ERROR: No BP column header.\n";
        exit;
}
$myBpColNum=$i;
$i=0;
while( ($Vals[$i] ne 'SNP') && ($i <= $#Vals) )
{
        $i++;
}
if ($Vals[$i] ne 'SNP')
{
        print "ERROR: No SNP column header.\n";
        exit;
}
$mySnpColNum=$i;
$JustPosIndex_caseEnrichDel=0;
$JustPosIndex_controlEnrichDel=0;
#$c="rm CNVR_caseEnrichDel.txt CNVR_controlEnrichDel.txt CNVR_caseEnrichDup.txt CNVR_controlEnrichDup.txt";$o=`$c`;print"$o\n";
open(CNVR_caseEnrichDel,">temp/$out"."caseEnrichDel.txt");
open(CNVR_controlEnrichDel,">temp/$out"."controlEnrichDel.txt");
while($snpStat=<DEL>)
{
	chomp($snpStat);
	$snpStat=~s/\r\n//;
	@stat=split(/\s+/,$snpStat);
        #if($stat[$myPColNum]>$CNVRMaxP)
        #{
        #        $CNVRMaxP=$stat[$myPColNum];
        #}
        #if($stat[$myPColNum]>$CNVRMaxPCon)
        #{
        #        $CNVRMaxPCon=$stat[$myPColNum];
        #}
	if($stat[$myPColNum]<=$maxPInclusion)
	{
		if((($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T') && $stat[$myBetaTColNum] > 0) || ($myCaseEnrichStat eq 'OR' && ($stat[$myORColNum] > 1 || $stat[$myORColNum] eq 'NA')))
		{
			$caseEnrichDelSnps++;
			if($caseEnrichDelSnps==1)
			{
				$CNVRStart=$stat[$myBpColNum];
				$CNVRMaxP=$stat[$myPColNum];
                                #$CNVRP=$stat[$myPColNum];
			}
			if($stat[$myPColNum]<$CNVRP || $caseEnrichDelSnps==1)
                        {
                        	$CNVRP=$stat[$myPColNum];
				$TagSnp=$stat[$mySnpColNum];
				if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
				{
					$CNVR_BTO=$stat[$myBetaTColNum];
				}
				elsif($myCaseEnrichStat eq 'OR')
				{
					$CNVR_BTO=$stat[$myORColNum];
				}
				else
				{
					print "ERROR: BETA, T, or OR column must be present.\n";
				}
                        }
			if($stat[$myPColNum]>$CNVRMaxP)
			{
				$CNVRMaxP=$stat[$myPColNum];
			}
			#print "$stat[$myPColNum]\n";
			if($stat[$myPColNum] == 0)
			{
				#print "NA is now 1\n";
				$stat[$myPColNum]=1;
			}
			#print -(log($stat[$myPColNum])/(log(10)))."\n";
			if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist) && (-(log($stat[$myPColNum])/(log(10))) < ($lastNegLogP + $mergePVar))&&(-(log($stat[$myPColNum])/(log(10))) > ($lastNegLogP - $mergePVar)))
			{
				#print "Extend CNVR\n";
			}
			else
			{
				if($caseEnrichDelSnps>1)
				{
					#print "End CNVR\n";
					if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist))
                                        {
                                        }
                                        else
                                        {
                                                $JustPosIndex_caseEnrichDel++;
                                        }
					$CNVREnd=$lastPos;
					print CNVR_caseEnrichDel "$lastChr $CNVRStart $CNVREnd $lastCNVRP $lastCNVR_BTO $lastTagSnp $CNVRMaxP $JustPosIndex_caseEnrichDel\n"; #Cases Controls Gene SIDs RFs
					$CNVRStart=$stat[$myBpColNum];
					$CNVRP=$stat[$myPColNum];
					$CNVRMaxP=$stat[$myPColNum];
					$TagSnp=$stat[$mySnpColNum];
					if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                        {
                                                $CNVR_BTO=$stat[$myBetaTColNum];
                                        }
                                        elsif($myCaseEnrichStat eq 'OR')
                                        {
                                                $CNVR_BTO=$stat[$myORColNum];
                                        }
					#$lastChrSigCaseEnrich=$lastChr;
					#$lastPosSigCaseEnrich=$CNVREnd;
				}
			}
			$lastChr=$stat[$myChrColNum];
			$lastPos=$stat[$myBpColNum];
			$lastNegLogP=-(log($stat[$myPColNum])/(log(10)));
			$lastCNVRP=$CNVRP;
			$lastCNVR_BTO=$CNVR_BTO;
			$lastTagSnp=$TagSnp;
		}
		else
		{
			#Control Enrich
			$controlEnrichDelSnps++;
			if($controlEnrichDelSnps==1)
			{
                                $CNVRStartCon=$stat[$myBpColNum];
				$CNVRMaxPCon=$stat[$myPColNum];
                                #$CNVRP=$stat[$myPColNum];
                        }
                        if($stat[$myPColNum]<$CNVRPCon || $controlEnrichDelSnps==1)
                        {
                                $CNVRPCon=$stat[$myPColNum];
				$TagSnpCon=$stat[$mySnpColNum];
                                if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                {
                                        $CNVR_BTOCon=$stat[$myBetaTColNum];
                                }
                                elsif($myCaseEnrichStat eq 'OR')
                                {
                                        $CNVR_BTOCon=$stat[$myORColNum];
                                }
                                else
                                {
                                        print "ERROR: BETA, T, or OR column must be present.\n";
                                }
                        }
                        if($stat[$myPColNum]>$CNVRMaxPCon)
                        {
                                $CNVRMaxPCon=$stat[$myPColNum];
                        }
                        #print "$stat[$myPColNum]\n";
                        if($stat[$myPColNum] == 0)
                        {
                                #print "NA is now 1\n";
                                $stat[$myPColNum]=1;
                        }

                        if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist) && (-(log($stat[$myPColNum])/(log(10))) < ($lastNegLogPCon + $mergePVar))&&(-(log($stat[$myPColNum])/(log(10))) > ($lastNegLogPCon - $mergePVar)))
                        {
                                #print "Extend CNVR\n";
                        }
                        else
                        {
                                if($controlEnrichDelSnps>1)
                                {
                                        #print "End CNVR\n";
					if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist))
                                        {
                                        }
                                        else
                                        {
                                                $JustPosIndex_controlEnrichDel++;
                                        }
                                        $CNVREndCon=$lastPosCon;
                                        print CNVR_controlEnrichDel "$lastChrCon $CNVRStartCon $CNVREndCon $lastCNVRPCon $lastCNVR_BTOCon $lastTagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDel\n"; #Cases Controls Gene SIDs RFs
                                        $CNVRStartCon=$stat[$myBpColNum];
                                        $CNVRPCon=$stat[$myPColNum];
					$CNVRMaxPCon=$stat[$myPColNum];
					$TagSnpCon=$stat[$mySnpColNum];
					if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                        {
                                                $CNVR_BTOCon=$stat[$myBetaTColNum];
                                        }
                                        elsif($myCaseEnrichStat eq 'OR')
                                        {
                                                $CNVR_BTOCon=$stat[$myORColNum];
                                        }
					#$lastChrSigControlEnrich=$lastChrCon;
					#$lastPosSigControlEnrich=$CNVREndCon;
                                }
                        }
                        $lastChrCon=$stat[$myChrColNum];
                        $lastPosCon=$stat[$myBpColNum];	
			$lastNegLogPCon=-(log($stat[$myPColNum])/(log(10)));
			$lastCNVRPCon=$CNVRPCon;
			$lastCNVR_BTOCon=$CNVR_BTOCon;
			$lastTagSnpCon=$TagSnpCon;
		}
	}

}
#Exit Condition
if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist)) ## $stat[$myChrColNum] eq $lastChr &&
{
}
else
{
        $JustPosIndex_caseEnrichDel++;
}
$CNVREnd=$lastPos;
if($caseEnrichDelSnps eq "")
{
}
else
{
	print CNVR_caseEnrichDel "$lastChr $CNVRStart $CNVREnd $CNVRP $CNVR_BTO $TagSnp $CNVRMaxP $JustPosIndex_caseEnrichDel\n";#$CasesControls\n"; #Cases Controls Gene SIDs RFs
}
if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist))
{
}
else
{
        $JustPosIndex_controlEnrichDel++;
}
$CNVREndCon=$lastPosCon;
if($controlEnrichDelSnps eq "")
{
}
else
{
	print CNVR_controlEnrichDel "$lastChrCon $CNVRStartCon $CNVREndCon $CNVRPCon $CNVR_BTOCon $TagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDel\n"; #Cases Controls Gene SIDs RFs
}

#DUPLICATIONS
$JustPosIndex_caseEnrichDup=0;
$JustPosIndex_controlEnrichDup=0;
open(CNVR_caseEnrichDup,">temp/$out"."caseEnrichDup.txt");
open(CNVR_controlEnrichDup,">temp/$out"."controlEnrichDup.txt");
$snpStat=<DUP>; #Discard Header (assume del and dup assoc have same header)
while($snpStat=<DUP>)
{
        chomp($snpStat);
        $snpStat=~s/\r\n//;
        @stat=split(/\s+/,$snpStat);
	#if($stat[$myPColNum]>$CNVRMaxP)
       	#{
        #	$CNVRMaxP=$stat[$myPColNum];
        #}
        #if($stat[$myPColNum]>$CNVRMaxPCon)
        #{
        #        $CNVRMaxPCon=$stat[$myPColNum];
        #}
        if($stat[$myPColNum]<=$maxPInclusion)
        {
                if((($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T') && $stat[$myBetaTColNum] > 0) || ($myCaseEnrichStat eq 'OR' && ($stat[$myORColNum] > 1 || $stat[$myORColNum] eq 'NA')))
                {
                        $caseEnrichDupSnps++;
                        if($caseEnrichDupSnps==1)
                        {
                                $CNVRStart=$stat[$myBpColNum];
				$CNVRMaxP=$stat[$myPColNum];
                                #$CNVRP=$stat[$myPColNum];
                        }
                        if($stat[$myPColNum]<$CNVRP || $caseEnrichDupSnps==1)
                        {
                                $CNVRP=$stat[$myPColNum];
                                $TagSnp=$stat[$mySnpColNum];
				if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                {
                                        $CNVR_BTO=$stat[$myBetaTColNum];
                                }
                                elsif($myCaseEnrichStat eq 'OR')
                                {
                                        $CNVR_BTO=$stat[$myORColNum];
                                }
                                else
                                {
                                        print "ERROR: BETA, T, or OR column must be present.\n";
                                }
                        }
                        if($stat[$myPColNum]>$CNVRMaxP)
                        {
                                $CNVRMaxP=$stat[$myPColNum];
                        }
                       	#print "$stat[$myPColNum]\n";
                        if($stat[$myPColNum] == 0)
                        {
                                #print "NA is now 1\n";
                                $stat[$myPColNum]=1;
                        }

                        if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist) && (-(log($stat[$myPColNum])/(log(10))) < ($lastNegLogP + $mergePVar))&&(-(log($stat[$myPColNum])/(log(10))) > ($lastNegLogP - $mergePVar)))
                        {
                                #print "Extend CNVR\n";
			}
                        else
                        {
                                if($caseEnrichDupSnps>1)
                                {
                                        #print "End CNVR\n";
					if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist))
                                        {
                                        }
                                        else
                                        {
                                                $JustPosIndex_caseEnrichDup++;
                                        }
                                        $CNVREnd=$lastPos;
                                        print CNVR_caseEnrichDup "$lastChr $CNVRStart $CNVREnd $lastCNVRP $lastCNVR_BTO $lastTagSnp $lastCNVRMaxP $JustPosIndex_caseEnrichDup\n"; #Cases Controls Gene SIDs RFs
                                        $CNVRStart=$stat[$myBpColNum];
                                        $CNVRP=$stat[$myPColNum];
					$CNVRMaxP=$stat[$myPColNum];
					$TagSnp=$stat[$mySnpColNum];
					if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                        {
                                                $CNVR_BTO=$stat[$myBetaTColNum];
                                        }
                                        elsif($myCaseEnrichStat eq 'OR')
                                        {
                                                $CNVR_BTO=$stat[$myORColNum];
                                        }
					#$lastChrSigCaseEnrich=$lastChr;
					#$lastPosSigCaseEnrich=$CNVREnd;
                                }
                        }
                        $lastChr=$stat[$myChrColNum];
                        $lastPos=$stat[$myBpColNum];
                        $lastNegLogP=-(log($stat[$myPColNum])/(log(10)));
			$lastCNVRP=$CNVRP;
			$lastCNVR_BTO=$CNVR_BTO;
			$lastTagSnp=$TagSnp;
			$lastCNVRMaxP=$CNVRMaxP;
                }
                else
                {
                        #Control Enrich
                        $controlEnrichDupSnps++;
                        if($controlEnrichDupSnps==1)
                        {
                                $CNVRStartCon=$stat[$myBpColNum];
				$CNVRMaxPCon=$stat[$myPColNum];
                                #$CNVRP=$stat[$myPColNum];
                        }
                        if($stat[$myPColNum]<$CNVRPCon || $controlEnrichDupSnps==1)
                        {
				#print"OK\n";
                                $CNVRPCon=$stat[$myPColNum];
				$TagSnpCon=$stat[$mySnpColNum];
                                if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                {
                                        $CNVR_BTOCon=$stat[$myBetaTColNum];
                                }
                                elsif($myCaseEnrichStat eq 'OR')
                                {
                                        $CNVR_BTOCon=$stat[$myORColNum];
                                }
                                else
                              	{
                                        print "ERROR: BETA, T, or OR column must be present.\n";
                                }
                        }
                        if($stat[$myPColNum]>$CNVRMaxPCon)
                        {
                                $CNVRMaxPCon=$stat[$myPColNum];
                        }
                        #print "$stat[$myPColNum]\n";
                        if($stat[$myPColNum] == 0)
                        {
                                #print "NA is now 1\n";
                                $stat[$myPColNum]=1;
                        }

                        if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist) && (-(log($stat[$myPColNum])/(log(10))) < ($lastNegLogPCon + $mergePVar))&&(-(log($stat[$myPColNum])/(log(10))) > ($lastNegLogPCon - $mergePVar)))
                        {
                                #print "Extend CNVR\n";
                        }
                        else
                        {
                                if($controlEnrichDupSnps>1)
                                {
					#print "End CNVR\n";
					if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist))
                                        {
                                        }
                                        else
                                        {
                                                $JustPosIndex_controlEnrichDup++;
                                        }
                                        $CNVREndCon=$lastPosCon;
                                        print CNVR_controlEnrichDup "$lastChrCon $CNVRStartCon $CNVREndCon $lastCNVRPCon $lastCNVR_BTOCon $lastTagSnpCon $lastCNVRMaxPCon $JustPosIndex_controlEnrichDup\n"; #Cases Controls Gene SIDs RFs
                                        $CNVRStartCon=$stat[$myBpColNum];
                                        $CNVRPCon=$stat[$myPColNum];
					$CNVRMaxPCon=$stat[$myPColNum];
					$TagSnpCon=$stat[$mySnpColNum];
					if($myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T')
                                        {
                                                $CNVR_BTOCon=$stat[$myBetaTColNum];
                                        }
                                        elsif($myCaseEnrichStat eq 'OR')
                                        {
                                                $CNVR_BTOCon=$stat[$myORColNum];
                                        }
					#$lastChrSigControlEnrich=$lastChrCon;
					#$lastPosSigControlEnrich=$CNVREndCon;
                                }
                        }
                        $lastChrCon=$stat[$myChrColNum];
                        $lastPosCon=$stat[$myBpColNum];
			#print "$stat[$myPColNum]\n";
                        $lastNegLogPCon=-(log($stat[$myPColNum])/(log(10)));
			$lastCNVRPCon=$CNVRPCon;
			$lastCNVR_BTOCon=$CNVR_BTOCon;
			$lastTagSnpCon=$TagSnpCon;
			$lastCNVRMaxPCon=$CNVRMaxPCon;
                }
        }

}
#Exit Condition
if($stat[$myChrColNum] eq $lastChr && $stat[$myBpColNum] < ($lastPos + $mergeDist))
{
}
else
{
        $JustPosIndex_caseEnrichDup++;
}
$CNVREnd=$lastPos;
if($caseEnrichDupSnps eq "")
{
}
else
{
	print CNVR_caseEnrichDup "$lastChr $CNVRStart $CNVREnd $CNVRP $CNVR_BTO $TagSnp $CNVRMaxP $JustPosIndex_caseEnrichDup\n"; #Cases Controls Gene SIDs RFs
}
if($stat[$myChrColNum] eq $lastChrCon && $stat[$myBpColNum] < ($lastPosCon + $mergeDist))
{
}
else
{
        $JustPosIndex_controlEnrichDup++;
}
$CNVREndCon=$lastPosCon;
if($controlEnrichDupSnps eq "")
{
}
else
{
	print CNVR_controlEnrichDup "$lastChrCon $CNVRStartCon $CNVREndCon $CNVRPCon $CNVR_BTOCon $TagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDup\n"; #Cases Controls Gene SIDs RFs
}

$c="grep . temp/$out"."caseEnrichDel.txt temp/$out"."controlEnrichDel.txt temp/$out"."caseEnrichDup.txt temp/$out"."controlEnrichDup.txt | sed '/:[\\ ]*\$/d' | sed 's/temp\\///' | sed 's/$out//' | sed 's/D/\td/' | sed 's/:/\t/' | sed 's/.txt//' | sed 's/Enrich//' | sort -k6,6g | awk '!_[\$1\$2\$10]++' | awk '{print \$3,\$4,\$5,\$6,\$7,\$8,\$9,\$1,\$2}' | sed 's/\ /\t/g' > temp/$out"."Report.txt";$o=`$c`;

open(REPORT,"temp/$out"."Report.txt");
open(REPORT2,">temp/$out"."Report2.txt");
open(ChrPosRanges,">temp/$out"."ChrPosRanges");
while($line=<REPORT>)
{
	chomp($line);
	@a=split(/\t/,$line);
	#print "$a[5]\n";
	$c="echo $a[5] > temp/$out"."TagSnp; $MyDirectoryPathPrefix"."PerlModules/./plink --bfile temp/$out$a[8] --extract temp/$out"."TagSnp --recode --allow-no-sex --out temp/$out"."TagSnp";$o=`$c`;#print "$c\n$o\n";
	$c="grep -v \"2 2\$\" temp/$out"."TagSnp.ped | awk 'BEGIN{cases=0}{if(\$6==2)cases++}END{ORS=\"\";print cases}'";$Cases=`$c`;
	$c="grep -v \"2 2\$\" temp/$out"."TagSnp.ped | awk 'BEGIN{controls=0}{if(\$6==1)controls++}END{ORS=\"\";print controls}'";$Controls=`$c`;
	$c="grep -v \"2 2\$\" temp/$out"."TagSnp.ped | awk '{print \$2\"\\t\"\$6}' | awk '{if(\$2==2)printf \$1\",\"}' | sed 's/,\$//'";$CaseSIDs=`$c`;
	$c="grep -v \"2 2\$\" temp/$out"."TagSnp.ped | awk '{print \$2\"\\t\"\$6}' | awk '{if(\$2==1)printf \$1\",\"}' | sed 's/,\$//'";$ControlSIDs=`$c`;
	if($Cases==0)
	{$CasesSIDs="NA";}
	if($Controls==0)
	{$ControlSIDs="NA";}
	print REPORT2 "$line\t$Cases\t$Controls\t$CaseSIDs\t$ControlSIDs\n";
	print ChrPosRanges "$a[0]:$a[1]-$a[2]\n";
}
$c="perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_knownGene.txt -knowngene "."-kgxref ".$MyDirectoryPathPrefix."GeneRef/".$build."_kgXref.txt -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/'| sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_SigCallsGene_Col23.txt";$o=`$c`;
#$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/".$out."ChrPosRanges_DEL.txt ".$MyDirectoryPathPrefix."GeneRef/".$build."_knownGene.txt -knowngene "."-kgxref ".$MyDirectoryPathPrefix."GeneRef/".$build."_kgXref.txt -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*\\t\\t//' | sed 's/^chr[^\\t]*\\t//' > temp/".$out."SigCallsGene_DEL_Col23.txt";

#RFs: SegDups (count, max, avg) DgvEntries      TeloCentro      AvgGC   AvgProbes       Recurrent       PopFreq PenMaxP_Freq_HighFreq           FreqInflated    Sparse  ABFreq  AvgConf AvgLength
$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_genomicSuperDups_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDup_Col23.txt";$o=`$command`;
$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_dgv_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk -F\"\\t\" '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDupDgv_Col23.txt";$o=`$command`;
$defFile= $MyDirectoryPathPrefix."GeneRef/".$build."_cytoBand.txt";
#Do once before as part of GeneRef setup #$c="join <(awk '{if(\$2==0)print}' GeneRef/hg19_cytoBand.txt | sort -k1,1) <(grep acen GeneRef/hg19_cytoBand.txt | sort -u -k1,1 | sort -k1,1) | join - <(grep acen GeneRef/hg19_cytoBand.txt | sort -k3,3nr | sort -u -k1,1 | sort -k1,1) | join - <(sort -k3,3nr GeneRef/hg19_cytoBand.txt | sort -u -k1,1 | sort -k1,1) | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttelomere\\n\"\$1\"\\t\"\$6\"\\t\"\$7\"\\tcentromere\\n\"\$1\"\\t\"\$10\"\\t\"\$11\"\\tcentromere\\n\"\$1\"\\t\"\$14\"\\t\"\$15\"\\ttelomere\"}' | sed 's/^chr//' > ".$MyDirectoryPathPrefix."GeneRef/".$build."_cytoBand_TeloCentro.bed";$o=`$c`;
$c="chmod +x ".$MyDirectoryPathPrefix."PerlModules/bedtools; ".$MyDirectoryPathPrefix."PerlModules/bedtools intersect -a temp/$out"."Report2.txt -b ".$MyDirectoryPathPrefix."GeneRef/".$build."_cytoBand_TeloCentro.bed -wo | awk '{print \$1\":\"\$2\"-\"\$3\"\\t\"\$(NF-1)}' > temp/".$out."CNVR_TeloCentro.txt";$o=`$c`;
$c="awk '{print \$1\":\"\$2\"-\"\$3}' temp/$out"."Report2.txt > a;perl $MyDirectoryPathPrefix"."PerlModules/Vlookup.pl temp/".$out."CNVR_TeloCentro.txt a | awk '{print \$2}' > b; mv b temp/".$out."CNVR_TeloCentro.txt";$o=`$c`;


$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_gc5Base_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDupDgvTcGc_Col23.txt";$o=`$command`;
$c="wc -l temp/$out"."fam | awk '{print \$1}'";$totalSamples=`$c`;
$c="awk -F\"\\t\" '{print (\$(NF-3)+\$(NF-2))/$totalSamples}' temp/$out"."Report2.txt > temp/".$out."CNVR_PopFreq.txt";$o=`$c`;
#peninsula
#inflated
$c="awk -F\"\\t\" '{print \$12}' temp/$out"."Report2.txt | tr ',' '\n' | sort | grep -v ^\$ | uniq -c > temp/$out"."CaseAndControlSIDs.txt; awk -F\"\\t\" '{print \$13}' temp/$out"."Report2.txt | tr ',' '\n'  | sort | grep -v ^\$ | uniq -c >> temp/$out"."CaseAndControlSIDs.txt; sort -k1,1nr temp/$out"."CaseAndControlSIDs.txt | head | awk '{print $2}' > temp/$out"."CaseAndControlSIDs_10MostInflated.txt";$o=`$c`;
#cytoband
$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_cytoBand_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDupDgvTcGcPerPenInfCyt_Col23.txt";$o=`$command`;
#recurrent
$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_somaticRearrangement_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDupDgvTcGcPerPenInfCytSom_Col23.txt";$o=`$command`;
#exon
$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/$out"."ChrPosRanges ".$MyDirectoryPathPrefix."GeneRef/".$build."_knownGene_Exons_SimFormat_AllCol_UniqueIDs -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/".$out."CNVR_ALL_CDupDgvTcGcPerPenInfCytSomExo_Col23.txt";$o=`$command`;

#PercentSamples
#$c="perl OutputUtilities/PercentSamples.pl ID_Category.txt SNP_SampleIDs.txt";$o=`$c`;

$c="paste temp/$out"."Report2.txt temp/$out"."CNVR* > temp/$out"."Report3.txt";$o=`$c`;
$c="sort -k4,4g temp/$out"."Report3.txt > temp/$out"."Report4.txt";$o=`$c`;
$c="cat temp/$out"."Report4.txt > temp/$out"."Report5.txt";$o=`$c`; #was cat header but now later

#printf( "%-10s", "hello" );
#printf( "%10d", 777);
#FAIL=CNVRLength/AvgLength<1000, DgvFreq>0.01/DgvEntries>10, MaxP>0.5, PopFreq>0.01, SegDups>10, Recurrent, FreqInflated>0.5, AvgConf<10
#Sort RF <=3, delP<5x10^-4 and ORDup>1 OR dupP<5x10^-4 and ORDup>1, (on exon)
open(INFLATED,"temp/$out"."CaseAndControlSIDs_10MostInflated.txt");
while($line=<INFLATED>)
{
	chomp($line);
        @a=split(/\s+/,$line);
	$InflatedSamples{$a[2]}=1;
}
open(REPORT,"temp/$out"."Report5.txt");
open(REPORT2,">temp/$out"."Report6.txt");
open(SNP_IDS,">temp/$out"."Snp_IDs.txt");
open(REPORT_BRIEF,">$out"."Report.txt");
printf REPORT_BRIEF ("%3s %11s %10s %10s %6s %8s %8s %6s %5s %5s\n","chr","start($build)","stop","p",$myCaseEnrichStat,"cases","controls","filter","type","tag");
$DgvEntries=0;
$SegDupsEntries=0;
while($line=<REPORT>)
{
	chomp($line);
	@a=split(/\t/,$line);
	$lengthTotal=$h_lengthTotal{$a[5]};
	if($a[9]+$a[10]>0)
	{
		#print"OK\n";
		$AvgLength=$lengthTotal/($a[9]+$a[10]);
	}
	else
	{
		$AvgLength=-9;
	}

	#DGV
	if(!($a[14]=~/-closest/))
	{
		@Dgv=split(/,/,$a[14]);
		$DgvEntries=$#Dgv+1;
	}
	if(!($a[13]=~/-closest/))
        {
                @SegDups=split(/,/,$a[13]);
                $SegDupsEntries=$#SegDups+1;
        }
	if($a[17] eq "NA-closest")
	{	
		$recurrent="N";
	}
	else
	{
		$recurrent="Y";
	}
	@caseIDs=split(/,/,$a[11]);
	@controlIDs=split(/,/,$a[12]);
	for($i=0;$i<=$#caseIDs;$i++)
	{
		if(exists($InflatedSamples{$caseIDs[$i]}))
		{
			$InflatedSamplesCount++;
		}
	}
	for($i=0;$i<=$#controlIDs;$i++)
        {
                if(exists($InflatedSamples{$controlIDs[$i]}))
                {
                        $InflatedSamplesCount++;
                }
        }
	if($a[9]+$a[10]>0)
	{
		$InflatedSamplesFreq=$InflatedSamplesCount/($a[9]+$a[10]);
	}
	else
	{
		$InflatedSamplesFreq=-9;
	}
	$confidenceTotal=$h_confidenceTotal{$a[5]};
        if($a[9]+$a[10]>0)
        {
                #print"OK\n";
                $AvgConfidence=$confidenceTotal/($a[9]+$a[10]);
        }
        else
        {
                $AvgConfidence=-9;
        }
	if($AvgLength<1000){$RF++;}
	if($DgvEntries>10){$RF++;}
	if($a[6]>0.5){$RF++;}#MaxP
	if($a[19]>0.01){$RF++;}#PopFreq
	if($SegDupsEntries>10){$RF++;}
	if($recurrent eq "Y"){$RF++;}
	if($InflatedSamplesFreq>0.5){$RF++;}
	if($AvgConfidence<10){$RF++;}
	#if($a[3]<0.0005 && $a[7] eq "case"){$RF--;}#p and direction green flag
	if(!($a[18]=~"-closest")){$RF--;}#exon green flag
	
	if($RF>2){$RF_PassFail="FAIL";$overallFail++;}
	else{$RF_PassFail="PASS";}	
	
	print REPORT2 "$line\t$AvgLength\t$DgvEntries\t$a[6]\t$a[19]\t$SegDupsEntries\t$recurrent\t$InflatedSamplesFreq\t$AvgConfidence\t$RF\t$RF_PassFail\n";
	
	printf REPORT_BRIEF ("%3d %11d %10d %10s %6s %8d %8d %6s %5s %5s\n",$a[0],$a[1],$a[2],$a[3],$a[4],$a[9],$a[10],$RF_PassFail,$a[8],$a[5]);
	
	$InflatedSamplesCount=0;
	$countLines++;
	$a[11]=~s/,/\ /g;
	$a[12]=~s/,/\ /g;
		print SNP_IDS "$a[5]\t";
		if($a[11] ne "")
		{
			print SNP_IDS " $a[11]";
		}
		if($a[12] ne "")
		{
			print SNP_IDS " $a[12]";
		}
		print SNP_IDS "\n";
	$DgvEntries=0;
	$SegDupsEntries=0;
	$RF=0;
}
#$countLines--;
$c="sort -k1,1 -k2,2n temp/$out$inputNoPath.rawcnv2 | sed 's/\\t\$//' > temp/$out$inputNoPath.rawcnv2Sorted";$o=`$c`;
$c="awk '{print \$10\"\\t\"\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9}' $out"."Report.txt | sed 's/_/\\t/' | awk '{print \$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11}' | sort -k1,1 -k2,2n | grep -v tag > $out"."Report.txtSorted";$o=`$c`;
$c=$MyDirectoryPathPrefix."PerlModules/bedtools intersect -a temp/$out$inputNoPath.rawcnv2Sorted -b $out"."Report.txtSorted -sorted -wo | awk '{if((\$4<2&&\$(NF-1)==\"del\")||(\$4>=2&&\$(NF-1)==\"dup\"))print}' > $out"."Report_ContributingCalls.txt";$o=`$c`;
if($countLines>0)
{
	$FailRate=$overallFail/$countLines;
}
if($FailRate>0.8)
{
	print "WARNING: High fail rate of ".sprintf("%.3f", $FailRate)." but some results may still be usable\n";
}

#printf( "%-10s", "hello" );
##printf( "%10d", 777);
##FAIL=CNVRLength/AvgLength<1000, DgvFreq>0.01/DgvEntries>10, MaxP>0.5, PopFreq>0.01, SegDups>10, Recurrent, FreqInflated>0.5, AvgConf<10
##Sort RF <=3, delP<5x10^-4 and ORDup>1 OR dupP<5x10^-4 and ORDup>1, (on exon)
close(REPORT2);
close(SNP_IDS);

if($batch ne "")
{
	$c="perl OutputUtilities/PercentSamples.pl $batch temp/$out"."Snp_IDs.txt";$o=`$c`;
	$c="cat ".$MyDirectoryPathPrefix."header temp/$out"."Report6.txt > temp/$out"."Report_Verbose.txt; mv temp/$out"."Report_Verbose.txt temp/$out"."Report6.txt";$o=`$c`;
	$c="cat temp/$out"."Snp_IDs.txt_CohortCounts.txt | cut -f 2- > temp/$out"."Snp_IDs.txt_CohortCounts_Cols.txt; paste temp/$out"."Report6.txt temp/$out"."Snp_IDs.txt_CohortCounts_Cols.txt > $out"."Report_Verbose.txt";$o=`$c`;
	
}
else
{
	$c="cat ".$MyDirectoryPathPrefix."header temp/$out"."Report6.txt > $out"."Report_Verbose.txt";$o=`$c`;
}

print"Output written to $out"."Report.txt\n";
print LOG "Output written to $out"."Report.txt\n";

#$GIFDel = median ( @DelSnps_x2 ) / 0.456 ; ####### http://en.wikipedia.org/wiki/Population_stratification 1-d.f. x2 distribution independence assumption! so somehow need to prune down to all CNVRs for all significance levels
#$GIFDup = median ( @DupSnps_x2 ) / 0.456 ; ####### Always 0 because median x2=0 (fisher p=1) from a majority of SNPs not having any CNVs need to prune out

open(ALLIDS,"temp/$out"."AllIDs.txt");
open(VCF,">$out"."Report_Verbose.vcf");
print VCF "##fileformat=VCFv4.3\n";
$DATE=`date +%Y%m%d`;
print VCF "##fileDate=$DATE";
print VCF "##reference=$build\n";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
while($line=<ALLIDS>)
{
	chomp($line);
	$AllIDs[$countIDs]=$line;
	print VCF "\t$line";
	$countIDs++;
}
print VCF "\n";
open(FILE,"temp/$out"."Report6.txt");
$line=<FILE>;#Header
while($line=<FILE>)
{
	chomp($line);
	@a=split(/\t/,$line);
	print VCF "$a[0]\t$a[1]\t$a[0]:$a[1]-$a[2]\t.\t<$a[8]>\t.\t$a[-1]\tEND=$a[2]\tGT";
	@cases=split(/,/,$a[11]);
	@controls=split(/,/,$a[12]);
	%VariantIDs=();
	for($i=0;$i<=$#cases;$i++)
	{
		$VariantIDs{$cases[$i]}=1;
	}
	for($i=0;$i<=$#controls;$i++)
        {
                $VariantIDs{$controls[$i]}=1;
        }
	for($i=0;$i<$countIDs;$i++)
	{
		if(exists($VariantIDs{$AllIDs[$i]}))
		{
			print VCF "\t0/1";
		}
		else
		{
			print VCF "\t0/0";
		}
	}
	print VCF "\n";
}

%MonthPastStart = ("0","Jan","1","Feb","2","Mar","3","Apr","4","May","5","Jun","6","Jul","7","Aug","8","Sep","9","Oct","10","Nov","11","Dec");
@timeData = localtime(time);
$Year = $timeData[5] + 1900;
if ($timeData[2] > 12)
{
        $timeData[2] -= 12;
        $AmOrPm = "PM";
}
elsif ($timeData[2] eq 12)
{
        $AmOrPm = "PM";
}
else
{
        $AmOrPm = "AM";
}
if($timeData[2]<10) {$timeData[2]="0".$timeData[2];}
if($timeData[1]<10) {$timeData[1]="0".$timeData[1];}
if($timeData[0]<10) {$timeData[0]="0".$timeData[0];}
print "Run Ended ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";
print LOG "Run Ended ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";

