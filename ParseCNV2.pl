#!/usr/bin/perl
##Joseph Glessner PhD
##Center for Applied Genomics, Children's Hospital of Philadelphia
##University of Pennsylvania Perelman School of Medicine
##Developed under advising of Hakon Hakonarson
##Copyright (C) GNU General Public License 2019 Joseph Glessner parsecnv@gmail.com
##Linux bash environment is best, but Windows Cygwin can also work.
if (@ARGV < 1 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--h" || $ARGV[0] eq "--help")
{
print "\nUSAGE MESSAGE:\n
perl ParseCNV2.pl -i CNV_Calls[.vcf/.rawcnv] [-c Cases.txt/-q Samples_QuantitativeTrait.txt] [arguments]

Required arguments:

-i	input[vcf/rawcnv/txt]
-c	cases
-b	build (if vcf can be inferred)

Optional arguments:

-o	output
-p	merge p variation
-d	merge distance
-t	transmission disequilibrium test
-m	max p inclusion
-bfile	plink bed file of SNP genotypes
-qc	quality control upfront
-log	penncnv log
-q	quantitative trait
-b	batch
-h	help

";
exit;
}
$myCommandLine = "Running Command: perl $0";
for ($i=0; $i <= $#ARGV;$i++)
{
        $myCommandLine .= " $ARGV[$i]";
}
$myCommandLine .= "\n";
use Getopt::Long;
GetOptions('i=s'=>\$input,'o=s'=>\$out,'b=s'=>\$build,'p=f'=>\$mergePVar,'d=i'=>\$mergeDist,'t'=>\$tdt,'m=f'=>\$maxPInclusion,'bfile=s'=>\$bfile,'qc'=>\$qc,'log=s'=>\$log,'c=s'=>\$cases,'q=s'=>\$quantitativeTrait,'batch=s'=>\$batch);
BEGIN {($_=$0)=~s{[^\\\/]+$}{};$_||="./"}  ##In case running elsewhere
$MyDirectoryPathPrefix = $_;
#print "MyDirectoryPathPrefix: $MyDirectoryPathPrefix\n";
if (not $out)
{
	$out="ParseCNV_";
}
if (not $input)
{
	print "ERROR: --i is a required input VCF or PennCNV rawcnv file or list of such files.\n";
	exit;
}
$inputNoPath=$input;
$inputNoPath=~s/.*\///;#get rid of potential path prefix
$dir = 'temp';
unless(-d $dir)
{
    mkdir $dir or die;
}
if($input =~ /.txt/)
{
	$c="rm temp/$out$inputNoPath.rawcnv2";$o=`$c`;
	open(LISTFILE,$input);
	while($line=<LISTFILE>)
	{
		chomp($line);
		
	if($line =~ /.vcf/)
	{       ##TODO get build from vcf
        	$c="perl InputFormatConversion/ConvertVCFtoPennCNV.pl $line | sort -k5,5 | grep -v '^#' | grep -v NOTICE >> temp/$out$inputNoPath.rawcnv2";$o=`$c`;
	}
	elsif($line =~ /.rawcnv/)
	{
        	$c="awk '{print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$8}' $line | sed 's/^chr//' | sed 's/:/\\t/' | sed 's/-/\\t/' | sed 's/\\tstate.,cn=/\\t/' | sed 's/\\tconf=/\\t/' | sort -k5,5 >> temp/$out$inputNoPath.rawcnv2";$o=`$c`;
	}
	else
	{
        	print "ERROR: --i must have a .vcf or .rawcnv extension in the .txt list file\n";
	}

	}
}
elsif($input =~ /.vcf/)
{	##TODO get build from vcf
	$c="perl InputFormatConversion/ConvertVCFtoPennCNV.pl $input | sort -k5,5 | grep -v '^#' | grep -v NOTICE > temp/$out$inputNoPath.rawcnv2";$o=`$c`;
}
elsif($input =~ /.rawcnv/)
{
	$c="awk '{print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$8}' $input | sed 's/^chr//' | sed 's/:/\\t/' | sed 's/-/\\t/' | sed 's/\\tstate.,cn=/\\t/' | sed 's/\\tconf=/\\t/' | sort -k5,5 > temp/$out$inputNoPath.rawcnv2";$o=`$c`;
}
else
{
	print "ERROR: --i must have a .vcf or .rawcnv extension\n";
}

$resfile = $out."ParseCNV.log";
open (LOG, ">$resfile");

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
print "Run Started ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";
print LOG "Run Started ".$timeData[2].":".$timeData[1].":".$timeData[0]." $AmOrPm on ".$MonthPastStart{$timeData[4]}." ".$timeData[3]." ".$Year."\n";
print $myCommandLine;
print LOG $myCommandLine;
print "MyDirectoryPathPrefix: $MyDirectoryPathPrefix\n";
if(!(-e $MyDirectoryPathPrefix."PerlModules/plink"))
{
	$c="unzip ".$MyDirectoryPathPrefix."PerlModules/plink.zip -d ".$MyDirectoryPathPrefix."PerlModules";$o=`$c`;
	$c="unzip ".$MyDirectoryPathPrefix."GeneRef/hg19_gc5Base_SimFormat_AllCol.sorted.zip -d ".$MyDirectoryPathPrefix."GeneRef";$o=`$c`;
	$c="unzip ".$MyDirectoryPathPrefix."GeneRef/hg19_knownGene_Exons_SimFormat_AllCol_UniqueIDs.zip -d ".$MyDirectoryPathPrefix."GeneRef";$o=`$c`;
}
#Check presence and version of dependencies
$cmd="which bash";$o=`$cmd`; print LOG "$o"; if($o=~" no "){print "ERROR: bash not found!\n";}
$cmd="which R";$o=`$cmd`; print LOG "$o"; if($o=~" no "){print "ERROR: R not found!\n";}
$cmd="which ".$MyDirectoryPathPrefix."PerlModules/plink";$o=`$cmd`; print LOG "$o"; if($o=~" no "){print "ERROR: plink not found!\n";}
#$c="awk '/MemAvailable/ {printf \$2/1000}' /proc/meminfo";$o=`$c`;print $o." MB RAM Available\n";print LOG $o." MB RAM Available\n";
$c="free -m | awk '{if(NR==2)printf\$4}'";$o=`$c`;print $o." MB RAM Available\n";print LOG $o." MB RAM Available\n";
if(-e plink.sexcheck)
{
	$c="rm plink.sexcheck";$o=`$c`;
}
if($qc)
{
	$c="perl InputUtilities/ParseCNV_QC.pl";
	if($input =~ /.rawcnv/)
	{
		$c.=" --rawcnv $input";
	}
	if($log)
	{
		$c.=" --log $log";
	}
	if($bfile)
	{
		#$cp="PerlModules/./plink --bfile $bfile --missing";$o=`$cp`;
		$cp="PerlModules/./plink --bfile $bfile --pca";$o=`$cp`;
		#$cp="PerlModules/./plink --bfile $bfile --genome";$o=`$cp`;#--min 0.2 to save disk space
		$cp="PerlModules/./plink --bfile $bfile --check-sex";$o=`$cp`;
		#$c.=" --callrate plink.imiss --popstrat plink.eigenvec --related plink.genome";
	}
	$o=`$c`;
	print "$o\n";
	$c="awk '{print \$1\"\\t\"\$4\"\\t\"\$5}' Cases.rawcnv_remove_QC_RemoveIDs.txt_wMaceQualityScore_remove_QC_RemoveCalls.txt_Indexes | sed 's/^chr//' | sed 's/:/\\t/' | sed 's/-/\\t/' | sed 's/\\tstate.,cn=/\\t/' | sort -k5,5 > temp/$out$inputNoPath.rawcnv2";$o=`$c`;
}
$c="awk '{print \$1\"_\"\$2\"\\n\"\$1\"_\"\$3}' temp/$out$inputNoPath.rawcnv2 | sort -u -k1,1 -k2,2n -t_ > temp/$out"."map";$o=`$c`;
%h_state=();
%chr_posIndex=();
open(MAP,"temp/$out"."map");
while($line=<MAP>)
{
	chomp($line);
	$h_state{$line}="-9";
	$chr_posIndex{$line}=$lineNum;
	#$h_id{$line}="";
	@chr_pos=split("_",$line);
	$diff[$lineNum]=$chr_pos[1]-$last_pos;
	$last_pos=$chr_pos[1];
	$lineNum++;
}
$diff[$lineNum]=1;#Exit at end

$c="awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u > temp/$out"."AllIDs.txt";$o=`$c`;

if(-e "plink.sexcheck")
{
	$c="awk '{if(NR>1)print \$2\"\\t\"\$4}' plink.sexcheck | sort > plink.sexcheck2; awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u > filteredSamples; join plink.sexcheck2 filteredSamples | awk '{print \$2}' > plink.sexcheck3";$o=`$c`;
	$c="awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u | cat $cases - | sort | uniq -c | paste - plink.sexcheck3 | awk '{print \"0 \"\$2\" 0 0 \"\$3\" \"\$1}' > temp/$out"."fam";$o=`$c`;
}
else
{
	if($cases ne "")
	{
		$c="awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u | fgrep -vf - $cases > temp/$out"."casesNotInCNV.txt";$o=`$c`;
		$c="wc -l < temp/$out"."casesNotInCNV.txt";$o=`$c`;chomp($o);
		if($o > 0)
		{print "WARNING: $o specified cases are not present in the CNV file as listed in temp/$out"."casesNotInCNV.txt.\n";}
		
		$c="awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u >temp/$out"."uniq_ids";$o=`$c`;
		$c="sed 's/\r//' $cases | sort -u | fgrep -vwf temp/$out"."casesNotInCNV.txt > temp/$out"."cases_NoCR_NoDup_Exists";$o=`$c`;
		$c="cat temp/$out"."cases_NoCR_NoDup_Exists temp/$out"."uniq_ids | sort | uniq -c | awk '{print \"0 \"\$2\" 0 0 0 \"\$1}' | fgrep -wf temp/$out"."uniq_ids > temp/$out"."fam";$o=`$c`;
		$c="awk '{print \$6}' temp/$out"."fam | sort -nr | uniq -c | awk '{if(NR==1){printf \$1\" cases, \"}else{printf \$1\" controls\"}}'";$o=`$c`;print"$o\n";print LOG "$o\n";
		
	}
	elsif($quantitativeTrait ne "")
	{
		$c="awk '{print \$5}' temp/$out$inputNoPath.rawcnv2 | sort -u >temp/$out"."uniqIDs; sort $quantitativeTrait >temp/$out"."qt.sort; join temp/$out"."uniqIDs temp/$out"."qt.sort | awk '{print \"0 \"\$1\" 0 0 0 \"\$2}' > temp/$out"."fam";$o=`$c`;
		$c="wc -l temp/$out"."uniqIDs | awk '{ORS=\"\";print \$1}'";$samples=`$c`;
		$c="wc -l $quantitativeTrait | awk '{ORS=\"\";print \$1}'";$traits=`$c`;
		$c="wc -l temp/$out"."fam | awk '{ORS=\"\";print \$1}'";$samplesWTrait=`$c`;
		print "$samples samples, $traits traits, and $samplesWTrait samples with trait\n";
	}
	
	if($tdt ne "")
	{
		if($bfile ne "")
		{
			$c="join temp/$out"."uniq_ids $bfile".".fam -1 1 -2 2 -o2.1,2.2,2.3,2.4,2.5,2.6 | sort -k2,2 > temp/$out"."fam";$o=`$c`;
			$c="wc -l < temp/$out"."uniq_ids";$o1=`$c`;chomp($o1);
			$c="wc -l < $bfile".".fam";$o2=`$c`;chomp($o2);
                	if($o2 > $o1)
                	{print "WARNING: ".($o2-$o1)." fam specified cases are not present in the CNV file.\n";}
		}
		else
		{
			print "ERROR: bfile must be specified with tdt option to provide family relationships\n";
			exit;
		}
	}
}
#$, = "|";
open(my $BEDDEL, '>:raw', 'temp/'.$out.'plinkDel.bed') or die "Unable to open: $!";
print $BEDDEL pack('B8','01101100') ; #Plink bed magic number
print $BEDDEL pack('B8','00011011') ; #Plink bed magic number
print $BEDDEL pack('B8','00000000') ; #Plink bed mode individual major
open(my $BEDDUP, '>:raw', 'temp/'.$out.'plinkDup.bed') or die "Unable to open: $!";
print $BEDDUP pack('B8','01101100') ; #Plink bed magic number
print $BEDDUP pack('B8','00011011') ; #Plink bed magic number
print $BEDDUP pack('B8','00000000') ; #Plink bed mode individual major
open(my $BEDDELDUP, '>:raw', 'temp/'.$out.'plinkDelDup.bed') or die "Unable to open: $!";
print $BEDDELDUP pack('B8','01101100') ; #Plink bed magic number
print $BEDDELDUP pack('B8','00011011') ; #Plink bed magic number
print $BEDDELDUP pack('B8','00000000') ; #Plink bed mode individual major
close($BEDDEL);
close($BEDDUP);
close($BEDDELDUP);
open(my $BEDDEL, '>>:raw', 'temp/'.$out.'plinkDel.bed') or die "Unable to open: $!";
open(my $BEDDUP, '>>:raw', 'temp/'.$out.'plinkDup.bed') or die "Unable to open: $!";
open(my $BEDDELDUP, '>>:raw', 'temp/'.$out.'plinkDelDup.bed') or die "Unable to open: $!";

#$| = 1; #Turn off buffering

open(FILE,"temp/$out$inputNoPath.rawcnv2");
$cnvLineNum=0;
$LinesOfCalls=`wc -l < temp/$out$inputNoPath.rawcnv2`;;
print "\rProgress:","0","\%","\r";
while($chrStaStoStaId=<FILE>)
{
	@a=split(/\t/,$chrStaStoStaId);
	$chr=$a[0];
	$start=$a[1];
	$stop=$a[2];
	$state=$a[3];
	$id=$a[4];
	$confidence=$a[5];
	$posIndex=$chr_posIndex{$chr."_".$start};
	#print $id."\n";
	#print "$posIndex\n";
	if($lastID eq $id || $cnvLineNum == 0)
	{
		#print "same id\n";
		#DONE: Only increment to next observed base rather than every one	
		#if($id eq "1749718487_A.int.csv"){print $diff[$posIndex]."_".$chr."_".$i."\n";}
		for($i=$start;$i<=$stop&$diff[$posIndex]>=0|$i==$start;$i+=$diff[$posIndex])
		{
			#if($id eq "1749718487_A.int.csv"){print $chr."_".$i."\n";}
			if(exists($h_state{$chr."_".$i}))
			{
				$h_state{$chr."_".$i}=$state;
				$h_lengthTotal{$chr."_".$i}+=$stop-$start+1;
				$h_confidenceTotal{$chr."_".$i}+=$confidence;
				#$h_countSamples{$chr."_".$i}++;
        			#$h_id{$chr."_".$i}.=$id;
        			$posIndex++;
			}
			else
			{
				print "h_state \$chr._.\$i Does Not Exist!\n";
			}
		}
	}
	else
	{
		$byte="";
		while (($key, $value) = each(%h_state))
		{
     			#print $key.",".$value."\n";
			#was [00,11,12,22] but for binary need [10,00,01,11]
			if($value=="-9" || $value>1){$byte.="11";}
			elsif($value=="0"){$byte.="00";}
			elsif($value=="1"){$byte.="01";}
			else{$byte.="10";}
			if(length($byte)==8)
                     	{
				#print "$byte\n";
                               	$packedForBedDel .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
				$byte="";
			}
		}
		print $BEDDEL $packedForBedDel;
		$packedForBedDel = "";
		if(length($byte)>0)
		{
			#print "$byte\n";
			print $BEDDEL pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
			$byte="";
		}
		
		while (($key, $value) = each(%h_state))
                {
                        #print $key.",".$value."\n";
                        if($value=="-9" || $value<2){$byte.="11";}
                        elsif($value=="4"){$byte.="00";}
                        elsif($value=="3"){$byte.="01";}
                        else{$byte.="10";}
                        if(length($byte)==8)
                        {
                                $packedForBedDup .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                                $byte="";
                        }
                }
		print $BEDDUP $packedForBedDup;
		$packedForBedDup= "";
                if(length($byte)>0)
                {
                        print $BEDDUP pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                        $byte="";
                }

		while (($key, $value) = each(%h_state))
                {
                       	#print $key.",".$value."\n";
                        if($value=="-9"){$byte.="11";}
                        elsif($value=="0" || $value=="4" || $value>4){$byte.="00";}
                        elsif($value=="1" || $value=="3" || $value=="2"){$byte.="01";}
                        else{$byte.="10";}
                        if(length($byte)==8)
                        {
                                $packedForBedDelDup .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                                $byte="";
                        }
			$h_state{$key}="-9";#Reinitialize for next iteration with new sample id
                	#$h_length{$key}=0;#Reinitialize for next iteration with new sample id
		}
		print $BEDDELDUP $packedForBedDelDup;
		$packedForBedDelDup = "";
                if(length($byte)>0)
                {
                        print $BEDDELDUP pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                        $byte="";
                }

		#if($id eq "1749718487_A.int.csv"){print $diff[$posIndex]."_".$chr."_".$i."\n";}
		for($i=$start;$i<=$stop&$diff[$posIndex]>=0|$i==$start;$i+=$diff[$posIndex])
                {
			#print $chr."_".$i."\n";
			#if($id eq "1749718487_A.int.csv"){print $chr."_".$i."\n";}
                	if(exists($h_state{$chr."_".$i}))
                        {
                                $h_state{$chr."_".$i}=$state;
                        	$h_lengthTotal{$chr."_".$i}+=$stop-$start+1;
                                $h_confidenceTotal{$chr."_".$i}+=$confidence;
			}
			$posIndex++;
                }
	}
	$lastID=$id;
	$cnvLineNum++;
	$Progress = ($cnvLineNum/$LinesOfCalls) *100;
	$Progress =~ s/\..*//;
	print("\rProgress:$Progress\%");
}
#print "EXIT\n";
#Exit last sample condition
open(BIM,">temp/$out$inputNoPath.rawcnv2.bim");
while (($key, $value) = each(%h_state))
{
	@Chr_Pos=split(/_/,$key);
	#TODO: Take into account non-human builds
	if($Chr_Pos[0]=="X"){$Chr_Pos[0]=23;}
	if($Chr_Pos[0]=="Y"){$Chr_Pos[0]=24;}
	if($Chr_Pos[0]=="XY"){$Chr_Pos[0]=25;}
	if($Chr_Pos[0]=="M" || $Chr_Pos[0]=="MT"){$Chr_Pos[0]=26;}
	print BIM "$Chr_Pos[0]\t$key\t0\t$Chr_Pos[1]\t1\t2\n";
}
                while (($key, $value) = each(%h_state))
                {
                        #print $key.",".$value."\n";
                        #was [00,11,12,22] but for binary need [10,00,01,11]
                        if($value=="-9" || $value>1){$byte.="11";}
                        elsif($value=="0"){$byte.="00";}
                        elsif($value=="1"){$byte.="01";}
                        else{$byte.="10";}
                        if(length($byte)==8)
                        {
                                #print "$byte\n";
				$packedForBedDel .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                                $byte="";
                        }
                }
		print $BEDDEL $packedForBedDel;
		$packedForBedDel = "";
                if(length($byte)>0)
                {
			#print "$byte\n";
                        print $BEDDEL pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                        $byte="";
                }

                while (($key, $value) = each(%h_state))
                {
                        #print $key.",".$value."\n";
                        if($value=="-9" || $value<2){$byte.="11";}
                        elsif($value=="4"){$byte.="00";}
                        elsif($value=="3"){$byte.="01";}
                        else{$byte.="10";}
                        if(length($byte)==8)
                        {
                                $packedForBedDup .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                                $byte="";
                        }
                }
		print $BEDDUP $packedForBedDup;
		$packedForBedDup = "";
                if(length($byte)>0)
                {
                        print $BEDDUP pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                        $byte="";
                }

                while (($key, $value) = each(%h_state))
                {
                        #print $key.",".$value."\n";
                        if($value=="-9"){$byte.="11";}
                        elsif($value=="0" || $value=="4" || $value>4){$byte.="00";}
                        elsif($value=="1" || $value=="3" || $value=="2"){$byte.="01";}
                        else{$byte.="10";}
                        if(length($byte)==8)
                        {
                                $packedForBedDelDup .= pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                                $byte="";
                        }
                        $h_state{$key}="-9";#Reinitialize for next iteration with new sample id
                }
		print $BEDDELDUP $packedForBedDelDup;
		$packedForBedDelDup = "";
                if(length($byte)>0)
                {
                        print $BEDDELDUP pack('b8',$byte) ;  ### small b to put into file backwards as in plink bed spec
                        $byte="";
                }
		close($BEDDEL);
		close($BEDDUP);
		close($BEDDELDUP);
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bed temp/$out"."plinkDel.bed --bim temp/$out$inputNoPath.rawcnv2.bim --fam temp/$out"."fam --make-bed --allow-no-sex --out temp/$out"."del";$o=`$c`;#print"$c\n";print"$o\n";
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bed temp/$out"."plinkDup.bed --bim temp/$out$inputNoPath.rawcnv2.bim --fam temp/$out"."fam --make-bed --allow-no-sex --out temp/$out"."dup";$o=`$c`;#print"$c\n";print"$o\n";
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bed temp/$out"."plinkDelDup.bed --bim temp/$out$inputNoPath.rawcnv2.bim --fam temp/$out"."fam --make-bed --allow-no-sex --out temp/$out"."deldup";$o=`$c`;
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."del --assoc fisher --allow-no-sex --out temp/$out"."del";$o=`$c`;
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."dup --assoc fisher --allow-no-sex --out temp/$out"."dup";$o=`$c`;
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."deldup --assoc fisher --allow-no-sex --out temp/$out"."deldup";$o=`$c`;

		#Get matching p-values to original ParseCNV with OR from assoc fisher
		$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."del --model fisher --allow-no-sex --out temp/$out"."del";$o=`$c`;
                $c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."dup --model fisher --allow-no-sex --out temp/$out"."dup";$o=`$c`;
                $c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."deldup --model fisher --allow-no-sex --out temp/$out"."deldup";$o=`$c`;
		
		$c="awk '{print \$2\"\\t\"\$0}' temp/$out"."del.assoc.fisher | sort -k1,1 > temp/$out"."del.assoc.fisher.forJoin";$o=`$c`;
		$c="awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"."del.model | sort -k1,1 > temp/$out"."del.model.forJoin";$o=`$c`;
		$c="echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"."del.assoc.fisher.model";$o=`$c`;
		$c="join temp/$out"."del.assoc.fisher.forJoin temp/$out"."del.model.forJoin | sort -k2,2 -k4,4n >> temp/$out"."del.assoc.fisher.model";$o=`$c`;
		$c="awk '{print \$2\"\\t\"\$0}' temp/$out"."dup.assoc.fisher | sort -k1,1 > temp/$out"."dup.assoc.fisher.forJoin";$o=`$c`;
                $c="awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"."dup.model | sort -k1,1 > temp/$out"."dup.model.forJoin";$o=`$c`;
                $c="echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"."dup.assoc.fisher.model";$o=`$c`;
                $c="join temp/$out"."dup.assoc.fisher.forJoin temp/$out"."dup.model.forJoin | sort -k2,2 -k4,4n >> temp/$out"."dup.assoc.fisher.model";$o=`$c`;
		$c="awk '{print \$2\"\\t\"\$0}' temp/$out"."deldup.assoc.fisher | sort -k1,1 > temp/$out"."deldup.assoc.fisher.forJoin";$o=`$c`;
                $c="awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"."deldup.model | sort -k1,1 > temp/$out"."deldup.model.forJoin";$o=`$c`;
                $c="echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"."deldup.assoc.fisher.model";$o=`$c`;
                $c="join temp/$out"."deldup.assoc.fisher.forJoin temp/$out"."deldup.model.forJoin | sort -k2,2 -k4,4n >> temp/$out"."deldup.assoc.fisher.model";$o=`$c`;
		if($tdt ne "")
        	{
			$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."del --tdt --allow-no-sex --out temp/$out"."del";$o=`$c`;
                	$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."dup --tdt --allow-no-sex --out temp/$out"."dup";$o=`$c`;
                	$c=$MyDirectoryPathPrefix."PerlModules/./plink --bfile temp/$out"."deldup --tdt --allow-no-sex --out temp/$out"."deldup";$o=`$c`;
		}	
		#$c="perl OutputUtilities/InsertPlinkPvalues.pl";$o=`$c`;
print("\rProgress:100\%\n");
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

