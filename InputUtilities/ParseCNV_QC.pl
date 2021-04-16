#!/usr/bin/perl
#Joseph Glessner
=head
Inputs:
PennCNV .log -> LRR_SD, GCWF
PennCNV .rawcnv -> CountCNV
Plink --missing .imiss / Illumina  LIMs Project Detail Report .csv / Illumina Genome Studio Samples Table / General format SampleID	CallRate -> CallRate
Eigenstrat smartPCA .pca.evec / Plink .mds -> PCA
Plink --genome .genome -> PI_HAT relatedness
make sure to include extremevalues_2.2.tar.gz http://cran.r-project.org/web/packages/extremevalues/index.html
For Illumina 550k data and related Illumina chip platforms, the key data quality metric thresholds we have observed are: call rate > 98%, SD LRR < 0.3, |GCWF| < 0.05, and count CNV< 100.
For Affymetrix 6.0 data, these measures include: call rate > 96%, SD LRR < 0.35, |GCWF| < 0.02, and count CNV < 80.
In addition, observations of quality metric modes from individual labs and sample sources are advisable to determine appropriate QC thresholds.
The distribution of these metric measures are constantly reviewed to include only those who fall within a linear mode of the quality metric outside exponential modes for any given genotyping platform. 
Run once with no thresholds specified and review automatically determined thresholds based on outliers in your dataset.
Then you may specify some adjusted threshold for a given quality metric based on your review of the plots.
=cut
BEGIN {($_=$0)=~s{[^\\\/]+$}{};$_||="./"}  ##In case running elsewhere
$MyDirPre = $_;
if(@ARGV<1)
{
print "USAGE: 
perl ParseCNV_QC.pl [arguments]

PLEASE RUN \"export R_LIBS=./\" first.

Arguments:
===Sample QC===
--log		PennCNV .log
--rawcnv	PennCNV .rawcnv
--callrate	Plink --missing .imiss OR Illumina  LIMs Project Detail Report .csv OR Illumina Genome Studio Samples Table OR General format SampleID	CallRate
--popstrat	Eigenstrat smartPCA .pca.evec OR Plink .mds
--related	Plink --genome .genome
--out		An output file prefix to label this run

Run once with no thresholds specified and review automatically determined thresholds based on outliers in your dataset.
Then you may specify some adjusted threshold for a given quality metric based on your review of the plots.
--qccallrate <float>        Call Rate threshold for inclusion of the sample
--qclrrsd <float>           LRR_SD threshold for inclusion of the sample
--qcgcwf <float>            GCWF threshold for inclusion of the sample
--qcnumcnv <int>            number of CNV calls for inclusion of the sample
--qcrelated <float>         sample relatedness for inclusion of the sample

===Call QC===
--numsnp <int>              minimum number of SNPs in CNV calls
--length <int>              minimum length of CNV calls
--confidence <float>        minimum confidence score of CNV calls
--quality <float>           minimum Mace quality score of CNV calls
--stratifyCN		    divide plots by CN state

--popstrat is just for plotting, no exclusion is currently supported based on popstrat.

A tab delimited QC Summary is provided: <myLog>_QC_Vals_Header.log


";
=head
=== Options under development ===
Other available PennCNV log quality metrics I do not use typically: LRR_mean (low/high), BAF_mean(low/high), BAF_SD(high), BAF_DRIFT(high), WF(low/high)
--qcbafdrift <float>        BAF_DRIFT threshold for inclusion of the sample
--qclrrmean <float>         LRR_mean threshold for inclusion of the sample
--qcbafmean <float>         BAF_mean threshold for inclusion of the sample
--qcbafsd <float>           BAF_SD threshold for inclusion of the sample
--qcbafdrift <float>        BAF_DRIFT threshold for inclusion of the sample
--qcwf <float>              WF threshold for inclusion of the sample
=cut
exit;
}
use Getopt::Long;
GetOptions('log=s'=>\$log,'rawcnv=s'=>\$rawcnv,'callrate=s'=>\$callrate,'popstrat=s'=>\$popstrat,'related=s'=>\$related,'out=s'=>\$out,'qccallrate=s'=>\$qccallrate,'qclrrsd=s'=>\$qclrrsd,'qcgcwf=s'=>\$qcgcwf,'qcnumcnv=s'=>\$qcnumcnv,'qcrelated=s'=>\$qcrelated,'numsnp=s'=>\$numsnp,'length=s'=>\$length,'confidence=s'=>\$confidence,'quality=s'=>\$quality,'stratifyCN'=>\$stratifyCN,'noSampleQC'=>\$noSampleQC,'noCallQC'=>\$noCallQC);

##perl ParseCNV_QC.pl --log PennCNV_Omni1.log --rawcnv file.rawcnv --callrate PCGC_Omni25-8v1_BED_Miss.imiss --popstrat file.pca.evec --related PCGC_Omni25-8v1_BED_GenomeAll.genome

if(!($log)){$log=1;print "WARNING:--log is strongly recommended to include since it contains crucial quality metric LRR_SD\n";} ## Allow missing arguments
if(!($rawcnv)){$rawcnv=1;}
if(!($callrate)){$callrate=1;print "WARNING:--callrate is strongly recommended to include since it contains crucial quality metric CallRate. BAF_SD may be best alternative\n";}
if(!($popstrat)){$popstrat=1;}
if(!($related)){$related=1;}

if(!($qccallrate) && !($qclrrsd) && !($qcgcwf) && !($qcnumcnv) && !($qcrelated)){print "WARNING: Review automatically determined thresholds based on outliers in your dataset. Then you may specify some adjusted threshold for a given quality metric based on your review of the plots using the qc* options.\n";}
if(!($qccallrate)){$qccallrate=-999;}
if(!($qclrrsd)){$qclrrsd=-999;}
if(!($qcgcwf)){$qcgcwf=-999;}
if(!($qcnumcnv)){$qcnumcnv=-999;}
if(!($qcrelated)){$qcrelated=.2;}

if(!($numsnp)){$numsnp=-999;}
if(!($length)){$length=-999;}
if(!($confidence)){$confidence=-999;}
if(!($quality)){$quality=-999;}

#Redirect screen output to files
open(NOT_STDOUT,">ParseCNV_QC_temp.log");
select(NOT_STDOUT);
open OLDERR, ">&",\*STDERR; 
open(STDERR, ">ParseCNV_QC_temp.err")  or print("Can't redirect stderr\n");

$c="grep summary $log | grep -v Xmean > $log"."_QC.log";print"$c\n";`$c`; #PennCNV_Omni1.log 

$c="head -1 $log"."_QC.log | sed 's/=[^ ]* / /g' | sed 's/=.*//' | sed 's/for.*:/for ChipID/' > Header_"."$log.txt";print"$c\n";`$c`;

$c="sed 's/ [^ ]*=/ /g' $log"."_QC.log | sed 's/\\://g' > $log"."_QC_Vals.log";print"$c\n";`$c`;

$c="cat Header_$log.txt $log"."_QC_Vals.log > $log"."_QC_Vals_Header.log";print"$c\n";`$c`;




$c="echo CountCNV ChipID> $rawcnv"."_CountCNVChipID.txt";print"$c\n";`$c`;
$c="awk '{print \$5}' $rawcnv | uniq -c | sed -e 's/^[ \\t]*//' >> $rawcnv"."_CountCNVChipID.txt";print"$c\n";`$c`;  # file.rawcnv   ../GeneBased/PQN_2+_8Cols.rawcnv_CLEANQC.rawcnv_CLEANQC.rawcnv_keep_PCGC_CaseIDs.txt



$c="perl ".$MyDirPre."/../PerlModules/GetCallRate.pl $callrate";print"$c\n";$output=`$c`;#### PCGC_Omni25-8v1_BED_Miss.imiss  
print"$output\n";

###perl GetCallRate.pl PCGC_PDR_7-19-12.csv

###perl GetCallRate.pl PCGC_Omni1_SamplesTable.txt


###file.pca.evec   file.mds
$c="head -1 $popstrat";print"$c\n";$output=`$c`;
#           #eigvals:   199.492    39.213     9.287     7.900     3.519     3.015     2.579     2.347     2.335     2.324
#		    FID                 IID    SOL           C1           C2           C3           C4
if($output=~/eigvals/)
{
	print "MATCH PCA\n";
	$c="cp $popstrat $popstrat"."_Standard.txt";print"$c\n";$output=`$c`;
	$PCAOrMDS="PCA";
}
elsif($output=~/FID/)
{
	print "MATCH MDS\n";
	$c="awk '{print \$2\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7}' $popstrat | sed '1d' > $popstrat"."_Standard.txt";print"$c\n";$output=`$c`;
	$PCAOrMDS="MDS";
}
else
{
	print "MATCH PCA\n";
        $c="awk '{print \$2\"\\t\"\$3\"\\t\"\$4}' $popstrat > $popstrat"."_Standard.txt";print"$c\n";$output=`$c`;
        $PCAOrMDS="PCA";
	print "No header on population stratification file so assuming PCA $popstrat\n";
}



###PCGC_Omni25-8v1_BED_GenomeAll.genome

###$qcrelated was .2
$c="awk '{if(\$10>$qcrelated){print \$2\"\\t\"\$4\"\\t\"\$10}}' $related > $related"."_3Col.genome";print"$c\n";$output=`$c`;  ###Added condition to lessen computational burden of displaying so many points

##Count IDs to determine minimal exclusions
$c="awk '{print \$2\"\\t\"\$4\"\\t\"\$10}' $related | awk '{if(\$3>$qcrelated){print \$1; print \$2}}' | sort | uniq -c | sed -e 's/^[ \\t]*//' | sort -n -r> $related"."_Gt2PihatCountIDs.txt";print"$c\n";$output=`$c`;

unless( open( RELATEDFILE, $related."_3Col.genome"))
{
	print "File not found!\n";
	exit;
}
unless( open( COUNTIDFILE, $related."_Gt2PihatCountIDs.txt"))  #### TO DO Start by removing IDs failing other QC metrics first
{
	print "File not found!\n";
	exit;
}
open(RELATEDQCREMOVEIDFILE,">$out"."QC_RemoveIDs.txt");
print RELATEDQCREMOVEIDFILE "ChipIDRemove Value Metric\n";
$PairsRemainingGt2Pihat=1;###Just initialize
#print"PairsResolvedByThisID:PairsRemainingGt2Pihat ";
for($i=0;$PairsRemainingGt2Pihat>0;$i++)
{
if($i>0)
{
	$c="mv RemainingGt2Pihat.txt Gt2Pihat.txt";`$c`;
	unless( open( RELATEDFILE, "Gt2Pihat.txt"))
	{
	print "File not found!\n";
	exit;
	}
}
$HighestCountIDLine=<COUNTIDFILE>; ### Next highest count ID taken each iteration of for loop
@Vals=split(/\s+/,$HighestCountIDLine);
$IDCount=$Vals[0];
$HighestCountID=$Vals[1];
open(RESFILE,">RemainingGt2Pihat.txt");
$PairsRemainingGt2Pihat=0;
$PairsResolvedByThisID=0;
while($line = <RELATEDFILE>)
{
	chomp($line);
	$line =~ s/[\r\n]+$//g;
	@Vals=split(/\s+/,$line);
	if($Vals[0] eq $HighestCountID || $Vals[1] eq $HighestCountID)
	{
		$PairsResolvedByThisID++;
	}
	else
	{
		print RESFILE "$line\n";
		$PairsRemainingGt2Pihat++;
	}
}
if($PairsResolvedByThisID>0 && $HighestCountID ne "IID1" && $HighestCountID ne "IID2")
{print RELATEDQCREMOVEIDFILE "$HighestCountID $IDCount"."PiHat>$qcrelated PiHat>$qcrelated"."HighestCountIDIterativelyRemoved\n";}#print"$PairsResolvedByThisID:$PairsRemainingGt2Pihat "}

close RELATEDFILE;
close RESFILE;
}
close RELATEDQCREMOVEIDFILE;
close COUNTIDFILE;
print"\n";

#$c="awk '{print \$1;print\$2}' $related"."_3Col.genome | sort |uniq --count >$related"."_AllIDCounts";print"$c\n";$output=`$c`;
#$c="R CMD INSTALL gWidgets_0.0-54.tar.gz -l ./";@output=`$c`;
#print "@output\n";
#$c="R CMD INSTALL digest_0.6.8.tar.gz -l ./";@output=`$c`;
#print "@output\n";
#$c="R CMD INSTALL gWidgetstcltk_0.0-55.tar.gz -l ./";@output=`$c`;
#print "@output\n";
$c="R CMD INSTALL ".$MyDirPre."../extremevalues_2.2.tar.gz -l ./";@output=`$c`;
print "@output\n";


#print "PLEASE RUN \"export R_LIBS=./\" if you get \"Error: cannot execute system command R CMD BATCH .AllRes.R\"\n";
#print "If still error, check .AllRes.Rout for specific error and submit to parsecnv@gmail.com if not resolvable.\n";
$c="export R_LIBS=./";`$c`; ##Somehow does not work

open (R, ">$out.AllRes.R");
print R <<END;
pdf(file=paste("$out","QC_Plot.pdf" ,sep=""), width=3,height=3, pointsize=6)

if("$log"!=1){a=read.table(paste("$log","_QC_Vals_Header.log",sep=""),header=TRUE, comment.char="")}
#one=read.table(paste("BafLrr/$out","AllRes.txt",sep=""),sep="\t")

if("$rawcnv"!=1){b=read.table(paste("$rawcnv","_CountCNVChipID.txt",sep=""),header=TRUE, comment.char="")}

if("$callrate"!=1){d=read.table(paste("$callrate","_CallRates.txt",sep=""),header=TRUE)}
if("$popstrat"!=1){e=read.table(paste("$popstrat","_Standard.txt",sep=""))}

if("$related"!=1){g=read.table(paste("$related","_3Col.genome",sep=""),header=TRUE)}

#install.packages("extremevalues_2.2.tgz",repos=FALSE)
#install.packages("extremevalues", dep=TRUE, lib=NULL)
par(mfrow=c(3,2),las=1,mar=c(4,4.2,1.5,1.3))
library()
library(extremevalues,lib.loc="./")

if("$callrate"!=1){

out <- getOutliers(sort(d\$CallRate),distribution="weibull")
x <- array(dim=c(1,length(d\$CallRate)))
x[out\$iLeft]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$CallRate)) }


if("$qccallrate"!=-999) { 
out\$limit[1] <- as.numeric("$qccallrate")
out\$iLeft = c(1:(sum(d\$CallRate < as.numeric("$qccallrate"))))
x <- array(dim=c(1,length(d\$CallRate)))
x[out\$iLeft]="blue"
out\$nOut[1]=sum(d\$CallRate < as.numeric("$qccallrate"))
}


x[max(out\$iLeft)+1:length(d\$CallRate)]="black"
#plot(sort(d$CallRate),ylim=c(.9,1),col=x)
plot(sort(d\$CallRate),ylim=c(.85,1),col=x,main=paste("CallRate >",round(out\$limit[1],4),"Fail=", out\$nOut[1]),xlab="Sample Index")#######out\$iLeft))  
report<-head(d[sort.list(d[,2]), ],max(out\$iLeft))
report[1:length(report\$CallRate),3]=paste("CallRate>",round(out\$limit[1],4),sep="")
###write.table("ChipIDRemove Value Metric", file = "QC_RemoveIDs.txt", append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE) ### append false for the first one to start new file if rerun (Now PI_HAT related above > clears)
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)



}


if("$log"!=1){
out <- getOutliers(sort(a\$LRR_SD),distribution="lognormal")
x <- array(dim=c(1,length(a\$LRR_SD)))
x[out\$iRight]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$LRR_SD)) }


if("$qclrrsd"!=-999) { 
out\$limit[2] <- as.numeric("$qclrrsd")
out\$iRight = c(((length(a\$LRR_SD)+1)-(sum(a\$LRR_SD > as.numeric("$qclrrsd")))):(length(a\$LRR_SD)))
x <- array(dim=c(1,length(a\$LRR_SD)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(a\$LRR_SD > as.numeric("$qclrrsd"))
}


x[1: min(out\$iRight)-1]="black"
#plot(sort(a\$LRR_SD), col=x)
plot(sort(a\$LRR_SD), col=x, main=paste("LRR_SD <",round(out\$limit[2],4),"Fail=", out\$nOut[2]),xlab="Sample Index")  ###sort(a\$LRR_SD)[min(out\$iRight)] instead of round(out\$limit[2],4)
q_Lrrsd<-tail(a[order(a\$LRR_SD), ],max(out\$iRight)-min(out\$iRight)+1)
report<-q_Lrrsd[,c("ChipID","LRR_SD")]
report[1:length(report\$ChipID),3]= paste("LRR_SD<", round(out\$limit[2],4) ,sep="")
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)


if(sum(!is.na(a\$GCWF))>2)
{
out <- getOutliers(sort(a\$GCWF),method="II")
x <- array(dim=c(1,length(a\$GCWF)))
x[out\$iRight]="blue"
x[out\$iLeft]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$GCWF)) }


if("$qcgcwf"!=-999) { 
out\$limit[2] <- as.numeric("$qcgcwf")
out\$iRight = c(((length(a\$GCWF)+1)-(sum(a\$GCWF > as.numeric("$qcgcwf")))):(length(a\$GCWF)))
x <- array(dim=c(1,length(a\$GCWF)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(a\$GCWF > as.numeric("$qcgcwf"))

out\$limit[1] <- as.numeric("-$qcgcwf")
out\$iLeft = c(1:(sum(a\$GCWF < as.numeric("-$qcgcwf"))))
###x <- array(dim=c(1,length(a\$GCWF)))
x[out\$iLeft]="blue"
out\$nOut[1]=sum(a\$GCWF < as.numeric("-$qcgcwf"))
}

x[(max(out\$iLeft)+1):(min(out\$iRight)-1)]="black"

#plot(sort(a\$GCWF), col=x,ylim=c(-.01,.01))
plot(sort(a\$GCWF), col=x, main=paste(round(out\$limit[1],4),"< GCWF <",round(out\$limit[2],4),"Fail=", out\$nOut[1]+ out\$nOut[2]),xlab="Sample Index") # out\$limit was not correct #sort(a\$GCWF)[max(out\$iLeft)] AND sort(a\$GCWF)[min(out\$iRight)] prints closest to specified threshold
q_Gcwf_Low <-head(a[order(a\$GCWF), ],max(out\$iLeft))
q_Gcwf_High<-tail(a[order(a\$GCWF), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_Gcwf_Low[,c("ChipID","GCWF")]
report[1:length(report\$ChipID),3]= paste(round(out\$limit[1],4),"<GCWF<",round(out\$limit[2],4),sep="")
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
report<-q_Gcwf_High[,c("ChipID","GCWF")]
report[1:length(report\$ChipID),3]= paste(round(out\$limit[1],4),"<GCWF<",round(out\$limit[2],4),sep="")
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}

}

if("$rawcnv"!=1){
out <- getOutliers(sort(b\$CountCNV),method="II")
x <- array(dim=c(1,length(b\$CountCNV)))
x[out\$iRight]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(b\$CountCNV)) }

if("$qcnumcnv"!=-999) { 
out\$limit[2] <- as.numeric("$qcnumcnv")
out\$iRight = c(((length(b\$CountCNV)+1)-(sum(b\$CountCNV > as.numeric("$qcnumcnv")))):(length(b\$CountCNV)))
x <- array(dim=c(1,length(b\$CountCNV)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(b\$CountCNV > as.numeric("$qcnumcnv"))

x[1: min(out\$iRight)-1]="black"
#plot(sort(b\$CountCNV), col=x)
plot(sort(b\$CountCNV), col=x, main=paste("CountCNV <",round(out\$limit[2],4) ,"Fail=", out\$nOut[2]),xlab="Sample Index")  ### was sort(b\$CountCNV)[min(out\$iRight)]
q_CountCNV <-tail(b[order(b\$CountCNV), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_CountCNV [,c("ChipID","CountCNV")]
report[1:length(report\$ChipID),3]= paste("CountCNV<", round(out\$limit[2],4) ,sep="") ### was sort(b\$CountCNV)[min(out\$iRight)]
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
else{
x[1: min(out\$iRight)-1]="black"
#plot(sort(b\$CountCNV), col=x)
plot(sort(b\$CountCNV), col=x, main=paste("CountCNV <",sort(b\$CountCNV)[min(out\$iRight)] ,"Fail=", out\$nOut[2]),xlab="Sample Index")
q_CountCNV <-tail(b[order(b\$CountCNV), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_CountCNV [,c("ChipID","CountCNV")]
report[1:length(report\$ChipID),3]= paste("CountCNV<", sort(b\$CountCNV)[min(out\$iRight)] ,sep="") 
write.table(report, file = paste("$out","QC_RemoveIDs.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
}

if("$popstrat"!=1){
mydata <- array(dim=c(length(e\$V2),2))
mydata[1:length(e\$V2),1]=e\$V2
mydata[1:length(e\$V2),2]=e\$V3
fit <- kmeans(mydata, 3)
library(cluster)
clusplot(mydata, fit\$cluster, color=TRUE,main="$PCAOrMDS Clusters")##, shade=TRUE, labels=2, lines=0 http://www.statmethods.net/advstats/cluster.html
}

if("$related"!=1){
out <- getOutliers(sort(g\$PI_HAT),method="II")
x <- array(dim=c(1,length(g\$PI_HAT)))
out\$iRight=1:length(g\$PI_HAT)  ###Just <0.2 keep condition $qcrelated
x[out\$iRight]="blue"
###x[1: min(out\$iRight)-1]="black" ###Just <0.2 keep condition $qcrelated
#plot(sort(g\$PI_HAT), col=x)
out\$nOut[2]=length(g\$PI_HAT)###Just <0.2 keep condition $qcrelated
plot(sort(g\$PI_HAT), col=x, main=paste("PI_HAT <", sort(g\$PI_HAT)[min(out\$iRight)],"Fail=", out\$nOut[2]),xlab="Sample Combination Index")
q_PI_HAT <-tail(g[order(g\$PI_HAT), ], max(out\$iRight)-min(out\$iRight)+1)
report<-tail(q_PI_HAT [,c("IID1","PI_HAT")], max(out\$iRight)-min(out\$iRight)+1)  ##"IID2", Just exclude first ID , better would be exclude highest count occurrence IDs, then remove resolved pairs with that ID, then exclude first ID when count ID is all one 
report[1:length(report\$IID1),3]=paste("PI_HAT<", sort(g\$PI_HAT)[min(out\$iRight)],sep="")
######write.table(report, file = "QC_RemoveIDs.txt", append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
######remove high count IDs iteratively above to minimize number of samples removed
}

AllQC<-read.table(paste("$out","QC_RemoveIDs.txt",sep="") ,header=TRUE,comment.char="")
#par(xpd=NA)
##text(0,-.08,paste(length(table(AllQC\$ChipIDRemove)),"unique fail any QC"),adj=0,cex=2)##Difficult because plot axis ranges unknown and applies to all plots, not just PI_HAT
mtext(paste(length(table(AllQC\$ChipIDRemove)),"unique samples fail any QC"), 1, padj=2.4)###put text below
#par(xpd=FALSE)






dev.off()
pdf(file=paste("$out","QC_Plot_2.pdf",sep=""), width=3,height=3, pointsize=6)

#other available PennCNV log quality metrics LRR_mean (low/high), BAF_mean(low/high), BAF_SD(high), BAF_DRIFT(high), WF(low/high)
par(mfrow=c(3,2),las=1,mar=c(4,4.2,1.5,1.3))
if("$log"!=1){
out <- getOutliers(sort(a\$LRR_mean),method="II")
x <- array(dim=c(1,length(a\$LRR_mean)))
x[out\$iRight]="blue"
x[out\$iLeft]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$LRR_mean)) }

x[(max(out\$iLeft)+1):(min(out\$iRight)-1)]="black"
#plot(sort(a\$LRR_mean), col=x,ylim=c(-.01,.01))
plot(sort(a\$LRR_mean), col=x, main=paste(sort(a\$LRR_mean)[max(out\$iLeft)+1],"<= LRR_mean <=",sort(a\$LRR_mean)[min(out\$iRight)-1],"Fail=", out\$nOut[1]+ out\$nOut[2]),xlab="Sample Index") # out\$limit was not correct
q_LRR_mean_Low <-head(a[order(a\$LRR_mean), ],max(out\$iLeft))
q_LRR_mean_High<-tail(a[order(a\$LRR_mean), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_LRR_mean_Low[,c("ChipID","LRR_mean")]
report[1:length(report\$ChipID),3]= paste(sort(a\$LRR_mean)[max(out\$iLeft)+1],"<=LRR_mean<=",sort(a\$LRR_mean)[min(out\$iRight)-1],sep="")
write.table("ChipIDRemove Value Metric", file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep="")  , append = FALSE,quote=FALSE,row.names=FALSE, col.names =FALSE) ### append false for the first one to start new file if rerun
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
report<-q_LRR_mean_High[,c("ChipID","LRR_mean")]
report[1:length(report\$ChipID),3]= paste(sort(a\$LRR_mean)[max(out\$iLeft)+1],"<=LRR_mean<=",sort(a\$LRR_mean)[min(out\$iRight)-1],sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)

#BAF_mean(low/high)
out <- getOutliers(sort(a\$BAF_mean),method="II")
x <- array(dim=c(1,length(a\$BAF_mean)))
x[out\$iRight]="blue"
x[out\$iLeft]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$BAF_mean)) }

x[(max(out\$iLeft)+1):(min(out\$iRight)-1)]="black"
#plot(sort(a\$BAF_mean), col=x,ylim=c(-.01,.01))
plot(sort(a\$BAF_mean), col=x, main=paste(sort(a\$BAF_mean)[max(out\$iLeft)+1],"<= BAF_mean <=",sort(a\$BAF_mean)[min(out\$iRight)-1],"Fail=", out\$nOut[1]+ out\$nOut[2]),xlab="Sample Index") # out\$limit was not correct
q_BAF_mean_Low <-head(a[order(a\$BAF_mean), ],max(out\$iLeft))
q_BAF_mean_High<-tail(a[order(a\$BAF_mean), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_BAF_mean_Low[,c("ChipID","BAF_mean")]
report[1:length(report\$ChipID),3]= paste(sort(a\$BAF_mean)[max(out\$iLeft)+1],"<=BAF_mean<=",sort(a\$BAF_mean)[min(out\$iRight)-1],sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
report<-q_BAF_mean_High[,c("ChipID","BAF_mean")]
report[1:length(report\$ChipID),3]= paste(sort(a\$BAF_mean)[max(out\$iLeft)+1],"<=BAF_mean<=",sort(a\$BAF_mean)[min(out\$iRight)-1],sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)

#BAF_SD(high)
out <- getOutliers(sort(a\$BAF_SD))
x <- array(dim=c(1,length(a\$BAF_SD)))
x[out\$iRight]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$BAF_SD)) }

x[1: min(out\$iRight)-1]="black"
#plot(sort(a\$BAF_SD), col=x)
plot(sort(a\$BAF_SD), col=x, main=paste("BAF_SD <", sort(a\$BAF_SD)[min(out\$iRight)],"Fail=", out\$nOut[2]),xlab="Sample Index")
q_Lrrsd<-tail(a[order(a\$BAF_SD), ],max(out\$iRight)-min(out\$iRight)+1)
report<-q_Lrrsd[,c("ChipID","BAF_SD")]
report[1:length(report\$ChipID),3]= paste("BAF_SD<", sort(a\$BAF_SD)[min(out\$iRight)] ,sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)

#BAF_DRIFT(high)
out <- getOutliers(sort(a\$BAF_DRIFT))
x <- array(dim=c(1,length(a\$BAF_DRIFT)))
x[out\$iRight]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$BAF_DRIFT)) }

x[1: min(out\$iRight)-1]="black"
#plot(sort(a\$BAF_DRIFT), col=x)
plot(sort(a\$BAF_DRIFT), col=x, main=paste("BAF_DRIFT <", sort(a\$BAF_DRIFT)[min(out\$iRight)],"Fail=", out\$nOut[2]),xlab="Sample Index")
q_Lrrsd<-tail(a[order(a\$BAF_DRIFT), ],max(out\$iRight)-min(out\$iRight)+1)
report<-q_Lrrsd[,c("ChipID","BAF_DRIFT")]
report[1:length(report\$ChipID),3]= paste("BAF_DRIFT<", sort(a\$BAF_DRIFT)[min(out\$iRight)] ,sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)

#WF(low/high)
out <- getOutliers(sort(a\$WF),method="II")
x <- array(dim=c(1,length(a\$WF)))
x[out\$iRight]="blue"
x[out\$iLeft]="blue"

if(max(out\$iLeft)==-Inf) { out\$iLeft = 0 }
if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(a\$WF)) }

x[(max(out\$iLeft)+1):(min(out\$iRight)-1)]="black"
#plot(sort(a\$WF), col=x,ylim=c(-.01,.01))
plot(sort(a\$WF), col=x, main=paste(sort(a\$WF)[max(out\$iLeft)+1],"<= WF <=",sort(a\$WF)[min(out\$iRight)-1],"Fail=", out\$nOut[1]+ out\$nOut[2]),xlab="Sample Index") # out\$limit was not correct
q_WF_Low <-head(a[order(a\$WF), ],max(out\$iLeft))
q_WF_High<-tail(a[order(a\$WF), ], max(out\$iRight)-min(out\$iRight)+1)
report<-q_WF_Low[,c("ChipID","WF")]
report[1:length(report\$ChipID),3]= paste(sort(a\$WF)[max(out\$iLeft)+1],"<=WF<=",sort(a\$WF)[min(out\$iRight)-1],sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
report<-q_WF_High[,c("ChipID","WF")]
report[1:length(report\$ChipID),3]= paste(sort(a\$WF)[max(out\$iLeft)+1],"<=WF<=",sort(a\$WF)[min(out\$iRight)-1],sep="")
write.table(report, file = paste("$out","QC_RemoveIDs_LRRmean_BAFmean_BAFSD_BAFDRIFT_WF.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}

END
print R "dev.off()\n";
close (R);
system ("R CMD BATCH $out.AllRes.R");
###Reopen STDOUT screen output for final results
select(STDOUT);
if ($? == 0) {
  #print "  - command executed successfully!\n";
  $c="head $out"."QC_RemoveIDs.txt";print"=======The Remove Sample IDs Report first 10 lines=======\n";$output=`$c`;if($noSampleQC){$output="NA: noSampleQC";$c="echo -n > $out"."QC_RemoveIDs.txt";`$c`};print"$output\n";  
}
else {
  print "Error: cannot execute system command R CMD BATCH $out.AllRes.R PLEASE TRY RUNNING \"export R_LIBS=./\" \n";
  print "If still error, check $out.AllRes.Rout (below) for specific error and if not resolveable submit to parsecnv\@gmail.com.\n";
  $c="tail $out.AllRes.Rout";print"=======The R Error log last 10 lines=======\n";$output=`$c`;print"$output\n";
}

#unlink ("$out.AllRes.R","$out.AllRes.Rout", ".RData");


$c="perl ".$MyDirPre."FilterCNV.pl $out"."QC_RemoveIDs.txt $rawcnv 5 remove";
`$c`;
if($noSampleQC){$c="cp $rawcnv $rawcnv"."_remove_".$out."QC_RemoveIDs.txt";`$c`;}
$c="sed 's/[^\ ]*=//g' $rawcnv"."_remove_".$out."QC_RemoveIDs.txt | sed 's/,//g' > a; awk '{print \$0\"\\t\"NR}' a > $rawcnv"."_remove_".$out."QC_RemoveIDs_JustNum.txt";
`$c`;
if($noSampleQC){}else{$c="awk '{print \$1}' $out"."QC_RemoveIDs.txt | sed '1d' | sort -u | wc -l"}
if($noSampleQC){$CountSamplesQCRemoved=0;}else{$CountSamplesQCRemoved=`$c`};
chomp($CountSamplesQCRemoved);

#TO DO Print Percent Samples and Calls Removed By QC

#TO DO Separate By CN State

$c="cat $rawcnv"."_remove_".$out."QC_RemoveIDs.txt | awk '{if(\$5~\"^/\"){prefix=\".\"}else{prefix=\"./\"};print \$1,\$2,\$3,\$4,prefix\$5,\$6,\$7,\$8}' | sed 's/chrX/chr23/' | sed 's/chrY/chr24/' | sed 's/state4,cn=2/state5,cn=3/' > $rawcnv"."_wPath";
`$c`;
$c="awk '{if(\$5~\"^/\"){prefix=\".\"}else{prefix=\"./\"};print \$1,\$2,\$3,\$4,prefix\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13,\$14}' $log > $log"."_wPath";
`$c`;
$c="perl ".$MyDirPre."../PerlModules/filter_cnv.pl $rawcnv"."_wPath -qclogfile $log"."_wPath -qcsumout QCsum.qcsum -out goodCNV.good.cnv -chroms 1-22";
@o=`$c`;print(join('',@o)."\n");
$c="mkdir tmp; Rscript ".$MyDirPre."R/Dependencies.R; Rscript ".$MyDirPre."R/R_script_convert_raw_cnv.R $rawcnv"."_wPath QCsum.qcsum ".$MyDirPre."./R";
@o=`$c`;print(join('',@o)."\n");
$c="Rscript ".$MyDirPre."R/R_script_calculate_quality_score.R ".$MyDirPre."./R";
@o=`$c`;print(join('',@o)."\n");
$c="awk '{print \$1\"_chr\"\$2\":\"\$3\"-\"\$4\"\t\"\$NF}' ".$MyDirPre."log_CNV_summary_dataframe.txt > atob";
`$c`;
$c="awk '{gsub(/.*\\//,\"\",\$5);print \$5\"_\"\$1\"\t\"\$0}' $rawcnv"."_remove_".$out."QC_RemoveIDs.txt > a";
`$c`;
$c="perl ".$MyDirPre."../PerlModules/Vlookup.pl atob a | sed 's/[^\\t]*\\t//' | awk '{print \$0\"\\t\"NR}' > $rawcnv"."_remove_".$out."QC_RemoveIDs.txt"."_wMaceQualityScore";
`$c`;

if($stratifyCN)
{	$c="awk '{print > \"$rawcnv"."_remove_"."$out"."QC_RemoveIDs_JustNum.txt_\"\$4}' $rawcnv"."_remove_"."$out"."QC_RemoveIDs_JustNum.txt";
	print($c);`$c`;
	$c="ls $rawcnv"."_remove_"."$out"."QC_RemoveIDs_JustNum.txt_*";
	@CallFiles=split('\n',`$c`);
	chomp(@CallFiles);
	print("@CallFiles\n");
	$CountCallFiles=$#CallFiles+1;
}
else
{
	$c="ls $rawcnv"."_remove_"."$out"."QC_RemoveIDs_JustNum.txt";
	@CallFiles=`$c`;
	chomp(@CallFiles);
	print("@CallFiles\n");
	$CountCallFiles=$#CallFiles+1;
}
open (R, ">$out.AllRes2.R");
print R <<END;

library(extremevalues,lib.loc="./")

write.table("ChipID Value Metric CallIndexRemove", file = paste("$out","QC_RemoveCalls.txt",sep="")  , append = FALSE,quote=FALSE,row.names=FALSE, col.names =FALSE)
i=1
print($CountCallFiles)
if($CountCallFiles>1)
{
filename <- system(paste("ls $rawcnv","_remove_","$out","QC_RemoveIDs_JustNum.txt_*",sep=""),intern=TRUE)
loopTimes=$CountCallFiles+1
}
if($CountCallFiles<=1)
{
filename <- system(paste("ls $rawcnv","_remove_","$out","QC_RemoveIDs_JustNum.txt",sep=""),intern=TRUE)
loopTimes=$CountCallFiles
}

while(i<= loopTimes)
{
if(i>2){CN=i}else{CN=i-1}
pdf(file=paste("$out","QC_Call_Plot_",CN,".pdf",sep=""), width=3,height=3, pointsize=6)

par(mfrow=c(2,2),las=1,mar=c(4,4.2,1.5,1.3))

if("$rawcnv"!=1){b=read.table(filename[i],header=FALSE, comment.char="")}

if("$rawcnv"!=1){
print("numsnp")
out <- getOutliers(sort(b\$V2),method="II")
x <- array(dim=c(1,length(b\$V2)))
x[1: min(out\$iRight)-1]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(b\$V2)) }

if("$numsnp"!=-999) { 
out\$limit[2] <- as.numeric("$numsnp")
out\$iRight = c(1:(sum(b\$V2 <= as.numeric("$numsnp"))))
x <- array(dim=c(1,length(b\$V2)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(b\$V2 <= as.numeric("$numsnp"))

x[max(out\$iRight):length(b\$V2)]="black"
#plot(sort(b\$V2), col=x)
plot(sort(b\$V2), col=x, main=paste("NumSNP <=",round(out\$limit[2],4) ,"Fail=", out\$nOut[2]),xlab="CNV Call Index",ylab="NumSNP")  ### was sort(b\$V2)[min(out\$iRight)]
q_CountCNV <-head(b[order(b\$V2), ], max(out\$iRight))
report<-q_CountCNV [,c("V5","V2","V2","V9")]
report[1:length(report\$V5),3]= paste("NumSNP<=", round(out\$limit[2],4) ,sep="") ### was sort(b\$V2)[min(out\$iRight)]
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
else{
x[out\$iRight]="black"
plot(sort(b\$V2), col=x, main=paste("NumSNP <=",sort(b\$V2)[min(out\$iRight)] ,"Fail=", min(out\$iRight)-1),xlab="CNV Call Index",ylab="NumSNP")
q_CountCNV <-head(b[order(b\$V2), ], min(out\$iRight)-1)
report<-q_CountCNV [,c("V5","V2","V2","V9")]
report[1:length(report\$V5),3]= paste("NumSNP<=", sort(b\$V2)[min(out\$iRight)] ,sep="") 
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}

print("length")
out <- getOutliers(sort(b\$V3),method="II")
x <- array(dim=c(1,length(b\$V3)))
x[1: min(out\$iRight)-1]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(b\$V3)) }

if("$length"!=-999) { 
out\$limit[2] <- as.numeric("$length")
out\$iRight = c(1:(sum(b\$V3 <= as.numeric("$length"))))
x <- array(dim=c(1,length(b\$V3)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(b\$V3 <= as.numeric("$length"))

x[max(out\$iRight):length(b\$V3)]="black"
plot(sort(b\$V3), col=x, main=paste("Length(bp) <=",round(out\$limit[2],4) ,"Fail=", out\$nOut[2]),xlab="CNV Call Index",ylab="Length(bp)")  ### was sort(b\$V2)[min(out\$iRight)]
q_CountCNV <-head(b[order(b\$V3), ], max(out\$iRight))
report<-q_CountCNV [,c("V5","V3","V3","V9")]
report[1:length(report\$V5),3]= paste("Length(bp)<=", round(out\$limit[2],4) ,sep="") ### was sort(b\$V2)[min(out\$iRight)]
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
else{
x[out\$iRight]="black"
#plot(sort(b\$V2), col=x)
plot(sort(b\$V3), col=x, main=paste("Length(bp) <=",sort(b\$V3)[min(out\$iRight)] ,"Fail=", min(out\$iRight)-1),xlab="CNV Call Index",ylab="Length(bp)")
q_CountCNV <-head(b[order(b\$V3), ], min(out\$iRight)-1)
report<-q_CountCNV [,c("V5","V3","V3","V9")]
report[1:length(report\$V5),3]= paste("Length(bp)<=", sort(b\$V3)[min(out\$iRight)] ,sep="") 
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
print("conf")
out <- getOutliers(sort(b\$V8),method="II")
x <- array(dim=c(1,length(b\$V8)))
x[1: min(out\$iRight)-1]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(b\$V8)) }

if("$confidence"!=-999) { 
out\$limit[2] <- as.numeric("$confidence")
out\$iRight = c(1:(sum(b\$V8 <= as.numeric("$confidence"))))
x <- array(dim=c(1,length(b\$V8)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(b\$V8 <= as.numeric("$confidence"))

x[max(out\$iRight):length(b\$V8)]="black"
plot(sort(b\$V8), col=x, main=paste("Confidence <=",round(out\$limit[2],4) ,"Fail=", out\$nOut[2]),xlab="CNV Call Index",ylab="Confidence")  ### was sort(b\$V2)[min(out\$iRight)]
q_CountCNV <-head(b[order(b\$V8), ], max(out\$iRight))
report<-q_CountCNV [,c("V5","V8","V8","V9")]
report[1:length(report\$V5),3]= paste("Confidence<=", round(out\$limit[2],4) ,sep="") ### was sort(b\$V2)[min(out\$iRight)]
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
else{
x[out\$iRight]="black"
plot(sort(b\$V8), col=x, main=paste("Confidence <=",sort(b\$V8)[min(out\$iRight)] ,"Fail=", min(out\$iRight)-1),xlab="CNV Call Index",ylab="Confidence")
q_CountCNV <-head(b[order(b\$V8), ], min(out\$iRight)-1)
report<-q_CountCNV [,c("V5","V8","V8","V9")]
report[1:length(report\$V5),3]= paste("Confidence<=", sort(b\$V8)[min(out\$iRight)] ,sep="") 
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}

print("quality")
b=read.table(paste("$rawcnv","_remove_","$out","QC_RemoveIDs.txt_wMaceQualityScore",sep=""),header=FALSE, comment.char="")
out <- getOutliers(sort(b\$V9),method="II")
x <- array(dim=c(1,length(b\$V9)))
#print(out)
x[out\$iLeft]="blue"

if(min(out\$iRight)==Inf) { out\$iRight = sum(!is.na(b\$V9)) }

if("$quality"!=-999) { 
out\$limit[2] <- as.numeric("$quality")
out\$iRight = c(1:(sum(b\$V9 <= as.numeric("$quality"))))
x <- array(dim=c(1,length(b\$V9)))
x[out\$iRight]="blue"
out\$nOut[2]=sum(b\$V9 <= as.numeric("$quality"))

x[max(out\$iRight):length(b\$V9)]="black"
plot(sort(b\$V9), col=x, main=paste("Quality <=",round(out\$limit[2],4) ,"Fail=", out\$nOut[2]),xlab="CNV Call Index",ylab="MaceQualityScore")  ### was sort(b\$V2)[min(out\$iRight)]
q_CountCNV <-head(b[order(b\$V9), ], max(out\$iRight))
report<-q_CountCNV [,c("V5","V9","V9","V10")]
report[1:length(report\$V5),3]= paste("Quality<=", round(out\$limit[2],4) ,sep="") ### was sort(b\$V2)[min(out\$iRight)]
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
else{
x[]="black" #x[(max(out\$iLeft)+1):length(b\$V9)]="black"
plot(sort(b\$V9), col=x, main=paste("MaceQualityScore <=",sort(b\$V9)[min(out\$iRight)] ,"Fail=", min(out\$iRight)-1),xlab="CNV Call Index",ylab="MaceQualityScore")
q_CountCNV <-head(b[order(b\$V9), ], max(out\$iLeft)-1)
report<-q_CountCNV [,c("V5","V9","V9","V10")]
report[1:length(report\$V5),3]= paste("MaceQualityScore<=", sort(b\$V9)[min(out\$iRight)] ,sep="") 
write.table(report, file = paste("$out","QC_RemoveCalls.txt",sep=""), append = TRUE,quote=FALSE,row.names=FALSE, col.names =FALSE)
}
}

AllQC2<-read.table(paste("$out","QC_RemoveCalls.txt",sep="") ,header=FALSE,comment.char="")

mtext(paste(length(table(AllQC2\$V4)),"unique calls fail any QC"), 1, padj=2.4)###put text below
i=i+1


END
print R "while (!is.null(dev.list()))  dev.off()\n}\n";
close (R);
system ("R CMD BATCH $out.AllRes2.R");
###Reopen STDOUT screen output for final results
select(STDOUT);
if ($? == 0) {
  #print "  - command executed successfully!\n";
  $c="head $out"."QC_RemoveCalls.txt";print"=======The Remove Call IDs Report first 10 lines=======\n";$output=`$c`;if($noCallQC){$output="NA: noCallQC"};print"$output\n";  
}
else {
  print "Error\n";
  print "Check $out.AllRes2.Rout (below) for specific error and if not resolveable submit to parsecnv\@gmail.com.\n";
  $c="tail $out.AllRes2.Rout";print"=======The R Error log last 10 lines=======\n";$output=`$c`;print"$output\n";
}
$c="awk '{print \$4}' $out"."QC_RemoveCalls.txt > $out"."QC_RemoveCalls.txt_Indexes";`$c`;
$c="perl ".$MyDirPre."FilterCNV.pl $out"."QC_RemoveCalls.txt_Indexes $rawcnv"."_remove_"."$out"."QC_RemoveIDs.txt_wMaceQualityScore"." 10 remove";`$c`;
if($noCallQC){$c="cp $rawcnv"."_remove_"."$out"."QC_RemoveIDs.txt_wMaceQualityScore $rawcnv"."_remove_"."$out"."QC_RemoveIDs.txt_wMaceQualityScore_remove_"."$out"."QC_RemoveCalls.txt_Indexes; echo -n > $out"."QC_RemoveCalls.txt";`$c`;}
$c="awk '{print \$4}' $out"."QC_RemoveCalls.txt | sed '1d' | sort -u | wc -l";
$CountCallsQCRemoved=`$c`;
chomp($CountCallsQCRemoved);
print "Sample $CountSamplesQCRemoved and Call $CountCallsQCRemoved QC Filtered File Written to: $rawcnv"."_remove_"."$out"."QC_RemoveIDs.txt_wMaceQualityScore_remove_".$out."QC_RemoveCalls.txt_Indexes\n";
$c="convert -quality 100 -density 600x600 $out"."QC_Call_Plot.pdf $out"."QC_Call_Plot.jpg";
#`$c`;
$c="convert -quality 100 -density 600x600 $out"."QC_Plot.pdf $out"."QC_Plot_2.pdf $out"."QC_Call_Plot.pdf -append $out"."QC_Sample+Call_Plot.jpg";
#`$c`;
