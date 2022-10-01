#!/usr/bin/perl
##Joseph Glessner PhD
##Center for Applied Genomics, Children's Hospital of Philadelphia
##University of Pennsylvania Perelman School of Medicine
##Developed under advising of Hakon Hakonarson
##Copyright (C) GNU General Public License 2022 Joseph Glessner parsecnv@gmail.com
##Linux bash environment is best, but Windows Cygwin can also work.
##To run with all bash commands printed: sed 's#;$o=`$c`;#;print("Run $c\\n");$o=`$c`;#' ParseCNV2_Plink2.pl > ParseCNV2_Plink2_PrintAllCmd.pl
if (   @ARGV < 1
    || "@ARGV" =~ /-h/
    || "@ARGV" =~ /-help/
    || "@ARGV" =~ /--h/
    || "@ARGV" =~ /--help/)
{
    print "\nUSAGE MESSAGE:\n
perl ParseCNV2.pl -i CNV_Calls[.vcf/.rawcnv] [-c Cases.txt/-q Samples_QuantitativeTrait.txt] [arguments]

Required arguments:

-i	input[vcf/rawcnv/txt]
-c	cases
-b	build (if vcf can be inferred)

Optional arguments:

-o	output
-covar	covariates file
-stat	statistic (fisher, logistic, or linear) or RvTests:
	Single variant: score, wald, exact, dominantExact, famLRT, famScore, famGrammarGamma, firth
	Burden: cmc, zeggini, mb, fp, exactCMC, cmcWald, rarecover, cmat, famcmc, famzeggini
	Variable threshold: price, analytic, famAnalytic
	Kernel: skat, skato, kbac, famSkat
	Meta-Analysis: score, dominant, recessive, cov, bolt, boltCov
-no_freq	No CNV frequency filtering, default is rare (minor allele frequency <=0.01) recurrent (minor allele count >= 2)
-p	merge p variation
-d	merge distance
-t	transmission disequilibrium test
-m	max p inclusion
-bfile	plink bed file of SNP genotypes
-qc	quality control upfront
-log	penncnv log
-q	quantitative trait
-batch	batch
-h	help

";
    exit;
}
$myCommandLine = "Running Command: perl $0";
for ( $i = 0 ; $i <= $#ARGV ; $i++ ) {
    $myCommandLine .= " $ARGV[$i]";
}
$myCommandLine .= "\n";

$MyDirectoryPathPrefix = $_;
print "MyDirectoryPathPrefix: $MyDirectoryPathPrefix\n";
use lib './PerlModules';

#use lib $MyDirectoryPathPrefix;
use lib $MyDirectoryPathPrefix. "PerlModules";

#use lib $_, $_."PerlModules";
use Getopt::Long;
GetOptions(
    'i=s'     => \$input,
    'o=s'     => \$out,
    'b=s'     => \$build,
    'p=f'     => \$mergePVar,
    'd=i'     => \$mergeDist,
    't'       => \$tdt,
    'm=f'     => \$maxPInclusion,
    'bfile=s' => \$bfile,
    'qc'      => \$qc,
    'log=s'   => \$log,
    'c=s'     => \$cases,
    'q=s'     => \$quantitativeTrait,
    'batch=s' => \$batch,
    'covar=s' => \$covar,
    'stat=s'  => \$statistic,
    'no_freq' => \$noFrequencyFilter
);
print "statistic=".$statistic."\n";
BEGIN { ( $_ = $0 ) =~ s{[^\\\/]+$}{}; $_ ||= "./" } ##In case running elsewhere
$MyDirectoryPathPrefix = $_;

#print "MyDirectoryPathPrefix: $MyDirectoryPathPrefix\n";
if ( not $out ) {
    $out = "ParseCNV_";
}
if ( not $input ) {
    print
"ERROR: --i is a required input VCF or PennCNV rawcnv file or list of such files.\n";
    exit;
}
$inputNoPath = $input;
$inputNoPath =~ s/.*\///;    #get rid of potential path prefix
$dir = 'temp';
unless ( -d $dir ) {
    mkdir $dir or die;
}
if ($noFrequencyFilter) {
    $FrequencyFilter = "";
}
else {
    $FrequencyFilter = "--mac 2 --max-maf 0.01";
}
#print "MyDirectoryPathPrefix: $MyDirectoryPathPrefix\n";
#Unzip larger files on github that had to be compressed to allow upload
if ( !( -e $MyDirectoryPathPrefix . "PerlModules/plink" ) ) {
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/plink.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/plink2.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules";
    $o = `$c`;
    $c = "mkdir ".$MyDirectoryPathPrefix."PerlModules/executable;" . "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/rvtest.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules/executable";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/vcfToBedpe.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/bcftools.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "PerlModules/tabix.zip -d "
      . $MyDirectoryPathPrefix
      . "PerlModules";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "GeneRef/$build"
      . "_gc5Base_SimFormat_AllCol.sorted.zip -d "
      . $MyDirectoryPathPrefix
      . "GeneRef";
    $o = `$c`;
    $c =
        "unzip "
      . $MyDirectoryPathPrefix
      . "GeneRef/$build"
      . "_knownGene_Exons_SimFormat_AllCol_UniqueIDs.zip -d "
      . $MyDirectoryPathPrefix
      . "GeneRef";
    $o = `$c`;
    if ( -e $MyDirectoryPathPrefix . "GeneRef/$build" . "_knownGene.txt.zip" ) {
     
	$c =
            "unzip "
          . $MyDirectoryPathPrefix
          . "GeneRef/*"
          . "_knownGene.txt.zip -d "
          . $MyDirectoryPathPrefix
          . "GeneRef";
        $o = `$c`;
        $c =
            "unzip "
          . $MyDirectoryPathPrefix
          . "GeneRef/*"
          . "_kgXref.txt.zip -d "
          . $MyDirectoryPathPrefix
          . "GeneRef";
        $o = `$c`;
    }
}

#Check presence and version of dependencies
$cmd = "which bash";
$o   = `$cmd`;
print LOG "$o";
if ( $o =~ " no " ) { print "ERROR: bash not found!\n"; }
$cmd = "which R";
$o   = `$cmd`;
print LOG "$o";
if ( $o =~ " no " ) { print "ERROR: R not found!\n"; }
$cmd = "which " . $MyDirectoryPathPrefix . "PerlModules/plink";
$o   = `$cmd`;
print LOG "$o";
if ( $o =~ " no " ) { print "ERROR: plink not found!\n"; }
$cmd = "which " . $MyDirectoryPathPrefix . "PerlModules/plink2";
$o   = `$cmd`;
print LOG "$o";
if ( $o =~ " no " ) { print "ERROR: plink2 not found!\n"; }

#$c="awk '/MemAvailable/ {printf \$2/1000}' /proc/meminfo";$o=`$c`;print $o." MB RAM Available\n";print LOG $o." MB RAM Available\n";
$c = "free -m | awk '{if(NR==2)printf\$4}'";
$o = `$c`;
print $o. " MB RAM Available\n";
print LOG $o . " MB RAM Available\n";
if ( -e "plink.sexcheck" ) {
    $c = "rm plink.sexcheck";
    $o = `$c`;
}

#sub function definition to reduce redundancy .txt vs. vcf/rawcnv input on the command line but have to be careful about scope of variables or pass references to variables like penncnv
#sub parse_vcf_vcfgz
#{
#	$input_vcf_vcfgz = $_[0];
#return ($res);
#}
#$sub_out = parse_vcf_vcfgz(6, 5);
#print $sub_out;

##VCF sub functions
sub vcf_check_build
{
	$vcf = $_[0];
	#Genome Build not always found in vcf header but sometimes so check to confirm
    $c = $cat_zcat . " $input | grep ^# | egrep 'NCBI34|hg16|b16|hs34' | wc -l";
    $o = `$c`;
    if ( $o > 0 ) { $BuildsInVcfHeader .= " NCBI34|hg16|b16|hs34"; }
	
    $c = $cat_zcat . " $input | grep ^# | egrep 'NCBI35|hg17|b17|hs35' | wc -l";
    $o = `$c`;
    if ( $o > 0 ) { $BuildsInVcfHeader .= " NCBI35|hg17|b17|hs35"; }
	
    $c = $cat_zcat . " $input | grep ^# | egrep 'NCBI36|hg18|b18|hs36' | wc -l";
    $o = `$c`;
    if ( $o > 0 ) { $BuildsInVcfHeader .= " NCBI36|hg18|b18|hs36"; }
	
    $c = $cat_zcat . " $input | grep ^# | egrep 'GRCh37|hg19|b19|hs37' | wc -l";
    $o = `$c`;
    if ( $o > 0 ) { $BuildsInVcfHeader .= " GRCh37|hg19|b19|hs37"; }
	
    $c = $cat_zcat . " $input | grep ^# | egrep 'GRCh38|hg38|b38|hs38' | wc -l";
    $o = `$c`;
    if ( $o > 0 ) { $BuildsInVcfHeader .= " GRCh38|hg38|b38|hs38"; }

    if ( $BuildsInVcfHeader ne "" ) {
        return ("Found Genome Builds in VCF Header:" . $BuildsInVcfHeader . "\n" );
    }
	else{
		return ("Found No Genome Builds matches in VCF Header\n" );
	}
}

sub vcf_extract_del_dup_IDs
{
	#vcf -> del and dup ID (Probe ID) lists
	$c =
        "$cat_zcat $input | awk '{if(\$0~/^#/){print > \"temp/$out"
      . "$inputNoPath"
      . "_del.vcf\";print > \"temp/$out"
      . "$inputNoPath"
      . "_dup.vcf\"}else{print}}' | awk 'BEGIN{IGNORE_CASE=1}{ALT=\$5;ALT_copy=\$5;count_del_matches=gsub(\/DEL|0|1\/, \"\",ALT); count_dip_matches=gsub(\/2\/, \"\",ALT_copy); if(length(count_dip_matches)>0 ){count_dip_matches_total++}; if( length(count_del_matches)>0 && count_del_matches>0 ){print \$3 > \"temp/$out"
      . "$inputNoPath"
      . "_DelSnpsFromVcf\";dels++;print >> \"temp/$out"
      . "$inputNoPath"
      . "_del.vcf\"}else{print \$3 > \"temp/$out"
      . "$inputNoPath"
      . "_DupSnpsFromVcf\";dups++;print >> \"temp/$out"
      . "$inputNoPath"
      . "_dup.vcf\"}}END{print \"dels=\"dels,\"dups=\"dups}'";
    $o = `$c`;
    print "c=" . $c . "\n";
    print $o. "\n";
}
sub vcf_to_pfile
{
	#vcf -> pfile (pgen,pvar,psam)
        $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink2 --vcf $input --make-pgen --allow-extra-chr --out temp/$out"
      . "$inputNoPath.del --extract temp/$out"
      . "$inputNoPath"
      . "_DelSnpsFromVcf";
      print "Running $c\n";
    $o = `$c`;    #All Deletion CN=0,1 pgen del file
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink2 --vcf $input --make-pgen --allow-extra-chr --out temp/$out"
      . "$inputNoPath.dup --extract temp/$out"
      . "$inputNoPath"
      . "_DupSnpsFromVcf";
    $o = `$c`;    #All Duplication CN>2 pgen dup file
	$c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink2 --vcf $input --make-pgen --allow-extra-chr --out temp/$out"
      . "$inputNoPath";
    $o = `$c`;    #All together pgen deldup file
    $c =
    "grep -v ^# temp/$out"
  . "$inputNoPath.pvar | awk -F\"\t\" '{if(\$5~/,/)print\$3}' > temp/$out". "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath --make-bed --allow-extra-chr --out temp/$out"
  . " --exclude temp/$out". "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath.del --make-bed --allow-extra-chr --out temp/$out"
  . "del --exclude temp/$out". "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath.dup --make-bed --allow-extra-chr --out temp/$out"
  . "dup --exclude temp/$out". "ExcludeMultiallelic";
$o = `$c`;

#$c = $MyDirectoryPathPrefix . "PerlModules/./plink --bfile temp/$out --impute-sex --make-bed --out temp/$out";$o = `$c`;
#$c = $MyDirectoryPathPrefix . "PerlModules/./plink --bfile temp/$out" . "del --impute-sex --make-bed --out temp/$out" . "del";$o = `$c`;
#$c = $MyDirectoryPathPrefix . "PerlModules/./plink --bfile temp/$out" . "dup --impute-sex --make-bed --out temp/$out" . "dup";$o = `$c`;
}

sub pvar_chr_pos_link
{
	#Define CHR_POS to ID (Probe ID) links for tracking in later reporting
	#ID -> CHR_POS
	#CHR_POS -> ID
	#CHR_POS -> REF
	#CHR_POS -> ALT
	#CHR_POS -> ALTtwo (Only the first ALT if comma separated list (multiallelic)
	my %PlinkVarID_to_CHR_POS_pvar;
    my %PlinkCHR_POS_to_VarID_pvar;
    my %PlinkCHR_POS_to_REF;
    my %PlinkCHR_POS_to_ALT;
    my %PlinkCHR_POS_to_ALTtwo;
    $myPvar = "temp/$out" . "$inputNoPath" . ".pvar";
    open( my $input_pvar, "<$myPvar" ) or die "cannot open $myPvar: $!\n";
    while ( my $line_pvar = <$input_pvar> ) {
        chomp($line_pvar);
        if ( !( $line_pvar =~ /^\#/ ) ) {
            my @vals = split( /\t/, $line_pvar );
            $PlinkVarID_to_CHR_POS_pvar{ $vals[2] } = $vals[0] . "_" . $vals[1];
            $ChrPos                                 = $vals[0] . "_" . $vals[1];
            $PlinkVarID_to_CHR_POS_pvar{ $vals[2] } = $ChrPos;
            $PlinkCHR_POS_to_VarID_pvar{$ChrPos}    = $vals[2];
            $PlinkCHR_POS_to_REF{$ChrPos}           = $vals[3];
            $PlinkCHR_POS_to_ALT{$ChrPos}           = $vals[4];
            if ( $vals[4] =~ /,/ ) { $vals[4] =~ s/,*//; }
            $PlinkCHR_POS_to_ALTtwo{$ChrPos} = $vals[4];
        }
    }
    close($input_pvar);
    #print "\$ChrPos="
    #  . $ChrPos
    #  . " \$PlinkCHR_POS_to_REF{\$ChrPos}="
    #  . $PlinkCHR_POS_to_REF{$ChrPos} . "\n";
    #$test_CHR_POS_to_VarID_pvar = "1_812283";
    #$VCF_VarID = $PlinkCHR_POS_to_VarID_pvar{$test_CHR_POS_to_VarID_pvar};
    #print "test_CHR_POS_to_VarID_pvar="
    #  . "$test_CHR_POS_to_VarID_pvar\t$VCF_VarID\n";
}

sub vcf_to_bedpe
{

	#Extract vcf lines with END= in INFO field for vcfToBedpe input
	#Might not be needed anymore
    $c =
        $cat_zcat
      . " $input | grep -v ^# | grep \"END=\" > temp/$out"
      . "$inputNoPath" . "_wEND";
    $o = `$c`;
    print "$o\n";

#;END= and ;SVLEN= would be more specific matching but miss if it was the first entry
#1       947113  UW_VH_22703     G       <CN0>   100     PASS    SVTYPE=DEL;CIEND=0,82;CIPOS=-179,0;END=948003;CS=DEL_union;MC=BI_GS_DEL1_B5_P0001_156;AC=4;AF=0.00079872;NS=2504;AN=5008;EAS_AF=0.0;EUR_AF=0.0;AFR_AF=0.003;AMR_AF=0.0;SAS_AF=0.0 GT      0|0

    $c =
        "python2 "
      . $MyDirectoryPathPrefix
      . "PerlModules/./vcfToBedpe -i temp/$out"
      . "$inputNoPath"
      . "_wEND > temp/$out"
      . "$inputNoPath"
      . "_wEND_vcfToBedpe";
    $o = `$c`;

    #average 2 3  5 6 since CIPOS used in bedpe
    open( OUT_vcfToBedpe, "temp/$out" . "$inputNoPath" . "_wEND_vcfToBedpe" );
    open( OUT_BedpeToMap, ">temp/$out" . ".map_presorting" );
    open( OUT_BedpeToRawcnv2,
        ">temp/$out" . "$inputNoPath" . "_wEND_vcfToBedpe.rawcnv2" );
    open( NoCNmatchedOUT, ">NoCNmatched.out" );
    $NoCNmatched              = 0;
    $numberMultipleALTs_total = 0;
    while (<OUT_vcfToBedpe>) {
        chomp;
        @_ = split( /\t/, $_ );
        my $rounded_start_pos = int( ( $_[1] + $_[2] ) / 2 + 0.5 );
        my $rounded_end_pos   = int( ( $_[4] + $_[5] ) / 2 + 0.5 );
        $CHR_POS = $_[0] . "_" . $rounded_start_pos;
        $CHR_POS_to_END{$CHR_POS} = $rounded_end_pos;
        $PlinkVarID_to_CHR_POS_vcfToBedpe{ $_[6] } =
          $_[1] . "_" . $rounded_start_pos;

        print OUT_BedpeToMap $_[0] . "_" . $rounded_start_pos . "\n";
        print OUT_BedpeToMap $_[0] . "_" . $rounded_end_pos . "\n";
    }
    close(OUT_BedpeToRawcnv2);
    close(OUT_BedpeToMap);
}

sub extract_info_vcf
{
	#Idea: awk tab info then split by ;
	#Using all variants here, not just those with an "END=" value
    #CHR START(POS) END SNP_ID REF ALT SITEPOST
    $c =
        $cat_zcat
      . " $input | grep -v ^## | awk '{if(NR==1){for(j=10;j<=NF;j++){SIDs[j]=\$j}}else{count=split(\$8,a,/;/);for(i=1;i<=count;i++){if(a[i]~/^END=/){gsub(/END\=/,\"\",a[i]);END_FOUND=1;myEND=a[i];if(myEND==\$2){myEND++}}};for(i=1;i<=count;i++){if(a[i]~/^SITEPOST=/){gsub(/SITEPOST\=/,\"\",a[i]);SITEPOST_FOUND=1;mySITEPOST=a[i]}};if(!(mySITEPOST)){mySITEPOST=1;};for(i=1;i<=count;i++){if(a[i]~/^SVLEN=/){gsub(/SVLEN\=/,\"\",a[i]);SVLEN_FOUND=1;mySVLEN=a[i]}};if(END_FOUND==1){}else if(SVLEN_FOUND==1){myEND=\$2+mySVLEN}else{myEND=\$2+1};for(j=10;j<=NF;j++){if(\$j!~/^0\\\/0/){Calls++;SamplesCols_wCalls[Calls]=\$j}};for(k=1;k<=Calls;k++){print \$1\"\\t\"\$2\"\\t\"myEND\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"mySITEPOST\"\\t\"SIDs[SamplesCols_wCalls[k]] > \"temp/$out"
      . "$inputNoPath"
      . "_wEND_CHR_POS_or_END_SNPID_wENDCol_wSitePostCol.txt\"}END_FOUND=0;SVLEN_FOUND=0;Calls=0;delete SamplesCols_wCalls}}'";
    print "$c\n";
    $o = `$c`;
    print "$o\n";
}

sub cn_extract
{
		# Use the *_wEND_CHR_POS_or_END_SNPID_wENDCol_wSitePostCol.txt file to attempt copy number state extraction
    $NoCNmatched=0;
    open( VCF_POS_and_END,
            "temp/$out"
          . "$inputNoPath"
          . "_wEND_CHR_POS_or_END_SNPID_wENDCol_wSitePostCol.txt" );

    #CHR START(POS) END SNP_ID REF ALT SITEPOST
    open( OUT_BedpeToRawcnv2, ">temp/$out" . "$inputNoPath" . ".rawcnv2" );
	open( OUT_BedpeToRawcnv, ">temp/$out" . "$inputNoPath" . ".rawcnv" );
    open( OUT_DEL_VarIDs,
        ">temp/" . $out . $inputNoPath . "_DelSnpsFromVcf_new" );
    open( OUT_DUP_VarIDs,
        ">temp/" . $out . $inputNoPath . "_DupSnpsFromVcf_new" );
    %CHR_POS_to_END;
    %PlinkVarID_to_CHR_POS;
    while ( $original_line = <VCF_POS_and_END> ) {
        chomp($original_line);
        @a             = split( /\t/, $original_line );
        $myREF         = $a[4];
        $myALT         = $a[5];
        $cn_forRawcnv2 = -1;

        my $number0 = () = $myALT =~ /0/gi;
        if ( $number0 > 0 ) { $cn_forRawcnv2 = 0; }
        my $number1 = () = $myALT =~ /1/gi;
        if ( $number1 > 0 ) { $cn_forRawcnv2 = 1; }
        my $numberDEL = () = $myALT =~ /DEL/gi;
        if ( $number1 > 0 ) { $cn_forRawcnv2 = 1; }
        my @countDEL = $myALT =~ /DEL/gi;
        if ( ( scalar @countDEL ) > 0 ) { $cn_forRawcnv2 = 1; }
        my $number2 = () = $myALT =~ /2/gi;
        if ( $number2 > 0 ) { $cn_forRawcnv2 = 2; }
        my $number3 = () = $myALT =~ /3/gi;
        if ( $number3 > 0 ) { $cn_forRawcnv2 = 3; }
        my @countDUP = $myALT =~ /DUP/gi;
        if ( ( scalar @countDUP ) > 0 ) { $cn_forRawcnv2 = 3; }
        my $number4 = () = $myALT =~ /4/gi;
        if ( $number4 > 0 ) { $cn_forRawcnv2 = 4; }

#my @countCNV = $myALT =~ /CNV/gi;if((scalar @countCNV)>0){$ALT_CNV++;$cn_forRawcnv2=1;}#Multiallelic CNV (both del and dup observed)
#$c="echo $myALT | awk -F, '{printf NF}'";$numberMultipleALTs=`$c`;
#if($numberMultipleALTs>1){$numberMultipleALTs_total++;}
        if ( $cn_forRawcnv2 == -1 ) {
            my $numberGt4 = () = $myALT =~ /\d+/gi;
            if ( $number4 > 0 ) { $cn_forRawcnv2 = 4; }
        }
        if   ( length($myREF) > length($myALT) ) { $cn_forRawcnv2 = 1; }
        else                                     { $cn_forRawcnv2 = 3; }

		#Output rawcnv2 as well for compatibility until transition to vcf stablized
        if ( $cn_forRawcnv2 > -1 ) {
            print OUT_BedpeToRawcnv2 $a[0] . "\t"
              . $a[1] . "\t"
              . $a[2] . "\t"
              . $cn_forRawcnv2
              . "\t" . $a[7] . "\t"
              . $a[6] . "\n";
			  
			%CN_to_StateCN=(0 => "state1,cn=0", 1 => "state2,cn=1", 2 => "state4,cn=2", 3 => "state5,cn=3", 4 => "state6,cn=4");
			print OUT_BedpeToRawcnv $a[0] . ":"
              . $a[1] . "-"
              . $a[2] . " "
			  . "numsnp=" . ($a[2]-$a[1])/1000 . " "
			  . "length=" . ($a[2]-$a[1])+1 . " "
              . $CN_to_StateCN{$cn_forRawcnv2}
              . " " . $a[7] . " "  #could try pulling SID from VCF header, plink2 psam
			  . $a[0] . "_" . $a[1] . " "
			  . $a[0] . "_" . $a[2] . " "
              . "conf=" . $a[6] . "\n";
			#chr1:1-3  numsnp=3 length=3  state1,cn=0 1_A.baflrr startsnp=rs1 endsnp=rs3 conf=1
			
            $CHR_POS                        = $a[0] . "_" . $a[1];
            $CHR_POS_to_END{$CHR_POS}       = $a[2];
            $PlinkVarID_to_CHR_POS{ $_[3] } = $CHR_POS;
            if ( $cn_forRawcnv2 < 2 ) { print OUT_DEL_VarIDs $_[3]; }
            else { print OUT_DUP_VarIDs $_[3]; }    #01234 DelSnpsFromVcf_new
        }
        else {
            $NoCNmatched++;
            print NoCNmatchedOUT "$myALT in $original_line" . "\n";
        }

#vcfToBedpe output format (may not be needed) sample ids in columns but not needed for individual calls
#CHROM_A        START_A END_A   CHROM_B START_B END_B   ID      QUAL    STRAND_A      STRAND_B        TYPE    FILTER  INFO    FORMAT  HG00096 HG00097
    }
    close(NoCNmatchedOUT);
    close(OUT_DEL_VarIDs);
    close(OUT_DUP_VarIDs);
    print
"$NoCNmatched variants listed in NoCNmatched.out had no CN matched out of (DEL,0,1) or (DUP,2,3,4,\>4) so will be excluded.\n";
    $c =
"$cat_zcat $input | grep -v ^\# | awk '{if(\$5~/,/)numMultAlts++}END{printf numMultAlts+0}'";
    $o                        = `$c`;
    $numberMultipleALTs_total = $o;
    print
"$numberMultipleALTs_total variants have multiple ALT values called Multiallelic CNV (both del and dup observed) so these will have multiple association out lines\n";
    close(OUT_BedpeToMap);
}

sub make_map_sorted
{
    $c = "sort -u -k1,1 -k2,2n -t_ temp/$out". ".map_presorting > temp/$out" . ".map";
    $o            = `$c`;
    %h_state      = ();
    %chr_posIndex = ();
    open( MAP, "temp/$out" . ".map" );
    while ( $line = <MAP> ) {
        chomp($line);
        $h_state{$line}      = "-9";
        $chr_posIndex{$line} = $lineNum;

        #$h_id{$line}="";
        @chr_pos        = split( "_", $line );
        $diff[$lineNum] = $chr_pos[1] - $last_pos;
        $last_pos       = $chr_pos[1];
        $lineNum++;
    }
    $diff[$lineNum] = 1;    #Exit at end
}

##RAWCNV sub functions
sub verify_rawcnv
{
    open( RAWCNVFILE,      $input );
    open( RAWCNVCLEANFILE, ">temp/$out" . "$inputNoPath.rawcnv" );
    open( RAWCNVBADFILE,   ">temp/$out" . "$inputNoPath.badcnv" );
    while (<RAWCNVFILE>) {
        chomp;
        if (
m/^(?:chr)?(\w+):(\d+)-(\d+)\s+numsnp=(\d+)\s+length=(\S+)\s+state(\d),cn=(\d)\s+(.+?)\s+startsnp=(\S+)\s+endsnp=(\S+)/
          )
        {
            print RAWCNVCLEANFILE "$_\n";
        }
        else {
            print RAWCNVBADFILE "$_\n";
            $skipped_line++;
        }
    }
    $skipped_line
      and print
"WARNING: $skipped_line lines were skipped due to unrecognizable formats as listed in temp/$out"
      . "$inputNoPath.badcnv\n";
    close(RAWCNVCLEANFILE);
}

sub rawcnv_to_vcf
{
#Old rawcnv -> rawcnv2 and vcf -> rawcnv2 approach
#Now rawcnv -> vcf -> plink2 pgen and vcf -> plink2 pgen for faster computation, more statistical/filtering methods and vcf is very popular now
#$c="awk '{print \$1\"\\t\"\$4\"\\t\"\$5\"\\t\"\$8}' temp/$out"."$inputNoPath.rawcnv | sed 's/^chr//' | sed 's/:/\\t/' | sed 's/-/\\t/' | sed 's/\\tstate.,cn=/\\t/' | sed 's/\\tconf=/\\t/' | sort -k5,5 > temp/$out"."$inputNoPath.rawcnv2";$o=`$c`;
#Convert rawcnv to vcf to allow same plink stats to be run (since .cnv input not supported in plink1.9 or plink2 yet)
    $c = "gzip -fk temp/$out" . "$inputNoPath.rawcnv 
sed 's/chr//' temp/$out"
      . "$inputNoPath.rawcnv | sed 's/:/\t/' | sed 's/-/\t/' | sed 's/\ /\t/' | gzip > temp/$out"
      . "$inputNoPath.rawcnv.bed.gz
RAWCNV_GZ_INPUT=temp/$out" . "$inputNoPath.rawcnv.gz
RAWCNV_BED_GZ_INPUT=temp/$out" . "$inputNoPath.rawcnv.bed.gz
zcat \$RAWCNV_GZ_INPUT > temp/$out" . "$inputNoPath.Input1;
zcat \$RAWCNV_BED_GZ_INPUT > \$RAWCNV_BED_GZ_INPUT.txt;
sort -C -k1,1 -k2,2n \$RAWCNV_BED_GZ_INPUT.txt;
if [ \$? == 0 ] ; then zcat \$RAWCNV_BED_GZ_INPUT > temp/$out" . "$inputNoPath.Input2; else sort -k1,1 -k2,2n -T temp \$RAWCNV_BED_GZ_INPUT.txt > temp/$out" . "$inputNoPath.Input2; fi;
awk '{if(NR==FNR){if(!(_[\$5])){_[\$5]=1}}\
      else{if(FNR==1){printf \"\#\#fileformat=VCFv4.3\\n\#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\";for(key in _){printf \"\\t\"key};print\"\"}\
CN=substr(\$6,11,length(\$6)-10)+0;\
if(CN==0){CN=1}\
if(CN==4){CN=3}\
if(lastChrPosAlt&&\$1\":\"\$2\":\"CN!=lastChrPosAlt)\
        {printf CHR_POS_9Cols;for(ID in _){ if(!(my_gt[ID])){printf \"\\t0/0:.\"}else{printf \"\\t\"my_gt[ID]}}print\"\";CHR_POS_9Cols=\"\";split(\"\", my_gt)};\
if(\$10!=\"\"){gsub(/conf=/,\"\",\$10)}else{\$10=\".\"};\
if(CN==0){svtype=\"DEL\";\
NegOrPosLen=\"-\";\
gt=\"1/1\"}else if(CN==1){svtype=\"DEL\";\
NegOrPosLen=\"-\";\
gt=\"0/1\"}else if(CN==2){svtype=\"DIP\";\
NegOrPosLen=\"\";\
gt=\"0/0\"}else{svtype=\"DUP\";\
NegOrPosLen=\"\";\
if(CN==3){gt=\"0/1\"}else{gt=\"1/1\"}};\
if(!(CHR_POS_9Cols)){CHR_POS_9Cols=\$1\"\\t\"\$2\"\\t\"\$1\":\"\$2\"-\"\$3\"\\tN\\t<CN=\"CN\">\\t.\\tPASS\\tSVTYPE=\"svtype\";END=\"\$3\";LEN=\"NegOrPosLen\$3-\$2+1\"\\tGT:GQ\"};\
my_gt[\$7]=gt\":\"\$10;\
lastChrPosAlt=\$1\":\"\$2\":\"CN\
}\
}END{printf CHR_POS_9Cols;for(ID in _){ if(!(my_gt[ID])){printf \"\\t0/0:.\"}else{printf \"\\t\"my_gt[ID]}}print\"\"}' temp/$out"."$inputNoPath.Input1 temp/$out"."$inputNoPath.Input2 > temp/$out" . "$inputNoPath.vcf; gzip -f temp/$out" . "$inputNoPath.vcf";
#<(zcat \$RAWCNV_GZ_INPUT) <(sort -C -k1,1 -k2,2n <(awk '{print \$1\"\\t\"\$2}' <(zcat \$RAWCNV_BED_GZ_INPUT)); if [ \$? == 0 ] ; then zcat \$RAWCNV_BED_GZ_INPUT; else mkdir temp; sort -k1,1 -k2,2n -T temp <(zcat \$RAWCNV_BED_GZ_INPUT); fi) | gzip > \$RAWCNV_BED_GZ_INPUT.vcf.gz";
print "Running $c";$o=`$c`;#print $o."\n";
#$o=`which mv`;print "which mv=".$o."\n";

$c = "zcat temp/$out" . "$inputNoPath.vcf.gz | grep ^#CHROM | head -1 | cut -f 10- | sed 's/\\t/\\n/g' > temp/$out"
  . "AllVcfIDs.txt";
$o = `$c`;
$c = "echo -e 'FID\\tIID\\tAffected' > temp/$out" . "file.pheno";
$o = `$c`;
$c =
    "cat $cases temp/$out"
  . "AllVcfIDs.txt | sort | uniq -c | awk '{count=\$1;printf \"0\\t\";for(i=2;i<=NF;i++){printf \$i}print\"\\t\"count}' >> temp/$out"
  . "file.pheno";
$o = `$c`;
}

##DONE sub function defintions

if ( $input =~ /\.txt$/ ) {

    #A text listfile with .rawcnv, vcf, or vcf.gz files listed
    $INPUT = $input
      ; #Use $INPUT as the .txt with a list of .rawcnv, vcf, or vcf.gz and $input as the individual file being processed at a time to match single .rawcnv, vcf, or vcf.gz input
    $c = "rm temp/$out" . "$inputNoPath.rawcnv2";
    $o = `$c`;
    open( LISTFILE, $INPUT );
	open( MergeVcfList, ">temp/$out" . $inputNoPath . "MergeVcfList.txt" );
    while ( $input = <LISTFILE> ) {
        chomp($input);
        if ( $input =~ /\.vcf$/ || $input =~ /\.vcf\.gz$/ ) {
				    if ( $input =~ /\.vcf\.gz$/ ) {
						$cat_zcat = "zcat";
					}
					else {
						$cat_zcat = "cat";
					}
				$sub_out = vcf_check_build($input);
				print $sub_out;
				$sub_out = vcf_extract_del_dup_IDs($input);
				print $sub_out;
				$sub_out = vcf_to_pfile($input);
				$sub_out = pvar_chr_pos_link;
				$sub_out = extract_info_vcf;
				$sub_out = cn_extract;
				#$sub_out = make_map_sorted;
				
				
				if ( $input =~ /\.vcf$/ )
				{
					$c="bgzip -c $input > $input temp/$out" . $inputNoPath . ".gz";$o=`$c`;
					print MergeVcfList "temp/$out" . $inputNoPath . ".gz\n";
					
				}
				else
				{
					print MergeVcfList "$input\n";
				}
				$c=$MyDirectoryPathPrefix . "PerlModules/" . "tabix -p vcf $input temp/$out" . $inputNoPath . ".rawcnv.bed.gz" . ".vcf.gz";
				#temp/$out" . $inputNoPath . ".rawcnv.bed.gz" . ".vcf.gz";
		}
        elsif ( $input =~ /\.rawcnv$/ ) {
		
		#.rawcnv output from PennCNV detect_cnv.pl
		#chr1:1-3  numsnp=3 length=3  state1,cn=0 1_A.baflrr startsnp=rs1 endsnp=rs3 conf=1
		#Prior to 7/1/22 inputs were all converted to a minimal set of fields all CNV detection algorithms should be expected to provide
		#rawcnv2 format
		#CHR START STOP CopyNumber SAMPLE CONFIDENCE/QualityScore
		#Now inputs are all converted to vcf since that is now a popular standard format and allows interaction with more tools moving forward
		
		#Checking lines match the rawcnv format since users may open rawcnv in excel to edit and don't delimit it properly by only spaces
            print "IN RAWCNV PARSING CODE\n";
	    if ($qc) {
		$sub_out = run_qc;
	    }
			$sub_out = verify_rawcnv;
			$sub_out = rawcnv_to_vcf;
			
			print MergeVcfList "temp/$out" . $inputNoPath . ".rawcnv.bed.gz" . ".vcf.gz";
			$c=$MyDirectoryPathPrefix . "PerlModules/" . "tabix -p vcf merge.$i.vcf.gz temp/$out" . $inputNoPath . ".rawcnv.bed.gz" . ".vcf.gz";
        }
        else {
            print
"ERROR: --i must have a .vcf, .vcf.gz, or .rawcnv extension in the .txt list file for $input\n";
        }

    }
	
	#Combine all files in the list file input converted to vcf
	$c=$MyDirectoryPathPrefix . "PerlModules/" . "bcftools merge -l temp/$out" . $inputNoPath . "MergeVcfList.txt -Oz -o temp/$out" . $inputNoPath . "merge.vcf.gz";$o=`$c`;print "$o\n";
	$input = "temp/$out" . $inputNoPath . "merge.vcf.gz";
	##Just let the list file merged vcf run through the input command line being vcf
}
##END of LIST input


if ( $input =~ /\.vcf$/ || $input =~ /\.vcf\.gz$/ ) {
    if ( $input =~ /\.vcf\.gz$/ ) {
        $cat_zcat = "zcat";
    }
    else {
        $cat_zcat = "cat";
    }
				$sub_out = vcf_check_build($input);
				print $sub_out;
				$sub_out = vcf_extract_del_dup_IDs($input);
				print $sub_out;
				$sub_out = vcf_to_pfile($input);
				$sub_out = pvar_chr_pos_link;
				$sub_out = extract_info_vcf;
				$sub_out = cn_extract;
				#$sub_out = make_map_sorted;
}

elsif ( $input =~ /\.rawcnv$/ ) {
		$original_input = $input;	
	        print "IN RAWCNV PARSING CODE\n";
			$sub_out = verify_rawcnv;
			$sub_out = rawcnv_to_vcf;
			$cat_zcat = "zcat";
			$input="temp/$out"."$inputNoPath".".vcf.gz";
			$sub_out = vcf_extract_del_dup_IDs($input="temp/$out"."$inputNoPath".".vcf.gz");
			$sub_out = vcf_to_pfile($input="temp/$out"."$inputNoPath".".vcf.gz");
			$sub_out = pvar_chr_pos_link($input="temp/$out"."$inputNoPath".".vcf.gz");
                        $sub_out = extract_info_vcf($input="temp/$out"."$inputNoPath".".vcf.gz");
                        $sub_out = cn_extract($input="temp/$out"."$inputNoPath".".vcf.gz");
			#$sub_out = make_map_sorted;

}

else {
    print "ERROR: --i must have a .vcf, .vcf.gz or .rawcnv extension\n";
	exit;
}

$resfile = $out . "ParseCNV.log";
open( LOG, ">$resfile" );

%MonthPastStart = (
    "0", "Jan", "1", "Feb", "2",  "Mar", "3",  "Apr",
    "4", "May", "5", "Jun", "6",  "Jul", "7",  "Aug",
    "8", "Sep", "9", "Oct", "10", "Nov", "11", "Dec"
);
@timeData = localtime(time);
$Year     = $timeData[5] + 1900;
if ( $timeData[2] > 12 ) {
    $timeData[2] -= 12;
    $AmOrPm = "PM";
}
elsif ( $timeData[2] eq 12 ) {
    $AmOrPm = "PM";
}
else {
    $AmOrPm = "AM";
}
if ( $timeData[2] < 10 ) { $timeData[2] = "0" . $timeData[2]; }
if ( $timeData[1] < 10 ) { $timeData[1] = "0" . $timeData[1]; }
if ( $timeData[0] < 10 ) { $timeData[0] = "0" . $timeData[0]; }
print "Run Started "
  . $timeData[2] . ":"
  . $timeData[1] . ":"
  . $timeData[0]
  . " $AmOrPm on "
  . $MonthPastStart{ $timeData[4] } . " "
  . $timeData[3] . " "
  . $Year . "\n";
print LOG "Run Started "
  . $timeData[2] . ":"
  . $timeData[1] . ":"
  . $timeData[0]
  . " $AmOrPm on "
  . $MonthPastStart{ $timeData[4] } . " "
  . $timeData[3] . " "
  . $Year . "\n";
print $myCommandLine;
print LOG $myCommandLine;

#if ($qc) {
sub run_qc {
    $c = "perl InputUtilities/ParseCNV_QC.pl";
    if ( $original_input =~ /.rawcnv/ ) {
        $c .= " --rawcnv $original_input"; print "$c\n";
    }
	elsif ( $input =~ /\.vcf$/ || $input =~ /\.vcf\.gz$/ ) {
		if ( $input =~ /\.vcf\.gz$/ ) {
			$cat_zcat = "zcat";
		}
		else {
			$cat_zcat = "cat";
		}
				$sub_out = vcf_check_build($input);
				print $sub_out;
				$sub_out = vcf_extract_del_dup_IDs($input);
				print $sub_out;
				$sub_out = vcf_to_pfile($input);
				$sub_out = pvar_chr_pos_link;
				$sub_out = extract_info_vcf;
				$sub_out = cn_extract;
				#$sub_out = make_map_sorted;
				
				$QC_input_original=$input;
				$input = "temp/$out" . $inputNoPath . ".rawcnv";
				$c .= " --rawcnv $input";
	}
    if ($log) {
        $c .= " --log $log";
    }
    if ($bfile) {
		#Could see if these commands work with plink2 to improve runtime
        $cp=$MyDirectoryPathPrefix."PerlModules/./plink --allow-extra-chr --bfile $bfile --missing";$o=`$cp`;
        $cp = $MyDirectoryPathPrefix."PerlModules/./plink --allow-extra-chr --bfile $bfile --pca";$o  = `$cp`;
		$cp=$MyDirectoryPathPrefix."PerlModules/./plink --allow-extra-chr --bfile $bfile --genome --min 0.4";$o=`$cp`;#--min 0.4 to save disk space
        $cp = $MyDirectoryPathPrefix."PerlModules/./plink --allow-extra-chr --bfile $bfile --check-sex";$o=`$cp`;
		$c.=" --callrate plink.imiss --popstrat plink.eigenvec --related plink.genome";
    }
    print "QC Command=".$c."\n";
    $o = `$c`;
    print "QC Command Output=".$o."\n";
    print "$o\n";
    #print "1\n";
    $c =
"awk '{print \$1\"\\t\"\$4\"\\t\"\$5}' Cases.rawcnv_remove_QC_RemoveIDs.txt_wMaceQualityScore_remove_QC_RemoveCalls.txt_Indexes | sed 's/^chr//' | sed 's/:/\\t/' | sed 's/-/\\t/' | sed 's/\\tstate.,cn=/\\t/' | sort -k5,5 > temp/$out"
      . "$inputNoPath.rawcnv2";
    $o = `$c`;

    $c = "cut -f 1 Cases.rawcnv_remove_QC_RemoveIDs.txt_wMaceQualityScore_remove_QC_RemoveCalls.txt_Indexes > $out"."$inputNoPath.QC_Filtered.rawcnv";
	$o=`$c`;
	##Reassign input variable since QC_Filtered rawcnv version should now be used
	$input=$out."$inputNoPath.QC_Filtered.rawcnv";
	$inputNoPath = $input;
	$inputNoPath =~ s/.*\///;    #get rid of potential path prefix
}
#print "2\n";

#$c="awk '{print \$1\"_\"\$2\"\\n\"\$1\"_\"\$3}' temp/$out"."$inputNoPath.rawcnv2 | sort -u -k1,1 -k2,2n -t_ > temp/$out"."map";$o=`$c`;
%h_state      = ();
%chr_posIndex = ();
open( MAP, "temp/$out" . ".map" );
while ( $line = <MAP> ) {
    chomp($line);
    $h_state{$line}      = "-9";
    $chr_posIndex{$line} = $lineNum;

    #$h_id{$line}="";
    @chr_pos        = split( "_", $line );
    $diff[$lineNum] = $chr_pos[1] - $last_pos;
    $last_pos       = $chr_pos[1];
    $lineNum++;
}
$diff[$lineNum] = 1;    #Exit at end
#print "3\n";

#$c="awk '{print \$5}' temp/$out"."$inputNoPath.rawcnv2 | sort -u > temp/$out"."AllIDs.txt";$o=`$c`;

if ( -e "plink.sexcheck" ) {
    #print "4\n";
    $c =
"awk '{if(NR>1)print \$2\"\\t\"\$4}' plink.sexcheck | sort > plink.sexcheck2; awk '{print \$5}' temp/$out"
      . "$inputNoPath.rawcnv2 | sort -u > filteredSamples; join plink.sexcheck2 filteredSamples | awk '{print \$2}' > plink.sexcheck3";
    $o = `$c`;
    #print "5\n";
    $c =
        "awk '{print \$5}' temp/$out"
      . "$inputNoPath.rawcnv2 | sort -u | cat $cases - | sort | uniq -c | paste - plink.sexcheck3 | awk '{print \"0 \"\$2\" 0 0 \"\$3\" \"\$1}' > temp/$out"
      . "fam";
    $o = `$c`;
}
    if ( $cases ne "" ) {
		#Not needed since now VCF is the primary format used instead of rawcnv2
        #print "6\n";
#$c="awk '{print \$5}' temp/$out"."$inputNoPath.rawcnv2 | sort -u | fgrep -vf - $cases > temp/$out"."casesNotInCNV.txt";$o=`$c`;
#$c="wc -l < temp/$out"."casesNotInCNV.txt";$o=`$c`;chomp($o);
#if($o > 0)
#{print "WARNING: $o specified cases are not present in the CNV file as listed in temp/$out"."casesNotInCNV.txt.\n";}
        #print "7\n";
#$c="awk '{print \$5}' temp/$out"."$inputNoPath.rawcnv2 | sort -u >temp/$out"."uniq_ids";$o=`$c`;
#$c="sed 's/\r//' $cases | sort -u | fgrep -vwf temp/$out"."casesNotInCNV.txt > temp/$out"."cases_NoCR_NoDup_Exists";$o=`$c`;
#$c="cat temp/$out"."cases_NoCR_NoDup_Exists temp/$out"."uniq_ids | sort | uniq -c | awk '{print \"0 \"\$2\" 0 0 0 \"\$1}' | fgrep -wf temp/$out"."uniq_ids > temp/$out"."fam";$o=`$c`;
#$c="awk '{print \$6}' temp/$out"."fam | sort -nr | uniq -c | awk '{if(\$2==2){printf \$1\" cases, \"}else{printf \$1\" controls\"}}'";$o=`$c`;print"$o\n";print LOG "$o\n";

    }
    if ( $quantitativeTrait ne "" ) {
        $c =
            "awk '{print \$1}' $quantitativeTrait > temp/$quantitativeTrait"
          . "_IDs; awk '{print \$5}' temp/$out"
          . "$inputNoPath.rawcnv2 | sort -u | fgrep -vf - temp/$quantitativeTrait"
          . "_IDs > temp/$out"
          . "casesNotInCNV.txt";
        $o = `$c`;
        $c = "wc -l < temp/$out" . "casesNotInCNV.txt";
        $o = `$c`;
        chomp($o);
        if ( $o > 0 ) {
            print
"WARNING: $o specified samples are not present in the CNV file as listed in temp/$out"
              . "casesNotInCNV.txt.\n";
        }

	$c = "cp temp/$out" . "AllVcfIDs.txt temp/$out" . "uniqIDs";$o = `$c`;
	#$c =
          #"awk '{print \$5}' temp/$out"
	 #. "$inputNoPath.rawcnv2 | sort -u >temp/$out"
	 #. "uniqIDs";
	  #$o = `$c`;
        $c =
            "sed 's/\r//' temp/$quantitativeTrait"
          . "_IDs | sort -u | fgrep -vwf temp/$out"
          . "casesNotInCNV.txt > temp/$out"
          . "cases_NoCR_NoDup_Exists";
        $o = `$c`;
        $c =
            "sort $quantitativeTrait >temp/$out"
          . "qt.sort; join temp/$out"
          . "cases_NoCR_NoDup_Exists temp/$out"
          . "qt.sort | awk '{print \"0 \"\$1\" 0 0 0 \"\$2}' > temp/$out"
          . "fam";
        $o       = `$c`;
        $c       = "wc -l temp/$out" . "uniqIDs | awk '{ORS=\"\";print \$1}'";
        $samples = `$c`;
        $c       = "awk '{ORS=\"\";if(NR==1)print NF-1}' $quantitativeTrait";
        $traits  = `$c`;
        $c       = "wc -l temp/$out" . "fam | awk '{ORS=\"\";print \$1}'";
        $samplesWTrait = `$c`;
        print
"$samples samples, $traits traits, and $samplesWTrait samples with trait\n";
    }

    if ( $tdt ) {
        if ( $bfile ne "" ) {
                $c = "cat temp/$out" . "AllVcfIDs.txt | sort > temp/$out" . "uniq_ids";
        	$o = `$c`;
	    $c =
                "join temp/$out"
              . "uniq_ids $bfile"
              . ".fam -1 1 -2 2 -o2.1,2.2,2.3,2.4,2.5,2.6 | sort -k2,2 > temp/$out"
              . "fam";
            $o  = `$c`;print "join=".$c."\n";
            $c  = "wc -l < temp/$out" . "uniq_ids";
            $o1 = `$c`;
            chomp($o1);
            $c  = "wc -l < $bfile" . ".fam";
            $o2 = `$c`;
            chomp($o2);

            if ( $o2 > $o1 ) {
                print "WARNING: "
                  . ( $o2 - $o1 )
                  . " fam specified cases are not present in the CNV file.\n";
            }
        }
        else {
            print
"ERROR: bfile must be specified with tdt option to provide family relationships\n";
            exit;
        }
    }
    #}

#Major shift from converting all input call formats to VCF rather than RAWCNV2
#Plink2 PFILE rather than BFILE
#temp/$out"."$inputNoPath
#Split Plink2 Generated PFILE into DEL, DUP, and DELDUP
#(base) [glessner@reslnvvhpc094 ParseCNV2]$ grep -v ^# temp/1KG_SudmantEtAl_wRvTestsALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.pvar | awk -F"\t" '{if($5~/,/)print$3}' > ExcludeMultiallelic
#(base) [glessner@reslnvvhpc094 ParseCNV2]$ PerlModules/./plink2 --pfile temp/1KG_SudmantEtAl_wRvTestsALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf --make-bed --out temp/outdel --exclude ExcludeMultiallelic

$c =
    "grep -v ^# temp/$out"
  . "$inputNoPath.pvar | awk -F\"\t\" '{if(\$5~/,/)print\$3}' > temp/$out"
  . "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath --make-bed --allow-extra-chr --out temp/$out"
  . " --exclude temp/$out"
  . "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath.del --make-bed --allow-extra-chr --out temp/$out"
  . "del --exclude temp/$out"
  . "ExcludeMultiallelic";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink2 --pfile temp/$out"
  . "$inputNoPath.dup --make-bed --allow-extra-chr --out temp/$out"
  . "dup --exclude temp/$out"
  . "ExcludeMultiallelic";
$o = `$c`;
$c = "head -1 temp/$out" . "del.bim";
$o = `$c`;
print $o. "\n";
$c = "head -1 temp/$out" . "del.fam";
$o = `$c`;
print $o. "\n";

=pod
$c = "mv temp/$out" . "del.bim temp/$out" . "del.bim.bak";
$o = `$c`;
$c =
"awk '{print \$1\"\\t\"\$1\"_\"\$4\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' temp/$out"
  . "del.bim.bak > temp/$out"
  . "del.bim";
$o = `$c`;

$c = "mv temp/$out" . "dup.bim temp/$out" . "dup.bim.bak";
$o = `$c`;
$c =
"awk '{print \$1\"\\t\"\$1\"_\"\$4\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' temp/$out"
  . "dup.bim.bak > temp/$out"
  . "dup.bim";
$o = `$c`;

$c = "mv temp/$out" . ".bim temp/$out" . ".bim.bak";
$o = `$c`;
$c =
"awk '{print \$1\"\\t\"\$1\"_\"\$4\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6}' temp/$out"
  . ".bim.bak > temp/$out"
  . ".bim";
$o = `$c`;
=cut
#$c="mv temp/$out"."del.bim temp/$out"."del.bim.bak2";$o=`$c`;
$c = "head -1 temp/$out" . "del.bim";
$o = `$c`;
print "Step1\n" . $o . "\n";
$c = "head -1 temp/$out" . "del.fam";
$o = `$c`;
print $o. "\n";
if   ( $input =~ /\.vcf\.gz$/ ) { $c = "zcat $input | "; }
else                            { $c = "cat $input | "; }
$c .= "grep ^#CHROM | head -1 | cut -f 10- | sed 's/\\t/\\n/g' > temp/$out"
  . "AllVcfIDs.txt";
$o = `$c`;

#print "Running:".$c."\n";$o=`$c`;print "Output:".$o."\n";
$c = "echo -e 'FID\\tIID\\tAffected' > temp/$out" . "file.pheno";
$o = `$c`;
$c =
    "cat $cases temp/$out"
  . "AllVcfIDs.txt | sort | uniq -c | awk '{count=\$1;printf \"0\\t\";for(i=2;i<=NF;i++){printf \$i}print\"\\t\"count}' >> temp/$out"
  . "file.pheno";
$o = `$c`;

open( PHENO, "temp/$out" . "file.pheno" );
my %IID_to_CaseControl;
while (<PHENO>) {
    chomp;
    @array = split( /\s+/, $_ );

    #print "@array\n$array[1]\n";
    $array[1] = "" . $array[1];
    $IID_to_CaseControl{ $array[1] } = $array[2];
}

$c = "head -1 temp/$out" . "del.bim";
$o = `$c`;                              #print "Step2\n".$o."\n";
$c = "head -1 temp/$out" . "del.fam";
$o = `$c`;                              #print $o."\n";
$c = "mv temp/$out" . ".fam temp/$out" . ".fam.bak";
$o = `$c`;
open( FAM_from_plink2,              "temp/$out" . ".fam.bak" );
open( FAM_from_plink2_wCaseControl, ">temp/$out" . ".fam" );

while (<FAM_from_plink2>) {
    chomp;
    @array = split( /\s+/, $_ );

    $array[1] = "" . $array[1];
    $IID = $array[1];

    print FAM_from_plink2_wCaseControl $array[0] . "\t"
      . $array[1] . "\t"
      . $array[2] . "\t"
      . $array[3] . "\t"
      . $array[4] . "\t"
      . $IID_to_CaseControl{ $array[1] } . "\n";
}
close(FAM_from_plink2);
close(FAM_from_plink2_wCaseControl);

$c = "mv temp/$out" . "del.fam temp/$out" . "del.fam.bak";
$o = `$c`;
open( FAM_from_plink2,              "temp/$out" . "del.fam.bak" );
open( FAM_from_plink2_wCaseControl, ">temp/$out" . "del.fam" );

while (<FAM_from_plink2>) {
    chomp;
    @array = split( /\s+/, $_ );

    $array[1] = "" . $array[1];
    $IID = $array[1];

    print FAM_from_plink2_wCaseControl $array[0] . "\t"
      . $array[1] . "\t"
      . $array[2] . "\t"
      . $array[3] . "\t"
      . $array[4] . "\t"
      . $IID_to_CaseControl{ $array[1] } . "\n";
}
close(FAM_from_plink2);
close(FAM_from_plink2_wCaseControl);

$c = "mv temp/$out" . "dup.fam temp/$out" . "dup.fam.bak";
$o = `$c`;
open( FAM_from_plink2,              "temp/$out" . "dup.fam.bak" );
open( FAM_from_plink2_wCaseControl, ">temp/$out" . "dup.fam" );

while (<FAM_from_plink2>) {
    chomp;
    @array = split( /\s+/, $_ );

    $array[1] = "" . $array[1];
    $IID = $array[1];

    print FAM_from_plink2_wCaseControl $array[0] . "\t"
      . $array[1] . "\t"
      . $array[2] . "\t"
      . $array[3] . "\t"
      . $array[4] . "\t"
      . $IID_to_CaseControl{ $array[1] } . "\n";
}
close(FAM_from_plink2);
close(FAM_from_plink2_wCaseControl);

#exit;
$c = "head -1 temp/$out" . "del.bim";
$o = `$c`;                              #print "Step2\n".$o."\n";
$c = "head -1 temp/$out" . "del.fam";
$o = `$c`;                              #print $o."\n";



#average 2 3  5 6 to get POS and END since CIPOS used in bedpe
open( OUT_vcfToBedpe, "temp/$out" . "$inputNoPath" . "_wEND_vcfToBedpe" );
open( OUT_BedpeToMap, ">temp/$out" . ".map_presorting" );
%CHR_POS_to_END;
while (<OUT_vcfToBedpe>) {

    #Remember Perl Arrays are 0 based not 1 based as in AWK
    chomp;
    @_ = split( /\t/, $_ );
    my $rounded_start_pos = int( ( $_[1] + $_[2] ) / 2 + 0.5 );
    my $rounded_end_pos   = int( ( $_[4] + $_[5] ) / 2 + 0.5 );
    $CHR_POS = $_[0] . "_" . $rounded_start_pos;
    $CHR_POS_to_END{$CHR_POS} = $rounded_end_pos;
    print OUT_BedpeToMap $_[0] . "\t"
      . $_[0] . "_"
      . $rounded_start_pos . "\t0\t"
      . $rounded_start_pos . "\n";
    print OUT_BedpeToMap $_[0] . "\t"
      . $_[0] . "_"
      . $rounded_end_pos . "\t0\t"
      . $rounded_end_pos . "\n";
}
close(OUT_BedpeToMap);
$c =
  "sort -u -k1,1 -k2,2n -t_ temp/$out" . ".map_presorting > temp/$out" . ".map";
$o = `$c`;

%h_state      = ();
%chr_posIndex = ();
open( MAP, "temp/$out" . ".map" );
while ( $line = <MAP> ) {
    chomp($line);
    $h_state{$line}      = "-9";
    $chr_posIndex{$line} = $lineNum;
    @chr_pos             = split( "_", $line );
    $diff[$lineNum]      = $chr_pos[1] - $last_pos;
    $last_pos            = $chr_pos[1];
    $lineNum++;
}
$diff[$lineNum] = 1;    #Exit at end
close(MAP);

if   ( $input =~ /\.vcf\.gz$/ ) { $c = "zcat $input | "; }
else                            { $c = "cat $input | "; }
$c .= "grep ^#CHROM | head -1 | cut -f 10- | sed 's/\\t/\\n/g' > temp/$out"
  . "AllVcfIDs.txt";
$o = `$c`;
$c = "echo -e 'FID\\tIID\\tAffected' > temp/$out" . "file.pheno";
$o = `$c`;
$c =
    "cat $cases temp/$out"
  . "AllVcfIDs.txt | sort | uniq -c | awk '{count=\$1;printf \"0\\t\";for(i=2;i<=NF;i++){printf \$i}print\"\\t\"count}' >> temp/$out"
  . "file.pheno";
$o = `$c`;

$c = "head -1 temp/$out" . "file.pheno";
$o = `$c`;                                 #print "Pheno ".$o."\n";

$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink --bfile temp/$out"
  . "del --assoc fisher --pheno temp/$out"
  . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out" . "del";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink --bfile temp/$out"
  . "dup --assoc fisher --pheno temp/$out"
  . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out" . "dup";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/./plink --bfile temp/$out"
  . " --assoc fisher --pheno temp/$out"
  . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out";
$o = `$c`;

#CALCULATE Variant ID based rare recurrent CNV association statistics
if ( $cases ne "" ) {

    if ( $statistic eq "" || $statistic eq "fisher") {
	print "stat=" . $statistic . "\n";
        #Get matching p-values to original ParseCNV with OR from assoc fisher
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink --bfile temp/$out"
          . "del --model fisher --pheno temp/$out"
          . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out" . "del";
        $o = `$c`;
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink --bfile temp/$out"
          . "dup --model fisher --pheno temp/$out"
          . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out" . "dup";
        $o = `$c`;
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink --bfile temp/$out"
          . " --model fisher --pheno temp/$out"
          . "file.pheno --allow-no-sex --allow-extra-chr $FrequencyFilter --out temp/$out";
        $o = `$c`;

        $c =
            "awk '{print \$2\"\\t\"\$0}' temp/$out"
          . "del.assoc.fisher | sort -k1,1 > temp/$out"
          . "del.assoc.fisher.forJoin";
        $o = `$c`;
        $c =
"awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"
          . "del.model | sort -k1,1 > temp/$out"
          . "del.model.forJoin";
        $o = `$c`;
        $c =
          "echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"
          . "del.assoc.fisher.model";
        $o = `$c`;
        $c =
            "join temp/$out"
          . "del.assoc.fisher.forJoin temp/$out"
          . "del.model.forJoin | sort -k2,2 -k4,4n >> temp/$out"
          . "del.assoc.fisher.model";
        $o = `$c`;
        $c =
            "awk '{print \$2\"\\t\"\$0}' temp/$out"
          . "dup.assoc.fisher | sort -k1,1 > temp/$out"
          . "dup.assoc.fisher.forJoin";
        $o = `$c`;
        $c =
"awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"
          . "dup.model | sort -k1,1 > temp/$out"
          . "dup.model.forJoin";
        $o = `$c`;
        $c =
          "echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"
          . "dup.assoc.fisher.model";
        $o = `$c`;
        $c =
            "join temp/$out"
          . "dup.assoc.fisher.forJoin temp/$out"
          . "dup.model.forJoin | sort -k2,2 -k4,4n >> temp/$out"
          . "dup.assoc.fisher.model";
        $o = `$c`;
        $c =
            "awk '{print \$2\"\\t\"\$0}' temp/$out"
          . ".assoc.fisher | sort -k1,1 > temp/$out"
          . ".assoc.fisher.forJoin";
        $o = `$c`;
        $c =
"awk '{if(\$5==\"GENO\")print \$2\"\\t\"\$NF\"\\t\"\$(NF-2)\"\\t\"\$(NF-1)}' temp/$out"
          . ".model | sort -k1,1 > temp/$out"
          . ".model.forJoin";
        $o = `$c`;
        $c =
          "echo SNP CHR SNP BP A1 F_A F_U A2 P_ASSOC OR P AFF UNAFF > temp/$out"
          . ".assoc.fisher.model";
        $o = `$c`;
        $c =
            "join temp/$out"
          . ".assoc.fisher.forJoin temp/$out"
          . ".model.forJoin | sort -k2,2 -k4,4n >> temp/$out"
          . ".assoc.fisher.model";
        $o = `$c`;
    }
    else {
        if ( $statistic eq "logistic" || $statistic eq "linear" ) {
		#Plink2 uses Generalized Linear Model (GLM) now
#Separate lines in glm association result for each of the multiple ALT alleles in multiallelic variants
#Plink2 does not have segment-based variants (extend POS to END, POS+LENGTH, or take difference of REF and ALT lengths) or “cnv” function as in the plink1.07
#"Rare Recurrent" filter with --mac 2 --max-maf 0.01

            if ( $covar ne "" ) {
                $glm_covar = "--glm --covar $covar";
            }
            else {
                $glm_covar = "--glm allow-no-covars";
            }
            if ( $input =~ /\.vcf$/ || $input =~ /\.vcf\.gz$/ || $input =~ /\.rawcnv$/) {

                #$cases
                if   ( $input =~ /\.vcf\.gz$/ ) { $c = "zcat $input | "; }
                elsif( $input =~ /\.vcf$/ ) { $c = "cat $input | "; }
		else {$input="temp/$out"."$inputNoPath".".vcf.gz";$c = "zcat $input | ";}
                $c .=
"grep ^#CHROM | head -1 | cut -f 10- | sed 's/\\t/\\n/g' > temp/$out"
                  . "AllVcfIDs.txt";
                $o = `$c`;
                $c = "echo -e 'IID\\tAffected' > temp/$out" . "file.pheno";
                $o = `$c`;
                $c =
                    "cat $cases temp/$out"
                  . "AllVcfIDs.txt | sort | uniq -c | awk '{count=\$1;for(i=2;i<=NF;i++){printf \$i\" \"}print\"\\t\"count}' >> temp/$out"
                  . "file.pheno";
                $o = `$c`;

                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/./plink2 --pfile temp/$out"
                  . "$inputNoPath $glm_covar --pheno temp/$out"
                  . "file.pheno --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
                  . "GLM_wCovar";
                $o = `$c`;
                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/./plink2 --pfile temp/$out"
                  . "$inputNoPath.del $glm_covar --pheno temp/$out"
                  . "file.pheno --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
                  . "GLM_wCovar.del";
                $o = `$c`;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/./plink2 --pfile temp/$out"
                  . "$inputNoPath.dup $glm_covar --pheno temp/$out"
                  . "file.pheno --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
                  . "GLM_wCovar.dup";
                $o = `$c`;
            }
        }
	### END STAT = LOGISTIC
	else{ ##RVTESTS
=pod
#RvTests http://zhanxw.github.io/rvtests/ Integration of Single Variant tests: wald, score and Groupwise tests: burden cmc, variable threshold price, and kernel skat, kbac.
#Xiaowei Zhan, Youna Hu, Bingshan Li, Goncalo R. Abecasis, and Dajiang J. Liu. RVTESTS: An Efficient and Comprehensive Tool for Rare Variant Association Analysis Using Sequence Data. Bioinformatics 2016 32: 1423-1426. doi:10.1093/bioinformatics/btw079
#wget https://github.com/zhanxw/rvtests/releases/download/v2.1.0/rvtests_linux64.tar.gz
#/mnt/isilon/cag_ngs/hiseq/glessner/ParseCNV2_EJHG_Revisions/ParseCNV2/PerlModules/executable/rvtest
#TO DO
#https://zhanxw.github.io/rvtests/#specify-groups-eg-burden-unit
#WORKS PerlModules/executable/./rvtest --inVcf ../1KG/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf --pheno 1KG_SudmantEtAldel.ped.rvtest.pheno --out rvtest_wald_output --single wald

#Make rvtest specific and plink similar pheno .ped file (really a plink .fam with header and needing y1 for the Affected column
=cut

        $c =
            "cat $cases temp/$out"
          . "AllVcfIDs.txt | sort | uniq -c | awk 'BEGIN{print\"fid iid fatid matid sex y1\"}{count=\$1;ID=\$2;print \"0\",ID,\"0\\t0\\t0\\t\"count}' > temp/$out"
          . "del.fam.rvtest.pheno";
        print "Running:" . $c . "\n";
        $o = `$c`;
        print "Output:" . $o . "\n";
        print "stat=" . $statistic . "\n";

        %SingleVariant_RvTests = (
            "score"           => 1,
            "wald"            => 1,
            "exact"           => 1,
            "dominantExact"   => 1,
            "famLRT"          => 1,
            "famScore"        => 1,
            "famGrammarGamma" => 1,
            "firth"           => 1
        );
        %Burden_RvTests = (
            "cmc"        => 1,
            "cmc"        => 1,
            "zeggini"    => 1,
            "mb"         => 1,
            "fp"         => 1,
            "exactCMC"   => 1,
            "cmcWald"    => 1,
            "rarecover"  => 1,
            "cmat"       => 1,
            "famcmc"     => 1,
            "famzeggini" => 1
        );
        %VariableThreshold_RvTests = (
            "price"       => 1,
            "analytic"    => 1,
            "famAnalytic" => 1
        );
        %Kernel_RvTests = (
            "skat"    => 1,
            "skato"   => 1,
            "kbac"    => 1,
            "famSkat" => 1
        );
        %MetaAnalysis_RvTests = (
            "score"     => 1,
            "dominant"  => 1,
            "recessive" => 1,
            "cov"       => 1,
            "bolt"      => 1,
            "boltCov"   => 1
        );

        #Needed for some rvtest statistics functions
        if( $input =~ /\.vcf$/ ){$c = "bgzip -c $input > $input.gz";$o = `$c`;}
        if( $input =~ /\.vcf$/ ){$my_input=$input.gz;}
        else{$my_input=$input;}
	$c =
            "bgzip -c temp/$out"
          . "$inputNoPath"
          . "_del.vcf > temp/$out"
          . "$inputNoPath"
          . "_del.vcf.gz";
        $o = `$c`;
        $c =
            "bgzip -c temp/$out"
          . "$inputNoPath"
          . "_dup.vcf > temp/$out"
          . "$inputNoPath"
          . "_dup.vcf.gz";
        $o = `$c`;

        if( $input =~ /\.vcf$/){$c = $MyDirectoryPathPrefix . "PerlModules/" . "tabix -p vcf $input.gz";$o = `$c`;}
        $c = $MyDirectoryPathPrefix . "PerlModules/" . "tabix -p vcf temp/$out" . "$inputNoPath" . "_del.vcf.gz";
        $o = `$c`;
        $c = $MyDirectoryPathPrefix . "PerlModules/" . "tabix -p vcf temp/$out" . "$inputNoPath" . "_dup.vcf.gz";
        $o = `$c`;

        if($statistic =~ /fam/)
	{
		$c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./vcf2kinship --inVcf $my_input --bn --out temp/$out". "_vcf2kinship";
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";
	}
	print "statistic $statistic\n";
	if ( $SingleVariant_RvTests{$statistic} == 1 ) {
        print "In SingleVariant\n";
		
	     $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf $my_input --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1 --single $statistic";
            if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
	    print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";
            $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf temp/$out"
              . "$inputNoPath"
              . "_del.vcf.gz --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1.del --single $statistic";
	      if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";
            $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf temp/$out"
              . "$inputNoPath"
              . "_dup.vcf.gz --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1.dup --single $statistic";
	      if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";
        }
        else {
#rvtest --inVcf input.vcf --pheno phenotype.ped --out output --geneFile refFlat_hg19.txt.gz --burden cmc --vt price --kernel skat,kbac
            if ( $Burden_RvTests{$statistic} == 1 ) {
                $RvTestsStatFlag = "--burden";
            }
            elsif ( $VariableThreshold_RvTests{$statistic} == 1 ) {
                $RvTestsStatFlag = "--vt";
            }
            elsif ( $Kernel_RvTests{$statistic} == 1 ) {
                $RvTestsStatFlag = "--kernel";
            }
            elsif ( $MetaAnalysis_RvTests{$statistic} == 1 ) {
                if ($covar) {
                    $RvTestsStatFlag =
"--covar $covar --inverseNormal --useResidualAsPhenotype --meta";
                }
                else {
                    $RvTestsStatFlag =
                      "--inverseNormal --useResidualAsPhenotype --meta";
                }
            }
            else {
                print "No Matching RvTest Found for $statistic\n";
            }

#Meta-Analysis: rvtest --inVcf input.vcf --pheno phenotype.ped --covar example.covar --covar-name age,bmi --inverseNormal --useResidualAsPhenotype  --meta score,cov --out output

            #Separate DEL DUP and DELDUP files
            #Create DEL DUP filtered VCFs
            $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf $input.gz --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1 --geneFile $MyDirectoryPathPrefix"
              . "GeneRef/"
              . $build
              . "_refFlat.txt.gz $RvTestsStatFlag $statistic";
	      if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";

            $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf temp/$out"
              . "$inputNoPath"
              . "_del.vcf.gz --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1.del --geneFile $MyDirectoryPathPrefix"
              . "GeneRef/"
              . $build
              . "_refFlat.txt.gz $RvTestsStatFlag $statistic";
	      if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";
            $c =
                $MyDirectoryPathPrefix
              . "PerlModules/executable/./rvtest --inVcf temp/$out"
              . "$inputNoPath"
              . "_dup.vcf.gz --pheno temp/$out"
              . "del.fam.rvtest.pheno --out temp/$out"
              . "_rvtest_"
              . $statistic
              . "_output1.dup --geneFile $MyDirectoryPathPrefix"
              . "GeneRef/"
              . $build
              . "_refFlat.txt.gz $RvTestsStatFlag $statistic";
	      if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
            print "Running:" . $c . "\n";
            $o = `$c`;
            print "Output:" . $o . "\n";

            if ( $MetaAnalysis_RvTests{$statistic} == 1 ) {
                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/executable/./rvtest --inVcf $input --pheno temp/$out"
                  . "del.fam.rvtest.pheno --out temp/$out"
                  . "_rvtest_"
                  . $statistic
                  . "_output1 $RvTestsStatFlag $statistic";
		  if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
                print "Running:" . $c . "\n";
                $o = `$c`;
                print "Output:" . $o . "\n";
                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/executable/./rvtest --inVcf temp/$out"
                  . $inputNoPath
                  . "_del.vcf.gz --pheno temp/$out"
                  . "del.fam.rvtest.pheno --out temp/$out"
                  . "_rvtest_"
                  . $statistic
                  . "_output1.del $RvTestsStatFlag $statistic";
		  if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
                print "Running:" . $c . "\n";
                $o = `$c`;
                print "Output:" . $o . "\n";

                $c =
                    $MyDirectoryPathPrefix
                  . "PerlModules/executable/./rvtest --inVcf temp/$out"
                  . $inputNoPath
                  . "_dup.vcf.gz --pheno temp/$out"
                  . "del.fam.rvtest.pheno --out temp/$out"
                  . "_rvtest_"
                  . $statistic
                  . "_output1.dup $RvTestsStatFlag $statistic";
		  if($statistic =~ /fam/){$c.=" --kinship temp/$out". "_vcf2kinship.kinship";}
                print "Running:" . $c . "\n";
                $o = `$c`;
                print "Output:" . $o . "\n";
            }

        }

        $c =
"awk -F\"\\t\" 'BEGIN{OFS=\"\\t\";RvTestsHeaderToMyStandard[\"Test\"]=\"SNP\";RvTestsHeaderToMyStandard[\"CHROM\"]=\"CHR\";RvTestsHeaderToMyStandard[\"POS\"]=\"BP\";RvTestsHeaderToMyStandard[\"Beta\"]=\"BETA\";RvTestsHeaderToMyStandard[\"Pvalue\"]=\"P\";} {test=\"Test\";if(test in RvTestsHeaderToMyStandard) {if(NR==1){ count=split(\$0,a,/\\t/); for(i=1;i<=NF;i++){ myTranslation=\"\"; if(\$i in RvTestsHeaderToMyStandard){ myTranslation=RvTestsHeaderToMyStandard[\$i]; if(myTranslation==\"SNP\"){mySNPColWithColon=i};ALLmyTranslations++;ALLmyTranslationsFields=ALLmyTranslationsFields myTranslation\" \"}; if(i<NF){if(myTranslation!=\"\"){printf myTranslation\"\\t\"}else{printf \$i\"\\t\"}}else{if(myTranslation!=\"\"){print myTranslation}else{print \$NF};}}} else{ for(i=1;i<NF;i++){ if(i==mySNPColWithColon){ gsub(/:/,\"_\",\$i)} ;printf \$i\"\\t\"} ;if(NF==mySNPColWithColon){ gsub(/:/,\"_\",\$NF)} ;print \$NF}\}\}'";

        #P, SNP, BETA/T/OR, CHR, BP required fields

        #combined del dup
        $cmd =
            "ls temp/$out"
          . "_rvtest_"
          . $statistic
          . "_output1*[^.log] | grep -v output1.del | grep -v output1.dup | grep -v .kinship\$ | tr -d '\\n'";
        $myRvTestOutput1File = `$cmd`;
        $myRvTestOutput2File =
          "temp/$out" . "_rvtest_" . $statistic . "_output2";
        $o = `$c $myRvTestOutput1File > $myRvTestOutput2File`;
        print "COMMAND=" . $c . "\n" . "OUTPUT=" . $o . "\\n";

        #del
        $cmd =
            "ls temp/$out"
          . "_rvtest_"
          . $statistic
          . "_output1.del*[^.log] | grep -v .kinship\$ | tr -d '\\n'";
        $myRvTestOutput1File = `$cmd`;
        $myRvTestOutput2File =
          "temp/$out" . "_rvtest_" . $statistic . "_output2.del";
        $o = `$c $myRvTestOutput1File > $myRvTestOutput2File`;
        print "COMMAND=" . $c . "\\n" . "OUTPUT=" . $o . "\\n";

        #dup
        $cmd =
            "ls temp/$out"
          . "_rvtest_"
          . $statistic
          . "_output1.dup*[^.log] | grep -v .kinship\$ | tr -d '\\n'";
        $myRvTestOutput1File = `$cmd`;
        $myRvTestOutput2File =
          "temp/$out" . "_rvtest_" . $statistic . "_output2.dup";
        $o = `$c $myRvTestOutput1File > $myRvTestOutput2File`;
        print "COMMAND=" . $c . "\\n" . "OUTPUT=" . $o . "\\n";
    }    # end of else for statement: if statistic eq "" else
    }
    #Star since RvTests adds to the end of --out based on the test

    #Efficiently added del and dup runs

}

#../plink2 --vcf /mnt/isilon/cag_ngs/hiseq/glessner/IMPUTE_CNV/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.vcf.gz --make-pgen --out 1KGP3_VCF_GZ_to_PlinkFiles_viaPlink2
#../plink2 --pfile 1KGP3_VCF_GZ_to_PlinkFiles_viaPlink2 \
#--glm (allow-no-covars) \
#--covar 1KGP3_VCF_GZ_to_PlinkFiles_viaPlink2.psam.4cols.covar \
#--pheno 1KGP3_VCF_GZ_to_PlinkFiles_viaPlink2.psam.pheno \
#--out 1KGP3_VCF_GZ_to_PlinkFiles_viaPlink2_GLM_wCovar
#(--glm allow-no-covars not recommended, find age/eth/sex for covar input)
#Separate lines in glm association result for each of the multiple ALT alleles in multiallelic variants
#--glm Basic association analysis on quantitative and/or case/control phenotypes. For each variant, a linear (for quantitative traits) or logistic (for case/control) regression is run with the phenotype as the dependent variable, and nonmajor allele dosage(s) and a constant-1 column as predictors.

if ( $quantitativeTrait ne "" ) {
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . "del --assoc --pheno $quantitativeTrait --allow-no-sex --allow-extra-chr --out temp/$out"
      . "del";
    if ( $covar ne "" ) { $c = $c . " --covar $covar"; }
    ;
    $o = `$c`;
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . "dup --assoc --pheno $quantitativeTrait --allow-no-sex --allow-extra-chr --out temp/$out"
      . "dup";
    if ( $covar ne "" ) { $c = $c . " --covar $covar"; }
    ;
    $o = `$c`;
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . " --assoc --pheno $quantitativeTrait --allow-no-sex --allow-extra-chr --out temp/$out";
    if ( $covar ne "" ) { $c = $c . " --covar $covar"; }
    ;
    $o = `$c`;

    if ( $statistic eq "linear" ) {
        if ( $covar ne "" ) {
            $glm_covar = "--glm hide-covar --covar $covar";
        }
        else {
            $glm_covar = "--glm allow-no-covars";
        }

        #Already set $FrequencyFilter at beginning
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink2 --pfile temp/$out"
          . "$inputNoPath $glm_covar --pheno $quantitativeTrait --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
          . "GLM_wCovar";
        $o = `$c`;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
	##TO DO: Extract AffectedTrait text from .pheno or output .logistic.hybrid
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink2 --pfile temp/$out"
          . "$inputNoPath.del $glm_covar --pheno $quantitativeTrait --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
          . "GLM_wCovar.del";
        $o = `$c`;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
        $c =
            $MyDirectoryPathPrefix
          . "PerlModules/./plink2 --pfile temp/$out"
          . "$inputNoPath.dup $glm_covar --pheno $quantitativeTrait --adjust $FrequencyFilter --allow-extra-chr --out temp/$out"
          . "GLM_wCovar.dup";
        $o = `$c`;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
    
	
     }

}
if ( $tdt ) {
     $c = "cp temp/$out" . "fam temp/$out" . "del.fam";$o = `$c`;
     $c = "cp temp/$out" . "fam temp/$out" . "dup.fam";$o = `$c`;
     $c = "cp temp/$out" . "fam temp/$out" . ".fam";$o = `$c`;
     $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . "del --tdt --pheno temp/$out"
      . "file.pheno --allow-no-sex --allow-extra-chr --out temp/$out" . "del";
    $o = `$c`;
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . "dup --tdt --pheno temp/$out"
      . "file.pheno --allow-no-sex --allow-extra-chr --out temp/$out" . "dup";
    $o = `$c`;
    $c =
        $MyDirectoryPathPrefix
      . "PerlModules/./plink --bfile temp/$out"
      . " --tdt --pheno temp/$out"
      . "file.pheno --allow-no-sex --allow-extra-chr --out temp/$out";
    $o = `$c`;

#--tdt normally computes parenTDT, transmission disequilibrium test, and combined test statistics, writing results to plink.tdt.
#--dfam extended version of the sib-TDT (see Spielman RS, Ewens WJ (1998) A sibship test for linkage in the presence of association: the sib transmission/disequilibrium test) which includes clusters of unrelated individuals, writing results to plink.dfam.
#--qfam quantitative trait family

}
    $c         = "wc -l <temp/$out" . "$inputNoPath.pvar";
    $countSNPs = `$c`;
    chomp($countSNPs);
#$c="perl OutputUtilities/InsertPlinkPvalues.pl";$o=`$c`;
print("\nDone Loading Input Files\n");
print "$countSNPs probes with observed CNV breakpoints\n";
if ( not $mergePVar ) {
    $mergePVar = 1;
}
if ( not $mergeDist ) {
    $mergeDist = 1000000;
}
if ( not $maxPInclusion ) {
    $maxPInclusion = 0.05;
}

if ( not $lowHighRisk ) {
    $lowHighRisk = "high";

    #print "NOTICE: --lowHighRisk set to 'high' by default\n";
    #print LOG "NOTICE: --lowHighRisk set to 'high' by default\n";
}
if ( $cases ne "" ) {
    if ( $statistic eq "fisher" ) {
        open( $DEL, "temp/$out" . "del.assoc.fisher.model" );
        open( $DUP, "temp/$out" . "dup.assoc.fisher.model" );
    }
    if ( $statistic eq "logistic" || $statistic eq "linear") {
        open( $DEL,
            "temp/$out" . "GLM_wCovar.del.Affected.glm.logistic.hybrid" )
          ;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
        open( $DUP,
            "temp/$out" . "GLM_wCovar.dup.Affected.glm.logistic.hybrid" )
          ;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
	  #$line = <DEL>;
	  #print "temp/$out" . "GLM_wCovar.del.*.glm.logistic.hybrid\n";
	  #print "$line\n";
    }
    if (   $SingleVariant_RvTests{$statistic} == 1
        || $Burden_RvTests{$statistic} == 1
        || $VariableThreshold_RvTests{$statistic} == 1
        || $Kernel_RvTests{$statistic} == 1
        || $MetaAnalysis_RvTests{$statistic} == 1 )
    {
        open( $DEL, "temp/$out" . "_rvtest_" . $statistic . "_output2.del" );
        open( $DUP, "temp/$out" . "_rvtest_" . $statistic . "_output2.dup" );
    }
}
elsif ( $quantitativeTrait ne "" ) {
    if ( $statistic ne "linear" ) {
        open( $DEL, "temp/$out" . "del.qassoc" );
        open( $DUP, "temp/$out" . "dup.qassoc" );
    }
    else {
        open( $DEL, "temp/$out" . "del.GLM_wCovar.*.glm.linear" )
          ;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
        open( $DUP, "temp/$out" . "dup.GLM_wCovar.*.glm.linear" )
          ;    #--> GLM_wCovar.AffectedTrait.glm.logistic.hybrid
    }

}
if ( $tdt ) {
    open( $DEL, "temp/$out" . "del.tdt" );
    open( $DUP, "temp/$out" . "dup.tdt" );
}
$line = <$DEL>;
#print "temp/$out" . "GLM_wCovar.del.Affected.glm.logistic.hybrid\n";
#print "$line\n";
@Vals = split( /\s+/, $line );
$i    = 0;
print "Statistics Header Values=" . "@Vals\n";
while ( ( uc( $Vals[$i] ) ne 'P' ) && ( uc( $Vals[$i] ) ne 'PVALUE' ) && ( uc( $Vals[$i] ) ne 'P-VALUE' ) && ( uc( $Vals[$i] ) ne 'PVALUETWOSIDE' ) && $i <= $#Vals ) {
    $i++;
}
if ( uc( $Vals[$i] ) ne 'P' && ( uc( $Vals[$i] ) ne 'PVALUE' ) && ( uc( $Vals[$i] ) ne 'P-VALUE' ) && ( uc( $Vals[$i] ) ne 'PVALUETWOSIDE' ) ) {
    print "ERROR: No P column header.\n";
    exit;
}
$myPColNum = $i;

$i = 0;
while ( ( uc( $Vals[$i] ) ne 'CHR' ) && ( uc( $Vals[$i] ) ne '#CHROM' ) && ( $i <= $#Vals ) ) {
    $i++;
}
if ( uc( $Vals[$i] ) ne 'CHR' && uc( $Vals[$i] ) ne '#CHROM') {
    print "ERROR: No CHR column header.\n";
    exit;
}
$myChrColNum = $i;
$i           = 0;
while ( ( uc( $Vals[$i] ) ne 'BP' ) && ( uc( $Vals[$i] ) ne 'POS' ) && ( $i <= $#Vals ) ) {
    $i++;
}
if ( uc( $Vals[$i] ) ne 'BP' && uc( $Vals[$i] ) ne 'POS') {
    print "ERROR: No BP column header.\n";
    exit;
}
$myBpColNum = $i;

$i         = 0;
while ( ( uc( $Vals[$i] ) ne 'SNP') && ( uc( $Vals[$i] ) ne 'ID' ) && ( $i <= $#Vals ) ) {
    $i++;
}
if ( uc( $Vals[$i] ) ne 'SNP' && uc( $Vals[$i] ) ne 'ID') {
    print "No SNP or ID column header so attempt to add CHR_POS column as SNP ID\n";
	my $fd = fileno ${$DEL};
	$Stats_Filename = readlink("/proc/$$/fd/$fd");
	$c="awk '{print \$0\"\\t\"\$($myChrColNum+1)\"_\"\$($myBpColNum+1)}' $Stats_Filename > $Stats_Filename.SNPColAdded";
    	print "Run $c\n";
	$o=`$c`;
	close($DEL);
	open($DEL, "$Stats_Filename.SNPColAdded");
	$i=$#Vals+1;
	$line=<$DEL>;
	@Vals = split( /\s+/, $line );

	my $fd = fileno ${$DUP};
        $Stats_Filename = readlink("/proc/$$/fd/$fd");
        $c="awk '{print \$0\"\\t\"\$($myChrColNum+1)\"_\"\$($myBpColNum+1)}' $Stats_Filename > $Stats_Filename.SNPColAdded";
        $o=`$c`;
	close($DUP);
        open($DUP, "$Stats_Filename.SNPColAdded");
	#$line=<$DUP>;
	#print "ERROR: No SNP column header.\n";
    	#exit;
}
$mySNPColNum   = $i;

$i             = 0;
$myBetaTColNum = 0;
$myORColNum    = 0;
while (uc( $Vals[$i] ) ne 'BETA'
    && uc( $Vals[$i] ) ne 'EFFECT'
    && ( uc( $Vals[$i] ) ne 'T' || $tdt )
    && uc( $Vals[$i] ) ne 'OR'
    && ( $i <= $#Vals ) )
{
    $i++;
}
if (   uc( $Vals[$i] ) ne 'BETA'
    && uc( $Vals[$i] ) ne 'EFFECT'
    && uc( $Vals[$i] ) ne 'T'
    && uc( $Vals[$i] ) ne 'OR' )
{
    print "No direction of effect field detected. Trying to use other indicators.\n";
    my $fd = fileno ${$DEL};
        $Stats_Filename = readlink("/proc/$$/fd/$fd");
	$c="awk '{if(NR==1){for(i=1;i<=NF;i++){if(\$i==\"CaseAF\"){CaseAFColNum=i}if(\$i==\"CtrlAF\"){CtrlAFColNum=i}}if(CaseAFColNum && CtrlAFColNum){print \$0\"\\tCaseAF-CtrlAF\"}else{print \"NoMatch CaseAF CtrlAF\";exit}}else {print \$0\"\\t\"\$CaseAFColNum-\$CtrlAFColNum}}' $Stats_Filename > $Stats_Filename.EFFECTColAdded";
        print "Run $c\n";
	$o=`$c`;
        close($DEL);
        open($DEL, "$Stats_Filename.EFFECTColAdded");
        $i=$#Vals+1;
        $line=<$DEL>;
	@Vals = split( /\s+/, $line );

	my $fd = fileno ${$DUP};
        $Stats_Filename = readlink("/proc/$$/fd/$fd");
        $c="awk '{if(NR==1){for(i=1;i<=NF;i++){if(\$i==\"CaseAF\"){CaseAFColNum=i}if(\$i==\"CtrlAF\"){CtrlAFColNum=i}}if(CaseAFColNum && CtrlAFColNum){print \$0\"\\tCaseAF-CtrlAF\"}else{print \"NoMatch CaseAF CtrlAF\";exit}}else {print \$0\"\\t\"\$CaseAFColNum-\$CtrlAFColNum}}' $Stats_Filename > $Stats_Filename.EFFECTColAdded";
	$o=`$c`;
        close($DUP);
        open($DUP, "$Stats_Filename.EFFECTColAdded");

	#print "ERROR: No BETA, EFFECT, T, or OR column header (at least 1 needed).\n";
    #exit;
}
if ( uc( $Vals[$i] ) eq 'BETA' || uc( $Vals[$i] ) eq 'T' || uc( $Vals[$i] ) eq 'EFFECT' || uc( $Vals[$i] ) eq 'CaseAF-CtrlAF') {
    $myBetaTColNum = $i;    #print"BETA or T matches\n";
}
if ( uc( $Vals[$i] ) eq 'OR' ) {
    $myORColNum = $i;       #print"OR matches\n";
}
if(uc( $Vals[$i] ) eq 'EFFECT')
{
	$Vals[$i]='BETA';
}
$myCaseEnrichStat = uc( $Vals[$i] );
print "myCaseEnrichStat=" . $myCaseEnrichStat . "\n";


print "myPColNum mySNPColNum myBetaTColNum myORColNum myChrColNum myBpColNum\n";
print $myPColNum. "\t"
  . $mySNPColNum . "\t"
  . $myBetaTColNum . "\t"
  . $myORColNum . "\t"
  . $myChrColNum . "\t"
  . $myBpColNum . "\n";

$JustPosIndex_caseEnrichDel    = 0;
$JustPosIndex_controlEnrichDel = 0;

#$c="rm CNVR_caseEnrichDel.txt CNVR_controlEnrichDel.txt CNVR_caseEnrichDup.txt CNVR_controlEnrichDup.txt";$o=`$c`;print"$o\n";
open( CNVR_caseEnrichDel,    ">temp/$out" . "caseEnrichDel.txt" );
open( CNVR_controlEnrichDel, ">temp/$out" . "controlEnrichDel.txt" );
while ( $snpStat = <$DEL> ) {
    chomp($snpStat);
    $snpStat =~ s/\r\n//;
    @stat = split( /\s+/, $snpStat );

    #if($stat[$myPColNum]>$CNVRMaxP)
    #{
    #        $CNVRMaxP=$stat[$myPColNum];
    #}
    #if($stat[$myPColNum]>$CNVRMaxPCon)
    #{
    #        $CNVRMaxPCon=$stat[$myPColNum];
    #}
    if ( $stat[$myPColNum] <= $maxPInclusion and $stat[$myPColNum] ne "NA" ) {
        if (
            (
                ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' )
                && $stat[$myBetaTColNum] > 0
            )
            || ( $myCaseEnrichStat eq 'OR'
                && ( $stat[$myORColNum] > 1 || $stat[$myORColNum] eq 'NA' ) )
          )
        {
            $caseEnrichDelSnps++;
            if ( $caseEnrichDelSnps == 1 ) {
                $CNVRStart = $stat[$myBpColNum];
                $CNVRMaxP  = $stat[$myPColNum];

                #$CNVRP=$stat[$myPColNum];
            }
            if ( $stat[$myPColNum] < $CNVRP || $caseEnrichDelSnps == 1 ) {
                $CNVRP  = $stat[$myPColNum];
                $TagSnp = $stat[$mySnpColNum];
                if ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' ) {
                    $CNVR_BTO = $stat[$myBetaTColNum];
                }
                elsif ( $myCaseEnrichStat eq 'OR' ) {
                    $CNVR_BTO = $stat[$myORColNum];
                }
                else {
                    print "ERROR: BETA, T, or OR column must be present.\n";
                }
            }
            if ( $stat[$myPColNum] > $CNVRMaxP ) {
                $CNVRMaxP = $stat[$myPColNum];
            }

            #print "$stat[$myPColNum]\n";
            if ( $stat[$myPColNum] == 0 ) {

                #print "NA is now 1\n";
                $stat[$myPColNum] = 1;
            }

            #print -(log($stat[$myPColNum])/(log(10)))."\n";
            if (
                   $stat[$myChrColNum] eq $lastChr
                && $stat[$myBpColNum] < ( $lastPos + $mergeDist )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) <
                    ( $lastNegLogP + $mergePVar ) )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) >
                    ( $lastNegLogP - $mergePVar ) )
              )
            {
                #print "Extend CNVR\n";
            }
            else {
                if ( $caseEnrichDelSnps > 1 ) {

                    #print "End CNVR\n";
                    if (   $stat[$myChrColNum] eq $lastChr
                        && $stat[$myBpColNum] < ( $lastPos + $mergeDist ) )
                    {
                    }
                    else {
                        $JustPosIndex_caseEnrichDel++;
                    }
                    $CNVREnd = $lastPos;
                    print CNVR_caseEnrichDel
"$lastChr $CNVRStart $CNVREnd $lastCNVRP $lastCNVR_BTO $lastTagSnp $CNVRMaxP $JustPosIndex_caseEnrichDel\n"
                      ;    #Cases Controls Gene SIDs RFs
                    $CNVRStart = $stat[$myBpColNum];
                    $CNVRP     = $stat[$myPColNum];
                    $CNVRMaxP  = $stat[$myPColNum];
                    $TagSnp    = $stat[$mySnpColNum];
                    if (   $myCaseEnrichStat eq 'BETA'
                        || $myCaseEnrichStat eq 'T' )
                    {
                        $CNVR_BTO = $stat[$myBetaTColNum];
                    }
                    elsif ( $myCaseEnrichStat eq 'OR' ) {
                        $CNVR_BTO = $stat[$myORColNum];
                    }

                    #$lastChrSigCaseEnrich=$lastChr;
                    #$lastPosSigCaseEnrich=$CNVREnd;
                }
            }
            $lastChr      = $stat[$myChrColNum];
            $lastPos      = $stat[$myBpColNum];
            $lastNegLogP  = -( log( $stat[$myPColNum] ) / ( log(10) ) );
            $lastCNVRP    = $CNVRP;
            $lastCNVR_BTO = $CNVR_BTO;
            $lastTagSnp   = $TagSnp;
        }
        else {
            #Control Enrich
            $controlEnrichDelSnps++;
            if ( $controlEnrichDelSnps == 1 ) {
                $CNVRStartCon = $stat[$myBpColNum];
                $CNVRMaxPCon  = $stat[$myPColNum];

                #$CNVRP=$stat[$myPColNum];
            }
            if ( $stat[$myPColNum] < $CNVRPCon || $controlEnrichDelSnps == 1 ) {
                $CNVRPCon  = $stat[$myPColNum];
                $TagSnpCon = $stat[$mySnpColNum];
                if ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' ) {
                    $CNVR_BTOCon = $stat[$myBetaTColNum];
                }
                elsif ( $myCaseEnrichStat eq 'OR' ) {
                    $CNVR_BTOCon = $stat[$myORColNum];
                }
                else {
                    print "ERROR: BETA, T, or OR column must be present.\n";
                }
            }
            if ( $stat[$myPColNum] > $CNVRMaxPCon ) {
                $CNVRMaxPCon = $stat[$myPColNum];
            }

            #print "$stat[$myPColNum]\n";
            if ( $stat[$myPColNum] == 0 ) {

                #print "NA is now 1\n";
                $stat[$myPColNum] = 1;
            }

            if (
                   $stat[$myChrColNum] eq $lastChrCon
                && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) <
                    ( $lastNegLogPCon + $mergePVar ) )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) >
                    ( $lastNegLogPCon - $mergePVar ) )
              )
            {
                #print "Extend CNVR\n";
            }
            else {
                if ( $controlEnrichDelSnps > 1 ) {

                    #print "End CNVR\n";
                    if (   $stat[$myChrColNum] eq $lastChrCon
                        && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist ) )
                    {
                    }
                    else {
                        $JustPosIndex_controlEnrichDel++;
                    }
                    $CNVREndCon = $lastPosCon;
                    print CNVR_controlEnrichDel
"$lastChrCon $CNVRStartCon $CNVREndCon $lastCNVRPCon $lastCNVR_BTOCon $lastTagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDel\n"
                      ;    #Cases Controls Gene SIDs RFs
                    $CNVRStartCon = $stat[$myBpColNum];
                    $CNVRPCon     = $stat[$myPColNum];
                    $CNVRMaxPCon  = $stat[$myPColNum];
                    $TagSnpCon    = $stat[$mySnpColNum];
                    if (   $myCaseEnrichStat eq 'BETA'
                        || $myCaseEnrichStat eq 'T' )
                    {
                        $CNVR_BTOCon = $stat[$myBetaTColNum];
                    }
                    elsif ( $myCaseEnrichStat eq 'OR' ) {
                        $CNVR_BTOCon = $stat[$myORColNum];
                    }

                    #$lastChrSigControlEnrich=$lastChrCon;
                    #$lastPosSigControlEnrich=$CNVREndCon;
                }
            }
            $lastChrCon      = $stat[$myChrColNum];
            $lastPosCon      = $stat[$myBpColNum];
            $lastNegLogPCon  = -( log( $stat[$myPColNum] ) / ( log(10) ) );
            $lastCNVRPCon    = $CNVRPCon;
            $lastCNVR_BTOCon = $CNVR_BTOCon;
            $lastTagSnpCon   = $TagSnpCon;
        }
    }

}

#Exit Condition
if (   $stat[$myChrColNum] eq $lastChr
    && $stat[$myBpColNum] <
    ( $lastPos + $mergeDist ) )    ## $stat[$myChrColNum] eq $lastChr &&
{
}
else {
    $JustPosIndex_caseEnrichDel++;
}
$CNVREnd = $lastPos;
if ( $caseEnrichDelSnps eq "" ) {
}
else {
    print CNVR_caseEnrichDel
"$lastChr $CNVRStart $CNVREnd $CNVRP $CNVR_BTO $TagSnp $CNVRMaxP $JustPosIndex_caseEnrichDel\n"
      ;    #$CasesControls\n"; #Cases Controls Gene SIDs RFs
}
if (   $stat[$myChrColNum] eq $lastChrCon
    && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist ) )
{
}
else {
    $JustPosIndex_controlEnrichDel++;
}
$CNVREndCon = $lastPosCon;
if ( $controlEnrichDelSnps eq "" ) {
}
else {
    print CNVR_controlEnrichDel
"$lastChrCon $CNVRStartCon $CNVREndCon $CNVRPCon $CNVR_BTOCon $TagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDel\n"
      ;    #Cases Controls Gene SIDs RFs
}

#DUPLICATIONS
$JustPosIndex_caseEnrichDup    = 0;
$JustPosIndex_controlEnrichDup = 0;
open( CNVR_caseEnrichDup,    ">temp/$out" . "caseEnrichDup.txt" );
open( CNVR_controlEnrichDup, ">temp/$out" . "controlEnrichDup.txt" );
$snpStat = <$DUP>;    #Discard Header (assume del and dup assoc have same header)
while ( $snpStat = <$DUP> ) {
    chomp($snpStat);
    $snpStat =~ s/\r\n//;
    @stat = split( /\s+/, $snpStat );

    #if($stat[$myPColNum]>$CNVRMaxP)
    #{
    #	$CNVRMaxP=$stat[$myPColNum];
    #}
    #if($stat[$myPColNum]>$CNVRMaxPCon)
    #{
    #        $CNVRMaxPCon=$stat[$myPColNum];
    #}
    if ( $stat[$myPColNum] <= $maxPInclusion and $stat[$myPColNum] ne "NA" ) {
        if (
            (
                ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' )
                && $stat[$myBetaTColNum] > 0
            )
            || ( $myCaseEnrichStat eq 'OR'
                && ( $stat[$myORColNum] > 1 || $stat[$myORColNum] eq 'NA' ) )
          )
        {
            $caseEnrichDupSnps++;
            if ( $caseEnrichDupSnps == 1 ) {
                $CNVRStart = $stat[$myBpColNum];
                $CNVRMaxP  = $stat[$myPColNum];

                #$CNVRP=$stat[$myPColNum];
            }
            if ( $stat[$myPColNum] < $CNVRP || $caseEnrichDupSnps == 1 ) {
                $CNVRP  = $stat[$myPColNum];
                $TagSnp = $stat[$mySnpColNum];
                if ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' ) {
                    $CNVR_BTO = $stat[$myBetaTColNum];
                }
                elsif ( $myCaseEnrichStat eq 'OR' ) {
                    $CNVR_BTO = $stat[$myORColNum];
                }
                else {
                    print "ERROR: BETA, T, or OR column must be present.\n";
                }
            }
            if ( $stat[$myPColNum] > $CNVRMaxP ) {
                $CNVRMaxP = $stat[$myPColNum];
            }

            #print "$stat[$myPColNum]\n";
            if ( $stat[$myPColNum] == 0 ) {

                #print "NA is now 1\n";
                $stat[$myPColNum] = 1;
            }

            if (
                   $stat[$myChrColNum] eq $lastChr
                && $stat[$myBpColNum] < ( $lastPos + $mergeDist )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) <
                    ( $lastNegLogP + $mergePVar ) )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) >
                    ( $lastNegLogP - $mergePVar ) )
              )
            {
                #print "Extend CNVR\n";
            }
            else {
                if ( $caseEnrichDupSnps > 1 ) {

                    #print "End CNVR\n";
                    if (   $stat[$myChrColNum] eq $lastChr
                        && $stat[$myBpColNum] < ( $lastPos + $mergeDist ) )
                    {
                    }
                    else {
                        $JustPosIndex_caseEnrichDup++;
                    }
                    $CNVREnd = $lastPos;
                    print CNVR_caseEnrichDup
"$lastChr $CNVRStart $CNVREnd $lastCNVRP $lastCNVR_BTO $lastTagSnp $lastCNVRMaxP $JustPosIndex_caseEnrichDup\n"
                      ;    #Cases Controls Gene SIDs RFs
                    $CNVRStart = $stat[$myBpColNum];
                    $CNVRP     = $stat[$myPColNum];
                    $CNVRMaxP  = $stat[$myPColNum];
                    $TagSnp    = $stat[$mySnpColNum];
                    if (   $myCaseEnrichStat eq 'BETA'
                        || $myCaseEnrichStat eq 'T' )
                    {
                        $CNVR_BTO = $stat[$myBetaTColNum];
                    }
                    elsif ( $myCaseEnrichStat eq 'OR' ) {
                        $CNVR_BTO = $stat[$myORColNum];
                    }

                    #$lastChrSigCaseEnrich=$lastChr;
                    #$lastPosSigCaseEnrich=$CNVREnd;
                }
            }
            $lastChr      = $stat[$myChrColNum];
            $lastPos      = $stat[$myBpColNum];
            $lastNegLogP  = -( log( $stat[$myPColNum] ) / ( log(10) ) );
            $lastCNVRP    = $CNVRP;
            $lastCNVR_BTO = $CNVR_BTO;
            $lastTagSnp   = $TagSnp;
            $lastCNVRMaxP = $CNVRMaxP;
        }
        else {
            #Control Enrich
            $controlEnrichDupSnps++;
            if ( $controlEnrichDupSnps == 1 ) {
                $CNVRStartCon = $stat[$myBpColNum];
                $CNVRMaxPCon  = $stat[$myPColNum];

                #$CNVRP=$stat[$myPColNum];
            }
            if ( $stat[$myPColNum] < $CNVRPCon || $controlEnrichDupSnps == 1 ) {

                #print"OK\n";
                $CNVRPCon  = $stat[$myPColNum];
                $TagSnpCon = $stat[$mySnpColNum];
                if ( $myCaseEnrichStat eq 'BETA' || $myCaseEnrichStat eq 'T' ) {
                    $CNVR_BTOCon = $stat[$myBetaTColNum];
                }
                elsif ( $myCaseEnrichStat eq 'OR' ) {
                    $CNVR_BTOCon = $stat[$myORColNum];
                }
                else {
                    print "ERROR: BETA, T, or OR column must be present.\n";
                }
            }
            if ( $stat[$myPColNum] > $CNVRMaxPCon ) {
                $CNVRMaxPCon = $stat[$myPColNum];
            }

            #print "$stat[$myPColNum]\n";
            if ( $stat[$myPColNum] == 0 ) {

                #print "NA is now 1\n";
                $stat[$myPColNum] = 1;
            }

            if (
                   $stat[$myChrColNum] eq $lastChrCon
                && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) <
                    ( $lastNegLogPCon + $mergePVar ) )
                && ( -( log( $stat[$myPColNum] ) / ( log(10) ) ) >
                    ( $lastNegLogPCon - $mergePVar ) )
              )
            {
                #print "Extend CNVR\n";
            }
            else {
                if ( $controlEnrichDupSnps > 1 ) {

                    #print "End CNVR\n";
                    if (   $stat[$myChrColNum] eq $lastChrCon
                        && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist ) )
                    {
                    }
                    else {
                        $JustPosIndex_controlEnrichDup++;
                    }
                    $CNVREndCon = $lastPosCon;
                    print CNVR_controlEnrichDup
"$lastChrCon $CNVRStartCon $CNVREndCon $lastCNVRPCon $lastCNVR_BTOCon $lastTagSnpCon $lastCNVRMaxPCon $JustPosIndex_controlEnrichDup\n"
                      ;    #Cases Controls Gene SIDs RFs
                    $CNVRStartCon = $stat[$myBpColNum];
                    $CNVRPCon     = $stat[$myPColNum];
                    $CNVRMaxPCon  = $stat[$myPColNum];
                    $TagSnpCon    = $stat[$mySnpColNum];
                    if (   $myCaseEnrichStat eq 'BETA'
                        || $myCaseEnrichStat eq 'T' )
                    {
                        $CNVR_BTOCon = $stat[$myBetaTColNum];
                    }
                    elsif ( $myCaseEnrichStat eq 'OR' ) {
                        $CNVR_BTOCon = $stat[$myORColNum];
                    }

                    #$lastChrSigControlEnrich=$lastChrCon;
                    #$lastPosSigControlEnrich=$CNVREndCon;
                }
            }
            $lastChrCon = $stat[$myChrColNum];
            $lastPosCon = $stat[$myBpColNum];

            #print "$stat[$myPColNum]\n";
            $lastNegLogPCon  = -( log( $stat[$myPColNum] ) / ( log(10) ) );
            $lastCNVRPCon    = $CNVRPCon;
            $lastCNVR_BTOCon = $CNVR_BTOCon;
            $lastTagSnpCon   = $TagSnpCon;
            $lastCNVRMaxPCon = $CNVRMaxPCon;
        }
    }

}

#Exit Condition
if (   $stat[$myChrColNum] eq $lastChr
    && $stat[$myBpColNum] < ( $lastPos + $mergeDist ) )
{
}
else {
    $JustPosIndex_caseEnrichDup++;
}
$CNVREnd = $lastPos;
if ( $caseEnrichDupSnps eq "" ) {
}
else {
    print CNVR_caseEnrichDup
"$lastChr $CNVRStart $CNVREnd $CNVRP $CNVR_BTO $TagSnp $CNVRMaxP $JustPosIndex_caseEnrichDup\n"
      ;    #Cases Controls Gene SIDs RFs
}
if (   $stat[$myChrColNum] eq $lastChrCon
    && $stat[$myBpColNum] < ( $lastPosCon + $mergeDist ) )
{
}
else {
    $JustPosIndex_controlEnrichDup++;
}
$CNVREndCon = $lastPosCon;
if ( $controlEnrichDupSnps eq "" ) {
}
else {
    print CNVR_controlEnrichDup
"$lastChrCon $CNVRStartCon $CNVREndCon $CNVRPCon $CNVR_BTOCon $TagSnpCon $CNVRMaxPCon $JustPosIndex_controlEnrichDup\n"
      ;    #Cases Controls Gene SIDs RFs
}

$c =
    "grep . temp/$out"
  . "caseEnrichDel.txt temp/$out"
  . "controlEnrichDel.txt temp/$out"
  . "caseEnrichDup.txt temp/$out"
  . "controlEnrichDup.txt | sed '/:[\\ ]*\$/d' | sed 's/temp\\///' | sed 's/$out//' | sed 's/D/\td/' | sed 's/:/\t/' | sed 's/.txt//' | sed 's/Enrich//' | sort -k6,6g | awk '!_[\$1\$2\$10]++' | awk '{print \$3,\$4,\$5,\$6,\$7,\$8,\$9,\$1,\$2}' | sed 's/\ /\t/g' > temp/$out"
  . "Report.txt";
$o = `$c`;

#2       1005300 1005300 0.00025375      -0.437555       2_1005300       0.0370378     control del
open( REPORT,       "temp/$out" . "Report.txt" );
open( REPORT2,      ">temp/$out" . "Report2.txt" );
open( ChrPosRanges, ">temp/$out" . "ChrPosRanges" );
$test_CHR_POS_to_VarID_pvar = "1_812283";

#print "\$ChrPos=".$ChrPos." \$PlinkCHR_POS_to_REF{\$ChrPos}=".$PlinkCHR_POS_to_REF{$ChrPos}."\n";
$VCF_VarID = $PlinkCHR_POS_to_VarID_pvar{$test_CHR_POS_to_VarID_pvar};

#CHR START(POS) END SNP_ID REF ALT SITEPOST
open( CHR_POS_or_END_SNP,
        "temp/$out"
      . "$inputNoPath"
      . "_wEND_CHR_POS_or_END_SNPID_wENDCol_wSitePostCol.txt" );
%CHR_POS_to_VarID_vcf;
%PlinkCHR_POS_to_REF;
%PlinkCHR_POS_to_ALT;

#123 5 ChrPosEnd Ref
while (<CHR_POS_or_END_SNP>) {
    chomp;
    @a                                = split( /\t/, $_ );
    $myCHR_POS                        = $a[0] . "_" . $a[1];
    $myCHR_POS_END                    = $a[0] . "_" . $a[1] . "_" . $a[2];;
    $CHR_POS_to_VarID_vcf{$myCHR_POS_END} = $a[3];
    $PlinkCHR_POS_to_REF{$myCHR_POS_END}  = $a[4];
    $PlinkCHR_POS_to_ALT{$myCHR_POS_END}  = $a[5];
}
close(CHR_POS_or_END_SNP);
print "test_CHR_POS_to_VarID_pvar="
  . "$test_CHR_POS_to_VarID_pvar\t$VCF_VarID\n";
$VCF_VarID = $CHR_POS_to_VarID_vcf{$test_CHR_POS_to_VarID_pvar};
print "test_CHR_POS_or_END_SNP=" . "$test_CHR_POS_to_VarID_pvar\t$VCF_VarID\n";

while ( $line = <REPORT> ) {
    chomp($line);
    @a = split( /\t/, $line );
    print "Here is array a contents initially @a\n";

    #print "$a[5]\n";
    #Could switch out bed for pgen but plink2 may not have recode
    #Update plink2 has --export ped
    #temp/$out"."$inputNoPath
    $CHR_POS_ID = $a[5]; # . "_" . $a[2];
    "$CHR_POS_ID";    # stringified
                      #$VCF_VarID=$PlinkCHR_POS_to_VarID_pvar{$CHR_POS_ID};
    $CHR_POS_ID =~ s/:|-/_/g;
    $VCF_VarID = $CHR_POS_to_VarID_vcf{$CHR_POS_ID};

#if($VCF_VarID) may need to add back in case of any variant exclusion
#Trying to reformulate the multiallelic to allow ped conversion when getting contributing samples and calls
#Error: .pgen file contains multiallelic variants, while .pvar does not.
#$c="awk -F\"\\t\" 'BEGIN{OFS=\"\\t\"}{gsub(/,.*/,\"\",\$5);print \$0}' temp/$out"."$inputNoPath".".pvar > temp/$out"."$inputNoPath".".pvar2; cp temp/$out"."$inputNoPath".".pvar temp/$out"."$inputNoPath".".pvar.bak; mv temp/$out"."$inputNoPath".".pvar2 temp/$out"."$inputNoPath".".pvar";$o=`$c`;
    print "CHR_POS_ID=" . $CHR_POS_ID . " VCF_VarID=" . $VCF_VarID . "\n";
    $myChrPos = $a[5]; #$a[0] . "_" . $a[1];
    $myChrPos =~ s/:|-/_/g;
    if ( $PlinkCHR_POS_to_ALT{$myChrPos} =~ /,/ ) {
	    #temp/1KG_SudmantEtAl_wRvTestsVcfHeader_and_multiallelic 
	    #$c=$cat_zcat . " " . $input . " | grep ^#CHROM > temp/$out" . "VcfHeader_and_multiallelic;" . $cat_zcat . " " . $input . " | grep -P \"$a[0]\\t$a[1]\\t\" >> temp/$out" . "VcfHeader_and_multiallelic";
	    #$c="awk '{for(i=10;i<=NF;i++){if(NR==1){IDs[i]=\$i}else{gsub(/:.*/,\"\",\$i);if(\$i!=\"0/0\"&&\$i!=\"0|0\"){print IDs[i]}}}}' temp/$out" . "VcfHeader_and_multiallelic | sort > temp/$out" . "_multiallelic_NonRefRef";$o=`$c`;
	    $c=$cat_zcat . " " . $input . " | grep ^\#CHROM > temp/$out" . "VcfHeader_and_multiallelic;" . $cat_zcat . " " . $input . " | grep -P \"$a[0]\\t$a[1]\\t\" >> temp/$out" . "VcfHeader_and_multiallelic";print "$c\n";$o=`$c`;print "$o\n";
            $c="awk '{for(i=10;i<=NF;i++){if(NR==1){IDs[i]=\$i}else{gsub(/:.*/,\"\",\$i);if(\$i!=\"0\/0\"&&\$i!=\"0\|0\"){print IDs[i]}}}}' temp/$out" . "VcfHeader_and_multiallelic > temp/$out" . "_multiallelic_NonRefRef";print "$c\n";$o=`$c`;print "$o\n";

	    #$PlinkCHR_POS_to_ALTtwo{$ChrPos} = $vals[4];
	    $c="awk '{if(NR>1){print \$2\"\\t\"\$3}}' temp/$out"."file.pheno | sort > temp/$out"."file.pheno2col";$o=`$c`;
	    $c="join temp/$out" . "_multiallelic_NonRefRef temp/$out"."file.pheno2col | awk 'BEGIN{cases=0}{if(\$2==2)cases++}END{ORS=\"\";print cases}'";
	    $Cases = `$c`;
	    $c="join temp/$out" . "_multiallelic_NonRefRef temp/$out"."file.pheno2col | awk 'BEGIN{controls=0}{if(\$2==1)controls++}END{ORS=\"\";print controls}'";
	    $Controls = `$c`;
	    $c="join temp/$out" . "_multiallelic_NonRefRef temp/$out"."file.pheno2col | awk '{if(\$2==2)printf \$1\",\"}' | sed 's/,\$//'";
    $CaseSIDs = `$c`;
    		$c="join temp/$out" . "_multiallelic_NonRefRef temp/$out"."file.pheno2col | awk '{if(\$2==1)printf \$1\",\"}' | sed 's/,\$//'";
    $ControlSIDs = `$c`;
	}
	else{	
    $c =
        "echo $VCF_VarID > temp/$out"
      . "TagSnp; $MyDirectoryPathPrefix"
      . "PerlModules/./plink2 --pfile temp/$out"
      . "$inputNoPath --snp $VCF_VarID --export ped --pheno temp/$out"
      . "file.pheno --allow-no-sex --allow-extra-chr --out temp/$out"
      . "TagSnp";
    $o = `$c`;

#$c="echo $a[5] > temp/$out"."TagSnp; $MyDirectoryPathPrefix"."PerlModules/./plink2 --bfile temp/$out$a[8] --extract temp/$out"."TagSnp --export ped --allow-no-sex --out temp/$out"."TagSnp";$o=`$c`;#print "$c\n$o\n";
    $myChrPos = $a[5]; #$a[0] . "_" . $a[1];
    $myChrPos_ = $a[5];
    $myChrPos =~ s/:|-/_/g;
    $myREF    = $PlinkCHR_POS_to_REF{$myChrPos};
    print "\$myREF= " . $myREF . "\n";
    $c = "grep -vP \"$myREF\\t$myREF\$\" temp/$out"
      . "TagSnp.ped | awk 'BEGIN{cases=0}{if(\$6==2)cases++}END{ORS=\"\";print cases}'";
    $Cases = `$c`;
    $c     = "grep -vP \"$myREF\\t$myREF\$\" temp/$out"
      . "TagSnp.ped | awk 'BEGIN{controls=0}{if(\$6==1)controls++}END{ORS=\"\";print controls}'";
    $Controls = `$c`;
    $c        = "grep -vP \"$myREF\\t$myREF\$\" temp/$out"
      . "TagSnp.ped | awk '{print \$2\"\\t\"\$6}' | awk '{if(\$2==2)printf \$1\",\"}' | sed 's/,\$//'";
    $CaseSIDs = `$c`;
    $c        = "grep -vP \"$myREF\\t$myREF\$\" temp/$out"
      . "TagSnp.ped | awk '{print \$2\"\\t\"\$6}' | awk '{if(\$2==1)printf \$1\",\"}' | sed 's/,\$//'";
    $ControlSIDs = `$c`;
   }
    if ( $Cases == 0 )    { $CasesSIDs   = "NA"; }
    if ( $Controls == 0 ) { $ControlSIDs = "NA"; }
    @chr_pos_end = split( /_/, $myChrPos );
    print "Here is array a contents before print @a\n";
    #Here is array a contents initially 17 9724255 9724255 0.03418 0.1348 17:9724255-9725522 0.03418 control del
    #Here is array a contents before print 17 9724255 9724255 0.03418 0.1348 17:9724255-9725522 0.03418 control del
    @chr_pos_end = split( /_/, $CHR_POS_ID );
    print "Here is array chr_pos_end contents before print @chr_pos_end\n";
    if($a[4]=="NA"){$a[4]="infinity"};
    if($a[4]=="0"){$a[4]="-infinity"};
	print REPORT2 "$chr_pos_end[0]\t$chr_pos_end[1]\t$chr_pos_end[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$Cases\t$Controls\t$CaseSIDs\t$ControlSIDs\n";

    #exit; #TESTING
    print ChrPosRanges "chr".$a[5]."\n"; #"$myChrPos\n"; #"$a[0]:$a[1]-$a[2]\n";
}
close(ChrPosRanges);
$c = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_knownGene.txt -knowngene "
  . "-kgxref "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_kgXref.txt -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/'| sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_SigCallsGene_Col23.txt";
$o = `$c`;

#$command = "perl ".$MyDirectoryPathPrefix."PerlModules/scan_region.pl temp/".$out."ChrPosRanges_DEL.txt ".$MyDirectoryPathPrefix."GeneRef/".$build."_knownGene.txt -knowngene "."-kgxref ".$MyDirectoryPathPrefix."GeneRef/".$build."_kgXref.txt -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*\\t\\t//' | sed 's/^chr[^\\t]*\\t//' > temp/".$out."SigCallsGene_DEL_Col23.txt";

#RFs: SegDups (count, max, avg) DgvEntries      TeloCentro      AvgGC   AvgProbes       Recurrent       PopFreq PenMaxP_Freq_HighFreq           FreqInflated    Sparse  ABFreq  AvgConf AvgLength
$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_genomicSuperDups_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDup_Col23.txt";
$o       = `$command`;
$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_dgv_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk -F\"\\t\" '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDupDgv_Col23.txt";
$o       = `$command`;
$defFile = $MyDirectoryPathPrefix . "GeneRef/" . $build . "_cytoBand.txt";

#Do once before as part of GeneRef setup #$c="join <(awk '{if(\$2==0)print}' GeneRef/hg19_cytoBand.txt | sort -k1,1) <(grep acen GeneRef/hg19_cytoBand.txt | sort -u -k1,1 | sort -k1,1) | join - <(grep acen GeneRef/hg19_cytoBand.txt | sort -k3,3nr | sort -u -k1,1 | sort -k1,1) | join - <(sort -k3,3nr GeneRef/hg19_cytoBand.txt | sort -u -k1,1 | sort -k1,1) | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\ttelomere\\n\"\$1\"\\t\"\$6\"\\t\"\$7\"\\tcentromere\\n\"\$1\"\\t\"\$10\"\\t\"\$11\"\\tcentromere\\n\"\$1\"\\t\"\$14\"\\t\"\$15\"\\ttelomere\"}' | sed 's/^chr//' > ".$MyDirectoryPathPrefix."GeneRef/".$build."_cytoBand_TeloCentro.bed";$o=`$c`;
$c =
    "chmod +x "
  . $MyDirectoryPathPrefix
  . "PerlModules/bedtools; "
  . $MyDirectoryPathPrefix
  . "PerlModules/bedtools intersect -a temp/$out"
  . "Report2.txt -b "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_cytoBand_TeloCentro.bed -wo | awk '{print \$1\":\"\$2\"-\"\$3\"\\t\"\$(NF-1)}' > temp/"
  . $out
  . "CNVR_TeloCentro.txt";
$o = `$c`;
$c =
    "awk '{print \$1\":\"\$2\"-\"\$3}' temp/$out"
  . "Report2.txt > a;perl $MyDirectoryPathPrefix"
  . "PerlModules/Vlookup.pl temp/"
  . $out
  . "CNVR_TeloCentro.txt a | awk '{print \$2}' > b; mv b temp/"
  . $out
  . "CNVR_TeloCentro.txt";
$o = `$c`;

$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_gc5Base_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDupDgvTcGc_Col23.txt";
$o = `$command`;

#$c="wc -l temp/$out"."fam | awk '{print \$1}'";$totalSamples=`$c`;
$c            = "wc -l temp/$out". "$inputNoPath.psam | awk '{print \$1-1}'";
$totalSamples = `$c`;
$c =
    "awk -F\"\\t\" '{print (\$(NF-3)+\$(NF-2))/$totalSamples}' temp/$out"
  . "Report2.txt > temp/"
  . $out
  . "CNVR_PopFreq.txt";
$o = `$c`;

#peninsula
#inflated
$c =
    "awk -F\"\\t\" '{print \$12}' temp/$out"
  . "Report2.txt | tr ',' '\n' | sort | grep -v ^\$ | uniq -c > temp/$out"
  . "CaseAndControlSIDs.txt; awk -F\"\\t\" '{print \$13}' temp/$out"
  . "Report2.txt | tr ',' '\n'  | sort | grep -v ^\$ | uniq -c >> temp/$out"
  . "CaseAndControlSIDs.txt; sort -k1,1nr temp/$out"
  . "CaseAndControlSIDs.txt | head | awk '{print $2}' > temp/$out"
  . "CaseAndControlSIDs_10MostInflated.txt";
$o = `$c`;

#cytoband
$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_cytoBand_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDupDgvTcGcPerPenInfCyt_Col23.txt";
$o = `$command`;

#recurrent
$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_somaticRearrangement_SimFormat_AllCol.sorted -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDupDgvTcGcPerPenInfCytSom_Col23.txt";
$o = `$command`;

#exon
$command = "perl "
  . $MyDirectoryPathPrefix
  . "PerlModules/scan_region.pl temp/$out"
  . "ChrPosRanges "
  . $MyDirectoryPathPrefix
  . "GeneRef/"
  . $build
  . "_knownGene_Exons_SimFormat_AllCol_UniqueIDs -knowngene -expandmax 5m | sed 's/^NOT_FOUND\\t0\\t[^\\t]*/NA\tNA/' | sed 's/^chr[^\\t]*\\t//' | awk '{if(\$2>0){print\$1\"-closest\"}else{print\$1}}' > temp/"
  . $out
  . "CNVR_ALL_CDupDgvTcGcPerPenInfCytSomExo_Col23.txt";
$o = `$command`;

#PercentSamples
#$c="perl OutputUtilities/PercentSamples.pl ID_Category.txt SNP_SampleIDs.txt";$o=`$c`;

$c =
    "paste temp/$out"
  . "Report2.txt temp/$out"
  . "CNVR* > temp/$out"
  . "Report3.txt";
$o = `$c`;
$c = "sort -k4,4g temp/$out" . "Report3.txt > temp/$out" . "Report4.txt";
$o = `$c`;
$c = "cat temp/$out" . "Report4.txt > temp/$out" . "Report5.txt";
$o = `$c`;    #was cat header but now later

#printf( "%-10s", "hello" );
#printf( "%10d", 777);
#FAIL=CNVRLength/AvgLength<1000, DgvFreq>0.01/DgvEntries>10, MaxP>0.5, PopFreq>0.01, SegDups>10, Recurrent, FreqInflated>0.5, AvgConf<10
#Sort RF <=3, delP<5x10^-4 and ORDup>1 OR dupP<5x10^-4 and ORDup>1, (on exon)
open( INFLATED, "temp/$out" . "CaseAndControlSIDs_10MostInflated.txt" );
while ( $line = <INFLATED> ) {
    chomp($line);
    @a = split( /\s+/, $line );
    $InflatedSamples{ $a[2] } = 1;
}

open( REPORT,       "temp/$out" . "Report5.txt" );
open( REPORT2,      ">temp/$out" . "Report6.txt" );
open( SNP_IDS,      ">temp/$out" . "Snp_IDs.txt" );
open( REPORT_BRIEF, ">$out" . "Report.txt" );
printf REPORT_BRIEF (
    "%3s %11s %10s %10s %6s %8s %8s %6s %5s %5s\n", "chr",
    "start($build)",                                "stop",
    "p",                                            $myCaseEnrichStat,
    "cases",                                        "controls",
    "filter",                                       "type",
    "tag"
);
$DgvEntries     = 0;
$SegDupsEntries = 0;

while ( $line = <REPORT> ) {
    chomp($line);
    @a           = split( /\t/, $line );
    $lengthTotal = $h_lengthTotal{ $a[5] };
    if ( $a[9] + $a[10] > 0 ) {

        #print"OK\n";
        $AvgLength = $lengthTotal / ( $a[9] + $a[10] );
    }
    else {
        $AvgLength = -9;
    }

    #DGV
    if ( !( $a[14] =~ /-closest/ ) ) {
        @Dgv        = split( /,/, $a[14] );
        $DgvEntries = $#Dgv + 1;
    }
    if ( !( $a[13] =~ /-closest/ ) ) {
        @SegDups        = split( /,/, $a[13] );
        $SegDupsEntries = $#SegDups + 1;
    }
    if ( $a[17] eq "NA-closest" ) {
        $recurrent = "N";
    }
    else {
        $recurrent = "Y";
    }
    @caseIDs    = split( /,/, $a[11] );
    @controlIDs = split( /,/, $a[12] );
    for ( $i = 0 ; $i <= $#caseIDs ; $i++ ) {
        if ( exists( $InflatedSamples{ $caseIDs[$i] } ) ) {
            $InflatedSamplesCount++;
        }
    }
    for ( $i = 0 ; $i <= $#controlIDs ; $i++ ) {
        if ( exists( $InflatedSamples{ $controlIDs[$i] } ) ) {
            $InflatedSamplesCount++;
        }
    }
    if ( $a[9] + $a[10] > 0 ) {
        $InflatedSamplesFreq = $InflatedSamplesCount / ( $a[9] + $a[10] );
    }
    else {
        $InflatedSamplesFreq = -9;
    }
    $confidenceTotal = $h_confidenceTotal{ $a[5] };
    if ( $a[9] + $a[10] > 0 ) {

        #print"OK\n";
        $AvgConfidence = $confidenceTotal / ( $a[9] + $a[10] );
    }
    else {
        $AvgConfidence = -9;
    }
    if ( $AvgLength < 1000 )          { $RF++; }
    if ( $DgvEntries > 10 )           { $RF++; }
    if ( $a[6] > 0.5 )                { $RF++; }    #MaxP
    if ( $a[19] > 0.01 )              { $RF++; }    #PopFreq
    if ( $SegDupsEntries > 10 )       { $RF++; }
    if ( $recurrent eq "Y" )          { $RF++; }
    if ( $InflatedSamplesFreq > 0.5 ) { $RF++; }
    if ( $AvgConfidence < 10 )        { $RF++; }

    #if($a[3]<0.0005 && $a[7] eq "case"){$RF--;}#p and direction green flag
    if ( !( $a[18] =~ "-closest" ) ) { $RF--; }     #exon green flag

    if ( $RF > 2 ) { $RF_PassFail = "FAIL"; $overallFail++; }
    else           { $RF_PassFail = "PASS"; }

    print REPORT2
"$line\t$AvgLength\t$DgvEntries\t$a[6]\t$a[19]\t$SegDupsEntries\t$recurrent\t$InflatedSamplesFreq\t$AvgConfidence\t$RF\t$RF_PassFail\n";

    printf REPORT_BRIEF (
        "%3d %11d %10d %10s %6s %8d %8d %6s %5s %5s\n",
        $a[0], $a[1],  $a[2],        $a[3], $a[4],
        $a[9], $a[10], $RF_PassFail, $a[8], $a[5]
    );

    $InflatedSamplesCount = 0;
    $countLines++;
    $a[11] =~ s/,/\ /g;
    $a[12] =~ s/,/\ /g;
    print SNP_IDS "$a[5]\t";
    if ( $a[11] ne "" ) {
        print SNP_IDS " $a[11]";
    }
    if ( $a[12] ne "" ) {
        print SNP_IDS " $a[12]";
    }
    print SNP_IDS "\n";
    $DgvEntries     = 0;
    $SegDupsEntries = 0;
    $RF             = 0;
}

#$countLines--;
$c =
    "sort -k1,1 -k2,2n temp/$out"
  . "$inputNoPath.rawcnv2 | sed 's/\\t\$//' > temp/$out"
  . "$inputNoPath.rawcnv2Sorted";
$o = `$c`;
$c =
"awk '{print \$10\"\\t\"\$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9}' $out"
  . "Report.txt | sed 's/_/\\t/' | awk '{print \$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5\"\\t\"\$6\"\\t\"\$7\"\\t\"\$8\"\\t\"\$9\"\\t\"\$10\"\\t\"\$11}' | sort -k1,1 -k2,2n | grep -v tag > $out"
  . "Report.txtSorted";
$o = `$c`;
$c =
    $MyDirectoryPathPrefix
  . "PerlModules/bedtools intersect -a temp/$out"
  . "$inputNoPath.rawcnv2Sorted -b $out"
  . "Report.txtSorted -sorted -wo | awk '{if((\$4<2&&\$(NF-1)==\"del\")||(\$4>=2&&\$(NF-1)==\"dup\"))print}' > $out"
  . "Report_ContributingCalls.txt";
$o = `$c`;
if ( $countLines > 0 ) {
    $FailRate = $overallFail / $countLines;
}
if ( $FailRate > 0.8 ) {
    print "WARNING: High fail rate of "
      . sprintf( "%.3f", $FailRate )
      . " but some results may still be usable\n";
}

#printf( "%-10s", "hello" );
##printf( "%10d", 777);
##FAIL=CNVRLength/AvgLength<1000, DgvFreq>0.01/DgvEntries>10, MaxP>0.5, PopFreq>0.01, SegDups>10, Recurrent, FreqInflated>0.5, AvgConf<10
##Sort RF <=3, delP<5x10^-4 and ORDup>1 OR dupP<5x10^-4 and ORDup>1, (on exon)
close(REPORT2);
close(SNP_IDS);

if ( $batch ne "" ) {
    $c =
      "perl OutputUtilities/PercentSamples.pl $batch temp/$out" . "Snp_IDs.txt";
    $o = `$c`;
    $c = "cat "
      . $MyDirectoryPathPrefix
      . "header temp/$out"
      . "Report6.txt > temp/$out"
      . "Report_Verbose.txt; mv temp/$out"
      . "Report_Verbose.txt temp/$out"
      . "Report6.txt";
    $o = `$c`;
    $c =
        "cat temp/$out"
      . "Snp_IDs.txt_CohortCounts.txt | cut -f 2- > temp/$out"
      . "Snp_IDs.txt_CohortCounts_Cols.txt; paste temp/$out"
      . "Report6.txt temp/$out"
      . "Snp_IDs.txt_CohortCounts_Cols.txt > $out"
      . "Report_Verbose.txt";
    $o = `$c`;

}
else {
    $c = "cat "
      . $MyDirectoryPathPrefix
      . "header temp/$out"
      . "Report6.txt > $out"
      . "Report_Verbose.txt";
    $o = `$c`;
}

print "Output written to $out" . "Report.txt\n";
print LOG "Output written to $out" . "Report.txt\n";

#$GIFDel = median ( @DelSnps_x2 ) / 0.456 ; ####### http://en.wikipedia.org/wiki/Population_stratification 1-d.f. x2 distribution independence assumption! so somehow need to prune down to all CNVRs for all significance levels
#$GIFDup = median ( @DupSnps_x2 ) / 0.456 ; ####### Always 0 because median x2=0 (fisher p=1) from a majority of SNPs not having any CNVs need to prune out

open( ALLIDS, "temp/$out" . "AllIDs.txt" );
open( VCF,    ">$out" . "Report_Verbose.vcf" );
print VCF "##fileformat=VCFv4.3\n";
$DATE = `date +%Y%m%d`;
print VCF "##fileDate=$DATE";
print VCF "##reference=$build\n";
print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
while ( $line = <ALLIDS> ) {
    chomp($line);
    $AllIDs[$countIDs] = $line;
    print VCF "\t$line";
    $countIDs++;
}
print VCF "\n";
open( FILE, "temp/$out" . "Report6.txt" );
$line = <FILE>;    #Header
while ( $line = <FILE> ) {
    chomp($line);
    @a = split( /\t/, $line );
    print VCF
      "$a[0]\t$a[1]\t$a[0]:$a[1]-$a[2]\t.\t<$a[8]>\t.\t$a[-1]\tEND=$a[2]\tGT";
    @cases      = split( /,/, $a[11] );
    @controls   = split( /,/, $a[12] );
    %VariantIDs = ();
    for ( $i = 0 ; $i <= $#cases ; $i++ ) {
        $VariantIDs{ $cases[$i] } = 1;
    }
    for ( $i = 0 ; $i <= $#controls ; $i++ ) {
        $VariantIDs{ $controls[$i] } = 1;
    }
    for ( $i = 0 ; $i < $countIDs ; $i++ ) {
        if ( exists( $VariantIDs{ $AllIDs[$i] } ) ) {
            print VCF "\t0/1";
        }
        else {
            print VCF "\t0/0";
        }
    }
    print VCF "\n";
}

%MonthPastStart = (
    "0", "Jan", "1", "Feb", "2",  "Mar", "3",  "Apr",
    "4", "May", "5", "Jun", "6",  "Jul", "7",  "Aug",
    "8", "Sep", "9", "Oct", "10", "Nov", "11", "Dec"
);
@timeData = localtime(time);
$Year     = $timeData[5] + 1900;
if ( $timeData[2] > 12 ) {
    $timeData[2] -= 12;
    $AmOrPm = "PM";
}
elsif ( $timeData[2] eq 12 ) {
    $AmOrPm = "PM";
}
else {
    $AmOrPm = "AM";
}
if ( $timeData[2] < 10 ) { $timeData[2] = "0" . $timeData[2]; }
if ( $timeData[1] < 10 ) { $timeData[1] = "0" . $timeData[1]; }
if ( $timeData[0] < 10 ) { $timeData[0] = "0" . $timeData[0]; }
print "Run Ended "
  . $timeData[2] . ":"
  . $timeData[1] . ":"
  . $timeData[0]
  . " $AmOrPm on "
  . $MonthPastStart{ $timeData[4] } . " "
  . $timeData[3] . " "
  . $Year . "\n";
print LOG "Run Ended "
  . $timeData[2] . ":"
  . $timeData[1] . ":"
  . $timeData[0]
  . " $AmOrPm on "
  . $MonthPastStart{ $timeData[4] } . " "
  . $timeData[3] . " "
  . $Year . "\n";

