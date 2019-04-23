##fileformat=VCFv4.3
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001
#1 2827694 rs2376870 CGTGGATGCGGGGAC C . PASS SVTYPE=DEL;END=2827708;HOMLEN=1;HOMSEQ=G;SVLEN=-14 GT:GQ 1/1:14
#2 321682 . T <DEL> 6 PASS SVTYPE=DEL;END=321887;SVLEN=-205;CIPOS=-56,20;CIEND=-10,62 GT:GQ 0/1:12
#2 14477084 . C <DEL:ME:ALU> 12 PASS SVTYPE=DEL;END=14477381;SVLEN=-297;CIPOS=-22,18;CIEND=-12,32 GT:GQ 0/1:12
#3 9425916 . C <INS:ME:L1> 23 PASS SVTYPE=INS;END=9425916;SVLEN=6027;CIPOS=-16,22 GT:GQ 1/1:15
#3 12665100 . A <DUP> 14 PASS SVTYPE=DUP;END=12686200;SVLEN=21100;CIPOS=-500,500;CIEND=-500,500 GT:GQ:CN:CNQ ./.:0:3:16.2
#4 18665128 . T <DUP:TANDEM> 11 PASS SVTYPE=DUP;END=18665204;SVLEN=76;CIPOS=-10,10;CIEND=-10,10 GT:GQ:CN:CNQ ./.:0:5:8.3
unless( open( SAMPSHEET, $ARGV[0]))
{
    print "NOTICE: VCF File Not Found! $ARGV[0]\n";
    exit;
}
$doneHeader='FALSE';
print "#chr\tstart\tstop\tcn\tsample\tquality\n";
while (defined ($inputLine = <SAMPSHEET>))
{
        chomp($inputLine);
        $inputLine =~ s/[\r\n]+$//;
        if($inputLine=~'^##reference=')
        {
		$inputLine =~ s/##reference=//;
		if($inputLine=~/19/ || $inputLine=~/37/)
		{
			$build="hg19";print "NOTICE: build=$build\n";
		}
		elsif($inputLine=~/18/ || $inputLine=~/36/)
		{
			$build="hg18";print "NOTICE: build=$build\n";
		}
		elsif($inputLine=~/38/)
		{
			$build="hg38";print "NOTICE: build=$build\n";
		}
		else
		{
			print "NOTICE: $inputLine found as reference does not match hg18/GRCh36, hg19/GRCh37, or hg38/GRCh38.\n";
		}
	}
	if($inputLine=~'^#CHR')
	{
		$header=$inputLine;
		@HeaderVals=split(/\s+/,$header);
		for($iSample=9;$iSample<=$#HeaderVals;$iSample++)
        	{
        		#print "$HeaderVals[$iSample]\n";
       		}
		$doneHeader='TRUE';
		print "NOTICE: DONE HEADER!\n";
	}
	if($doneHeader eq "TRUE" && $inputLine!~'^#CHR')
	{
		@Vals=split(/\s+/,$inputLine);
		$Vals[0]=~s/chr//;
		$chr=$Vals[0];
		$start=$Vals[1];
		$id=$Vals[2];
		$ref=$Vals[3];
		$ref_len=length($ref);
		$alt=$Vals[4];
		$alt_len=length($alt);
		$alt=~s/\<CN//g;
		$alt=~s/\>//g;
		@alts=split(/\,/,$alt);
		$qual=$Vals[5];
		$filter=$Vals[6];
		$info=$Vals[7];
		@info_semi=split(/;/,$info);
		$times_type=0;
		$times_end=0;
		$times_len=0;
		for($i=0;$i<=$#info_semi;$i++)
		{
			@info_semi_equal=split(/=/,$info_semi[$i]);
			if($info_semi[$i]=~'TYPE')
			{
				if($times_type==0)
				{
					$type=$info_semi_equal[1];
				}
				$times_type++;
			}
			if($info_semi[$i]=~'^END')
			{
				if($times_end==0)
				{
					$end=$info_semi_equal[1];
				}
				$times_end++;
			}
			if($info_semi[$i]=~'SVLEN' || $info_semi[$i]=~'CNVLEN' || $info_semi[$i]=~'LENGTH')
			{
				if($times_len==0)
				{
					$len=$info_semi_equal[1];
					if($len<0 && $type eq "")
					{
						$type="DEL";
					}
					elsif($len>0 && $type eq "")
					{
						$type="DUP";
					}
				}
				$times_len++;
			}
		}
		if($end eq "")
		{
			if($len eq "")
			{
				if($ref_len == $alt_len)
				{
					print "NOTICE: Cannot figure out variant end!\n";
				}
				else
				{
					$end=$start+abs($ref_len-$alt_len)+1;
				}
			}
			else
			{
				$end=$start+abs($len)+1;
			}	
		}
		$format=$Vals[8];
		@format_colon=split(/:/,$format);
		for($i=0;$i<=$#format_colon;$i++)
                {
			if($format_colon[$i] eq "GT")
			{
				$iGT=$i;
			}
			elsif($format_colon[$i] eq "GQ")
			{
				$iGQ=$i;
			}
		}
		for($iSample=9;$iSample<=$#HeaderVals;$iSample++)
                {
			#print $iSample."\n";
			@GT=split(/:/,$Vals[$iSample]);
			if(($GT[$iGT] eq '0/1' || $GT[$iGT] eq '1/0' || $GT[$iGT] eq '1' || $GT[$iGT] eq '2' || $GT[$iGT] eq '0|1' || $GT[$iGT] eq '1|0') && $#alts<1)
			{
				if($type=~"DEL"||$type=~"Del"||$type=~"del"||$GT[$iGT] eq '1')
				{
					$cn=1;
				}
				else
				{
					$cn=3;
				}
				print $chr."\t".$start."\t".$end."\t".$cn."\t".$HeaderVals[$iSample]."\t".$GT[$iGQ]."\n";
			}
			elsif(($GT[$iGT] eq '1/1' || $GT[$iGT] eq '1|1') && $#alts<1)
			{
				if($type=~"DEL"||$type=~"Del"||$type=~"del")
                                {
                                        $cn=0;
                                }
                                else
                                {
                                        $cn=4;
                                }
				print $chr."\t".$start."\t".$end."\t".$cn."\t".$HeaderVals[$iSample]."\t".$GT[$iGQ]."\n";
			}
			elsif($GT[$iGT] eq '0/0' || $GT[$iGT] eq '0|0')
			{
				#Ref/Ref is not reported diploid cn=2
				#print "Ref/Ref";
			}
			else
			{
				#@GTs=();
				#print"$#GTs ";
				if($GT[$iGT]=~/\//)
				{@GTs=split(/\//,$GT[$iGT]);}
				else
				{@GTs=split(/\|/,$GT[$iGT]);}
				#print "$#alts";
				if($GTs[0]==0)
                                {$cn=1;}
				elsif($alts[$GTs[0]-1] > -1)
				{$cn=$alts[$GTs[0]-1];}#print"START=$cn ";}
				#elsif($GTs[0]==$alts[1])
				#{$cn=$alts[1];}
				else
				{print "NOTICE: Could not find $GTs[0] in alts!\n";}
				if($GTs[1]==0 && $#GTs>0)
                                {$cn+=1;}
                                elsif($alts[$GTs[1]-1] > -1 && $#GTs>0)
                                {$cn+=$alts[$GTs[1]-1];}#print"END+=$alts[$GTs[1]-1] ";}
				#elsif($GTs[1]==0)
				#{$cn+=1;}
				else
				{print "NOTICE: Could not find $GTs[1] in alts!\n";}

                                if($type=~"DEL"||$type=~"Del"||$type=~"del"||$GT[$iGT] eq '1')
                                {
                                        $cn=1;
                                }
                                else
                                {
                                        $cn=3;
                                }

				print $chr."\t".$start."\t".$end."\t".$cn."\t".$HeaderVals[$iSample]."\t".$GT[$iGQ]."\n";
			}
		}
	}
}
