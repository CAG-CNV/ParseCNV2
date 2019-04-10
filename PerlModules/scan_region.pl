#!/usr/bin/perl
use warnings;
use strict;
use Carp;
use Pod::Usage;
use Getopt::Long;

our $VERSION = 			'$Revision: 94 $';
our $LAST_CHANGED_DATE =	'$LastChangedDate: 2007-09-20 01:07:23 -0400 (Thu, 20 Sep 2007) $';

our ($verbose, $help, $man);
our ($queryfile, $dbfile);
our ($snp_flag, $score_threshold, $normscore_threshold, $overlap, $dbregion, $append, $condense_query, $mce_flag, $phastcons_flag, $evofold_flag, $tfbs_flag, 
	$wgrna_flag, $refgene_flag, $segdup_flag, $knowngene_flag, $test_flag,
	$autoexpand, $expandleft, $expandright, $expandmax, $bothside, $strand, $minoverlap, $kgxref, $reflink, $name2_flag, $quiet, $maxoverlap, 
	$maxquerydbratio, $minquerydbratio);
GetOptions('verbose'=>\$verbose, 'help'=>\$help, 'man|m'=>\$man, 'snp_flag'=>\$snp_flag, 'score_threshold=f'=>\$score_threshold, 
	'normscore_threshold=i'=>\$normscore_threshold, 'overlap'=>\$overlap, 'dbregion'=>\$dbregion,
	'append'=>\$append, 'condense_query'=>\$condense_query, 
	'mce_flag'=>\$mce_flag, 'phastcons_flag'=>\$phastcons_flag, 'evofold_flag'=>\$evofold_flag, 'tfbs_flag'=>\$tfbs_flag, 'wgrna_flag'=>\$wgrna_flag, 
	'refgene_flag'=>\$refgene_flag, 'segdup_flag'=>\$segdup_flag, 'knowngene_flag'=>\$knowngene_flag, 'test_flag'=>\$test_flag,
	'expandleft=s'=>\$expandleft, 'expandright=s'=>\$expandright, 'expandmax=s'=>\$expandmax, 'bothside'=>\$bothside,
	'strand=s'=>\$strand, 'minoverlap=f'=>\$minoverlap, 'kgxref=s'=>\$kgxref, 'reflink=s'=>\$reflink, 'name2_flag'=>\$name2_flag, 'quiet'=>\$quiet,
	'maxoverlap=f'=>\$maxoverlap, 'maxquerydbratio=f'=>\$maxquerydbratio, 'minquerydbratio=f'=>\$minquerydbratio) or pod2usage ();

$help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
$man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
@ARGV or pod2usage (-verbose=>0, -exitval=>1, -output=>\*STDOUT);
@ARGV == 2 or pod2usage ("Syntax error");

($queryfile, $dbfile) = @ARGV;

#preparing default parameters
$snp_flag and pod2usage ("Error in argument: the --snp_flag argument is in development");
$strand and $strand =~ m/^(plus|minus|forward|reverse)$/ || pod2usage ("Error in argument: the --strand argument can be set as only plus (forward) or minus (reverse)");
$minoverlap ||= 0;
$expandleft ||= 0; $expandright ||= 0; $expandmax ||= 0;
if ($expandmax) {
	$expandmax =~ s/k$/000/i; $expandmax =~ s/m$/000000/i;
	$expandmax =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandmax argument should be a positive integer (suffix of k and m are acceptable)");
}
if ($expandleft) {
	$expandleft =~ s/k$/000/i; $expandleft =~ s/m$/000000/i;
	$expandleft =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandleft argument should be a positive integer (suffix of k and m are acceptable)");
}
if ($expandright) {
	$expandright =~ s/k$/000/i; $expandright =~ s/m$/000000/i;
	$expandright =~ m/^\d+$/ or pod2usage ("Error in argument: the --expandright argument should be a positive integer (suffix of k and m are acceptable)");
}

if ($mce_flag or $evofold_flag or $tfbs_flag or $wgrna_flag or $segdup_flag) {
	my ($chr_loc, $count_query_bp) = readChrLoc ($queryfile);
	$condense_query and ($chr_loc, $count_query_bp) = condenseOverlap ($chr_loc);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in query to scan\n" and exit(0);	#there is no query regions in the query-location-file
	
	if ($mce_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'mce');
	} elsif ($evofold_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'evofold');
	} elsif ($tfbs_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'tfbs');
	} elsif ($wgrna_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'wgrna');
	} elsif ($segdup_flag) {
		scanUCSCRegion ($chr_loc, $dbfile, $count_query_bp, 'segdup');
	}
} elsif ($knowngene_flag or $refgene_flag) {
	if ($knowngene_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'knowngene', $kgxref, $name2_flag);
	} elsif ($refgene_flag) {
		scanUCSCGene ($queryfile, $dbfile, 0, 'refgene', $kgxref, $name2_flag);
	}
} elsif ($phastcons_flag) {							#special treatment for phsatcons file
	my ($chr_loc, $count_query_bp) = readChrLoc ($queryfile);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in query to scan\n" and exit(0);	#there is no query regions in the query-location-file
	print STDERR "NOTICE: the --phastcons flag requires condensing query regions into a set of non-overlapping regions\n";
	($chr_loc, $count_query_bp) = condenseOverlap ($chr_loc);
	scanPhastCons ($chr_loc, $dbfile, $count_query_bp);
} elsif ($test_flag) {
	scanRegion ($queryfile, $dbfile);
} else {
	scanRegion ($queryfile, $dbfile);
}





sub readChrLoc {
	my ($queryfile) = @_;
	my ($count_query_bp, $count_query_region, %chr_loc);	
	if ($queryfile eq 'stdin') {
		*QUERY = *STDIN;
	} else {
		open (QUERY, $queryfile) or confess "Error: cannot read from query file $queryfile: $!";		#!
	}
	while (<QUERY>) {
		/\S/ or next;												#skip empty lines
		s/[\r\n]+$//;												#get rid of newline characters
		#if ($_ =~ /chr:-/) {confess "No Results Here!\n";}	#JOE
		m/^(chr)?(\w+):(\d+)(\-(\d+))?(.*)/ or confess "Error: invalid record in chromosome file: <$_>";	#JOE	#might contain only one locatoin (SNP) or a start-end (for chromosome regions)
		my ($chr, $loc_start, $loc_end, $other_info) = ($2, $3, $5, $6);
		my $switch = 0;												#switch start and end position
		
		defined $loc_end or $loc_end = $loc_start;
		defined $other_info or $other_info = '';
		if ($loc_start > $loc_end) {										#in-del mutation such as rs34899723 in human
			$switch = 1;
			($loc_start, $loc_end) = ($loc_end, $loc_start);
		}
		push @{$chr_loc{$chr}}, [$loc_start, $loc_end, $other_info, $switch];
	
		$count_query_bp += ($loc_end-$loc_start+1);
		$count_query_region++;
	}
	for my $chr (keys %chr_loc) {
		@{$chr_loc{$chr}} = sort {$a->[0] <=> $b->[0]} @{$chr_loc{$chr}};				#sort the chr_location by the start position for each chr
	}
	$verbose and print STDERR "NOTICE: A total of $count_query_bp (possibly overlapped) base pairs found in $queryfile in ${\(scalar keys %chr_loc)} chromosomes: " , join (" ", sort keys %chr_loc), "\n";
	return (\%chr_loc, $count_query_bp, $count_query_region);
}

sub scanPhastCons {
	my ($chr_loc, $dbfile, $count_query_bp) = @_;
	my ($chr, $start, $current_loc, $step, $qloc, $qpointer);
	my ($count_db_bp, $count_db_bp_found) = (0, 0);
	my (%chr_finish, %chr_pointer);
	
	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		if (m/^fixedStep chrom=chr(\w+) start=(\d+) step=(\d+)$/) {
			#$verbose and print STDERR "NOTICE: next header $_";
			($chr, $start, $current_loc, $step) = ($1, $2-1, $2-1, $3);	#the next line is the real start (so first minus 1 to get the real start)
			$qloc = $chr_loc->{$chr} || undef;
			$qpointer = $chr_pointer{$chr} || 0;
			$step == 1 or confess "Error: multiple steps found in phastcons file";
		} else {
			length ($_) == 6 or die "Error: invalid record in db file: <$_>";
			$current_loc++;
			$count_db_bp++;

			$verbose and $count_db_bp =~ m/000000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_bp base of dbfile $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
			$qloc or next;					#the $qloc is the query location array for this chromosome
			$chr_finish{$chr} and next;

			for my $i ($qpointer .. @$qloc-1) {
				my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};
				if ($qend < $current_loc) {
					$chr_pointer{$chr} = $i;
					$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
					next;
				} elsif ($qstart > $current_loc) {
					last;
				} else {
					$count_db_bp_found++;
					print $overlap?"$chr:$current_loc-$current_loc$qinfo":"$chr:$qstart-$qend$qinfo";
					print $append?"\t$_":"\n";
				}
			}
		}
	}
}

sub scanRegion {
	my ($queryfile, $dbfile) = @_;
	my ($db_chr_loc, $count_db_bp, $count_db_region) = readChrLoc ($dbfile);
	($db_chr_loc, $count_db_bp, $count_db_region) = condenseOverlap ($db_chr_loc);

	my ($chr_loc, $count_query_bp, $count_query_region) = readChrLoc ($queryfile);
	$count_query_bp or print STDERR "WARNING: there is NOTHING in query to scan\n" and exit(0);	#there is no query regions in the query-location-file
	$condense_query and ($chr_loc, $count_query_bp, $count_query_region) = condenseOverlap ($chr_loc);

	my ($count_db_bp_found, %chr_finish, %chr_pointer) = (0);
	for my $chr (sort sortChr keys %$db_chr_loc) {
		for my $nextregion (@{$db_chr_loc->{$chr}}) {
			my ($start, $end) = @$nextregion;
			
			my $qloc = $chr_loc->{$chr} or next;				#skip this chromosome because query does not have anything on it
			$chr_finish{$chr} and next;					#skip this chromosome because it has been finished processing by query
	
			my $pointer = $chr_pointer{$chr} || 0;				#current pointer for this chromosome
			for my $i ($pointer .. @$qloc-1) {
				my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};
	
				if ($maxquerydbratio and ($qend-$qstart+1)/($end-$start+1) > $maxquerydbratio) {
					next;			#maximum allowable query/db size ratio difference (query is too big to claim overlap/similarity with db)
				}
				if ($minquerydbratio and ($qend-$qstart+1)/($end-$start+1) < $minquerydbratio) {
					next;			#minimum allowable query/db size ratio difference (query is too small to claim overlap/similarity with db)
				}
	
				if ($qend < $start) {					#query end is before the start position of this db region
					#db:            <------------------------->
					#query: <--->
					$chr_pointer{$chr} = $i;			#move the pointer only when qend<start (I thought this for a long time)
					$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
					next;
				} elsif ($qend <= $end) {				#move pointer to this region? (probably NOT, especially when one query overlap with another)
					if ($qstart >= $start) {			#query contained completely within db region
						#db:      <-------------------------->
						#query:       <------------------>
						$count_db_bp_found += ($qend-$qstart+1);
						if ($minoverlap) {
							1;				#query is contained within hit
						}

						if ($dbregion) {
							print "chr$chr:$start-$end";	#there is no qinfo for db region
						} else {
							print "chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					} else {					#query overlap but upstream of db region
						#db:       <------------------------->
						#query: <---------------------->
						$count_db_bp_found += ($qend-$start+1);
						if ($minoverlap) {
							if (($qend-$start+1)/($qend-$qstart+1) < $minoverlap and ($qend-$start+1)/($end-$start+1) < $minoverlap) {
								next;
							}
						}
						
						if ($dbregion) {
							print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$start-$end";
						} else {
							print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					}
				} elsif ($qstart <= $end) {
					if ($qstart >= $start) {			#query overlap but downstream of db region
						#db:      <------------------------>
						#query:        <----------------------->
						$count_db_bp_found += ($end-$qstart+1);
						if ($minoverlap) {
							if (($end-$qstart+1)/($qend-$qstart+1) < $minoverlap and ($end-$qstart+1)/($end-$start+1) < $minoverlap) {
								next;
							}
						}
						if ($dbregion ) {
							print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$start-$end";
						} else {
							print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					} else {					#db region completely contained within query
						#db:      <------------------------->
						#query: <------------------------------>
						$count_db_bp_found += ($end-$start+1);
						if ($minoverlap) {
							1;
						}
						if ($dbregion) {
							print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$start-$end";
						} else {
							print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
						}
						print "\n";
					}
				} else {
					#db:      <---------->
					#query:       	         <----------------------->
					last;						#should examine next db region
				}
			}
		}
		$verbose and $count_db_region =~ m/00000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_region line of template file $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
	}
	$quiet or print STDERR "NOTICE: Total of $count_db_region template chr_regions (total of $count_db_bp) examined and found $count_db_bp_found / $count_query_bp (${\($count_db_bp_found / $count_query_bp)}) interested base pairs that belong to these template regions\n";

}

#this subroutine is used to scan various files downloaded from UCSC table browser (usually tab-delimited files). Many of them have the same or similar file formats
sub scanUCSCRegion {
	my ($chr_loc, $dbfile, $count_query_bp, $dbtype) = @_;
	my (%chr_warning, %chr_finish, %chr_pointer);			#warning for no chr information in query file, chr in query has been finished processing, pointer to current query chr info
	my ($count_db_bp, $count_db_region, $count_db_bp_found, $count_db_region_found) = (0, 0, 0, 0);
	my @record;
	my ($chr, $start, $end, $score, $normscore);			#score is defined as LOD for MCE, Z for TFBS and conf for Evofold; normscore is between 0-1000
	
	
	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		s/[\r\n]+$//;							#deleting the newline characters
		@record = split (/\t/, $_);
		$record[0] eq '#bin' and next;					#comment line
		if ($dbtype eq 'mce') {
			@record == 6 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[4], $record[5]);
			$score =~ s/^lod=// or confess "Error: invalid lod score designation (no 'lod=' found) in dbfile $dbfile: <$_>";
		} elsif ($dbtype eq 'evofold') {
			@record == 10 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[5], $record[5]);
		} elsif ($dbtype eq 'tfbs') {
			@record == 8 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[7], $record[5]);
		} elsif ($dbtype eq 'wgrna') {
			@record == 10 or confess "Error: invalid record in dbfile $dbfile (expecting 10 fields in evofold file): <$_>";
			($chr, $start, $end, $score, $normscore) = ($record[1], $record[2], $record[3], $record[5], $record[5]);
		} elsif ($dbtype eq 'segdup') {
			@record == 30 or confess "Error; invalid record in dbfile $dbfile (expecting 30 fields in segdup file): <$_>";
			($chr, $start, $end, $score, $normscore) = @record[1, 2, 3, 5, 5];
		} else {
			confess "Error: the dbtype $dbtype cannot be handled by the current program";
		}

		if ($dbtype =~ m/^(mce|evofold|tfbs|wgrna|segdup)$/) {
			$score_threshold and $score < $score_threshold and next;			#if --score_threshold is set, the low scoring segment will be skipped
			$normscore_threshold and $normscore < $normscore_threshold and next;		#if --normscore_threshold is set, the low scoring segment will be skipped
			$start++;									#due to the zero-opening coordinate system in UCSC
		}

		$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>";
		$count_db_region++;
		$count_db_bp += ($end-$start+1);
		
		my $qloc = $chr_loc->{$chr} or next;				#skip this chromosome because query does not have anything on it
		$chr_finish{$chr} and next;					#skip this chromosome because it has been finished processing by query

		my $pointer = $chr_pointer{$chr} || 0;				#current pointer for this chromosome
		for my $i ($pointer .. @$qloc-1) {
			my ($qstart, $qend, $qinfo, $qswitch) = @{$qloc->[$i]};

			if ($qend < $start) {					#query end is before the start position of this db region
				#db:            <------------------------->
				#query: <--->
				$chr_pointer{$chr} = $i;			#move the pointer only when qend<start (I thought this for a long time)
				$i == @$qloc-1 and $chr_finish{$chr} = 1;	#this chromosome is DONE
				next;
			} elsif ($qend <= $end) {				#move pointer to this region? (probably NOT, especially when one query overlap with another)
				if ($qstart >= $start) {			#query contained completely within db region
					#db:      <-------------------------->
					#query:       <------------------>
					$count_db_bp_found += ($qend-$qstart+1);
					if ($dbregion) {
						print "chr$chr:$start-$end";	#there is no qinfo for db region
					} else {
						print "chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				} else {					#query overlap but upstream of db region
					#db:       <------------------------->
					#query: <---------------------->
					$count_db_bp_found += ($qend-$start+1);
					if ($minoverlap) {
						if (($qend-$start+1)/($qend-$qstart+1) < $minoverlap and ($qend-$start+1)/($end-$start+1) < $minoverlap) {
							next;
						}
					}
					if ($dbregion) {
						print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$start-$qend$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				}
			} elsif ($qstart <= $end) {
				if ($qstart >= $start) {			#query overlap but downstream of db region
					#db:      <------------------------>
					#query:        <----------------------->
					$count_db_bp_found += ($end-$qstart+1);
					if ($minoverlap) {
						if (($end-$qstart+1)/($qend-$qstart+1) < $minoverlap and ($end-$qstart+1)/($end-$start+1) < $minoverlap) {
							next;
						}
					}
					if ($dbregion ) {
						print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$qstart-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				} else {					#db region completely contained within query
					#db:      <------------------------->
					#query: <------------------------------>
					$count_db_bp_found += ($end-$start+1);
					if ($dbregion) {
						print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$start-$end";
					} else {
						print $overlap?"chr$chr:$start-$end$qinfo":"chr$chr:$qstart-$qend$qinfo";
					}
					$append and print "\t$score\t$normscore";
					print "\n";
				}
			} else {
				last;						#should examine next db region
			}
		}
		$verbose and $count_db_region =~ m/00000$/ and print STDERR "NOTICE: Processing chr $chr in $count_db_region line of template file $dbfile with $count_db_bp_found/$count_query_bp interested base pairs found in overlapping regions\n";
	}
	$quiet or print STDERR "NOTICE: Total of $count_db_region template chr_regions (total of $count_db_bp) examined and found $count_db_bp_found / $count_query_bp (${\($count_db_bp_found / $count_query_bp)}) interested base pairs that belong to these template regions\n";

}

#this subroutine scan a queryfile against a knownGene file or refGene file from UCSC to identify the nearby genes for each query_region in queryfile
#instead of scanning each line in the template file, I decided to just read everything in memory, then process each query sequentially. this is conceptually different and opposite the scanDBRegion subroutine
sub scanUCSCGene {
	my ($queryfile, $dbfile, $count_query_bp, $dbtype, $kgxref, $name2_flag) = @_;
	my (%chr_warning, %chr_finish, %chr_pointer);				#warning for no chr information in query file, chr in query has been finished processing, pointer to current query chr info
	my ($count_db_bp, $count_db_region) = (0, 0);
	my (%db_chr_loc);
	my ($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend);
	my ($gene_xref);

	open (DB, $dbfile) or confess "Error: cannot read from database file $dbfile: $!";
	while (<DB>) {
		m/^#/ and next;							#comment line
		s/[\r\n]+$//;							#deleting the newline characters
		my @record = split (/\t/, $_);

		if ($dbtype eq 'knowngene') {
			@record == 12 or confess "Error: invalid record in dbfile $dbfile (expecting 12 fields in knownGene file): <$_>";
			($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend) = @record[0..6];				#human hg17
		} elsif ($dbtype eq 'refgene') {
			if (@record == 16) {
				($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend) = @record[1..7];		#human hg18, mouse
				$name2_flag and $name = $record[12];
			} elsif (@record == 10) {
				($name, $chr, $dbstrand, $start, $end, $cdsstart, $cdsend) = @record[0..6];			#rat, human hg17
				$name2_flag and pod2usage ("Error in argument: the --name2_flag argument can be used only for refGene table with name2 annotations (such as those in hg18 and mm8 database)");
			} else {
				confess "Error: invalid record in template-location-file $dbfile (expecting 16 or 10 tab-delimited fields in refGene file): <$_>";
			}
		}
		
		if (defined $strand) {						#only consider a particular strand (by default both strand will be considered)
			$strand eq 'plus' || $strand eq 'forward' and $dbstrand eq '+' || next;
			$strand eq 'minus' || $strand eq 'reverse' and $dbstrand eq '-' || next;
		}

		$start++;							#due to the zero-opening coordinate system in UCSC
		$chr =~ s/^chr// or confess "Error: invalid chromosome designation (no 'chr' found) in dbfile $dbfile: <$_>";
		$count_db_region++;
		$count_db_bp += ($end-$start+1);
		push @{$db_chr_loc{$chr}}, [$start, $end, $name];
	}
	for my $chr (sort keys %db_chr_loc) {					#sort dbregion to make sure that smaller regions occur before bigger ones
		@{$db_chr_loc{$chr}} = sort {$a->[0] <=> $b->[0] or $a->[1] <=> $b->[1]} @{$db_chr_loc{$chr}};
	}
	
	#now scan each query_loc against all sorted db entry in the template file
	#by default, any refgene overlapping with query will be printed out
	#if no refgene overlap, then the closest gene (as specified by --expandmax) will be printed out
	#db: ---------          -------------------                  -----        ---------------
	#        --------    -------------------------    ---          -------------                ----------
	#query:    <------->                     <-------------------------->
	#                       <--------->                                             <-->
	#                                                       <-->

	if ($kgxref) {		#in the output, rather than using gene identifier (name), use the corresponding gene symbol in the known gene cross-reference file
		$gene_xref = readKgXref ($kgxref);
	} elsif ($reflink) {
		$gene_xref = readRefLink ($reflink);
	}
		
	if ($queryfile eq 'stdin') {
		*QUERY = *STDIN;
	} else {
		open (QUERY, $queryfile) or confess "Error: cannot read from queryfile $queryfile: $!";
	}
	while (<QUERY>) {
		/\S/ or next;												#skip empty lines
		s/[\r\n]+$//;												#get rid of newline characters
		#if ($_ =~ /chr:-/) {print "No Results Here!\n";}	#JOE	
		m/^(chr)?(\w+):(\d+)(\-(\d+))?(.*)/ or  confess "Error: invalid record in chromosome file: <$_>";	#might contain only one locatoin (SNP) or a start-end (for chromosome regions)
		my ($chr, $qstart, $qend, $qinfo, $qswitch) = ($2, $3, $5||$3, $6||'');
		if ($qstart > $qend) {
			$qswitch = 1;
			($qstart, $qend) = ($qend, $qstart);
		}
		
		my (@closest_gene, $closest_dist, @overlap);
		my ($qstart_expand, $qend_expand) = ($qstart - $expandmax, $qend + $expandmax);
		$expandleft and $qstart_expand = $qstart - $expandleft;
		$expandright and $qend_expand = $qend + $expandright;
		
		my $dbloc = $db_chr_loc{$chr};
		if (not $dbloc) {								#for example, when chr is 22_random, there is NO genes there
			print "NOT_FOUND\t0\t$chr:$qstart-$qend\t$qinfo\n";
			next;
		}
		for my $j (0 .. @$dbloc-1) {
			my ($dbstart, $dbend, $dbname) =@{$dbloc->[$j]};
			$qstart_expand > $dbend and next;					#the dbend has not reached expanded query yet
			if ($qend >= $dbstart and $qstart <= $dbend) {				#overlap found between query and db!
				push @overlap, $dbname;
			}
			if (not @overlap and $expandleft || $expandright || $expandmax) {
				if ($qend_expand >= $dbstart and $qstart_expand <= $dbend) {	#overlap found between expanded query and db
					my $dist;
					if ($qend <= $dbstart) {
						$dist = $dbstart-$qend;
					} else {
						$dist = $qstart-$dbend;
					}

					if ($closest_dist and $dist < $closest_dist) {
						$closest_dist = $dist;
						@closest_gene = $dbname;
					} elsif (not defined $closest_dist or $dist == $closest_dist) {
						$closest_dist = $dist;
						push @closest_gene, $dbname;

					}
				}
			}
			$qend_expand < $dbstart and last;	#the dbstart has passed expanded query
		}
		
		if ($kgxref or $reflink) {
			@overlap = map {$gene_xref->{$_} || $_} @overlap;
			@closest_gene = map {$gene_xref->{$_} || $_} @closest_gene;
		}
			
		
		my %seen;
		@overlap = sort grep {!$seen{$_}++} @overlap;
		@closest_gene = sort grep {!$seen{$_}++} @closest_gene;
		
		if (@overlap) {
			print "chr$chr:$qstart-$qend$qinfo\t", join(",", @overlap), "\t0\n";
		} elsif (@closest_gene) {
			print "chr$chr:$qstart-$qend$qinfo\t", join(",", @closest_gene), "\t$closest_dist\n";
		} else {
			print "chr$chr:$qstart-$qend$qinfo\tNOT_FOUND\tNOT_FOUND\n";
		}
	}
}

#sometimes two exons enclose with each other, sometimes two regions locate in two strands so have overlaps, so it is important to eliminate overlap positions
#this subroutine reduce overlaps, when feeded with a series of chromosome regions in two arrays, one for start position and one for end position
sub condenseOverlap {
	my ($chr_loc) = @_;
	my ($length, $newlength, $newcount) = (0, 0, 0);
	for my $chr (keys %$chr_loc) {
		my @loc = @{$chr_loc->{$chr}};
		my @newloc;
		my $pointer = 0;
		my ($prestart, $preend) = ($loc[0]->[0], $loc[0]->[1]);
		$length += $loc[0]->[1]-$loc[0]->[0]+1;
		
		while (1) {
			$pointer++;					#process next segment
			$pointer == @loc and last;			#finish processing the loc array
			$length += ($loc[$pointer]->[1]-$loc[$pointer]->[0]+1);
			if ($loc[$pointer]->[0] <= $preend) {		#start of next element less than end of previous element
				if ($loc[$pointer]->[1] >= $preend) {	#make sure the next element is not contained within the previous element
					$preend = $loc[$pointer]->[1];
				}
			} else {
				push @newloc, [$prestart, $preend, '', ''];
				$newlength += ($preend-$prestart+1);
				($prestart, $preend) = ($loc[$pointer]->[0], $loc[$pointer]->[1]);
			}
		}
		push @newloc, [$prestart, $preend, '', ''];		#process the last element
		$newlength += ($preend-$prestart+1);
		$newcount += @newloc;

		$chr_loc->{$chr} = \@newloc;
	}
	$quiet or print STDERR "NOTICE: After cleaning query, the length of query changes from $length to $newlength\n";
	return ($chr_loc, $newlength, $newcount);
}

sub sortChr {
	my ($a, $b) = @_;		#this line is required for functioning properly and I do not know why
	if ($a =~ m/^(\d+)$/) {
		if ($b =~ m/^(\d+)$/) {
			return $a<=>$b;
		} else {
			return -1;
		}
	} elsif ($b =~ m/^(\d+)$/) {
		return 1;
	} else {
		return $a cmp $b;
	}
}

sub readKgXref {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (XREF, $inputfile) or confess "Error: cannot read from kgxref file $inputfile: $!";
	while (<XREF>) {
		m/^#/ and next;
		my @record = split (/\t/, $_);
		if ($gene_xref{$record[0]}) {			#BC003168 occur twice in kgxref file (OSBPL10, BC003168)
			if ($gene_xref{$record[0]} =~ m/^(BC|AK)\d+$/) {
				$gene_xref{$record[0]} = $record[4];
			}
		} else {
			$gene_xref{$record[0]} = $record[4];
		}
	}
	close (XREF);
	return (\%gene_xref);
}

sub readRefLink {
	my ($inputfile) = @_;
	my (%gene_xref);
	open (REFLINK, $inputfile) or confess "Error: cannot read from refLink file $inputfile: $!";
	while (<REFLINK>) {
		m/^#/ and next;
		my @record = split (/\t/, $_);
		$gene_xref{$record[2]} = $record[0];
	}
	close (REFLINK);
	return (\%gene_xref);
}


=head1 SYNOPSIS

 scan_region.pl [arguments] <query-location-file> <template-location-file>

 Optional arguments:
        -h, --help                      print help message
        -m, --man                       print complete documentation
        -v, --verbose                   use verbose output
            --snp_flag			query is SNP (same start and stop position)
            --overlap			print overlapped portion of region only
            --dbregion			print db region (default is to print query region)
            --append			append extra information from annotation file to output
            --condense_query		condense and eliminate overlapping regions in query
            --minoverlap <float>	minimum portion of overlap to decide concordance of query and db
            --score_threshold <float>	score threshold for inclusion in output
            --normscore_threshold <float>	normalized score threshold for inclusion in output

            --mce_flag			dbfile is UCSC MCE annotation file
            --evofold_flag		dbfile is UCSC EvoFold annotation file
            --tfbs_flag			dbfile is UCSC TFBS annotation file
            --wgrna_flag		dbfile is UCSC wgRna annotation file
            --knowngene_flag		dbfile is UCSC knownGene annotation file
            --refgene_flag		dbfile is UCSC refGene annotation file
            --phastcons_flag		dbfile is UCSC phastcons conservation score file

            --expandleft <int>		expand left side of query regions (overwrite --expandmax)
            --expandright <int>		expand right side of query regions (overwrite --expandmax)
            --expandmax <int>		size of maximum expansion for query region to find overlap
            --bothside			output genes in both side of expansion (not implemented yet!)

            --kgxref <file>		use gene symbol in known gene xref file in output
            --reflink <file>		use gene symbol in refGene link file in output
            --name2_flag		use name2 annotation in refGene file in output
            --quiet			suppress printing progress messages

 Function: scan genomic regions in a query-location-file against a template-
 location-file (dbfile), which contain various sequence features for genomic regions

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--verbose>

use verbose output.

=item B<--snp_flag>

THIS ARGUMENT IS NOT IMPLEMENTED YET. specify that the query file contains 
position information for SNPs that contain only one base pair. This option 
significantly speeds up the program, since only one point in genome, rather than 
a region in genome, will be probed for each input line in query file.

=item B<--overlap>

instead of printing the query region, only print the overlapped portion of the 
query region and template region.

=item B<--dbregion>

print the region in tempalte file, rather than query file, when an overlapped 
hit is found

=item B<--append>

append the score and normscore for the overlapped template region to the output 
for template files downloaded as UCSC tables.

=item B<--condense_query>

condense overlapped regions in the query file into non-overlapped regions. When 
this argument is set, the annotation for each query (the strings after the 
chromosome location in each line of the query file) will not in the output.

=item B<--minoverlap>

the minimum percentage of overlap between query and template chromosome region 
to infer an overlapping event. By default, even a single base pair overlap will 
be considered as overlap, but a lot of times people prefer to use something like 
0.5 (50% overlap) to make sure that the query and template regions have high 
concordance.

=item B<--score_threshold>

specify the score threshold in the template file to include in the search for 
overlaps. This argument is file format dependent.

=item B<--normscore_threshold>

specify the normalized score threshold in the template file to include in the 
search for overlaps. This argument is file format dependent.

=item B<--mce_flag>

specify that the template file is in MCE format from UCSC genome browser

=item B<--evofold_flag>

specify that the template file is in EvoFold format from UCSC genome browser

=item B<--tfbs_flag>

specify that the template file is in TFBS format from UCSC genome browser

=item B<--wgrna_flag>

specify that the template file is in WGRNA format from UCSC genome browser

=item B<--knowngene_flag>

specify that the template file is in knownGene format from UCSC genome browser

=item B<--refgene_flag>

specify that the template file is in refGene format from UCSC genome browser

=item B<--phastcons_flag>

specify that the template file is in phastcons format from UCSC genome browser

=item B<--expandleft>

expand the query region on the left side (5' in forward strand, 3' in reverse 
strand) to find overlap (used in conjunction with --refgene or --knowngene 
argument)

=item B<--expandright>

expand the query region on the right side (3' in forward strand, 5' in reverse 
strand) to find overlap (used in conjunction with --refgene or --knowngene 
argument)

=item B<--expandmax>

maximum expansion size of the query region on both side to find at least one 
overlap (used in junction with --refgene or --knowngene argument)

=item B<--bothside>

THIS ARGUMENT IS NOT IMPLEMENTED YET. process and display gene information on 
both side of the expansion. (By default only the closet gene from either side is 
in output). The genes from two side are separated by ";" in tab-delimited 
output.

=item B<--kgxref>

specify a cross-reference file for the knownGene track in UCSC genome browser, 
so that in the output, the gene identifier (gene name or refseq id) are replaced 
by the gene symbol specified in the kgxref file. (If not found in the kgxref 
file, the gene identifiers are still used)

=item B<--name2_flag>

this argument is used in conjunction with the --refgene argument, to specify 
that the alternative gene symbol in the "name2" field in the refGene file be 
printed in the output.

=back

=head1 DESCRIPTION

This program is used to retrieve shared chromosome region from a given location 
file based on another template location file that has been sorted by chromosome 
and start location. Both the query-location-file and the template-location-file 
contains one location per line. THE TEMPLATE FILE MUST BE SORTED BY CHROMOSOME 
AND START POSITIONS!

For example, a sample query-location-file is shown below (The "chr" prefix is 
optional for each line):

	chr3:2000-349990
	chr19:32333-52333

If the query-location-file is for SNPs, then it is okay to use only the start 
location in the line.

A sample template-location-file can be the same format as the query-location-
file, or can be a tab-delimted file (with three columns: chr, start, end), or it 
can be generated by various tracks from UCSC genome browser. Each of the file 
format is described in detail below.

=over 8

=item * B<plain format>

This is the simplest format: either in the form of "chr3:2000-349990", or in the 
form of "chr3\t2000\t349000", where "\t" indicates a tab character. Each region 
occupy one single line.

=item * B<General overview of UCSC genome browser tables>

When processing template files as UCSC tables, this program works by first 
reading all information from the query file and store them in memory, then scan 
each line of the second file (the template file) to find overlaps. (I use the 
name "template file" to refer to the various annotation tables in the UCSC 
genome browser.) THEREFORE, THE TEMPLATE FILE CANNOT CONTAIN TWO LOCATIONS THAT 
SHARE OVERLAPS. IN ADDITION, THE TEMPLATE FILE MUST BE SORTED FIRST BY 
CHROMOSOME THEN BY START LOCATION. The chromosome can be sorted numerically or 
alphabetically (since there will usually be X, Y and MT chromosomes), but start 
location should be always sorted numerically.

Generally speaking, I recommended saving the template file with the same file 
name as the original database file, but prefixing with the version and organism 
of the corresponding genome. For example, for human genome, the template file 
name should start with hg17_ or hg18_, to eliminate confusions when dealing with 
different genome builds. For example, after downloading the file evofold.txt.gz 
for human May 2004 assembly from the UCSC genome browser database, you can use 
the sort command to rename the file as hg17_evofold.sorted (see details below).

=item * B<MCE format>

The --mce_flag argument specifies that template file is in MCE format from the 
UCSC genome browser. There are two ways to download the file: You can download a 
MCE template file using the UCSC Table Browser (http://www.genome.ucsc.edu/cgi-
bin/hgTables). Choose group as "comparative genomics", track as "most 
conserved", table as "phastConsElements17way" (table names may differ for other 
species other than human). After getting the output file (save it as something 
like mm8_phastcons), then sort the file by commands such as 

	sort -k 2,2 -k 3,3n mm8_phastcons > mm8_phastcons.sorted

Alternatively, you can also directly download the MCE annotation file as 
http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/phastConsElements17way.t 
xt.gz. Note that when directly downloading the annotation file, there is no 
header line (the first line in the file) that specify the name of each tab-
delimited column in the file. You can then use

	gunzip phastConsElements17way.txt.gz

to unzip the file and generate a text file phastConsElements17way.txt, then run

	sort -k 2,2 -k 3,3n phastConsElements17way.txt > hg17_phastcons.sorted

I personally prefer using the Table Browser to download data, because it is much 
easier to use the "group", "track" and "table" option to find the desired 
content in Table Browser. However, the "wget" command is a much faster way if 
you know exactly what table you want.

A sample sorted MCE file is shown below:

	585     chr1    1865    1943    lod=31  307
	585     chr1    2039    2096    lod=104 444
	585     chr1    2473    2566    lod=179 505
	585     chr1    2873    2917    lod=107 447
	585     chr1    3081    3135    lod=58  378

The --score_threshold argument operates on the fifth field in the line (such as 
"lod=31"), and the -- normscore_threshold argumet operates on the sixth field 
(such as 307).

To scan the queryfile against a MCE file:

	scan_region.pl queryfile hg18_phastcons.sorted -mce

The output will be multiple identical lines, each representing one hit in query with the 
MCE file. The first a few lines are:

	19:32333-52333
	19:32333-52333
	19:32333-52333
	19:32333-52333
	19:32333-52333

If only the region in MCE file that overlap with query is desired, you can use:

	scan_region.pl queryfile hg18_phastcons.sorted -mce --dbregion

The first a few lines in output file are:

	19:35715-35751
	19:37371-37466
	19:38216-38406
	19:38437-38504
	19:38551-38678

If only the overlapped region is desired, you can use:

	scan_region.pl queryfile hg18_phastcons.sorted -mce -overlap

The output will be the overlapped chromosome regions between query and MCE file. 
The first a few output lines are:

	19:35715-35751
	19:37371-37466
	19:38216-38406
	19:38437-38504
	19:38551-38678

In this case, the query region is very big, so that the overlapped regions are 
usually identical to the dbregion (region in the MCE file).

=item * B<EvoFold format>

The EvoFold file from UCSC genome browser contains genomic segments that are 
predicted to retain stable RNA structural folds in different species. You can 
download a EvoFold template file using the UCSC Table Brower, by choosing group 
as "Gene and gene prediction track", track as "EvoFold", and table as "evofold". 
After downloading the file and save it as something like hg18_evofold, use the 
command to sort the chromosome positions:

	sort -k 2,2 -k 3,3n hg18_evofold > hg18_evofold.sorted

Alternatively, you can directly download the EvoFold annotation file using the 
command below, and then sort the resulting uncompressed file:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/evofold.txt.gz

The first a few lines of a sorted EvoFold table is shown below:

	591     chr1    886773  886798  608_0_-_96      96      -       25      (((((.....(((...))).)))))       0.97,0.98,0.99,0.99,0.9,0.78,0.87,0.99,0.98,0.98,0.21,0.22,0.11,0.99,0.95,0.91,0.11,0.22,0.21,0.98,0.9,0.99,0.99,0.98,0.97
	591     chr1    888417  888435  617_0_+_156     156     +       18      ((((((......))))))      0.97,0.99,0.99,1.0,0.99,0.92,1.0,1.0,0.99,1.0,1.0,1.0,0.92,0.99,1.0,0.99,0.99,0.97

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "96").

=item * B<TFBS format>

The TFBS file from UCSC genome browser contains predicted transcription factor 
binding sites (TFBS), based on positional weight matrices (PWM) from the TRANFAC 
database. You can download a TFBS template file using the UCSC Table Browser, by 
choosing group as "Expression and Regulation", track as "TFBS Conserved" and 
table as "tfbsConsSites". After downloading the file and save it as something 
like hg18_tfbs, use the following command to sort the chromosome positions:

	sort -k 2,2 -k 3,3n tfbsConsSites > hg18_tfbsConsSites.sorted

Alternatively, you can download the annotation file directly using the 
command below, and then sort the resulting uncompressed file:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/tfbsConsSites.txt.gz

The first a few lines of a sorted TFBS table is shown below:

	585     chr1    1832    1848    V$ARP1_01       818     -       2.01
	585     chr1    3211    3232    V$NRSF_01       749     -       2.19
	585     chr1    3603    3623    V$YY1_02        790     -       1.93
	585     chr1    3949    3956    V$NKX25_01      1000    -       2.48
	585     chr1    3996    4014    V$CART1_01      798     -       1.75

The --score_threshold argument operates on the eighth field in the line (such as 
"2.01"), and the --normscore_threshold argumet operates on the sixth field 
(such as 818).

=item * B<wgRna format>

The wgRna table from UCSC genome browser contains microRNA and small nucleolar 
RNA information. You can download a wgRna template file using the UCSC Table 
Browser, by choosing group as "Gene and gene prediction tracks", track as 
"sno/miRNA", and table as "wgRna". After downloading the file and save it as 
something like hg18_wgrna, use the following command to sort the chromosome positions:

	sort -k 2,2 -k 3,3n hg18_wgrna > hg18_wgrna.sorted


Alternatively, you can download the annotation file directly using the 
command below, and then sort the resulting uncompressed file:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/database/wgRna.txt.gz

The first a few lines of a sorted wgRna table is shown below:

	593     chr1    1092346 1092441 hsa-mir-200b    960     +       1092402 1092425 miRna
	593     chr1    1093105 1093195 hsa-mir-200a    960     +       1093120 1093142 miRna
	593     chr1    1093105 1093195 hsa-mir-200a    960     +       1093158 1093180 miRna
	593     chr1    1094247 1094330 hsa-mir-429     960     +       1094297 1094319 miRna
	611     chr1    3467118 3467214 hsa-mir-551a    480     -       3467133 3467154 miRna

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "960").

=item * B<segdup format>

The segDup table from UCSC genome browser contains chromosome regions with 
segmental duplications, as well as the "target" regions that match these 
duplications. You can download a segdup template file using the UCSC Table 
Browser, by choosing group as "Variation and Repeats", track as "Segmental 
Dups", and table as "genomicSuperDups". After downloading the file and save it 
as something like hg17_segdup, use the following command to sort the chromosome 
locations:

	sort -k 2,2 -k 3,3n hg17_segdup > hg17_segdup.sorted

Alternatively, you can download the annotation file directly using the 
command below, and then sort the resulting uncompressed file:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/genomicSuperDups.txt.gz

The first a few lines of a sorted segdup table is shown below:

	585     chr1    465     30596   No.1139,chr2:114046528  570211  -       chr2    114046528       114076216       243018229       1139    1000    N/A     Filtered        N/A     N/A     build35/align_both/0012//both064929     30176   43      531     29645   29282   363     128     235     0.987755        0.986324        0.012346        0.0123574
	585     chr1    486     30596   No.2251,chr9:844        582086  +       chr9    844     30515   138429268       2251    1000    N/A     Filtered        N/A     N/A     build35/align_both/0013//both065185     30137   14      491     29646   29474   172     64      108     0.994198        0.993729        0.00582435      0.00582657

The --score_threshold and --normscore_threshold arguments both operate on the 
sixth field in the line (such as "570211").

=item * B<knownGene format>

The knownGene table from UCSC genome browser contains gene identifier (name), 
location, exon count and location, Swiss-Prot identifiers. It does not contain 
the gene symbol information, so an additional step is necessary to interrogate 
the kgXref table to find the gene symbol for each gene identifier. You can 
download a knownGene template file using the UCSC Table Browser, by choosing 
group as "Gene and gene prediction tracks", track as "Known Genes", table as 
"knownGene", and save the file as something like hg17_knowngene. You can 
download the kgXref file using the table "kgXref" and save it as something like 
hg17_kgxref. There is no need to sort the knownGene file or the kgXref file.

Alternatively, you can directly download the knownGene annotation file by:

	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg17/database/knownGene.txt.gz

Unless the --condense_query argument is set, the output order is the same as the 
query locations in the query file. Those query that does not match any knownGene 
will still be printed, although "NOT_FOUND" will show as the first 2 columns in 
the output.

To scan a query file against the knownGene table, use:

	scan_region.pl --knowngene queryfile hg17_knowngene

The output is:

	NM_006614       0       3:2000-349990
	NOT_FOUND       NOT_FOUND       19:32333-52333

The generated output line contains gene identifiers and their distances (should 
be ZERO) that overlap with query location.

To use gene symbols in the generated output, use:

	scan_region.pl --knowngene --kgxref hg17_kgxref queryfile hg17_knowngene

The output is:

	CHL1    0       3:2000-349990
	NOT_FOUND       NOT_FOUND       19:32333-52333

To identify surrounding genes to the query location, use:

	scan_region.pl --knowngene -expandmax 1m queryfile hg17_knowngene

The output is:

	NM_006614       0       3:2000-349990
	AF346307        3643    19:32333-52333

This command will expand the query location by 1m (1 million base pairs) in both 
sides (use --expandleft and --expandright if you only want to find overlapping 
genes to the left or right of the query location), and print out closest 
overlapping genes as well as their distance to the query.

=item * B<refGene format>

The refGene table from UCSC Genome Browser contains genome location of RefSeq 
mRNA (usually with NM_ suffix in its identifier). You can download a refGene 
template file using the UCSC Table Browser, by choosing group as "Gene and Gene 
Prediction Tracks", track as "RefSeq Genes", table as "refGene", and save the 
file as something like hg17_refgene. Additionally, you can download the "kgXref"
table and save it as something like hg17_kgxref.

The syntax for scanning refGene table is similar to that scanning the knownGene 
table. However, some refgene table (such as hg18) contains more information than 
other refgene table from earlier genome assembly (such as hg17). In that case, 
the --name2 argument can be used to replace the gene name (the "name" field in 
the refGene table) in output by the alternative gene name ("name2" field in the 
refGene table). This eliminates the need to use kgXref table and the --kgxref 
argument, but I am not sure about the reliability of the "name2" field.

To scan a query file against the refGene table:

	scan_region.pl --refgene --name2 queryfile hg18_refgene

The output is:

	CHL1    0       3:2000-349990
	NOT_FOUND       NOT_FOUND       19:32333-52333

To identify surroudning genes as well:

	scan_region.pl --refgene --name2 -expandright 1m queryfile hg18_refgene

The output is:

	CHL1    0       3:2000-349990
	OR4F17  9346    19:32333-52333

Generally speaking, refGene table contains less genes than the knownGene table, 
but the genes in refGene is better annotated and usually have gene symbol 
associated with them.

=item * B<phastCons format>

Despite its name, do not confuse it with the "most conserved" track of the 
genome browser. The phastcons file merely contains conservation scores for each 
genome position that can be aligned to other species (in comparison, the MCE 
file contains the conservation scores for the top 5% most conserved genomic 
regions in the genome).

Since the file contains a score for each genomic position, the file size is 
extremely large (typically >10GB), and the scanning takes a lot of time (up to a 
whole day). In most circumstances, I recommended using the MCE file for sequence 
conservation analysis, since the phastcons score is less meaningful than the MCE 
score, which use a hidden Markov model to identify and summarize a genomic 
region that is conserved during evolution.

=back

=cut                                       