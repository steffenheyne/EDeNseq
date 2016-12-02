#!/usr/bin/perl

use strict;

#use Statistics::Descriptive;

my %BED         = ();
my %seq2feat    = ();
my %idx2feature = ();
my %feature2idx = ();
my %feature2BED = ();

my $bed_file = $ARGV[0];    ## EDeNseq index bed file
my $res_file = $ARGV[1];    ## EDeNseq results file

my $line = "";

## read in required info from result file header, stop at end of header
open( RES, "zcat $res_file |" );
while ( $line = <RES> ) {
	chomp $line;
	if ( $line =~ /^##/ ) {
		next;
	} elsif ( $line =~ /^#HIST_IDX/ ) {
		my @tmp = split( " ", $line );
		$idx2feature{ $tmp[1] } = $tmp[3];
		$feature2idx{ $tmp[3] } = $tmp[1];

		#print "RES HIST IDX\t", $tmp[1], " -> ", $tmp[3], "\n";
	} elsif ( $line =~ /#SEQ/ ) {
		last;
	}
}

## read in used index bed file
open( BED, $bed_file );
while ( my $line = <BED> ) {
	chomp $line;
	my @tmp = split( '\t', $line );
	$BED{ $tmp[0] }         = \@tmp;
	$feature2BED{ $tmp[3] } = \@tmp;

	#print $tmp[0], " SEQ ", $BED{$tmp[0]}->[3], "\n";
}
close(BED);

print "keys in bed", scalar keys %feature2BED, "\n";

my %counts;
my $num = 0;

my $numNotEval = 0;

## continue with results file
{
	while ( $line = <RES> ) {
		chomp $line;
		my @tmp = split( "\t", $line );
		my @id  = split( "-",  $tmp[0] );
		if ($tmp[0] =~ /reference_seq/){
		  @id = join("_",(split( "_",  $tmp[0] ))[0..2]);
		  #print "test: ".$id[0]."\n";
		}
		my @max = split( ",",  $tmp[9] );
		my @all = split( ",",  $tmp[7] );

		if ( $num % 1000000 == 0 ) {
			print $num. " analyzed seqs\n";
		}

		my $targetFeature = "";
		if ( exists $BED{ $id[0] } ) {
			$targetFeature = $BED{ $id[0] }->[3];
		} else {

			#print "cannot eval " . $id[0] . " - no BED entry\n";
			#$targetFeature = -1;
			$numNotEval++;
		}

		if ( $targetFeature != "" ) {

			my $targetIdx = $feature2idx{$targetFeature};

      if (!(exists $counts{$targetFeature})){
        $counts{$targetFeature}->{"SEQS"}  = "-";
        $counts{$targetFeature}->{"NOCLASS"} = "-";
        $counts{$targetFeature}->{"UNIQ"} = "-";
		    $counts{$targetFeature}->{"HIT"} = "-";
		    $counts{$targetFeature}->{"NUNIQ"} = "-";
        $counts{$targetFeature}->{"NOMAX"} = "-";
        $counts{$targetFeature}->{"NOMATCH"} = "-";
      }

			$num++;
			$counts{$targetFeature}->{"SEQS"}++;

			if ( @tmp <= 7 ) {

				#print "Not CLASSIFIED","\n";
				$counts{$targetFeature}->{"NOCLASS"}++;

				# print $id[0]."-".$id[1]."\tNOCLASS\n";
				#next;
			} elsif ( @max == 1 && $targetFeature == $idx2feature{ $max[0] } ) {

				#print "Unique  MATCH max $max[0]","\n";
				$counts{$targetFeature}->{"UNIQ"}++;
				$counts{$targetFeature}->{"HIT"}++;
			} 	elsif (	( @max > 1 ) && ( ( $tmp[9] =~ /\,$targetIdx\,/ ) || ( $tmp[9] =~ /^$targetIdx\,/ ) ) ) {

				#print "NON Unique match #max=".@max." ".join(":",@max)."\n";
				$counts{$targetFeature}->{"NUNIQ"}++;
				$counts{$targetFeature}->{"HIT"}++;
			}
			elsif (	( @all >= 1 ) && (    ( $tmp[7] =~ /\,$targetIdx\,/ ) || ( $tmp[7] =~ /^$targetIdx\,/ ) ) )
			{

				#print "NON MAX match\n";
				$counts{$targetFeature}->{"NOMAX"}++;

				#  print $id[0]."-".$id[1]."\tNOMAX\n";
			} else {

				#print "No MATCH #max=".@max." ".$tmp[6],"\n";
				$counts{$targetFeature}->{"NOMATCH"}++;

				#   print $id[0]."-".$id[1]."\tNOMATCH\n";
			}

			if ( $num % 7000000 == 0 ) {
				#	 last;
			}
		}
	}
}
close(RES);

#print "keys:".join("\n",keys %counts)."\n";
print $num."\n";
my %sum_avg;
foreach my $feature (sort { $feature2BED{$a}->[7] cmp $feature2BED{$b}->[7] } keys %counts )
{
	print "RESULTS\t", $feature, "\t",
	  sprintf( "%25.25s", $feature2BED{$feature}->[7] ), "\t", "SEQS", "\t",
	  $counts{$feature}->{"SEQS"}, "\t";
	foreach my $key ( ( "UNIQ", "NUNIQ", "HIT", "NOMAX", "NOMATCH", "NOCLASS" ) ) 
	{
		print $key. "\t" . $counts{$feature}->{$key}, "\t" . sprintf( "%.3f",abs($counts{$feature}->{$key} / $counts{$feature}->{"SEQS"} )) . "\t";
		$sum_avg{ $key . "_SUM" } += $counts{$feature}->{$key};
		$sum_avg{ $key . "_AVG" } += $counts{$feature}->{$key};
	}
	print "\n";
	$num -= $counts{$feature}->{"SEQS"};
	$sum_avg{"SEQS_SUM"} += $counts{$feature}->{"SEQS"};
}

print "\nRESULTS_ALL\tALL\t", sprintf( "%20.20s", "ALL" ), "\t", "SEQS", "\t", $sum_avg{"SEQS_SUM"}, "\t";
foreach my $feature ( "UNIQ", "NUNIQ", "HIT", "NOMAX", "NOMATCH", "NOCLASS" ) {
	print $feature. "\t" . $sum_avg{ $feature . "_SUM" }, "\t" . sprintf( "%.3f",$sum_avg{ $feature . "_AVG" } / ( $sum_avg{"SEQS_SUM"} ) ) ."\t";
}

print "\n\n                   seqs left: ", $num,        "\n";
print "Not evaluated (unknown) seqs: ",     $numNotEval, "\n\n";
