#!/usr/bin/perl

use strict;

#use Statistics::Descriptive;

my %BED         = ();
my %seq2feat    = ();
my %idx2feature = ();
my %feature2idx = ();
my %feature2BED = ();



my $bed_file = $ARGV[0]; ## EDeNseq index bed file
my $res_file = $ARGV[1]; ## EDeNseq results file



my $line = "";

## read in required info from result file header, stop at end of header
open( RES, "zcat $res_file |" );
while ( $line = <RES> ) {
  chomp $line;
  if ( $line =~ /^##/ ) {
    next 
  } elsif ( $line =~ /^#HIST_IDX/ ) {
    my @tmp = split( " ", $line );
    $idx2feature{ $tmp[1] } = $tmp[3];
    $feature2idx{ $tmp[3] } = $tmp[1];
    #print "RES HIST IDX\t", $tmp[1], " -> ", $tmp[3], "\n";
  } elsif ( $line !~ /^#/ ) {
    last;
  }
}


## read in used index bed file
open( BED, $bed_file );
while ( my $line = <BED> ) {
  chomp $line;
  my @tmp = split( '\t', $line );
  $BED{ $tmp[0] } = \@tmp;
  $feature2BED{ $tmp[3] } = \@tmp;
  #print $tmp[0], " SEQ ", $tmp[3], "\n";
}
close(BED);

my %counts;
my $num = 0;

## continue with resulkts file
do {
  chomp $line;
  my @tmp = split("\t",$line);
  my @id = split("-",$tmp[0]);
  my @max = split(",",$tmp[9]);
  my @all = split(",",$tmp[7]);
  
  $num++;
  if ($num%1000000 == 0){
    print $num." analyzed seqs\n";
  }

  my $targetFeature = "";
  if (!exists $BED{$id[0]}){
    #print "cannot eval ".$id[0]." - no BED entry\n";
    #$targetFeature = -1;
    #next;
  } else {
    $targetFeature = $BED{$id[0]}->[3];
  }
  my $targetIdx = $feature2idx{$targetFeature};
  
  if (!exists $counts{$targetFeature}){
    #$counts{$targetFeature} = {};
  }
  
  $counts{$targetFeature}->{"SEQS"}++;
  
  if (@tmp <=7){
    #print "Not CLASSIFIED","\n";
    $counts{$targetFeature}->{"NOCLASS"}++;
    #next;
  } elsif (@max == 1 && $targetFeature == $idx2feature{$max[0]}) {
   #print "Unique  MATCH max $max[0]","\n";
   $counts{$targetFeature}->{"UNIQ"}++;
   $counts{$targetFeature}->{"HIT"}++;
  } elsif ((@max > 1) && (($tmp[9] =~ /\,$targetIdx\,/) || ($tmp[9] =~ /^$targetIdx\,/)) ) {
    #print "NON Unique match #max=".@max." ".join(":",@max)."\n";
    $counts{$targetFeature}->{"NUNIQ"}++;
    $counts{$targetFeature}->{"HIT"}++;
  } elsif ((@all >= 1) && (($tmp[7] =~ /\,$targetIdx\,/) || ($tmp[7] =~ /^$targetIdx\,/)) ) {
    #print "NON MAX match\n";
    $counts{$targetFeature}->{"NOMAX"}++;
  } else {
    #print "No MATCH #max=".@max." ".$tmp[6],"\n";
    $counts{$targetFeature}->{"NOMATCH"}++;
  }
  
} while ( $line = <RES> );

close(RES);

#print "keys:".join(":",keys %counts)."\n";
my %sum_avg;
foreach my $feature (sort {$feature2BED{$a}->[7] cmp $feature2BED{$b}->[7]} keys %counts){
  print "RESULTS\t",$feature,"\t",sprintf("%25.25s",$feature2BED{$feature}->[7]),"\t","SEQS","\t",$counts{$feature}->{"SEQS"},"\t";
  foreach my $key (("UNIQ","NUNIQ","HIT","NOMAX","NOMATCH","NOCLASS")){
    print $key."\t".$counts{$feature}->{$key},"\t".sprintf("%.3f",$counts{$feature}->{$key}/$counts{$feature}->{"SEQS"})."\t";  
   $sum_avg{$key."_SUM"} += $counts{$feature}->{$key};
   $sum_avg{$key."_AVG"} += $counts{$feature}->{$key}; 
  }
  print "\n";
  $num -= $counts{$feature}->{"SEQS"};
  $sum_avg{"SEQS_SUM"} += $counts{$feature}->{"SEQS"};
}

print "\nRESULTS_ALL\t",sprintf("%25.25s","ALL"),"\t","SEQS","\t",$sum_avg{"SEQS_SUM"},"\t";
foreach my $feature ("UNIQ","NUNIQ","HIT","NOMAX","NOMATCH","NOCLASS"){
  print $feature."\t".$sum_avg{$feature."_SUM"},"\t".sprintf("%.3f",$sum_avg{$feature."_AVG"}/($sum_avg{"SEQS_SUM"}))."\t";
}


print "\n\nNot evaluated (unknown) seqs: ",$num,"\n\n";