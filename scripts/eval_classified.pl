#!/usr/bin/perl

use strict;

#use Statistics::Descriptive;

my %BED         = ();
my %seq2feat    = ();
my %idx2feature = ();
my %feature2idx = ();
my %feature2BED = ();



my $bed_file = $ARGV[0];
my $res_file = $ARGV[1];


open( RES, "zcat $res_file |" );

my $line = "";
while ( $line = <RES> ) {
  chomp $line;
  if ( $line =~ /^#HIST_IDX/ ) {
    my @tmp = split( " ", $line );
    $idx2feature{ $tmp[1] } = $tmp[3];
    $feature2idx{ $tmp[3] } = $tmp[1];
    #print "RES HIST IDX\t", $tmp[1], " -> ", $tmp[3], "\n";
  } else {
    last;
  }
}

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

while ( $line = <RES> ) {
  chomp $line;
  my @tmp = split("\t",$line);
  my @id = split("-",$tmp[0]);
  my @max = split(",",$tmp[6]);
  my @all = split(",",$tmp[4]);
  
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
  my $targetIdx = $feature2idx{$BED{$id[0]}->[3]};
  
  if (!exists $counts{$targetFeature}){
    #$counts{$targetFeature} = {};
  }
  
  $counts{$targetFeature}->{"SEQS"}++;
  
  if (@tmp <=4){
    #print "Not CLASSIFIED","\n";
    $counts{$targetFeature}->{"NOCLASS"}++;
    next;
  } elsif (@max == 1 && $targetFeature == $idx2feature{$max[0]}) {
   #print "Unique  MATCH max $max[0]","\n";
   $counts{$targetFeature}->{"UNIQ"}++;
   $counts{$targetFeature}->{"HIT"}++;
  } elsif ((@max > 1) && (($tmp[6] =~ /\,$targetIdx\,/) || ($tmp[6] =~ /^$targetIdx\,/)) ) {
    #print "NON Unique match #max=".@max." ".join(":",@max)."\n";
    $counts{$targetFeature}->{"NUNIQ"}++;
    $counts{$targetFeature}->{"HIT"}++;
  } elsif ((@all >= 1) && (($tmp[4] =~ /\,$targetIdx\,/) || ($tmp[4] =~ /^$targetIdx\,/)) ) {
    #print "NON MAX match\n";
    $counts{$targetFeature}->{"NOMAX"}++;
  } else {
    #print "No MATCH #max=".@max." ".$tmp[6],"\n";
    $counts{$targetFeature}->{"NOMATCH"}++;
  }
  
}

close(RES);

#print "keys:".join(":",keys %counts)."\n";

foreach my $feature (sort {$feature2BED{$a}->[7] cmp $feature2BED{$b}->[7]} keys %counts){
  print "RESULTS\t",$feature,"\t",sprintf("%25.25s",$feature2BED{$feature}->[7]),"\t","SEQS","\t",$counts{$feature}->{"SEQS"},"\t";
  foreach my $key (("UNIQ","NUNIQ","HIT","NOMAX","NOMATCH","NOCLASS")){
    print $key."\t".$counts{$feature}->{$key},"\t".sprintf("%.3f",$counts{$feature}->{$key}/$counts{$feature}->{"SEQS"})."\t";  
  }
  print "\n";
  $num -= $counts{$feature}->{"SEQS"}; 
}

print "\nNum Seqs left ",$num,"\n";