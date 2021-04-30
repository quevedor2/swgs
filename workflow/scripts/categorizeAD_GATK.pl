#!/usr/bin/perl
use List::Util qw/sum/;

my ($filename, $ref_col, $alt_col) = @ARGV;
open(FH, "<", $filename) or die $!;

while(<FH>){
  ## Isolate the samples Ref and Alt allelic depth
  my @spl=split(/\t/, $_);
  my $ref = $spl[$ref_col];
  my $alt = $spl[$alt_col];

  ## Categorize the depth per SNP
  my $sums = $ref + $alt;
  if($sums == 0){
    print "0\n"; ## NO COVERAGE
  } else {
    my $frac = $ref / $sums;
    if ($frac == 1) {
      print "1\n"; ## REF.HOMOZYGOUS
    } elsif ($frac == 0) {
      print "3\n"; ## ALT.HOMOZYGOUS
    } else {
      print "2\n"; ## HETEROZYGOUS
    }
  }
}
close(FH)
