##-------------------------------------------------------------------------##

#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my %mutType = (
                "CT" => 1,
                "TC" => 1,
                "AG" => 1,
                "GA" => 1,
                "GT" => 2,
                "TG" => 2,
                "GC" => 2,
                "CG" => 2,
                "CA" => 2,
                "AC" => 2,
                "AT" => 2,
                "TA" => 2,
                "ct" => 1,
                "tc" => 1,
                "ag" => 1,
                "ga" => 1,
                "gt" => 2,
                "tg" => 2,
                "gc" => 2,
                "cg" => 2,
                "ca" => 2,
                "ac" => 2,
                "at" => 2,
                "ta" => 2
                
);

#
# Well Characterized Alignment Pairings
#   A well characterized pairing is one between a fixed
#   consensus base ( not IUB codes ) and a fixed query
#   base ( not IUB codes ).
#
my %wellCharacterizedBases = (
                               'CG' => 1,
                               'CA' => 1,
                               'CC' => 1,
                               'CT' => 1,
                               'TA' => 1,
                               'TC' => 1,
                               'TG' => 1,
                               'TT' => 1,
                               'AA' => 1,
                               'AC' => 1,
                               'AG' => 1,
                               'AT' => 1,
                               'GA' => 1,
                               'GC' => 1,
                               'GG' => 1,
                               'GT' => 1,
                               'cg' => 1,
                               'ca' => 1,
                               'cc' => 1,
                               'ct' => 1,
                               'ta' => 1,
                               'tc' => 1,
                               'tg' => 1,
                               'tt' => 1,
                               'aa' => 1,
                               'ac' => 1,
                               'ag' => 1,
                               'at' => 1,
                               'ga' => 1,
                               'gc' => 1,
                               'gg' => 1,
                               'gt' => 1,
);


# Subroutine to calculate Kimura Divergence
sub calcKimuraDivergence {
    my ($subjSeq, $querySeq) = @_;

    my $divCpGMod = 1; # You may want to specify this value or modify it accordingly

    my $DEBUG = 0;

    my $deletions              = 0;
    my $insertions             = 0;
    my $transitions            = 0;
    my $transversions          = 0;
    my $CpGSites               = 0;
    my $wellCharacterizedBases = 0;
    my $pSBase                 = "";
    my $prevTrans              = 0;
    my @sBases                 = split //, $subjSeq;
    my @qBases                 = split //, $querySeq;

  foreach my $i ( 0 .. ( length( $subjSeq ) - 1 ) ) {

    if ( $sBases[ $i ] eq "-" ) {
        $insertions++;
        next;
    }
    if ( $qBases[ $i ] eq "-" ) {
        $deletions++;
        #next;

    }

    $wellCharacterizedBases++
        if ( $wellCharacterizedBases{ $qBases[ $i ] . $sBases[ $i ] } );

    if ( $divCpGMod && $pSBase eq "c" && $sBases[ $i ] eq "g" ) {

      # CpG
      $CpGSites++;
      my $mt = $mutType{ $qBases[ $i ] . $sBases[ $i ] } || 0;
      if ( $mt == 1 ) {
        $prevTrans++;
      }
      elsif ( $mt == 2 ) {
        $transversions++;
      }
      if ( $prevTrans == 2 ) {
        # CpG sites contains 2 transitions ( treat as 1 trans )$subjSeq
        $prevTrans = 1;
      }
      elsif ( $prevTrans == 1 ) {

        # CpG sites contains 1 transition ( treat as 1/10 trans )
        $prevTrans = 1 / 10;
      }
    }
    else {
      $transitions += $prevTrans;
      $prevTrans = 0;

      # Normal
      my $mt = $mutType{ $qBases[ $i ] . $sBases[ $i ] } || 0;
      if ( $mt == 1 ) {

        # Delay recording transition for CpG accounting
        $prevTrans = 1;
        #print "i\n"
      }
      elsif ( $mt == 2 ) {
        $transversions++;
        #print "t\n"
      }
      #print "\n"
    }
    $pSBase = $sBases[ $i ];
  }
  $transitions += $prevTrans;

  my $kimura = 100.00;
  if ( $wellCharacterizedBases >= 1 ) {
    my $p          = $transitions / $wellCharacterizedBases;
    my $q          = $transversions / $wellCharacterizedBases;
    my $logOperand = ( ( 1 - ( 2 * $p ) - $q ) * ( 1 - ( 2 * $q ) )**0.5 );
    if ( $logOperand > 0 ) {
      $kimura = ( abs( ( -0.5 * log( $logOperand ) ) ) * 100 );
    }
  }

  return ( $kimura, $transitions, $transversions, $wellCharacterizedBases, $CpGSites , $insertions, $deletions, length($subjSeq));
}

# Get a list of all files in the current directory


my @nameinfo = @ARGV;
my $naem = $nameinfo[0];
my $class = $nameinfo[1];
my $outname = "summary.$naem.tsv";

open(my $outfile, '>', $outname) or die "Could not open file '$outname' $!";    
print $outfile "sequence\torder\tk80adjust\ttrans\ttransv\tdefBases\tnumCpGs\tinsert\tdelet\talignLength\tfam\tidentifier(s)\n";

my @files = glob("*.aln");
if (@files == 0) {
   exit
} 
my @nameparts = split /\./, $files[0];
my $name = $nameparts[0];
# Iterate through each file
foreach my $file (@files) {
    # Open the input file
    open(my $fh, '<', $file) or die "Could not open file '$file' $!";

    my ($concencussequence, $querysequence) = ("", "");
    my $current_sequence = 1;
    my $identifiers = '';
    # Read the input file
    while (my $line = <$fh>) {
        chomp $line;
        if ($line =~ /^>/) {
            # If the line starts with ">", it's a header line
            # If it's not the first header line, switch to the second sequence
            if ($line =~ /\)\s*(.*)/) {
                $identifiers = $1; # Capture the identifiers
            }
	    if ($current_sequence == 1 && $concencussequence) {
                $current_sequence = 2;
                next;
            }
        } else {
            # Concatenate sequence lines for the current sequence
            if ($current_sequence == 1) {
                $concencussequence .= $line;
            } elsif ($current_sequence == 2) {
                $querysequence .= $line;
            }
        }
    }

    my ( $div, $transi, $transv, $wellCharBases, $numCpGs , $insertions, $deletions, $alignmentlength) = calcKimuraDivergence($concencussequence, $querysequence);
    
    print $outfile "$file\t$class\t$div\t$transi\t$transv\t$wellCharBases\t$numCpGs\t$insertions\t$deletions\t$alignmentlength\t$naem\t$identifiers\n";
    close($fh);

}
close($outfile);
