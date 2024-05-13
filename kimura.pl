##-------------------------------------------------------------------------##
sub calcKimuraDivergence {
  my $this            = shift;
  my %nameValueParams = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $divCpGMod = $nameValueParams{'divCpGMod'};

  my $DEBUG = 0;

  my $transitions            = 0;
  my $transversions          = 0;
  my $CpGSites               = 0;
  my $wellCharacterizedBases = 0;
  my $pSBase                 = "";
  my $prevTrans              = 0;
  my @sBases                 = split //, $this->{'subjSeq'};
  my @qBases                 = split //, $this->{'querySeq'};
  foreach my $i ( 0 .. ( length( $this->{'subjSeq'} ) - 1 ) ) {
    next if ( $sBases[ $i ] eq "-" );

    $wellCharacterizedBases++
        if ( $wellCharacterizedBases{ $qBases[ $i ] . $sBases[ $i ] } );

    if ( $divCpGMod && $pSBase eq "C" && $sBases[ $i ] eq "G" ) {

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
        # CpG sites contains 2 transitions ( treat as 1 trans )
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
      }
      elsif ( $mt == 2 ) {
        $transversions++;
      }
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

  return ( $kimura, $transitions, $transversions, $wellCharacterizedBases,
           $CpGSites );
}
