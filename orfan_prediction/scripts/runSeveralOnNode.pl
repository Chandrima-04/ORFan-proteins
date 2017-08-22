`module add easy`;
`module add heimdal`;

$command = $ARGV[0];

@commands = split( "&", $command );

foreach $aCommand (@commands) {

    chomp($aCommand);

    if ( !( $aCommand eq "" ) ) {
        system( $aCommand. "&" );
    }
}
