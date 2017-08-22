#!/usr/bin/perl
# $Id: compare_rnacode.pl 3106 2011-09-05 13:45:40Z dmessina $
# created by Dave Messina on 2010-06-14

use Modern::Perl;
use Getopt::Long;
use File::Spec;
use File::Basename;

my ( $dir1, $dir2 );
GetOptions(
    'dir1:s' => \$dir1,
    'dir2:s' => \$dir2,
);

my $usage = "
compare_rnacode.pl --dir1 <dir> --dir2 <dir>

where --dir1 and --dir2 are the two dirs of RNAcode output to compare

";
die $usage unless (defined $dir1 && defined $dir2);

# global counters
my $total   = 0;
my $longer  = 0;
my $shorter = 0;
my $betterp = 0;
my $worsep  = 0;

# open first dir and make list of files
my ( @files1, @files2 );
my $dirlookup = {
    'dir1' => $dir1,
    'dir2' => $dir2,
};

my $filelookup = {
    'dir1' => \@files1,
    'dir2' => \@files2,
};

foreach my $dirname ( 'dir1', 'dir2' ) {
    my $dir      = $dirlookup->{$dirname};
    my $filelist = $filelookup->{$dirname};

    next unless ( -d $dir );    # skip any non-directories
    next if $dir =~ /^\./;      # skip . .. and all other dirs beginning .

    opendir( my $dirhandle, $dir ) or die "couldn't open $dir:$!\n";

    while ( my $entry = readdir($dirhandle) ) {
        if ( $entry =~ /\.rnacode$/ ) {
            my $path     = $dir . '/' . $entry;
            my $fullpath = File::Spec->rel2abs($path);
            push @$filelist, $fullpath;
        }
    }
    closedir($dirhandle);
}

FILE: foreach my $file (@files1) {
    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );

    ## original, more complicated matching code
    # my ( $root, $cluster, $replicate ) =
      # $base =~ /((methjann-(?:non)*coding_\d+)_(\d+))/;
    # my $oth_file = $root . '_';

    # now I'm using same filenames, so should be easy
    my $oth_file = $base;

    # find matching file in other dir
    foreach my $item (@files2) {
        if ($item =~ $oth_file) {
            if (!-e $item) {
                say STDERR "no matching file for $base";
                next FILE;
            }

            open(my $fh1, '<', $file) or die "couldn't open $file: $!\n";
            my $string1 = <$fh1>;
            chomp $string1;
            
            open(my $fh2, '<', $item) or die "couldn't open $item: $!\n";
            my $string2 = <$fh2>;
            chomp $string2;

            next FILE unless (defined $string1 && defined $string2);
            say $base;
            print $string1, "\n";
            print $string2, "\n";

            my @fields1 = split /\t/, $string1;
            my @fields2 = split /\t/, $string2;

            # length
            if ($fields2[5] != $fields1[5]) {
                my $len_diff = $fields2[5] - $fields1[5];
                my $p_diff = $fields2[10] - $fields1[10];
                $len_diff > 0 ? $longer++ : $shorter++;
                say STDERR join "\t", '  length', $base, $len_diff, $p_diff; 
            }

            # p-val
            if ($fields2[10] != $fields1[10]) {
                my $p_diff = $fields2[10] - $fields1[10];
                    # say STDERR join "\t", 'any P', $base, $fields1[10], $fields2[10], $p_diff; 
                if ($p_diff < -0.1)    {
                    $betterp++;
                    say STDERR join "\t", 'better P', $base, $fields1[10], $fields2[10], $p_diff; 
                }
                elsif ($p_diff > 0.1 ) {
                    $worsep++;
                    say STDERR join "\t", ' worse P', $base, $fields1[10], $fields2[10], $p_diff; 
                }
            }
            $total++;
        }
    }
}

say STDERR '';
say STDERR '  longer', "\t", $longer;
say STDERR ' shorter', "\t", $shorter;
say STDERR 'better P', "\t", $betterp;
say STDERR ' worse P', "\t", $worsep;
say STDERR '   total', "\t", $total;
