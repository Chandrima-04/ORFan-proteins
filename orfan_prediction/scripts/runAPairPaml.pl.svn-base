use strict;
use warnings;
use File::Basename;

my $rid = shift @ARGV;

my $dir = "/scratch/davepaml$rid";
my $outdirpath =
  '/afs/pdc.kth.se/home/d/dmessina/dave3/virus/20100301/pairs_codeml';
my $bindir = '/afs/pdc.kth.se/home/d/dmessina/bin/paml44_clust/bin';

mkdir($dir) or die "couldn't mkdir $dir: $!\n";

foreach my $file (@ARGV) {

    my ( $base, $path, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
    my $seqfile  = $file;
    my $treefile = $path . '/' . $base . '.tre';

    my $seqonclust  = $dir . '/' . $base . $suffix;
    my $treeonclust = $dir . '/' . $base . '.tre';

    # extract cluster id from path to create outdir for that cluster
    my $outdir;
    if ( $path =~ /(cluster\d+)/ ) {
        $outdir = $outdirpath . '/' . $1;
        if ( !-e $outdir ) {
            mkdir($outdir) or die "couldn't make $outdir: $!\n";
        }
    }
    else { $outdir = $outdirpath; }

    system("cp $seqfile $dir")  && die "couldn't cp $seqfile $dir: $!\n";
    system("cp $treefile $dir") && die "couldn't cp $treefile $dir: $!\n";

    # first run of codeml
    open( FH, ">$dir/codeml.ctl" );
    print FH "seqfile = $seqonclust\n",
      "outfile = $dir/$base.out\n",
      "treefile = $treeonclust\n",
      "runmode = 0\n",
      "verbose = 0\n",
      "seqtype = 1\n",
      "CodonFreq = 2\n",
      "model = 0\n",
      "clock = 1\n",
      "cleandata = 1\n",
      "NSsites = 0 3\n";
    close(FH);
    system("cd $dir; $bindir/codeml");

    # copy files back to my disk
    system("cp -r $dir $outdir")
      && die "couldn't cp $dir to $outdir: $!\n";
}

system("rm -rf $dir") && die "couldn't rm $dir: $!\n";
