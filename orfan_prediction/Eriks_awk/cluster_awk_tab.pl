#!/usr/bin/perl
eval 'exec /usr/bin/perl -S $0 ${1+"$@"}'
    if $running_under_some_shell;
			# this emulates #! processing on NIH machines.
			# (remove #! line above if indigestible)

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches

open(NAMES_FH, $ARGV[0]) || die "Cannot open file $ARGV[0].";

# Make single-linkage clusters from a list of pairwise links (matches)
# Script was originally made for blast output (see Gerald/THE_NEW_VH_UNIVERSE)

# Requires file (first cmd-line argument) with all sequence names in current dir.
# Most recent version adapted to work with output from selfcomp.crossmatch

$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

# seq[seqname] = cluster number of 'seqname'

while (($_ = &Getline2('NAMES_FH'),$getline_ok)) {
    $seq{$Fld1} = ++$n;# pointer seq -> cluster
    $names{$Fld1} = $Fld1;# all seqs

    ;
}

while (<>) {
    next if /^#/;
    ($Fld1,$Fld2,$Fld3,$Fld4,$Fld5) = split(' ', $_, -1);

    $s1 = $Fld1;
    $s2 = $Fld2;

    # Join clusters of s1 and s2 by setting all s2 clusters to the s1 cluster

    $clus1 = $seq{$s1};
    $clus2 = $seq{$s2};

    foreach $i (keys %names) {
	if ($seq{$i} eq $clus2) {	#???
	    $seq{$i} = $clus1;
	}
    }
}

# Count cluster sizes
$clus = 0;
for ($j = 1; $j <= $n; $j++) {
    foreach $i (keys %names) {
	if ($seq{$i} == $j) {	#???
	    $clussize{$j}++;
	}
    }
}

$clus = 0;
for ($j = 1; $j <= $n; $j++) {
    $first = 1;
    foreach $i (keys %names) {
	if ($seq{$i} == $j && $clussize{$j} > 1) {	#???
	    if ($first) {
		print 'Cluster ' . ++$clus . ':';
		$first = 0;
	    }
	    print $i;
	}
    }
    if (!$first) {
	print '';
    }
}

sub Getline2 {
    ($fh) = @_;
    if ($getline_ok = (($_ = <$fh>) ne '')) {
	($Fld1,$Fld2,$Fld3,$Fld4,$Fld5) = split(' ', $_, -1);
    }
    $_;
}
