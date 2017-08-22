#!/usr/bin/perl
# $Id: make_slaves.pl 2968 2011-05-09 11:59:38Z dmessina $

use strict;
use warnings;

die if @ARGV != 2;
my ($start, $end) = @ARGV;

my $retval;
my $last_job = 0;
for (my $i = $start; $i <= $end; $i += 50) {
    my $first = $i;
    my $last  = ($i + 49 >= $end) ? $end : $i+49;

    print $first, "\t", $last, "\n";

    my $command = 'for f in {'. $first . '..' . $last . '} ;' .
	' do echo $f ; esubmit -J ' . $first . '_' . $last . ' -v -n 1 -t 2000 ';

    unless ( $i == $start ) {
	$command = $command . ' -F ' . $last_job . ' ';
    }

    $command .= ' ~/dave4/virus/20110303/blast_slave.sh $f a06c11n07.pdc.kth.se ; done';
#    $command .= " ~/dave6_nobak/virus/test.sh ; done";


    # run it
    print STDERR $command, "\n";
    $retval = `$command`;
    print STDERR $retval, "\n";
    if ($retval =~ /jid=(\d+)/g) {
	my @jobs = $1;
	$last_job = $jobs[-1];
    }
    defined ($last_job) or die "no last job\n";
}

# my $c = 'esubmit -v -n 1 -t 2000 echo ' . $start;
# print STDERR $c, "\n";
# my $retval = `$c`;
# my $last_job;
# if ($retval =~ /jid=(\d+)/) {
#    $last_job = $1;
# }
# defined ($last_job) or die "no last job\n";

# my $newc = 'esubmit -v -n 1 -t 2000 -F ' . $last_job . " echo \"last job was $last_job\" ";
# print STDERR $newc, "\n";
# my $newretval = `$newc`;
# print STDERR "DONE\n";
