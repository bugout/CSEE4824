#!/usr/bin/perl

use warnings;
use strict;

my @lines = <>;
my @out = ();
my $op;
my $num_args = $#ARGV + 1;
my $infile = $num_args ? $ARGV[0] : 'all.mat';
my $outfile = 'in.mat';

for (@lines) {
    push @out, $_;

    if (/# op: (\w+)/) {
	$op = $1;
    }

    if (/# --+/ and @out) {
	# print to file
	open FH, ">$outfile" or die "can't open '$outfile': $!";
	foreach (@out) {
	    print FH $_;
	}
	close FH or die "cannot close '$outfile': $!";
	@out = ();

	# invoke octave
	system("octave -q $op.m");
	if ($? == -1) {
	    print STDERR "failed to execute: $!\n";
	    exit 1;
	}
	elsif ($? & 127) {
	    printf STDERR "child died with signal %d, %s coredump\n",
	    ($? & 127),  ($? & 128) ? 'with' : 'without';
	    exit 1;
	}
	else {
	    my $errval = $? >> 8;
	    exit 1 if $errval;
	}
    }
}
