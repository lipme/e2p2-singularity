#!/usr/bin/perl -w

# Chuan Wang @ Rhee-lab, Department of Plant Biology, Carnegie Institution for Science
# Date: June, 2015

# Description: Generating jobscripts to run E2P2 on clusters.
# Input: fasta file folder, script folder, a sample of job script
# Output: Files in the script folder

my $usage = <<USAGE;

  Usage: $0 <fasta-folder> <script-folder> <sample-job-script> <full-path-runE2P2.py>

USAGE

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;
use File::Copy;

$|=1;   # Output prior
my $path = dirname(File::Spec->rel2abs(__FILE__));  # Path of this script
my $time = sprintf "%d/%02d/%02d-%02d:%02d:%02d", map { $$_[5]+1900, $$_[4]+1, $$_[3], $$_[2], $$_[1], $$_[0] } [localtime];
my $cmdline = "$0 ".join(' ',@ARGV);

print STDERR "\nStarts at [$time]\nCmd: $cmdline\n\n";

if (@ARGV<2) {
    die $usage;
}

# read status file
my ($fafolder, $jobfolder) = map {s!/$!!; $_} @ARGV;
`mkdir -p $jobfolder`;

foreach my $fa (`ls $fafolder/*`) {
    chomp $fa;
    next unless -f $fa;
    my $jobname = $fa;
    $jobname =~ s/.*\///;
    my $cmd = "$ARGV[3] -i $fa -o $fa.e2p2v3";
    open(SCRPT, ">$jobfolder/e2p2v3-$jobname") or die $!;
    if ($ARGV[2]) {
        open(SP, "<$ARGV[2]") or die $!;
        while (<SP>) {
            s/\$JOBNAME/blast-$jobname/;
            s/\$CMD/$cmd/;
            print SCRPT;
        }
        close SP;
    } else {
        print SCRPT "#!/bin/bash\n$cmd\n";
    }
    close SCRPT;
}
