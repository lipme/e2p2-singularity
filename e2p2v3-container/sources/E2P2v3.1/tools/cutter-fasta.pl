#!/usr/bin/perl

# Chuan Wang
# 20150617
# Plant Metabolic Network
# Department of Plant Biology
# Carnegie Institution for Science
# Stanford, CA 94305

# this script cuts FASTA files into specified number of smaller files

use strict;
use warnings;

my $nn=0;
my $i=0;
my ($file,$num)=@ARGV;
my $flength=-s "$file";
$flength=0 if (!$flength);
my $leng=$flength/$num;

open IN,"$file";
my $outfile="${file}_seg_"."$nn";
open OUT,">$outfile";
select OUT;
while (<IN>) {

	$i+=length;
	if ($i > $leng && /^>/) {
		close OUT;
		$nn++;
		my $outfile = sprintf "${file}_seg_%02d", $nn;
		open OUT,">$outfile";
		select OUT;
		$i=0;
	}
	print;
}
close OUT;
select STDOUT;
