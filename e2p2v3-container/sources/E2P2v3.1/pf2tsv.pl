#!/usr/bin/env perl

use strict;

my @a_in = `cat $ARGV[0]`;
chomp @a_in;

my %h_out = ();
my $id;
my $ec;
foreach my $line (@a_in)
{
	if ($line =~ /^ID\s+(\S+)/)
	{
		$id = $1;
		$h_out{$id} = [] unless defined $h_out{$id};
	}
	if ($line =~ /^EC\s+(\S+)/)
	{
		$ec = $1;
		push (@{$h_out{$id}}, $ec);
	}
}

foreach my $id (sort keys %h_out)
{
	print "$id\t" . join(';', @{$h_out{$id}}) . "\n" if (scalar  @{$h_out{$id}} > 0);
}

