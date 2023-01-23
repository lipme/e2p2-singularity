#!/usr/local/bin/perl -w 


# Chuan Wang
# 20150617
# Plant Metabolic Network
# Department of Plant Biology
# Carnegie Institution for Science
# Stanford, CA 94305

# this script converts EC numbers to offical MetaCyc reaction ids in the .pf file

use FindBin;

my $usage = <<USAGE;

  Usage: $0 <pf-file>
  
USAGE

# Check arguments
if (@ARGV<1) {
    die $usage;
}


# read metacyc general mapping file
my %map;
open(GMAP, "<$FindBin::Bin/data/metacyc-RXN-EC.mapping") or die $!;

while (<GMAP>) {
    next if (/^ID|^#/);
    chomp;
    my @r = split /\t/;
    $r[2] = '-' if !$r[2];
    $r[2] =~ s/ -- .*$//;
    $map{$r[2]}{$r[0]} ++;
}

close GMAP;


# read remove list
my %remove;
open(NO, "<$FindBin::Bin/data/to-remove-non-small-molecule-metabolism.txt") or die $!;

while (<NO>) {
    next if (/^ID|^#/);
    chomp;
    my @r = split /\t/;
    $remove{$r[0]} ++;
}

close NO;


# read sub-rxn file
my %subrxn;
open(R, "<$FindBin::Bin/data/metacyc-sub-reactions") or die $!;

while (<R>) {
    next if (/^ID|^#/);
    chomp;
    my @r = split /\t/;
    $subrxn{$r[0]}{$r[1]} ++;
}

close R;


# read metacyc superseded file
my %ecsuperseded;
open(S, "<$FindBin::Bin/data/EC-superseded") or die $!;

while (<S>) {
    next if (/^ID|^#/);
    chomp;
    my @r = split /\t/;
    map {s/^EC-//} @r;
    $ecsuperseded{$r[2]}{$r[0]} ++;
}

close S;


# read metacyc mapping file
my %ec2rid;
open(MAP, "<$FindBin::Bin/data/metacyc-RXN-official-EC.mapping") or die $!;

while (<MAP>) {
    next if (/^ID|^#/);
    chomp;
    my @r = split /\t/;
    my @s = split /\|/, $r[1];
    $ec2rid{$r[0]}{$_} ++ foreach @s;
}

close MAP;

# read pf file
my %rxn;
my $id = '';
my $protein = '';

open(my $IN, "<$ARGV[0]") or die $!;
my $outfile = $ARGV[0];
$outfile =~ s/\.pf$//;
open(my $out, '>', "$outfile.orxn.pf") or die $1;

while (<$IN>) {
    chomp;
    if (/^ID\t(.*)$/) {
        $id = $1;
        $protein = $_. "\n";
    } elsif (/^METACYC\t(.*)$/) {
        $protein .= $_. "\n" unless ($remove{$1} || $rxn{$id}{$1});
        $rxn{$id}{$1} ++;
    } elsif (/^EC\t(.*)$/) {
        my $ec = $1;
        next if $remove{$ec};
        if ($ecsuperseded{$ec}) {
            foreach my $c (keys %{$ecsuperseded{$ec}}) {
                next if $remove{$c};
                my $ecprinted = 0;
                foreach my $rx (keys %{$ec2rid{$c}}) {
                    next if $remove{$rx};
                    $protein .= "METACYC\t$rx\n" unless $rxn{$id}{$rx};
                    $rxn{$id}{$rx} ++;
                    $ecprinted ++;
                }
                if (!$ecprinted && scalar(keys %{$map{$c}}) > 0) {
                    foreach my $rx (keys %{$map{$c}}) {
                        next if $remove{$rx};
                        $protein .= "METACYC\t$rx\n#unofficial\n" unless $rxn{$id}{$rx};
                        $rxn{$id}{$rx} ++;
                    }
                }
            }
        } else {
            my $ecprinted = 0;
            foreach my $rx (keys %{$ec2rid{$ec}}) {
                next if $remove{$rx};
                $protein .= "METACYC\t$rx\n" unless $rxn{$id}{$rx};
                $rxn{$id}{$rx} ++;
                $ecprinted ++;
            }
            if (!$ecprinted && scalar(keys %{$map{$ec}}) > 0) {
                foreach my $rx (keys %{$map{$ec}}) {
                    next if $remove{$rx};
                    $protein .= "METACYC\t$rx\n#unofficial\n" unless $rxn{$id}{$rx};
                    $rxn{$id}{$rx} ++;
                }
            }
        }
    } elsif (/^\/\//) {
        $protein .= $_ . "\n";
        print $out $protein if (keys %{$rxn{$id}} > 0);
        $protein = '';
    } else {
        $protein .= $_ . "\n";
    }
}

close $IN;
close $out;
