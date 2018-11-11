#!/usr/bin/perl 
#===========================================================================
#
#         FILE: check_JuncEntropy.pl
#        USAGE: perl check_JuncEntropy.pl
#
#       AUTHOR: Wang Meng, mengwang55@gmail.com
#      VERSION: 1.0
#      CREATED: 03/19/2018 11:07:51 AM
#===========================================================================

use strict;
use warnings;
use List::Util qw/sum/;

my ($tmp,$sampn,$readl) = @ARGV;
my $inBed      = "$tmp/$sampn.bed12";
my $outJunc    = "$tmp/$sampn.Junc.bed";
my $outEntropy = "$tmp/$sampn.Entropy.bed";

my ( %Junc, %Entropy );
open IN, $inBed || die "No input Bed12 file provided";
while (<IN>) {
    chomp;
    my @row = split;
    next unless $row[9] > 1;
    my @blocks       = split /,/, $row[10];
    my @blocks_start = split /,/, $row[11];
    next if sum(@blocks) != $readl;
    for my $i ( 1 .. ( $row[9] - 1 ) ) {
        my $start = $row[1] + $blocks[ $i - 1 ] + $blocks_start[ $i - 1 ];
        my $end   = $row[1] + $blocks_start[$i];
        next if $blocks[ $i - 1 ] < 6 || $blocks[$i] < 6;
        $row[3] =~ s/\/[12]//g;
        my $read = $row[3];
        $Junc{ join( "\t", $row[0], $start, $end ) }{$read}++;
        $Entropy{ join( "\t", $row[0], $start, $end ) }{ $row[10] }++;
    }
}

open OUT_J, ">$outJunc"    || die "No output Junction file provided";
open OUT_E, ">$outEntropy" || die "No output Entropy file provided";

# Junctions
for my $k ( keys %Junc ) {
    my @row = split /\t/, $k;
    my $count = keys %{ $Junc{$k} };
    print OUT_J join( "\t", @row[ 0, 1, 2 ], $count ), "\n";
}

#entropy
for my $k ( keys %Entropy ) {
    map { print OUT_E join( "\t", $k, $_, $Entropy{$k}{$_} ), "\n" }
      keys %{ $Entropy{$k} };
}
