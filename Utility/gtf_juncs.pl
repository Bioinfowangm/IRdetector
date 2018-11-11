#!/usr/bin/perl 
#===========================================================================
#
#         FILE: gtf_juncs.pl
#        USAGE: perl gtf_juncs.pl
#
#       AUTHOR: Wang Meng, mengwang55@gmail.com
#      VERSION: 1.0
#      CREATED: 03/25/2018 04:37:19 PM
#===========================================================================

use strict;
use warnings;

my ( %strand, %Exon, %chr );
while (<>) {
    chomp;
    my @row = split;
    next unless $row[2] eq 'exon';
    my ($transcript_id) = /transcript_id "(.*?)";/;
    push @{ $Exon{$transcript_id} }, @row[ 3, 4 ];
    $strand{$transcript_id} = $row[6];
    $chr{$transcript_id}    = $row[0];
}

for my $k ( keys %Exon ) {
    my @coordinates = sort { $a <=> $b } @{ $Exon{$k} };
    for my $i ( 1 .. ( @coordinates / 2 - 1 ) ) {
        print join( "\t",
            $chr{$k},
            $coordinates[ $i * 2 - 1 ],
            $coordinates[ $i * 2 ]-1,
            $strand{$k} ),
          "\n";
    }
}
