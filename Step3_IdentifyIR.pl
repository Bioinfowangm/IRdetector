#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::Handle;
use Bio::SeqIO;
use List::Util qw/sum/;

# Parameter section
my $man   = 0;
my $help  = 0;
my $readl = 75;
GetOptions(
    'ReadLength=i' => \$readl,
    'help|?'       => \$help,
    man            => \$man
) or pod2usage(2);
pod2usage(0) if $help;
pod2usage( -exitval => 0, -verbose => 2 ) if $man;

# Creat folder Step3
if ( -e "Step3" ) {
    warn "File/Folder Step3 exist...\n";
}
else {
    mkdir "Step3";
}

# Check the coverage of every selected introns across samples
for my $bed (<./Step1/*Junc.bed>) {
    my ($SampName) = $bed =~ m#Step1/(.*).Junc.bed#;
    &CheckIntronCoverage($SampName);
}

# Filter introns: retain those with 100% coverage and at least mapped 10 reads
my ( %all, %retained );
open I_g, "<./Step2/potential_introns.bed";
while (<I_g>) {
    chomp;
    $all{$_} = 1;
}
for my $file (<./Step3/*candidateIR>) {
    open I_f, $file;
    while (<I_f>) {
        chomp;
        my @row = split;
        next unless $all{ join( "\t", @row[ 0 .. 3 ] ) };
        $retained{ join( "\t", @row[ 0 .. 3 ] ) } = 1
          if $row[4] >= 1 && $row[5] >= 10;
    }
}

# For each Candidate IR, calculate the PIR
for my $bed (<./Step1/*Junc.bed>) {
    my ($SampName) = $bed =~ m#Step1/(.*).Junc.bed#;
    &ObtainIR($SampName);
}

sub ObtainIR {
    my $samp = shift;
    my ( %intron, %junction, %events );
    open I_f, "./Step3/$samp.candidateIR";
    open O_f, ">./Step3/$samp.IR";
    while (<I_f>) {
        chomp;
        my @row = split;
        next unless $retained{ join( "\t", @row[ 0 .. 3 ] ) };
        $intron{ join( "\t", @row[ 0 .. 3 ] ) } = $row[5];
    }
    open I_f2, "./Step1/${samp}.Junc.bed";
    while (<I_f2>) {
        chomp;
        my @row = split;
        $junction{ join( "\t", @row[ 0, 1, 2 ] ) } = $row[3];
    }
    open I_events, "<./Step2/potential_introns_intersect.bed";
    while (<I_events>) {
        chomp;
        my @row = split;
        next unless $retained{ join( "\t", @row[ 0 .. 3 ] ) };
        $events{ join( "\t", @row[ 0 .. 3 ] ) }{ join( "\t", @row[ 4 .. 7 ] ) }
          = 1;
    }

    for my $j ( keys %events ) {
        my @pos         = split /\s+/, $j;
        my $length      = $pos[2] - $pos[1];
        my @tj          = keys %{ $events{$j} };
        my $intron_read = $intron{$j} ? $intron{$j} : 0;
        my $junction_read;
        map {
            my @row = split /\s+/, $_;
            $junction_read +=
                $junction{ join( "\t", @row[ 0, 1, 2 ] ) }
              ? $junction{ join( "\t", @row[ 0, 1, 2 ] ) }
              : 0;
        } @tj;
        my $intron_read_normalized = $intron_read / ( $length + $readl - 11 );
        my $junction_read_normalized = $junction_read / ( $readl - 11 );
        next unless ( $intron_read + $junction_read ) > 0;
        my $pir = $intron_read_normalized /
          ( $intron_read_normalized + $junction_read_normalized );
        print O_f join( "\t",
            $j, $intron_read, $junction_read, $intron_read_normalized,
            $junction_read_normalized, $pir ),
          "\n";
    }
}

sub CheckIntronCoverage {
    my $samp = shift;
    open I_F1,
"intersectBed -a ./Step1/${samp}.bed12 -b ./Step2/potential_introns.bed -wa -wb|"
      || die "Error input for CheckIntronCoverage subroutine";
    open O_F1, ">./Step3/${samp}.candidateIR";
    my ( %coverage, %readcount );
    while (<I_F1>) {
        chomp;
        my @row = split;

        #    $row[5] =~ tr/-+/+-/ if $row[3] =~ /\/1/;
        #    next unless $row[5] eq $row[15];
        next if $row[9] > 1;
        $row[3] =~ s/\/[12]//;
        my $read         = $row[3];
        my @blocks       = split /,/, $row[10];
        my @block_starts = split /,/, $row[11];

        my $start = $row[1] + 1;
        my $end   = $row[2];
        next
          if ( $row[2] < ( $row[13] + 6 ) || $row[14] < ( $row[1] + 6 ) )
          ;    # read should overlap at least 6 bps of intron region

        map {
            if ( $_ >= ( $row[13] + 1 ) && $_ <= $row[14] ) {
                $coverage{ join( "\t", @row[ 12 .. 15 ] ) }{$_}++;
                $readcount{ join( "\t", @row[ 12 .. 15 ] ) }{$read}++;
            }
        } $start .. $end;
    }

    for my $k ( keys %coverage ) {
        my @row = split /\s/, $k;
        my $c = ( keys %{ $coverage{$k} } ) / ( $row[2] - $row[1] );
        my $r = keys %{ $readcount{$k} };
        print O_F1 join( "\t", $k, $c, $r ), "\n";
    }
    O_F1->autoflush(1);
}

__END__
=head1 NAME

Step3_IdentifyIR.pl

=head1 SYNOPSIS

perl Step3_IdentifyIR.pl [options]

 Options:
   -ReadLength      the length of reads, required
   -help            brief help message
   -man             full documentation


=head1 OPTIONS

=over 8

=item B<-ReadLength>

Required parameter. The length of single-end read, or either end of paired-end reads (default is 75)

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This is step 3>, which identifies candidate IRs (100% coverage and at least 10 reads in at least 1 sample), and calculate the PIR value of these IRs

=cut
