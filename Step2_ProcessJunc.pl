#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use IO::Handle;
use Bio::SeqIO;
use List::Util qw/sum/;

# Parameter section
my $man    = 0;
my $help   = 0;
my $refgtf = "";
my $gmfa   = "";
GetOptions(
    'Refgtf=s'   => \$refgtf,
    'GenomeFa=s' => \$gmfa,
    'help|?'     => \$help,
    man          => \$man
) or pod2usage(2);
pod2usage(0) if $help;
pod2usage( -exitval => 0, -verbose => 2 ) if $man;

# Check Paramters
die "ERROR: No Genome assembly provided\n" unless $gmfa;
my %Junc_ref;
if ($refgtf) {
    print("Reference gtf is provided!\n");

    # gtf to bed
    &gtf_juncs;

    open I_rb, "./Step2/Junctions_ref.bed";
    while (<I_rb>) {
        chomp;
        my @row = split;
        $Junc_ref{ join( "\t", @row[ 0, 1, 2 ] ) } = $row[3];
    }
}
else {
    print("No reference junctions are provided!\n");
}

# Create folder Step2
if ( -e "Step2" ) {
    warn "File/Folder Step2 exist...\n";
}
else {
    mkdir "Step2";
}

# Obtain Entropy and Strand information of detected junctions;
# Filter the observed junctions, requiring them to have Entropy >=2, and already been annotated, if the "Overlap" option is specified
my $rEntropy = &Calc_Entropy;
my %Entropy  = %$rEntropy;
my $rStrand  = &Get_Strand;
my %Strand   = %$rStrand;

open O_I_filtered, ">./Step2/Junctions_all.bed";
for my $k ( keys %Entropy ) {
    next unless $Entropy{$k} >= 2;
    next unless $Strand{$k};
#    print "$k\t$Strand{$k}\n";
    if ($refgtf) {
        print O_I_filtered join( "\t", $k, $Strand{$k} ), "\n"
          if $Junc_ref{$k};
    }
    else {
        print O_I_filtered join( "\t", $k, $Strand{$k} ), "\n";
    }
}
O_I_filtered->autoflush(1);


# Selected from the junctions for those smallest intron units that are able to interrogate
my %junction;
open I_i,
"intersectBed -a ./Step2/Junctions_all.bed -b ./Step2/Junctions_all.bed -wa -wb|";
open O_i, ">./Step2/potential_introns_intersect.bed";
while (<I_i>) {
    chomp;
    my @row = split;
    next unless $row[3] eq $row[7];
    $junction{ join( "\t", @row[ 0 .. 3 ] ) }{ join( "\t", @row[ 4 .. 7 ] ) } =
      1;
}

my (%partI,%partE);
for my $j ( keys %junction ) {
    my @a = split /\s+/, $j;
    for my $tj ( keys %{ $junction{$j} } ) {
        my @b = split /\s+/, $tj;
        $partI{$j} = 1
          if ( $b[1] > $a[1] && $b[1] < $a[2] )
          || ( $b[2] > $a[1] && $b[2] < $a[2] );
    }
}

if ($refgtf) {
    my $rpartE = &Junc_filter_exon($refgtf);
    %partE = %$rpartE;
}

Line: for my $j ( keys %junction ) {
    for my $tj ( keys %{ $junction{$j} } ) {
#        print "$j\n";
        next Line if $partI{$j} || $partE{$j};
        print O_i "$j\t$tj\n" unless $partI{$j};
    }
}

system("cut -f 1-4 ./Step2/potential_introns_intersect.bed |sort|uniq >./Step2/potential_introns.bed");

sub Junc_filter_exon {

    # filter2, excluding ones partially overlapping other exons
    # moreover, remove ones containing other shorter exons
    my $refgtf = shift;
    my (%part);
    open I_i2, "intersectBed -a ./Step2/Junctions_all.bed -b $refgtf -wa -wb|";
    while (<I_i2>) {
        chomp;
        my @row = split;
        next unless $row[6] eq "exon";
        $part{ join( "\t", @row[ 0 .. 3 ] ) } = 1
          if ( $row[7] >= $row[1] && $row[7] <= $row[2] )
          || ( $row[8] >= $row[1] && $row[8] <= $row[2] )
          || ( $row[7] >= $row[1] && $row[8] <= $row[2] );
    }
    \%part;
}

sub gtf_juncs {
    system("perl ./Utility/gtf_juncs.pl $refgtf >./Step2/Junctions_ref.bed");
}

sub Calc_Entropy {
    my ( %forEntropy, %Entropy );
    open I_fE, "cat ./Step1/*Entropy.bed|"
      || die "ERROR: No files input for entropy calculation!\n";
    while (<I_fE>) {
        chomp;
        my @row = split;
        $forEntropy{ join( "\t", @row[ 0, 1, 2 ] ) }{ $row[3] } += $row[4];
    }

    open O_e, ">./Step2/Entropy.bed";
    for my $k ( keys %forEntropy ) {
        my @values = values %{ $forEntropy{$k} };
        my $E;
        map { my $p = $_ / sum(@values); $E += -$p * log($p) / log(2) } @values;
        print O_e join( "\t", $k, $E ), "\n";
        $Entropy{$k} = $E;
    }
    O_e->autoflush(1);
    \%Entropy;
}

sub Get_Strand {
    my ( %Junc, %Junc_strand );
    open I_fJ, "cat ./Step1/*.Junc.bed|"
      || die "ERROR: No files input for Junction info!\n";
    while (<I_fJ>) {
        chomp;
        my @row = split;
        $Junc{ join( "\t", @row[ 0, 1, 2 ] ) } += $row[3];
    }
    open O_j, ">./Step2/Junctions.tmp1.bed";
    for my $k ( keys %Junc ) {
        print O_j join( "\t", $k, $Junc{$k} ), "\n";
    }
    O_j->autoflush(1);

    system(
"seqtk subseq $gmfa ./Step2/Junctions.tmp1.bed >./Step2/Junctions_seq.tmp1.fasta"
    );

    my $in = Bio::SeqIO->new(
        -file   => "./Step2/Junctions_seq.tmp1.fasta",
        -format => "fasta"
    );
    while ( my $seqobj = $in->next_seq() ) {
        my $id     = $seqobj->id;
        my $seq    = $seqobj->seq;
        my $first2 = substr( $seq, 0, 2 );
        my $last2  = substr( $seq, length($seq) - 2, 2 );
        my $four   = "$first2$last2";
        my ( $chr, $start, $end ) = $id =~ /(.*?):(\d+)-(\d+)/;
        $start = $start - 1;
        my $strand;

        if ( $four eq 'GTAG' || $four eq 'GCAG' || $four eq 'ATAC' ) {
            $strand = "+";
        }
        elsif ( $four eq 'CTAC' || $four eq 'CTGC' || $four eq 'GTAT' ) {
            $strand = '-';
        }
        next unless $strand;

        #        print join("\t",$chr,$start,$end),"\n";
        $Junc_strand{ join( "\t", $chr, $start, $end ) } = $strand;
    }
    \%Junc_strand;
}

__END__
=head1 NAME

Step2_ProcessJunc.pl

=head1 SYNOPSIS

perl Step2_ProcessJunc.pl [options]

 Options:
   -Refgtf          The gtf file of a  reference annotation, optional but highly recommended
   -GenomeFa        The fasta sequence of a reference genome, required
   -help            brief help message
   -man             full documentation


=head1 OPTIONS

=over 8

=item B<-Refgtf>

The GTF file should have all exons of each transcript included. The chromosome/contig name of the GTF file and the fasta reference should match (i.e., both with 'chr' prefix, or without)
If this GTF is not provided, <All identified junctions> with Entropy >=2 will be considered. If yes, <Identified junctions that are annotated> with Entropy >=2 will be considered.
The GTF file is also used to filter introns that partially overlap an exon if provided

=item B<-GenomeFa>

The fasta sequence of ther reference is used to determine the strand of junctions, based on the splice site types: GT-AG,GC-AG, and AT-AC;

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This is step 2>, which filters junctions by Entropy, determines strand. And importantly, define the smallest intron from a set of overlapping introns, and removes those partially overlap annotated exons

=cut
