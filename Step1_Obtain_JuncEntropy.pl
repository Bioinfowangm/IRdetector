#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Parameter section
my $man = 0;
my $help = 0;
my $bam = "";
my $mapq = 50;
my $sampn = "";
my $readl = 75;
GetOptions(
    'BAM=s' => \$bam,
    'MapQual=i' => \$mapq,
    'ReadLength=i' => \$readl,
    'SampName=s' => \$sampn,
    'help|?' => \$help, 
    man => \$man) 
    or pod2usage(2);
pod2usage(0) if $help;
pod2usage(-exitval => 0, -verbose => 2) if $man;


# Create folder Step1
if(-e "Step1"){
    warn "File/Folder Step1 exist...\n";
}
else{
    mkdir "Step1";
}
system("samtools view -q $mapq -b $bam |bamToBed -bed12 -i -  > Step1/${sampn}.bed12");
system("perl ./Utility/check_JuncEntropy.pl Step1 $sampn $readl");



__END__
=head1 NAME

Step1_Obtain_JuncEntropy.pl

=head1 SYNOPSIS

perl Step1_Obtain_JuncEntropy.pl [options]

 Options:
   -BAM             input BAM file, required
   -SampName        name of the sample, required
   -ReadLength      the length of reads, required
   -MapQual         the least mapping quality specified [50]
   -help            brief help message
   -man             full documentation


=head1 OPTIONS

=over 8

=item B<-BAM>

Required parameter. The input BAM file to be processed

=item B<-SampName>

Required parameter. The name of the sample. This will be 
the prefix of output files

=item B<-ReadLength>

Required parameter. The length of single-end read, or either end of paired-end reads

=item B<-MapQual>

Reads with mapping quality >= this value are retained. 
Default is 50 for uniquely mapped reads from TopHat

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This is step 1>, which identifies the Junctions of a bam file, and prepares information to calculate entropy in the next step.
If there are multiple samples, they should be run one-by-one, or parallel

=cut
