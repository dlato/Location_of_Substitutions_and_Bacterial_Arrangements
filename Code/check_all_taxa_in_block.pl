#!/usr/bin/perl
########################################################################################################################
# xmfa_example.pl
# Last modified:	2014-06-03
# this will go through the xmfa file and print out in a new file
# "*_all_taxa_in_block" all the blocks that have all the taxa in it.
# this taxa number can be changed below from 5 to whatever need to
# input the xmfa file as an argument
########################################################################################################################
use strict;
use warnings;
use lib qw(/home/wilson/PerlMod);
use Bio::AlignIO;
my $file =$ARGV[0]; #"AK83chrom1_BL225C_chrom.xmfa";
my $tax_num =$ARGV[1]; #"AK83chrom1_BL225C_chrom.xmfa";
#my @taxa= ();#sets an empty array (which is really a column vector not
my %taxa = (); #here we use a hash because we want the taxa to be associated with their genomes. if we used an array it would just put them in an order which is not what we want, we need the sequence stored for each block to be specific to the genome
#a row vector
my @seq_lens= ();
my $in = Bio::AlignIO->new(-file => $file, -format => 'xmfa');


$file=~ s/\.xmfa//;
my $name= $file . "_all_taxa_in_block";

open (FILE, ">$name") || die "cant open file";
for ( my $aln_num = 1; my $aln = $in->next_aln; $aln_num += 1 ){ #bioperl can recognize each alignment block from the xmfa file and it gets each one and that is what the -> next_aln does
    my @seqs = $aln->each_seq();
#below if just says that we need to have more than 5  sequences in our
#alignment (because it is pointless to get the sequences that did not
#match with anything (they will be useless in determining a tree). we
#want all the seqs present in each block
	for ( my $seq_num = 0; $seq_num < scalar(@seqs); $seq_num += 1 ){
	}
    if ($tax_num == scalar(@seqs)){
	print FILE "mv Block$aln_num" . "_* GoodBlocks/\n";
    } # if
} # for

close (FILE);


