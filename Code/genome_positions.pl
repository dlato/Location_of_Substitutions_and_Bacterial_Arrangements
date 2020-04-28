#!/usr/bin/perl
use strict;
use warnings;
#########################################################
#########################################################
#this script will go through the mafft file to get the positions in
#the alignment and then store them.
#it will then go through the treefile and will make sure that the
#positions are printed out in the same order as the tree file so that
#the call ancestor script will work properly
#it will ultimatly print out commands for each of the positions in the
#alignment so it can go through call ancestor and reconstruct the
#ancestral positions for each of the positions in the alignment.
#
#to run the script:
#genome_positions.pl Block17.mafft treefile > out
#########################################################
#########################################################


#initalize variables
my $block_name = "";
my @seqName = ();
my @aln = ();
my $max = 0;
my $number = 0;
my @aln_chars = ();
my @branch = ();
my @from1 = ();
my @to1 = ();
my @site = ();
my @from2 = ();
my @to2 = ();
my @prob = ();
my @start_genome_pos = ();
my @reversed = ();
my $aln_len = 0;

#MAFFT FILE INFO
#grabbing the alignment from the mafft file as well as the position in
#the genome information.
open(FILE, '<:encoding(UTF-8)', $ARGV[0]) || die "can not open file: $!";
if($ARGV[0]=~/(\S+)\.\S+/) { $block_name = $1; } else { $block_name = $ARGV[0]; } 
my $i=-1;
while(<FILE>) {
    if(/^>/) {
	$aln[++$i]="";
	$seqName[$i] = $_;
    } else { chomp; $aln[$i] .=  $_; }
}
close(FILE);
$number = $i+1; 
for($max=0, $i=0; $i<$number; $i++) { 
    if(length($aln[$i]) > $max) { $max = length($aln[$i]); } 
}
for($i=0; $i<$number; $i++) { 
    while(length($aln[$i]) < $max) { $aln[$i] .= "-"; } 
}

#spliting up the alignment by character
for($i=0; $i<$number; $i++) { 
    my @temp = split //, $aln[$i];
    $aln_len = scalar(@temp);
    push (@aln_chars, \@temp); #pushes stuff into an array of array
}

#obtaining the start positions for each taxa and keeping track of if
#the seq is reversed
@start_genome_pos = @seqName;
my @aln_taxa_order = @start_genome_pos;
for (my $j=0; $j<scalar(@start_genome_pos); $j++) {
    if($start_genome_pos[$j] =~ /Reversed/) {
	$reversed[$j] = 1;
	$start_genome_pos[$j] =~ m/-/;
	$start_genome_pos[$j] =~ s/$`//;
	$start_genome_pos[$j] =~ s/-//;
	$start_genome_pos[$j] =~ m/:/;
	$start_genome_pos[$j] =~ s/$'//;
	$start_genome_pos[$j] =~ s/ //;
	$start_genome_pos[$j] =~ s/://;
	$start_genome_pos[$j] =~ s/\n//;
	$start_genome_pos[$j] = $start_genome_pos[$j] +1;
    } else {
	$reversed[$j] = 0;
	$start_genome_pos[$j] =~ m/-/;
	$start_genome_pos[$j] =~ s/$'//;
	$start_genome_pos[$j] =~ s/-//;
	$start_genome_pos[$j] =~ m/:/;
	$start_genome_pos[$j] =~ s/$`//;
	$start_genome_pos[$j] =~ s/://;
	$start_genome_pos[$j] =~ m/ /;
	$start_genome_pos[$j] =~ s/$`//;
	$start_genome_pos[$j] =~ s/ //;
	$start_genome_pos[$j] =~ s/\n//;
	$start_genome_pos[$j] = $start_genome_pos[$j] -1;
    }#else
}#for

#getting the order of the taxa in alignment
for (my $j=0; $j<scalar(@start_genome_pos); $j++) {
    if($reversed[$j] == 1) {
	$aln_taxa_order[$j] =~ m/:/; 
	$aln_taxa_order[$j] =~ s/$'//; 
	$aln_taxa_order[$j] =~ s/>//; 
	$aln_taxa_order[$j] =~ s/ //; 
	$aln_taxa_order[$j] =~ s/://; 
	$aln_taxa_order[$j] =~ s/\n//; 
    } else {
	$aln_taxa_order[$j] =~ m/:/; 
	$aln_taxa_order[$j] =~ s/$'//; 
	$aln_taxa_order[$j] =~ s/>//; 
	$aln_taxa_order[$j] =~ s/\n//; 
	$aln_taxa_order[$j] =~ s/://; 
    }#else
}#for


#read in treefile for the correct order of the taxa
my $treefile = "";
my @treefile_order = ();
open(FILE, '<:encoding(UTF-8)', $ARGV[1]) || die "can not open file: $!";
while(<FILE>) {
    $treefile = $_;
}#while
close(FILE);
@treefile_order = split /:/, $treefile;
for (my $a =0; $a < scalar(@treefile_order); $a++) {
    if ($treefile_order[$a] !~ /NC/) {
	$treefile_order[$a] = "";
    } else {  
	$treefile_order[$a] =~ s/^[^NC]*(?=NC)//;#removes everything before the NC
    }#else
}#for


#below is the part that prints out the actual positions and the call
#ancestor command
my @started = 0 x @aln_taxa_order;
#printing the starting positions for each site in the alignment
for (my $k = 0; $k < $aln_len; $k++) {
    print "/home/dlato/CODE/call_ancestor.pl ";
    for (my $b = 0; $b < scalar(@treefile_order); $b++) {
	if ($treefile_order[$b] ne "") {
	    for (my $c = 0; $c < scalar(@aln_taxa_order); $c++) {
		if ($treefile_order[$b] eq $aln_taxa_order[$c]) {
		    if ($aln_chars[$c][$k] ne "-") {#not a gap
			$started[$c]++;
			#altering the starting positions when seq
			#starts with a gap which affects the ending
			#positions (off by one)
			if ($k == 0) {
			    if ($reversed[$c] == 0) {
				$start_genome_pos[$c]++;
			    } else {
				$start_genome_pos[$c]--;
			    }#else
			}#if
			if ($started[$c] == 1) {
			    print "$start_genome_pos[$c] ";
			    goto OUT;
			}#if
			if ($reversed[$c] == 0 ) {
			    $start_genome_pos[$c]++;
			    print "$start_genome_pos[$c] ";
			} else {
			    $start_genome_pos[$c]--;
			    print "$start_genome_pos[$c] ";
			}#else
		    } else {#gap
			#altering the starting positions when seq
			#starts with a gap which affects the ending
			#positions (off by one)
			if ($k == 0) {
			    if ($reversed[$c] == 0) {
				$start_genome_pos[$c]++;
			    } else {
				$start_genome_pos[$c]--;
			    }#else
			}#if
			print "$start_genome_pos[$c] ";
		    }#else
		    OUT:
		} else {
		    #do nothing because the treefile order and aln
		    #order are different.
		}#else
	    }#for
	}#if
    }#for
    print ">> $block_name" . "_ancestor_out\n";
#    print ">> $block_name" . "_ancestor_out_" . count "\n";
}#for
