#!/usr/bin/perl
#
# reads a bunch of characters, then uses paml to reconstruct the
# ancestral characters.  Does this by translating the characters to
# fake amino acids and then translates them back again.
# (the characters/species must be in the same order as the treefile
# and must be not more than 20 in number).
# order of the positions entered in the command needs to be the same
# order that they appear in the treefile (when reading from left to
# right)
# 
#  syntax: call_anscestors.pl  <list_of_characters>
#          (requires files: treefile)
#          (creates  files: seqfile/codeml.ctl/rst/rst1/rub/mlc/lnf/orig.treefile)
# 

$i=0;
@aa = qw(a c d e f g h i k l m n p q r s t v w y); 
@characters = @ARGV; 

#makes an AA for each position and stores the info
foreach $x (sort @characters) {
    if(! $fake{$x} ) { $fake{ $x } = $aa[$i]; $i++; }
    if($i>20) { die "A maximum of 20 different characters is permited\n"; }
}

#opens the treefile and gets the taxa names and branch info and puts
#it in a new format/tree file that it uses later on
open(FILE,"< treefile") || die "Can't open file treefile \n";
while(<FILE>) { $string .= $_; }
close(FILE);
system("cp treefile orig.treefile");
$string =~ s/:0\.00000/:0.00001/g;
$string =~ s/\)[0-9.]+/)/g;
open(FILE,"> treefile") || die "Can't open file treefile \n";
print FILE $string;
close(FILE);

# gets the taxa names from the treefile to use later
$string =~ s/\d+\s+\d+//g;
$string =~ s/:[0-9.]*//g;
$string =~ s/[,;)(]/ /g;
@taxa = split /\s+/, $string;

#creates a seqfile which  has the taxa names and the corrisponding AA
#that was assigned to each taxa/position
$i=0;
open(FILE,"> seqfile") || die "Can't open file seqfile \n";
print FILE " $#taxa  1\n";
foreach $x (@taxa) {
    print FILE "$taxa[$i+1]\n$fake{$characters[$i]}\n";
    $i++;
}
close(FILE);

#uses codeml in paml (sub routine specified below) to do the
#reconstruction
&create_codeml;
#below specifies what version of PAML to use
system("/usr/local/paml/current/bin/codeml codeml.ctl > /dev/null");
# system("/usr/local/paml/current/bin/codeml codeml.ctl ");
%inverse = reverse %fake;

#creats rst file which is just prints out the results of the output
#from PAML but not the actual reconstructions. those are just spit out
#to STDERR
open(FILE,"< rst") || die "Can't open file rst \n";
while(<FILE>) {
    if($_ =~ /TreeView/) {
	while(<FILE>) { 
	    $labelled_tree .= $_;
	    if($_ =~ /;/) { last; }
	}
#	print STDERR "For the tree\n$labelled_tree\n";
#	print "#For the tree\n#$labelled_tree\n";
    }
    if($_ =~ /List of extant and reconstructed sequences/) {
	$_ = <FILE>;
	$_ = <FILE>;
	$_ = <FILE>;
	while(<FILE>) {
	    if($_ =~ /(.*)\s+([A-Y])/) { 
#		print "$inverse{lc $2}\t"; 
		print "$inverse{lc $2}\t"; 
	    }
	    if($_ =~ /^$/) { last; }
	}
	last;
    }
}
close(FILE);
print "\n";

#subroutine that creates the codeml file that is used during the PAML
#part of this code. 
sub create_codeml {
    open(FILE,"> codeml.ctl") || die "Can't create file codeml.ctl\n";
    print FILE <<EOF;
      seqfile = seqfile   * sequence data filename
     treefile = treefile  * tree structure file name
      outfile = mlc       * main result file name

        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1  * 0: concise; 1: detailed, 2: too much
      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 2  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*        ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 0
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
  fix_blength = 2 * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
EOF
close(FILE);
}
