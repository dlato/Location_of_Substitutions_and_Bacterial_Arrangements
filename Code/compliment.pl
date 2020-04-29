#!/usr/bin/perl

while(@ARGV) {
    $arg = shift @ARGV;
    if($arg eq '-h') { 
            print << "help_text";
This program reads fasta formated sequence files from either a
file or from stdin and writes to stdout the equivalent compliment
formated sequence file.  It can do multiple records per file.
    Syntax: compliment file
help_text
            exit(1); next; }
#
#   print STDOUT "The arg is $arg\n";
#   if($arg eq "") { open(FILE, STDIN) || die "Error: Can't open stdin\n"; }
#   else { open(FILE, $arg) || die "Error: Can't open file $arg\n"; }
   open(FILE, $arg) || die "Error: Can't open file $arg\n"; 
   $i=-1;
   while (<FILE>) {
     if(/^>/) { 
	 $i++; 
	 chop $_; 
	 $title[$i]=$_; 
	 $title[$i]=~ s/^> //; 
	 $title[$i]=~ s/^>//; 
#  Don't care about the length; let PHYLIP take care of itself
#         $len=length($title[$i]);
#         if($len <10) { for($j=0; $j < 9-$len; $j++) { $title[$i] .= " ";}}
#         else { $title[$i] = substr($title[$i],0,10);}
     }
     elsif ($_ =~ /\s*(\d+)\s+(\d+)/){ $line = $_;}
     else { $seq[$i] .= $_; # chop $seq[$i]; 
           }
   }
   close(FILE);
   for($j=0; $j<=$i; $j++) {
       $seq[$j] = reverse $seq[$j];
       $seq[$j] =~ s/^\s*//;
       $seq[$j] =~ s/a/%/ig;
       $seq[$j] =~ s/t/A/ig;
       $seq[$j] =~ s/%/T/ig;
       $seq[$j] =~ s/c/%/ig;
       $seq[$j] =~ s/g/C/ig;
       $seq[$j] =~ s/%/G/ig;
       print STDOUT "> $title[$j]\n";
       print STDOUT "$seq[$j]\n";
   }
}
