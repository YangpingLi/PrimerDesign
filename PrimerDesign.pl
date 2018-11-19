print "
					Gene-specific primers design
					
Gene-specific primers should be designed following the five simple guidelines:
a. The 3-terminal base should be A or T.
b. The five 3-terminal bases should include no more than two G or C.
c. The four 3-terminal bases should not find a perfect match within the primer being designed and within primers that are going to be used in PCR with this primer (adapterspecific primers and another gene-specific primer to amplify positive control, if applicable).
d. The length of the primer should be at least 20 bases.
e. The annealing temperature of the primer should be equal or higher than 60oC; calculated by the formula: 4x(G+C) + 2x(A+T) + 3
\n";
###########################################################################
sub parse_command_line {
    if(!@ARGV){usage();}
    else{
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-r$/)  { $genome   = shift @ARGV; }
	elsif ($_ =~ /^-s$/)  { $gff  = shift @ARGV; }
	elsif ($_ =~ /^-c$/)  { $cds  = shift @ARGV; }
	elsif ($_ =~ /^-g$/)  { $gene_list{shift @ARGV}++; }
        elsif ($_ =~ /^-sca$/){ $sca  = shift @ARGV; }
        elsif ($_ =~ /^-gl$/) { $gene_list  = shift @ARGV; }
        elsif ($_ =~ /^-b$/)  { $blast  = shift @ARGV; }
        elsif ($_ =~ /^-gs$/) { $geneseq  = shift @ARGV; }
        elsif ($_ =~ /^-plmin$/) { $pl_min  = shift @ARGV; }
        elsif ($_ =~ /^-plmax$/) { $pl_max  = shift @ARGV; }
	else {
	    usage();
	}
    }
    }
}

sub usage {
    print STDERR <<EOQ; 
USAGE:
    perl PrimerDesign.pl [-h]   
    r  : reference genome, in fasta format
    c  : reference cds sequences, in fasta format
    s  : genome stucture file, in gff3 format
    g  : gene name to perform Primer design
    sca: scaffold ID, generate a fasta file of this scaffold
    gl : gene list for primer design
    b  : if run blast for the primer designed [0]
    gs : gene sequence, in fasta format
    plmin: minmum length of PCR product  [400]
    plmax: maximum length of PCR product [800]
    h  : display the help information
	
EOQ
exit(0);
}

sub primer_design (@_){
   my $T;
   my $candi_p;
   my $n_G;
   my $n_C;
   my $seq=$_[0];
   my $sign=$_[1];
   my $length=length($seq);
   my @primer=();
   for(my $i=0;$i<$length;$i++){
      $T=0;
      for(my $l=$min_Plen;$T<$Tmax;$l++){
         $candi_p=substr($seq,$i,$l);
         if($i+$l>$length){last;}
         if($sign eq antisense){$candi_p=revcomp($candi_p);}
         $T=tmvalue($candi_p);
         if($T<$Tmin||$T>$Tmax){;}
         else{
            my $A=substr($candi_p,$l-1,);
            my $AAAAA=substr($candi_p,$l-5);
            my $AAAA=substr($candi_p,$l-4);
            if($A eq C || $A eq G){;}
            elsif(GC_content($candi_p)>0.6){;}
            elsif(GC_content($candi_p)<0.4){;}
            else{
               $n_G=($AAAAA=~s/G/G/g);
               $n_C=($AAAAA=~s/C/C/g);
               if($n_G+$n_C>2){;}
               else{
                  my $revcompAAAA=revcomp($AAAA);
                  if($candi_p=~/$revcompAAAA/){;}
                  elsif($candi_p=~/AAA/||$candi_p=~/TTT/||$candi_p=~/CCC/||$candi_p=~/GGG/){;}
                  else{
                     push (@primer,$candi_p);
                     $p_start{$candi_p}=$i+1;
                     $p_end{$candi_p}=$i+$l;
                     $T{$candi_p}=$T;
                     $GC{$candi_p}=GC_content($candi_p);
                  }
               }
            }
         }
      }
   }
   return @primer;
}


sub tmvalue(@_){
   my $seq=$_[0];
   my $G=($seq=~s/G/G/g);
   my $C=($seq=~s/C/C/g);
   my $A=($seq=~s/A/A/g);
   my $T=($seq=~s/T/T/g);
   4*($G+$C)+2*($A+$T)+3;
}

sub GC_content(@_){
   my $seq=$_[0];
   my $G=($seq=~s/G/G/g);
   my $C=($seq=~s/C/C/g);
   my $A=($seq=~s/A/A/g);
   my $T=($seq=~s/T/T/g);
   ($G+$C)/($G+$C+$A+$T);
}


sub revcomp(@_){
   my $dna = shift;
   my $revcomp = reverse($dna);
   $revcomp =~ tr/ACGTacgt/TGCAtgca/;
   return $revcomp;
}

sub cross_dimer(@_){
  my $s=$_[0];
  my $a=$_[1];
  $rev_sAAAA=revcomp(substr($s,(length($s)-3)));
  $rev_aAAAA=revcomp(substr($a,(length($a)-3)));
  (($s=~/$rev_aAAAA/)||($a=~/$rev_sAAAA/))
}

$min_Plen=20;
$Tmin=60;
$Tmax=75;
$blast=0;
$pl_min=100;
$pl_max=200;
parse_command_line();
########################################################## getseq #######
if($sca){
   print "Calling for $sca"."'s sequence...\n\n";
   open GE,"<$genome";
   open OUT,">$sca".".fasta";
   my $jud=0;
   while(<GE>){
      chomp;
      if(/>/){
         s/>//;
         if($_ eq $sca){
             print OUT ">$_\n";
             $jud=1;
         }
         else{$jud=0;}
      }
      elsif($jud){print OUT "$_\n";}
      else{;}
   }
   $seq{$name}=$seq;
   close GE;
   close OUT;
   print "Got $sca"."'s sequence in $sca.fasta\n\n";
   exit(0);
}
else{;}
############################################ Loading genome data #########
if($genome){
print "Genome file is $genome\n";
print "Loading genome data...\n";

my $seq="";
my $name="";
%seq=();

open GE,"<$genome";
while(<GE>){
   chomp;
   if(/>/){
      $seq{$name}=$seq;
	  $seq="";
      s/>//;
	  @_=split(/\s+/);
	  $name=$_[0];
   }
   else{$seq="$seq"."$_";}
}
$seq{$name}=$seq;
close GE;
print "Genome data was loaded!\n\n";
}
else{;}
######################################## Loading genome structrue data ####
if($gff){
open GFF,"<$gff";
print "Structure file is $gff\n";
print "Loading genome structrue data...\n";
while(<GFF>){
   chomp;
   if(/^#/){print "$_\n";}
   else{
      @_=split(/\t/);
      @info=split(/\;|\=/,$_[-1]);
      ${$info[2]}{$info[1]}=$info[-1];
      $genome{$info[1]}=$_[0];
      $start{$info[1]}=$_[3];
      $end{$info[1]}=$_[4];
      $dir{$info[1]}=$_[6];
      push(@{${$_[2]}{$info[-1]}},$info[1]);
      if(/mRNA/){$gene_s{$info[1]}=$_[3];$gene_e{$info[1]}=$_[4];}
      elsif(/CDS/){
         push(@{$cds_s{$info[3]}},$_[3]-1);
         push(@{$cds_e{$info[3]}},$_[4]-1);
         push(@{$cds_l{$info[3]}},abs($_[4]-$_[3])+1);
      }
   }
}
close GFF;
print "Genome structure data was loaded!\n\n";
}
else{;}
############################################## get mRNA sequence ###########
##  Need to be modified 


#################################################### Load cds data #########
if($cds){
   print "cds file is $cds\n";
   print "Loading cds data...\n";
   open CDS,"<$cds";
   while(<CDS>){
      chomp;
      if(/>/){
         $seq{$name}=$seq;
	     $seq="";
         s/>//;
	     $name=$_;
      }
      else{$seq="$seq"."$_";}
   }
   $seq{$name}=$seq;
   close CDS;
   print "CDS data was loaded!\n\n";
}
else{;}

##################################################### get gene list ########
if($gene_list){
   print "Get gene list for primer design...\n";
   open FIG,"<$gene_list";
   while(<FIG>){
      chomp;
      if($_ eq ""){;}
      else{$gene_list{$_}++;}
   }
   close FIG;
   print "Got gene list!\n\n"
}
else{;}
############################################################################
if($geneseq){
   print "gene sequence file is $geneseq\n";
   open FIG,"<$geneseq";
   while(<FIG>){
      s/\r//g;
      chomp;
      if(/>/){
         $_=~s/>//;
         $gene=$_;
         $seq{$gene}="";
         $gene_list{$gene}++;
      }
      else{$seq{$gene}="$seq{$gene}"."$_";}
   }
   system("mkdir PrimerDesign");
   foreach $gene(keys %gene_list){
      open OUT,">./PrimerDesign/$gene"."_PrimerPairs";
      open OUT1,">./PrimerDesign/$gene"."_PrimerPairs.fasta";
      print "Primer design for gene $gene...\n\n";
      print OUT "The CDS sequence of $gene is :\n$seq{$gene}\n\n";
      @sense=primer_design($seq{$gene},sense);
      @anti=primer_design($seq{$gene},antisense);
      print OUT "PrimerPairs\tsense\tp_start\tp_end\tTm_value\tGC_content\tantisense\tp_start\tp_end\tTm_value\tGC_content\tProduct_length\n";
      $cnt=1;
      foreach $s (@sense){
         foreach $a (@anti){
            $deta_T=abs($T{$s}-$T{$a});
            $p_l=$p_end{$a}-$p_start{$s};
            if($p_l>$pl_max){;}
            elsif($p_l<$pl_min){;}
            elsif(abs($T{$s}-$T{$a})>2){;}
            elsif(cross_dimer($s,$a)){;}
            else{
               printf OUT "PrimerPairs-$cnt\t$s\t$p_start{$s}\t$p_end{$s}\t$T{$s}\t%.2f\t$a\t$p_start{$a}\t$p_end{$a}\t$T{$a}\t%.2f\t$p_l\n",$GC{$s},$GC{$a};
               print OUT1">PrimerPairs-$cnt-sense-$p_start{$s}-$p_end{$s}\n$s\n>PrimerPairs-$cnt-antisense-$p_start{$a}-$p_end{$a}\n$a\n";
               $cnt++;
            }
         }
      }
      close OUT;
      close OUT1;
   }
}
else{
   ##################################################### Primer Design ########
   system("mkdir PrimerDesign");
   foreach $gene(keys %gene_list){
      $num_cds=@{$cds_e{$gene}};
      open OUT,">./PrimerDesign/$gene"."_PrimerPairs";
      open OUT1,">./PrimerDesign/$gene"."_PrimerPairs.fasta";

      print "Primer design for gene $gene...\n\n";
      print OUT "The CDS sequence of $gene is :\n$seq{$gene}\n\n";
      print OUT "The CDS sequence contains $num_cds exon:\n";

      for(my $i=0;$i<@{$cds_l{$gene}};$i++){
         print OUT "@{$cds_s{$gene}}[$i]\t@{$cds_e{$gene}}[$i]\t@{$cds_l{$gene}}[$i]\t$dir{$gene}\n";
      }
      print OUT "\n";
      
      if($dir{$gene} eq "-"){
         @{$cds_l{$gene}}=reverse(@{$cds_l{$gene}});
         @{$cds_s{$gene}}=reverse(@{$cds_s{$gene}});
         @{$cds_e{$gene}}=reverse(@{$cds_e{$gene}});
      }
      
      my $node=0;
      foreach (@{$cds_l{$gene}}){$node=$node+$_;push(@{$cds_n{$gene}},$node);}
      print OUT "The node site is @{$cds_n{$gene}}\n";
      
      my $cds_seq="";
      for(my $i=0;$i<@{$cds_s{$gene}};$i++){
         my $cds=substr($seq{$genome{$gene}},@{$cds_s{$gene}}[$i],@{$cds_l{$gene}}[$i]);
         $cds_seq="$cds_seq"."$cds";
      }
      if($dir{$gene} eq "-"){$cds_seq=revcomp($cds_seq);}
      print OUT "The cds sequence that I made myself is:\n>$gene\n$cds_seq\n";
      if($cds_seq eq $seq{$gene}){print "The sequence is OK!\n\n";}
      else{print "The sequence of CDS cannot match!\n\n";}

      @sense=primer_design($seq{$gene},sense);
      @anti=primer_design($seq{$gene},antisense);

      foreach my $s(@sense){$flag{$s}="N";}
      foreach my $a(@anti){$flag{$a}="N";}
      print OUT "PrimerPairs\tsense\tp_start\tp_end\tTm_value\tGC_content\tCross_exon\tantisense\tp_start\tp_end\tTm_value\tGC_content\tCross_exon\tProduct_length\n";
      my $cnt=1;
      foreach $s (@sense){
         foreach my $node (@{$cds_n{$gene}}){
            if($node>$p_start{$s}){
               if($node<$p_end{$s}){$flag{$s}="Y";}
            }
         }
         foreach $a (@anti){
            foreach my $node (@{$cds_n{$gene}}){
               if($node>$p_start{$a}){
                  if($node<$p_end{$a}){$flag{$a}="Y";}
               }
            }
            $deta_T=abs($T{$s}-$T{$a});
            $p_l=$p_end{$a}-$p_start{$s};
            if($p_l>$pl_max){;}
            elsif($p_l<$pl_min){;}
            elsif(abs($T{$s}-$T{$a})>2){;}
            elsif(cross_dimer($s,$a)){;}
            elsif($flag{$s}){
               printf OUT "PrimerPairs-$cnt\t$s\t$p_start{$s}\t$p_end{$s}\t$T{$s}\t%.2f\t$flag{$s}\t$a\t$p_start{$a}\t$p_end{$a}\t$T{$a}\t%.2f\t$flag{$a}\t$p_l\n",$GC{$s},$GC{$a};
               print OUT1">PrimerPairs-$cnt-sense-$p_start{$s}-$p_end{$s}\n$s\n>PrimerPairs-$cnt-antisense-$p_start{$a}-$p_end{$a}\n$a\n";
               $cnt++;
            }
            elsif($flag{$a} eq "Y"){
               printf OUT "PrimerPairs-$cnt\t$s\t$p_start{$s}\t$p_end{$s}\t$T{$s}\t%.2f\t$flag{$s}\t$a\t$p_start{$a}\t$p_end{$a}\t$T{$a}\t%.2f\t$flag{$a}\t$p_l\n",$GC{$s},$GC{$a};
               print OUT1">PrimerPairs-$cnt-sense-$p_start{$s}-$p_end{$s}\n$s\n>PrimerPairs-$cnt-antisense-$p_start{$a}-$p_end{$a}\n$a\n";
               $cnt++;
            }
            else{;}
         }
      }
      close OUT;
      close OUT1;
   }
}
################################################################### DO blast configuration #################
if($blast){
   print "Blastall for Primer designed...\n";
   system("formatdb -i $cds -p F -o T");
   opendir (DIR, "PrimerDesign");
   @dire = readdir DIR;
   for($i=0;$i<@dire;$i++){
      if(@dire[$i]=~/fasta/){
           @name=split(/_/,$dire[$i]);
           system("blastall -p blastn -d $cds -i ./PrimerDesign/$dire[$i] -o ./PrimerDesign/$name[0]-hit-$cds -e 1e-1 -m 8 -v 10 -b 10 -a 10");
      }
   }
}
############################################################################################################
print "\nSucceed running PrimerDesign!\n";
