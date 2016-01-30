#Created by T. Scott Chen, Keating Lab, MIT

my $inputfile = shift;
my $output = shift;

my $ilp_model_file = "library_design.mod"; 

# modify the path to the glpsol binary as necessary
my $glpsol = "/usr/local/bin/glpsol";

open (my $ifh, "<", $inputfile) or die "Cannot open file $inputfile";

my @positions;
my @codons;

my %hit;
my @sublist = ();
my $total_spec_size = 0; 
while (<$ifh>) {    #Read codons_output input file
  if (/^([A-Z][0-9]+)/) { #match position name
    my $site = $1;  #store position names in $site
    push(@positions, $site); #store sites in @positions
    if (scalar(@sublist) > 0) {  #if sublist is populated?
      my @tl = @sublist;  
      push(@codons, \@tl)  #push elements of temporary list (from sublist) into @codons
    }
    @sublist = (); #empty sublist
    %hit = ();
  } elsif (/^[A-Z][A-Z][A-Z] /) {  #matches codon lines
    chomp;
    my @lines = split; #@lines holds elements of codon lines
    my $codon = trim($lines[0]);
    my $aa_str = trim($lines[1]);
    my $dna_size = trim($lines[2]); #number of trinucleotides
    my $score = trim($lines[3]);  #number of preferred (non-disrupt) residues
    my $perc = trim($lines[4]);  #percent of codons encoding preferred residues
    my $h = {};
    $h->{dna_size} = $dna_size;  #put stats into $h with keys
    $h->{score} = $score;
    $h->{aa_str} = $aa_str;
    $h->{perc} = $perc;
    $h->{codon} = $codon;

    my $label_str = "$dna_size" . "_$score";
    if (!defined($hit{$label_str}))    #this deals with codons of same size and score
    {
        $hit{$label_str} = 1;        
        push(@sublist, $h);   #this populates sublist with the codons and stats from $h
    }
  }  
}

if (scalar(@sublist) > 0) {
  my @tl = @sublist;
  push(@codons, \@tl)  #add elements of sublist to @codons
}

close($ifh);

my $L = {};

open (my $ofh, ">", $output) or die "Cannot write into $output\n";

#####################Writing temporary data file with parameters defined:
my $ilp_data_file = "ILP_temp.dat"; #temporary data file
my $length = scalar(@positions);  #length = number of positions

open (my $dfh, ">", $ilp_data_file) or die "Cannot write into $ilp_data_file\n";

my $ns = 0;
printf $dfh "data;\n\n";   #constructing temp data file
printf $dfh "param num_posn := %d ;\n", $length; #number of positions
for (my $i = 0; $i < $length; $i++) {  #for each position
  $ns += scalar(@{$codons[$i]});  #count total number of codon choices = $ns 
}
printf $dfh "param num_nodes := %d ;\n", $ns; #total number of codons across all positions
printf $dfh "param posn_size :=\n";
for (my $i = 0; $i < $length; $i++) { #print position sizes = number of codons at that position
  printf $dfh "%d %d\n", $i+1, scalar(@{$codons[$i]}); 
}
printf $dfh ";\n\n";

printf $dfh "param costVTOT :=\n";  #variable for DNA size
my $k = 1;
for (my $i = 0; $i < $length; $i++) {
  foreach my $h (@{$codons[$i]}) {   #$h contains codon stats
    my $total_size = log($h->{dna_size})/log(10); #print codon number followed by log size
    printf $dfh "$k $total_size\n";
    $k++; 
  }
}
printf $dfh ";\n\n";

printf $dfh "param costVSCORE :=\n";  #variable for score (number of preferred residues)
$k = 1;
for (my $i = 0; $i < $length; $i++) {
  foreach my $h (@{$codons[$i]}) {  #print codon number followed by log score (num preferred)
    my $score = log($h->{score})/log(10);
    printf $dfh "$k $score\n";
    $k++;
  }
}
printf $dfh ";\n\n";

printf $dfh "end;\n";
close($dfh);                
################End of writing temporary data file

my $ilp_output = "ILP_temp.ilpout"; #temp ILP output file
#RUN ILP: 
system("$glpsol --math $ilp_model_file --data $ilp_data_file --output $ilp_output");
#model_file =.mod, data_file=temp data file written above

################Processing the ILP output
my @consnames = ('score', 'totalsize'); #Need to add name for MSAscore
my @tcnames = @consnames;
my $name = shift @tcnames;
foreach my $n (@tcnames) {
  $name .= "|$n";
}
$name = "(" . $name . ")";
my @X;
my %E; 
my $i = 0;
open (my $ilpfh, "<", $ilp_output) or die "Cannot open file $ilp_output"; 
while (<$ilpfh>) {last if ($_ =~ /-----/);}  #end when you get to dashes?
foreach my $line (<$ilpfh>) {
  if ($line =~ / $name /) {
    my @arr = split " ", $line;
    $E{$1} = $arr[2];   #get score and total size numbers
  }
  next if ($line !~ /X\[/);  #look for codon list lines
  my @arr = split " ", $line;
  if ($i >= $ns) {die "Too many node variables (only $ns expected)!\n";}
  $X[$i] = $arr[3] + 0.0;   #holds whether or not a codon was chosen (0/1)
  $i++;
}
close($ilpfh);
  
if ($E{score} == 0) {
  printf $ofh "Feasible solution not found\n";
}

my @sol;
$k = 0;
for (my $i = 0; $i < $length; $i++) {
  foreach my $h (@{$codons[$i]}) {
    if ($X[$k]) {   #true if codon was chosen
      push(@sol, $h);   #get codon identity and stats from codons list
    }
    $k++;
  }
}
foreach my $n (@consnames) {
  $L->{$n} = $E{$n};
} 
$L->{sols} = \@sol;

my $total_aa_size = 1;
my $useful_perc = 1.0;
my $score = 1;

for (my $i = 0; $i < $length; $i++) {
  my $pos = $positions[$i];
  printf $ofh "$pos ";
  my $h = $L->{sols}->[$i];
  my $codon = $h->{codon};
  my $aa_str = $h->{aa_str};
  my $perc = $h->{perc};
  my $aa_size = length $aa_str;
  $total_aa_size = $total_aa_size * $aa_size; #iteratively multiplies number of aa in each codon
  $useful_perc = $useful_perc * $perc; #useful percent=product(percent of trinucleotides encoding preferred aa)
  $score = $score * $h->{score}; #score=product(number of preferred aa ecoded)
  printf $ofh "$codon $aa_str\n";
}


my $expscore= $E{score};   #10^ of this score = "Score", i.e., number of sequences composed entirely of preferred residues
my $total_dna_score= $E{totalsize};

printf $ofh "Total size in DNA sequences (10^_) = $total_dna_score\n";
printf $ofh "Total size in protein sequences = $total_aa_size\n";
printf $ofh "Score = %d\n", $score;
printf $ofh "Useful fraction = %.4f\n", $useful_perc;

close($ofh);

##########
#system("rm $ilp_data_file $ilp_output");  #Include to delete temporary input and output files
##########

sub trim {
   my $string = shift;
   $string =~ s/^\s+//;
   $string =~ s/\s+$//;
   return $string;
}
