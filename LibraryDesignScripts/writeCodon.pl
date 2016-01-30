#Created by T. Scott Chen, Keating lab, MIT
# read in a list of amino acids that need to be included at each designed site (required - sp)
my $nfile = shift;
# read in a file that contains score for each amino acid (preferred - nd)
my $efile = shift;
# the output file
my $ofile = shift;

# read the codon file
my $codon_combo_file = 'codon_combos.txt';

# read in the codon combinations
my $C = {};
open (my $cfh, "<", $codon_combo_file) or die "Cannot open $codon_combo_file\n";
while (<$cfh>) {
  chomp;
  my @lines = split;
  $C->{$lines[0]}->{size} = trim($lines[1]); #codon size hash
  my $aa_str;
  my $count_hash = {};
  for (my $i = 2; $i < scalar(@lines); $i = $i + 2) {
    my $aa = trim($lines[$i]);
    $aa_str .= $aa;
    my $count = trim($lines[$i+1]);
    $count_hash->{$aa} = $count;
  }
  $C->{$lines[0]}->{str} = $aa_str;  #aas encoded by codon
  $C->{$lines[0]}->{count} = $count_hash;  #number of times each aa is encoded in codon
}
close($cfh);

# read in the file that contains information about required amino acids at each site
open (my $nfh, "<", $nfile) or die "Cannot open $nfile\n";
my $N = {};
while (<$nfh>) {
  chomp;
  if (/^([A-Z][0-9]+) ([A-Z]+)/) {  #match position letter, number, and aa string
    my $site = $1;   #site=position name
    my $aa_str = $2;  #aa_str=required aa at that position
     
    my @aas = split "", $aa_str;
    foreach my $aa (@aas) {
      $N->{$site}->{$aa} = 1;  #$N holds required residues, by position (site)
    }
  }
}
close($nfh);

# read in the preferred file that contains information about scores for point mutants 
open (my $efh, "<", $efile) or die "Cannot open $efile\n";
my $index = 0;
my $E = {};
my $site;
my $aa;
while (<$efh>) {
  chomp;
  if (/^([A-Z][0-9]+)/) {   #matches position name/number
    my $site = $1;
    my @lines = split;
    shift @lines;
    for (my $i = 0; $i < scalar(@lines); $i = $i + 2) {
      my $aa = $lines[$i];   #get preferred aa
      my $score = $lines[$i+1];  #get 'score' (0 or 1) for that aa
      $E->{$site}->{$aa} = $score;  #$E holds scores for preferred aa, organized by position
    } 
  }
}
close($efh);

open (my $ofh, ">", $ofile) or die "Cannot write into $ofile\n";
foreach my $site (sort resid_sort keys %$E) {    #sort position names
  printf $ofh "$site\n";

  my $sols = {};
  foreach my $cname (keys %$C) {  #cname = codon name
    my $aa_str = $C->{$cname}->{str};
    my $codon_size = $C->{$cname}->{size};  #number of trinucleotides
    my $count = $C->{$cname}->{count};   #number of times aa is encoded by codon
    my $include = 1;
    foreach my $aa (keys %{$N->{$site}}) {  #required aa
      if (!defined($count->{$aa})) {   #Only include codons that encode all required aa
        $include = 0;
        last; # break out of loop
      }
    }
    if ($include == 0) {
      next;
    }

    my $score = 0.0;   #score = number of preferred (non-disrputive) aa
    my $useful_size = 0.0;  #number of trinucleotides encoding preferred aa
    foreach my $aa (keys %{$E->{$site}}) {  #preferred residues
      if (defined($count->{$aa})) {   #if required?
        $score += $E->{$site}->{$aa};    #add score for preferred residue (0 or 1)
        if ($E->{$site}->{$aa} > 0) {  #if residue is non-disruptive (preferred) (score is 1)   
          $useful_size += $count->{$aa};   #add number of times aa is encoded by codon
        }  
      }
    }
#$sols holds codons that (meet criteria?)
    $sols->{$cname}->{perc} = $useful_size / $codon_size;  #percent of codons encoding preferred aa
    $sols->{$cname}->{codon} = $cname; #codon name
    $sols->{$cname}->{size} = $codon_size;  #number of trinucleotides encoded
    $sols->{$cname}->{aa_str} = $aa_str;  #list of aa encoded
    $sols->{$cname}->{score} = $score;  #number of preferred aa encoded
  }

  # Build output. Each possible codon must be compare to all previously accepted codons. 
  my $PO = {};
  NEWSOL: # codon not yet in PO
  foreach my $codon1 (keys %$sols) {
    my $h1 = $sols->{$codon1};
    OLDSOL: # codon already in PO, but being compared to possible new codon
    foreach my $codon2 (keys %$PO) {
      my $h2 = $sols->{$codon2};
      if (($h1->{aa_str} eq $h2->{aa_str}) and ($h1->{size} > $h2->{size})) {
        next NEWSOL; #if two codons encode the same aas and codon1 has more trincleotides, then do not include codon1
      }

      if (($h1->{size} == $h2->{size}) and ($h1->{score} == $h2->{score})) { #if codon size and number of preferred residues are equal for two codons
        if ($h1->{aa_str} =~ /Z/ and $h2->{aa_str} !~ /Z/) {  #if codon1 includes a stop codon and codon2 doesn't, then skip codon1
          next NEWSOL;
        } elsif ($h1->{aa_str} !~ /Z/ and $h2->{aa_str} =~ /Z/) { #if codon1 doesn't include a stop codon and codon2 does, then delete codon2 and add codon1 to list
          delete $PO->{$codon2};
          next OLDSOL;
        }
      } elsif (($h1->{size} >= $h2->{size}) and ($h1->{score} <= $h2->{score})) { #if codon1 size is larger than codon2 and codon1 score is smaller
        next NEWSOL;  #then do not include codon1
      } elsif (($h1->{size} <= $h2->{size}) and ($h1->{score} >= $h2->{score})) { #if codon1 is smaller than codon2 and codon1 has a larger score
        delete $PO->{$codon2};  #then remove codon2
        next OLDSOL;
      }

    }
    $PO->{$codon1} = 1;
  }

  foreach my $c_name (keys %$PO) {   #Print output to output file
    printf $ofh "%s %s %s %.2f %.2f\n", $c_name, $C->{$c_name}->{str}, $sols->{$c_name}->{size}, $sols->{$c_name}->{score}, $sols->{$c_name}->{perc};
  }  
}

close($ofh);

sub resid_sort {   #Position sorting function: sort by position number. Commented out code sorts by wt AA name. 
#    $chain_a = substr $a, 0, 1;
    $res_a = substr $a, 1;
#    $chain_b = substr $b, 0, 1;
    $res_b = substr $b, 1;

#    $chain_a cmp $chain_b
#      ||
    $res_a <=> $res_b
}

sub trim {   #line trimming function
   my $string = shift;
   $string =~ s/^\s+//;
   $string =~ s/\s+$//;
   return $string;
}    

