# This file was modified from version downloaded at at http://compbio.cs.princeton.edu/scplp/cks_lp/
# to conform to the NU MathProg modeling language, which is a subset of
# AMPL.

# The problem structure: the number of variable residues (num_posn), 
# total number of codons (num_nodes), and the set of positions (POSN). 
# and the set of codons (V).

param num_posn integer >= 2;
param num_nodes integer >= num_posn;
set POSN := 1..num_posn;
set V := 1..num_nodes;


# The structure of the positions: convert posn_size array into sets of 
# nodes C

param posn_size {i in POSN} integer > 0;
check: sum {i in POSN} posn_size[i] = num_nodes;

set C {i in POSN} := {u in V : sum {j in POSN : j < i} posn_size[j] < u and u <= sum {j in POSN : j <= i} posn_size[j]};
check: (setof{i in POSN, j in C[i]}(j) symdiff V) within {};


# The self-energies (costV)

param costVTOT {V} default 0.0;
param costVSCORE {V} default 0.0;


# The variables

var X {V} binary >= 0;  # node variables

# Objective

maximize score: (sum {v in V} (costVSCORE[v]) * X[v]);

subject to column {i in POSN}:
    sum {u in C[i]} X[u] = 1;    

# additional constrains

subject to totalsize: sum {v in V} costVTOT[v] * X[v], <= 7.0;

end;
