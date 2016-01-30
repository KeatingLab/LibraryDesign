# This file was modified by Gevorg Grigoryan from the downloaded version
# to conform to the NU MathProg modeling language, which is a subset of
# AMPL.
# Date: 05/04/05

# THIS SOFTWARE IS LICENSED UNDER THE GNU PUBLIC LICENSE VERSION 2.
# See LICENSE file available at http://compbio.cs.princeton.edu/scplp/cks_lp/
# or downloaded with this software. If you use this software please cite the
# paper mentioned at the above website.

#
# The problem structure: the number of variable residues (num_posn), 
# total number of codons (num_nodes), the set of positions (POSN) 
# and the set of rotamers (V).
#
param num_posn integer >= 2;
param num_nodes integer >= num_posn;
set POSN := 1..num_posn;
set V := 1..num_nodes;

#
# The structure of the positions: convert posn_size array into sets of 
# nodes C
#
param posn_size {i in POSN} integer > 0;
check: sum {i in POSN} posn_size[i] = num_nodes;

set C {i in POSN} := {u in V : sum {j in POSN : j < i} posn_size[j] < u and u <= sum {j in POSN : j <= i} posn_size[j]};
check: (setof{i in POSN, j in C[i]}(j) symdiff V) within {};

#
# The self-energies (costV)
#
param costVTOT {V} default 0.0;
param costVSCORE {V} default 0.0;
#param costVTRS {V} default 0.0;

#
# The variables
#
var X {V} binary >= 0;  # node variables

# minimization
#minimize energy: (sum {v in V} costV[v] * X[v]) +
#                 (sum {(u,v) in E} costE[u,v] * Y[u,v]);
#
#minimize energy: (sum {v in V} costVT[v] * X[v]);
#
maximize score: (sum {v in V} (costVSCORE[v]) * X[v]);


subject to column {i in POSN}:
    sum {u in C[i]} X[u] = 1;

# additional constrains

subject to totalsize: sum {v in V} costVTOT[v] * X[v], <= 7.0;

#subject to trs: sum {v in V} costVTRS[v] * X[v], <= 1000.0;

end;
