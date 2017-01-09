NAME
     lineage - Computes clonal lineages between sequences and groups them
                             together

SYNOPSIS
     perl lineage.pl <input file> <output file>

DESCRIPTION
     This script takes input from the merge.pl script and computes clonal
     lineages from the sequences presented. Clonal lineages are determined based on
     V, D and J similarity, allowing for one variation in the three, as well as
     CDR3 length and homology. CDR3 homology cut off levels can be adjusted in the
     lineage.inf settings file which should be contained in the same directory as
     lineage.pl.

     It contains the following settings:

        lineage_rep_size_cutoff:(0-X)
        cdr3_homology:(0-100)
        v_Variations:(y|n)
        d_Variations:(y|n)
        j_Variations:(y|n)
        useVDJ:(y|n)
        runIsotypes:(y|n)
        runWithIdenticals:(y|n)
        runWithoutIdenticals:(y|n)
        runV_genes:(y|n)

     Cdr3_homology can be a number ranging from 0 - 100 representing the percent
        homology of the cdr3 segments between sequences.

     Lineage_rep_size_cutoff can be a whole number ranging from 0 to the size of
        the largest lineage created by the program.  It controls what lineages are
        output to the '_LINEAGE_REPS.txt' file.
 
     vdjVariations can be either y or n indicating whether to allow variations in
        matching of v,d and j gene id's.
 
     useVDJ can be used to specify whether to consider v d and j gene assignments
        while clustering sequences.
 
     runIsotypes can be used to specify whether to run lineage on the individual
        isotypes as created by merge
 
     runWithIdenticals can be used to specify whether to run lineage on the
        identical output as created by merge
 
     runWithoutIdenticals can be used to specify whether to run lineage on the
        non-identical output as created by merge
 
     runV_genes can be used to specify whether to run lineage on the individual
        v_gene output as created by merge

AUTHOR
     Christopher Fucile & Alex Rosenberg
 
     Laboratory of Ignacio Sanz
     Division of Allergy, Immunology and Rheumatology
     University of Rochester Medical Center
     Copyright (c) 2013, All Rights Reserved
 
     Permission is hereby granted, free of charge, to any person obtaining a copy of this software
     and associated documentation files (the "Software"), to deal in the Software without restriction,
     including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
     and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
     subject to the following conditions:
 
     The above copyright notice and this permission notice shall be included in all copies or substantial
     portions of the Software.
 
     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
     NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

