NAME
     comparelineage - Compares lineages between populations and generates a list of links between them

SYNOPSIS
     perl comparelineage.pl <output name>

DESCRIPTION
     This script reads a directory consisting of the output from multiple lineage454.pl runs and compares
     the results.  It generates links between lineages within different populations when the sequences
     within both lineages have equal CDR3 length, XX% CDR3 homology and have matching IDs for at least two
     of the V, D and J segments.  Settings allow for the matching of lineage groups based solely on cdr3 length
     and homology as well.
 
     To say two lineages are matched is to say a specific percent of the smaller of the
     two groups meet the cut off value specified by lineageSimilarity in terms of the similarity measures used.
     For example, given the settings: cdr3_homology=90, lineageSimilarity=95, useVDJ=n; then two lineages
     are matched if 95% of the sequences between them have 90% cdr3_homology while ignoring VDJ gene names.

SETTINGS
     cdr3_homology:(0-100)
        percent identity between cdr3 fragments when comparing sequences within lineages
        calculated using hamming distance

     lineageSimilarity:(0-100)
        -percent identity between two lineages in regards to v, d and j gene id's and cdr3 from sequences.
        -determines the number of matches between two lineages divides the number by the smaller
            of the two lineages and multiply by 100 to give percent
 
     useVDJ:'y'|'n'
        allows use or ignores VDJ gene names for matching of lineages
 
     X_Variations:'y'|'n'
        specifies whether to allow specific variations for v, d or j.  useVDJ must be set to 'y' for these
        to have an affect.

     lineage_rep_size_cutoff:(0-X)
        cut off for tab delimited file output consisting of matched lineages in file LARGEST_GROUPED_LINEAGES.txt

OUTPUT
     comparelineage454_log.txt
        log file of the run

     GROUPED_LINEAGE_COUNTS.txt
        Includes Grouped Lineage ID, V-gene, Total number of lineages, Total number
        of sequences, Number of lineages per population and number of sequences per 
        population
 
     GROUPED_LINEAGES_###.FASTA
        FASTA format of sequences in grouped lineage ###

     GROUPED_LINEAGES_###.INFO
        merge454.pl formatted style of sequences in grouped lineage ###

     LARGEST_GROUPED_LINEAGES.txt
        tab delimited file with specific columns from .INFO files only lineages above
        the settings cutoff included

     POPULATION_MATCHES.txt
        tab delimited file of populations and their corresponding number of matches

     RUN_SETTINGS.txt
        settings that were used for the run
 
CAVEATS
     If sequences with NO_DGENE_ALIGNMENT are allowed to be matched into lineages though merge454.pl
     & lineage454.pl, and lineageSimilarity = 100 in settings for comparelineage454.pl while using VDJ comparisons,
     any lineage or singletons with sequences containing NO_DGENE_ALIGNMENT will not be matched.  This is due
     to the fact that sequences are not allowed to use NO_DGENE_ALIGNMENT to be matched upon at all.

     Lineages can not have identical prefix names on them

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

