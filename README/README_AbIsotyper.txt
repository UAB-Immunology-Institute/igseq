NAME
     AbIsotyper - Reads in FNA files which contain sequences in FASTA format and
     determines the Ab isotype that the sequence came from based on constant primer
     regions. Will search for exact matches to labels specified in AbIsotyper_labels.inf
     as well as truncated versions of these labels.  Also includes options for 
     sample reduction of sequences either percentage based or by a set number.

VERSION
     2013.03.29.1

DEPENDENCIES
     Text::LevenshteinXS
     Carp
     File::Path
     Parallel::ForkManager

SYNOPSIS
     perl AbIsotyper.pl <input directory> -diff:(1-X)
     perl AbIsotyper.pl <input directory> -diff:(1-X)  < -size:(1-X) OR -percent:(1-100) >
     perl AbIsotyper.pl <input directory> -diff:(1-X) -seed:(1-X)  < -size:(1-X) OR -percent:(1-100) >

DESCRIPTION
     This script reads a directory consisting of the FASTA formatted files and compares 
      the sequences to a list of labels in the AbIsotyper_labes.inf file. Matched labels 
      will add an identifier to the beginning of the FASTA comment line for that sequence.
      If no match is found a U will be used.  It will declare near matches up to the 
      specified cut off value.  That is if a sequence and a label are dissimilar by 1 bp 
      and the cut off is 2 it will label that sequence.

SETTINGS
     -diff:(1-X)
        Takes a number which will be used as a cutoff value for variance in similarity.
        Units is bp. (2 is suggested)

     -seed:(1-X)  OPTIONAL
            Seed for a random number generator. Set to a specific value for 
        reproducibility.
    
     -size:(1-X)  OPTIONAL
        Takes a number which will be used as a cutoff value for sample reduction.
    
     -percent:(1-100)  OPTIONAL
        Takes a percentage which will be used as a cutoff value for sample reduction.

OUTPUT
     Generates identical files with the same respective names and the input modified
      as outlined above.  As well as additional files based on options representing a
      portion of the sequences as specified in the settings.

CAVETS
     Don't try to specify a sequence label for U.

REVISIONS
     2012.10.22
      - Integrated ability for sample size reduction.
 
     2012.10.18
      - Fixed a bug causing truncated sequences not directly matching to be labeled as U 
       (~1% of sequences)
      - Parallelized files, should see significant speed up while doing multiple files
      - Memory utilization improved

AUTHOR
     Christopher Fucile
 
     Division of Allergy, Immunology and Rheumatology
     University of Rochester Medical Center
     Copyright (c) 2012, All Rights Reserved
 
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

