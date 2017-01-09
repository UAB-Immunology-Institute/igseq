NAME
     fastaSubSample.pl 
       - Randomly selects a number of sequences from a FASTA file.

VERSION
     2013.03.29.1

DEPENDENCIES
     Parallel::ForkManager

SYNOPSIS
     perl fastaSubSample.pl <input directory>  < -size:(1-X) OR -percent:(1-100) >
     perl fastaSubSample.pl <input directory>  < -size:(1-X) OR -percent:(1-100) > -seed:(1-X)

DESCRIPTION
     This script reads a directory consisting of the FASTA formatted files and will
     randomly select a subset

SETTINGS
     -seed:(1-X)  OPTIONAL
            Seed for a random number generator. Set to a specific value for 
        reproducibility.
    
     -size:(1-X)  OPTIONAL
        Takes a number which will be used as a cutoff value for sample reduction.
    
     -percent:(1-100)  OPTIONAL
        Takes a percentage which will be used as a cutoff value for sample reduction.

OUTPUT
     Generates identical files with the the number of sequences added to the name.

CAVETS
REVISIONS
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

