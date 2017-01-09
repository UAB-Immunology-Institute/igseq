NAME
     merge - merges primary data files from IHMMune and V-Quest for sequencing analysis

SYNOPSIS
     perl merge454.pl <IHMMune file> <Vquest 2 file> <Vquest 4 file> <Vquest 8 file> <output file>

DESCRIPTION
     This script takes files from IHMMune and V-Quest and combines their results into one file.  V-Quest files needed
     include file numbers 2, 4 and 8.  Files are combined based on identical sequencing IDs.  Output is a tab delimited
     file consisting of the IHMMune entry first then the V-Quest entry after for each sequence.  Options for including
     sequences that did not return a significant D gene alignment can be specified from the settings file as well as a
     cut off for the difference in the number of basepairs for sequences to not be included.  See below for more details.

     Settings:
        See README_merge file and merge.inf


     Notes on including alignments with no D gene match:
        - If you are planning on allowing variability of v,d and j segments when running the lineage script
            it is advised that alignments with no D gene match are included in the merged output.
            This allows these sequences to match based on all other criteria besides the D gene segment.

AUTHOR
     Alex Rosenberg & Christopher Fucile

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


