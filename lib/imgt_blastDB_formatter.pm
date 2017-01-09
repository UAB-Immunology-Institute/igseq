package imgt_blastDB_formatter;

################################################################################
#
# imgt_blastDB_formatter.pm
#
# Parses fasta files downloaded from IMGT for use by NCBI C++ Toolkit
# Requires the NCBI C++ Toolkit's BLAST components to be installed
# Output is used as input for makeblastdb and subsequently for igblast(n|p)
#
# Author: Christopher Fucile
#
# Division of Allergy, Immunology and Rheumatology
# University of Rochester Medical Center
#
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the "Software"), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# Copyright (c) 2013, All Rights Reserved
#
# Usage:
#
#   check_database()
#
################################################################################


use strict;
use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(check_database);
use LWP::Simple;
use HTTP::Status qw(is_error);
use Carp;
use Data::Dumper;
use File::Path qw(make_path remove_tree);
use Cwd;

our $VERSION = '2.0.1';


################################
#   Constants for URL addresses
################################

use constant IMGT_URL => 'http://imgt.org/download/GENE-DB/';

open my $INF, "<". "IMGT_DATABASES.txt" or die "Error opening databse file list. $!\n";
my @db_files = <$INF>;
chomp @db_files;
my ($AA_GAPLESS_FILE) 		= grep {/AA.*?without/i} @db_files;
my ($NT_ALL_GAPLESS_FILE) 	= grep {/nt.*?without.*?all/i} @db_files;
my ($NT_INFRAME_GAPLESS_FILE) = grep {/nt.*?without.*?inframe/i} @db_files;
print "$AA_GAPLESS_FILE\n$NT_ALL_GAPLESS_FILE\n$NT_INFRAME_GAPLESS_FILE\n";

use constant RELEASE => 'RELEASE';
use constant DB_DIR => getcwd.'/database/';
use constant CUSTOM_DB_DIR => getcwd.'/custom_database_fasta/';

################################
#           Subroutines
################################


################################
#   Check which release version
#   of the database exists,
#   update if necessary.

sub check_database{
    make_path(DB_DIR);

    open my $IFH_SPECIES, "< species.txt" or die "Error opening species file $!\n";
    my @species = <$IFH_SPECIES>;
    chomp @species;
    close $IFH_SPECIES;
    
    print "Checking local files to IMGT...\n";

    my ($nt_custom_db, $aa_custom_db) = check_custom_db(\@species);

    if (-r DB_DIR.RELEASE){
        # if there is a release file check the version
        
        open my $INF, "<". DB_DIR.RELEASE or die "Error checking database version $!\n";
        my $old_release = <$INF>;
        chomp($old_release);
        close $INF;
        croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.RELEASE, DB_DIR.RELEASE));
        open $INF, "<". DB_DIR.RELEASE;
        my $new_release = <$INF>;
        chomp($new_release);
        close $INF;
        
        opendir my $IND, DB_DIR or die "Error opening directory $!\n";
        my @databases = readdir $IND;
        closedir $IND;
        
        my $new_specie = 0;
        foreach my $specie (@species){
            $specie =~ s/\s/_/g;
            $new_specie = 1 if !grep { /$specie/ } @databases;
            print "New species found in list, re-making databases\n" if $new_specie;
            last if $new_specie;
        }
        
        # files are up to date, exit

        if ($new_release eq $old_release && !$new_specie ){
            return "Files are up to date.\n".
                    "No new species in list.\n";
        }

    }

#remove_tree(DB_DIR, {keep_root => 1});
    
    # if no release file download new versions of files
croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.RELEASE, DB_DIR.RELEASE));

    # no error downloading occured, files out of date, update the database

    #gapped files removed due to makeblastdb not being able to handle them.
    #croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.AA_GAP_FILE, DB_DIR.AA_GAP_FILE));
    #croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.NT_INFRAME_GAP_FILE, DB_DIR.NT_INFRAME_GAP_FILE));
    print "Downloading files from IMGT...\n";
croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.$AA_GAPLESS_FILE, DB_DIR.$AA_GAPLESS_FILE));
croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.$NT_ALL_GAPLESS_FILE, DB_DIR.$NT_ALL_GAPLESS_FILE));
croak "Error downloading file$!\n" if is_error(getstore(IMGT_URL.$NT_INFRAME_GAPLESS_FILE, DB_DIR.$NT_INFRAME_GAPLESS_FILE));

    #gapped files removed due to makeblastdb not being able to handle them.
    #makeblastdb(format_file(DB_DIR.AA_GAP_FILE));
    #makeblastdb(format_file(DB_DIR.NT_INFRAME_GAP_FILE));
    
    makeblastdb( format_file( DB_DIR.$AA_GAPLESS_FILE, \@species, $aa_custom_db));
    makeblastdb( format_file( DB_DIR.$NT_ALL_GAPLESS_FILE, \@species, $nt_custom_db));
    makeblastdb( format_file( DB_DIR.$NT_INFRAME_GAPLESS_FILE, \@species));
    
}


sub check_custom_db{
    my $species = shift;
    
    my @nt_results;
    my @aa_results;
    
    opendir my $ODH, CUSTOM_DB_DIR or die "Error opening custom database directory $!\n";
    my @custom_dbs = grep {/\.fasta$/} readdir $ODH;
    closedir $ODH;
    
    foreach my $db (@custom_dbs){
        my $db_file = $db;
        $db =~ s/\_nt.*//ig;
        $db =~ s/\_/ /ig;
        die "Error: Can't find species $db from custom database in list.\n".
            "Check spelling.\n" if !grep { /$db/ } @{$species};
        if ($db_file =~ /nt/i){
            open my $INF, "<". CUSTOM_DB_DIR.$db_file or die "Error opening custom db file $db_file $!\n";
            my @tmp = <$INF>;
            close $INF;
            push @nt_results, \@tmp;
        }elsif( $db_file =~ /AA/i){
            open my $INF, "<". CUSTOM_DB_DIR.$db_file or die "Error opening custom db file $db_file $!\n";
            my @tmp = <$INF>;
            close $INF;
            
            push @aa_results, \@tmp;
        }else{
            die "Error: type of custom database not defined properly.\n";
        }
    }
    return (\@nt_results, \@aa_results);

}


sub makeblastdb {
    my @db_files = @_;
    foreach my $file(@db_files){
        my $dbname = $file;
        $dbname =~ s/IMGTGENEDB-ReferenceSequences.fasta-//i;
        my $title = $dbname;
        $title =~ s/.*\///;
        my $type;
        $type = "nucl" if $title =~ /nt/i;
        $type = "prot" if $title =~ /aa/i;
        my $command = "makeblastdb -in $file -out $dbname -title $title -dbtype $type -hash_index -parse_seqids";
        print "\n$command\n\n";
        system($command);
    }
}

sub format_file {
    my @results;
    my $input_file = shift;
    my $species = shift;
    
    my $custom_db = shift;

    # open fasta file from imgt
    open my $IFH, "< $input_file" or die "Error opening $input_file $!\n";
    $input_file =~ s|.*?\/database\/||;

    die "Error in FASTA file\nCheck formatting\n" if <$IFH> !~ /^\>/;
    seek $IFH, 0, 0;
    my @fasta_file = <$IFH>;
    close $IFH;
    
    
    if (defined($custom_db)){
        if (scalar(@{$custom_db}) > 0){
            foreach my $db (@{$custom_db}){
                push @fasta_file, @{$db};
            }
        }
    }

    foreach my $type(@{$species}){
        my $type_name = $type;
        $type_name =~ s/\s/\_/;
        
        open my $OFH_v, ">".DB_DIR."imgt_$type_name\_v_$input_file" or die "Error opening $type_name\_v $!\n";
        open my $OFH_d, ">".DB_DIR."imgt_$type_name\_d_$input_file" or die "Error opening $type_name\_d $!\n";
        open my $OFH_j, ">".DB_DIR."imgt_$type_name\_j_$input_file" or die "Error opening $type_name\_j $!\n";

        push @results, (DB_DIR."imgt_$type_name\_v_$input_file",
        DB_DIR."imgt_$type_name\_d_$input_file",
        DB_DIR."imgt_$type_name\_j_$input_file");
        
        my $j = 0;
        while($j <= $#fasta_file){
            # iterate through fasta array and allow the sequence portion to be
            # in multiple subsequent cells.

            my $x = 1;# number of lines (cells) the sequence takes up
            if($fasta_file[$j] =~ /.*?$type.*V-REGION.*?/i ){
                $fasta_file[$j] =~ s/^\>(.*?)\|(.*?)\|.*?\|.*?\|.*?\|.*?\|.*?\|(.*?)\|.*/\>lcl\|$2_$1_$3/;

                print $OFH_v $fasta_file[$j];
                while($fasta_file[$j+$x] !~ /^\>/){ 
                    print $OFH_v $fasta_file[$j+$x];
                    $x++;
                }
            }elsif($fasta_file[$j] =~ /.*?$type.*D-REGION.*?/i ){
                $fasta_file[$j] =~ s/^\>(.*?)\|(.*?)\|.*?\|.*?\|.*?\|.*?\|.*?\|(.*?)\|.*/\>lcl\|$2_$1_$3/;

                print $OFH_d $fasta_file[$j];
                while($fasta_file[$j+$x] !~ /^\>/){
                    print $OFH_d $fasta_file[$j+$x];
                    $x++;
                }
            }elsif($fasta_file[$j] =~ /.*?$type.*J-REGION.*?/i ){
                $fasta_file[$j] =~ s/^\>(.*?)\|(.*?)\|.*?\|.*?\|.*?\|.*?\|.*?\|(.*?)\|.*/\>lcl\|$2_$1_$3/;

                print $OFH_j $fasta_file[$j];
                while($fasta_file[$j+$x] !~ /^\>/){
                    print $OFH_j $fasta_file[$j+$x];
                    $x++;
                }
            }
            $j += $x;
        }

        close $OFH_v;
        close $OFH_d;
        close $OFH_j;
        print "Created Databases:\n\t".
            DB_DIR."imgt_$type_name\_v_$input_file\n\t".
            DB_DIR."imgt_$type_name\_d_$input_file\n\t".
            DB_DIR."imgt_$type_name\_j_$input_file\n";
    }
    return @results;
}


sub usage {
    print "Takes the following arguments:\n";
    print "Name of fasta file to be formatted.\n";
}

=head1 NAME

 imgt_blastDB_formatter.pl - Parses fasta files downloaded from IMGT for use by
  NCBI C++ Toolkit.  Requires the NCBI C++ Toolkit's BLAST components to be 
  installed. Output is used as input for makeblastdb and subsequently for 
  igblast(n|p).

=head1 VERSION
 
 1.0.0
 
=head1 DEPENDENCIES
 
 LWP::Simple
 HTTP::Status
 Carp
 
=head1 SYNOPSIS
 
 perl imgt_blastDB_formatter.pl
 
=head1 DESCRIPTION
 
 Downloads files from http://imgt.org/download/GENE-DB/ and attempts to format
 them using makeblastdb which is apart of the NCBI C++ Toolkit.  All files are
 intended to be in /usr/local/database/  Only tested on OS X 10.8.2
 
=head1 SETTINGS
 
 N/A
 
=head1 OUTPUT
 
 Generates identical files with the same respective names and the input modified
 as outlined above.  As well as additional files based on options representing a
 portion of the sequences as specified in the settings.
 
=head1 CAVETS
 
 
 
=head1 REVISIONS
 
 
 
=head1 AUTHOR
 
 Christopher Fucile
 
 Division of Allergy, Immunology and Rheumatology
 University of Rochester Medical Center
 
 Copyright (c) 2012, All Rights Reserved
 
=cut