package circos_config;

###################################################################################################
#
# circos_config.pm
#   - package to generate output from compareLineage for Circos visualization
#
#
# by Christopher Fucile
########################
#
# Division of Allergy, Immunology and Rheumatology
# University of Rochester Medical Center
# Copyright (c) 2013, All Rights Reserved
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
#
# USAGE:
#
#    use circos_config;
#
###################################################################################################

use warnings;
use strict;
use Data::Dumper;

our $VERSION = '2013.03.29.1';

sub print_config {
	my $path		= shift;
    my $output_dir  = shift;
    my $popsizes    = shift;
    my $cutoffs		= shift;
    my $lin_cutoffs = shift;
    my $uniq_iso	= shift;

    my %settings;
    open my $input_settings, "< $path/config/circos.inf";
    while(<$input_settings>){
        chomp;
        my $tmp_line = $_;
        next if $tmp_line =~ /^\#|^\s|^\n|^$/;

        my @line = split( ":", $tmp_line);
        $settings{$line[0]} = $line[1];
    }

    print_ideogram(     $output_dir, $settings{'label_size'}, $settings{'label_radius'});
    print_lineage(      $output_dir, $popsizes, $settings{'numSeq'}, $settings{'image_radius'}, $lin_cutoffs);
    print_links(        $output_dir, $cutoffs, $settings{'connections_track'});
    print_mutation_axis($output_dir );
    print_plots(        $output_dir, $settings{'isotype_track'}, $settings{'connections_track'}, $settings{'v_gene_mut_track'}, $cutoffs, $uniq_iso);
    print_ticks(    	$output_dir, $settings{'tick_spacing'}, $settings{'tick_size'}, $settings{'tick_label_size'}, $settings{'tick_label_spacing'});
}


sub print_ticks {
    my $output_dir 			= shift;
    my $tick_spacing		= shift;
    my $tick_size  			= shift;
    my $tick_label_size     = shift;
    my $tick_label_spacing  = shift;
    
    open my $ticks, "> $output_dir/ticks.conf" or die "Error opening file $!\n";
    
    print $ticks    "show_ticks          = yes\n".
                    "show_tick_labels    = yes\n".
                    "<ticks>\n".
                    "radius           = 1.001r\n".
                    "color            = black\n".
                    "thickness        = 2p\n".
                    "multiplier       = 1\n".
                    "format           = %d\n".
                    "<tick>\n".
                    "spacing        = ".$tick_spacing."u\n".
                    "size           = ".$tick_size."p\n".
                    "</tick>\n".
                    "<tick>\n".
                    "spacing        = ".$tick_label_spacing."u\n".
                    "size           = 30p\n".
                    "show_label     = yes\n".
                    "label_size     = ".$tick_label_size."p\n".
                    "label_offset   = 10p\n".
                    "format         = %d\n".
                    "</tick>\n".
                    "</ticks>\n";
    close $ticks;
    

}



sub print_plots {
    my $output_dir          = shift;
    my $isotype_track       = shift;
    my $connections_track   = shift;
    my $v_gene_mut_track    = shift;
    my $cutoffs				= shift;
    my $uniq_iso			= shift;
    
    my @colors = ('blue','dorange','green','purple','yellow','red','black');
    my $counter = 0;
    my $colorstring = "";
    map {$colorstring .= $colors[$counter].","; $counter++; $counter = 0 if $counter == 6;} sort keys %{$uniq_iso};
    $colorstring =~ s/,$//;
    my $output =    "<plots>\n".
                    "<plot>\n".
                    "show = no\n".
                    "type = heatmap\n".
                    "stroke_thickness=1\n".
                    "stroke_color = black\n".
                    "file = circos_v_variability_data.txt\n".
                    "color = ylorrd-9-seq\n".
                    "r0=.975r\n".
                    "r1=.975r+25p\n".
                    "</plot>\n".
                    "<plot>\n".
                    "show = no\n".
                    "type = heatmap\n".
                    "stroke_thickness=1\n".
                    "stroke_color = black\n".
                    "file = circos_d_variability_data.txt\n".
                    "color = ylorrd-9-seq\n".
                    "r0=.95r\n".
                    "r1 = .95r+25p\n".
                    "</plot>\n".
                    "<plot>\n".
                    "show= no\n".
                    "type = heatmap\n".
                    "stroke_thickness=1\n".
                    "stroke_color = black\n".
                    "file = circos_j_variability_data.txt\n".
                    "color = ylorrd-9-seq\n".
                    "r0=.925r\n".
                    "r1=.925r+25p\n".
                    "</plot>\n";

    #isotype track
        $output .=  "<plot>\n";
        $output .=  $isotype_track eq 'y' ? "show = yes\n" : "show = no\n";
		$output .=  "type = histogram\n".
                    "file = circos_lineage_isotype_data.txt\n".
                    "r0 =1.25r\n".
                    "r1 =1.34r\n".
                    "min = 0\n".
                    "max = 100\n".
                    "color = vlgrey\n".
                    #"fill_color =  blue\n".
                    #			  IgA	,  IgG,   IgM , Unknown
                    "#			  ".join(",", sort keys %{$uniq_iso})."\n".
                    "fill_color = $colorstring\n".
                    "fill_under = yes\n".
                    "thickness = 1\n".
                    "sort_bin_values = no\n".
                    "extend_bin = no\n".
                    "</plot>\n";


    #connections track
        $output .=  "<plot>\n";
        $output .=	$connections_track eq 'y' ? "show = yes\n" : "show = no\n";
		$output .=	"type = heatmap\n".
                    "stroke_thickness=1\n".
                    "stroke_color = lgrey\n".
                    "file = circos_connections_data.txt\n".
                    "color = ylgnbu-9-seq\n".
                    "r0 = .9999r\n".
                    "r1 = .95r\n".
                    "</plot>\n";
    
    #v gene mutation track
        $output .=  "<plot>\n";
        $output .=	$v_gene_mut_track eq 'y' ? "show = yes\n" : "show = no\n";
		$output .=  "type = line\n".
                    "file = circos_V_mutations_data.txt\n".
                    "thickness = 1\n".
                    "orientation = out\n".
                    "color = vdgrey\n".
                    "r0 = 1.14r\n".
                    "r1 = 1.24r\n".
                    "z=1\n".
                    "min = 0\n".
                    "max = 30\n".
                    "background_color = vvlgrey\n".
                    "background_stroke_color = black\n".
                    "<<include mutations_axis.conf>>\n".
                    "</plot>\n".
                    "<plot>\n".
                    "type = line\n".
                    "file = circos_V_avg_mutation_data.txt\n";
       $output .=	$v_gene_mut_track eq 'y' ? "show = yes\n" : "show = no\n";
       $output .=	"thickness = 3\n".
                    "orientation = out\n".
                    "color = red\n".
                    "r0 = 1.14r\n".
                    "r1 = 1.24r\n".
                    "min = 0\n".
                    "max = 30\n".
                    "z=1\n".
                    "</plot>\n";
    
	#lineage band highlights
    $output .= "<plot>\n".
                "type = highlight\n".
                "file = circos_lineage_bands.txt\n".
                "r1   = dims(ideogram,radius_outer)\n".
                "r0   = dims(ideogram,radius_inner)\n".
                "stroke_thickness = 2\n".
                "stroke_color     = black\n".
                "<rules>\n".
                "<rule>\n".
                "condition  =";
                my $count = 1;
                foreach my $pop (keys %{$cutoffs}){
                	$output .= " (var(chr1) eq \"$pop\" && var(end1) > $cutoffs->{$pop}) ";
                	if ($count != scalar keys %{$cutoffs}){ $output .= " ||" }else{ $output .= "\n" }
                	$count++;
                }
	$output .= "fill_color = eval(sprintf(\"%s_a3\",var(fill_color)))\n".
                "flow       = continue\n".
                "</rule>\n".
                "<rule>\n".
                "condition = ";
				$count = 1;
                foreach my $pop (keys %{$cutoffs}){
                	$output .= " (var(chr1) eq \"$pop\" && var(end1) <= $cutoffs->{$pop}) ";
                	if ($count != scalar keys %{$cutoffs}){ $output .= " ||" }else{ $output .= "\n" } 
                	$count++;
                }
	$output .= "stroke_thickness =0\n".
                "fill_color = vlgrey\n".
                "show      = yes\n".
                "</rule>\n".
                "</rules>\n".
                "</plot>\n".
                "</plots>\n";
    
    open my $plots, "> $output_dir/plots.conf" or die "Error opening file $!\n";
    print $plots $output;
    close $plots;

}


sub print_mutation_axis {
    my $output_dir = shift;
    
    open my $mutation_axis, "> $output_dir/mutations_axis.conf" or die "Error opening file $!\n";
    print $mutation_axis    "<axes>\n".
                            "show = yes\n".
                            "thickness = 1\n".
                            "color = lgrey\n".
                            "<axis>\n".
                            "spacing = 2\n".
                            "color = grey\n".
                            "</axis>\n".
                            "<axis>\n".
                            "position = 0r\n".
                            "color = vdgrey\n".
                            "</axis>\n".
                            "<axis>\n".
                            "position = 1r\n".
                            "color = vdgrey\n".
                            "</axis>\n".
                            "</axes>\n";
    close $mutation_axis;
}



sub print_links {
    my $output_dir  	  = shift;
    my $cutoffs     	  = shift;
    my $connections_track = shift;
    
    my $line;
    if($cutoffs){
        $line = "<rule>\n".
                "condition  = ";
                my $count = 1;
                foreach my $pop (keys %{$cutoffs}){
                	$line .= " ((var(chr1) eq \"$pop\" && var(end1) > $cutoffs->{$pop}) ||";
					$line .= "  (var(chr2) eq \"$pop\" && var(end2) > $cutoffs->{$pop})) ";
                	if ($count != scalar keys %{$cutoffs}){ $line .= " ||" }else{ $line .= "\n" }
                	$count++;
                }
		$line .="show       = yes\n".
                "z          = 10\n".
                "thickness  = 5\n".
                "</rule>\n".
                "<rule>\n".
                "condition  = ";
                $count = 1;
                foreach my $pop (keys %{$cutoffs}){
                	$line .= " ((var(chr1) eq \"$pop\" && var(end1) <= $cutoffs->{$pop}) ||";
					$line .= "  (var(chr2) eq \"$pop\" && var(end2) <= $cutoffs->{$pop})) ";
                	if ($count != scalar keys %{$cutoffs}){ $line .= " ||" }else{ $line .= "\n" }
                	$count++;
                }
		$line .="show       = yes\n".
                "ribbon     = no\n".
                "color      = vlgrey\n".
                "z          = 1\n".
                "thickness  = 1\n".
                "</rule>\n";
    }else{
        $line = "\n";
    }
    
    open my $links_conf, "> $output_dir/links.conf" or die "Error opening file $!\n";
    
    print $links_conf   "<links>\n".
                        "<link>\n".
                        "file               = circos_lineage_links.txt\n";
    print $links_conf   "radius             = .95r\n" if $connections_track eq 'y';
    print $links_conf   "radius             = .99r\n"   if $connections_track ne 'y';
    print $links_conf   "bezier_radius      = 0.1r\n".
                        "color              = grey\n".
                        "thickness          = 5\n".
                        "ribbon             = yes\n".
                        "flat               = yes\n".
                        "stroke_color       = vdgrey\n".
                        "stroke_thickness   = 1\n".
                        "<rules>\n".
                        "<rule>\n".
                        "condition          = var(intrachr)\n".
                        "show               = no\n".
                        "</rule>\n".
                        $line.
                        "</rules>\n".
                        "</link>\n".
                        "</links>\n";
    close $links_conf;
    
}


sub print_lineage {
    my $output_dir   = shift;
    my $popsizes     = shift;
    my $lineage_size = shift;
    my $image_radius = shift;
    my $lin_cutoffs = shift;

    my $line;
    if($lineage_size){
    	if ($lineage_size =~ 'd'){
    		$line = "chromosomes_display_default    = no\n".
					"chromosomes = ";
			foreach my $pop (keys %{$lin_cutoffs}){
				$line .= "$pop:". ($lin_cutoffs->{$pop}) ."-);";
			}
    	}else{
			$line = "chromosomes_display_default    = no\n".
					"chromosomes = ";
			foreach my $pop (keys %{$popsizes}){
				$line .= "$pop:". ($popsizes->{$pop} - $lineage_size) ."-);";
			}
		}
    }else{
        $line = "chromosomes_display_default    = yes\n";
    }
    open my $lineage_conf, "> $output_dir/lineage.conf" or die "$!\n";
    
    print $lineage_conf "karyotype                      = circos_lineage_populations.txt\n".
                        "chromosomes_units              = 1\n".
                        $line."\n".
                        "<<include ideogram.conf>>\n".
                        "<<include ticks.conf>>\n".
                        "<<include links.conf>>\n".
                        "<<include plots.conf>>\n".
                        "<colors>\n".
                        "chr1* = red\n".
                        "chr2* = orange\n".
                        "chr3* = green\n".
                        "chr4* = blue\n".
                        "chr5* = yellow\n".
                        "<<include etc/colors.conf>>\n".
                        "</colors>\n".
                        "<image>\n".
                        "<<include etc/image.conf>>\n".
                        "24bit = yes\n".
                        "radius*=".$image_radius."p\n".
                        "</image>\n".
                        "<<include etc/colors_fonts_patterns.conf>>\n".
                        "<<include etc/housekeeping.conf>>\n";
    close $lineage_conf;
}


sub print_ideogram {
    my $output_dir   = shift;
    my $label_size   = shift;
    my $label_radius = shift;
    
    open my $output_file, "> $output_dir/ideogram.conf" or die "Error opening file $!\n";
    print $output_file  "<ideogram>\n".
                        "<spacing>\n".
                        "default                = 0.02r\n".
                        "break                  = 0.1u\n".
                        "axis_break_at_edge     = yes\n".
                        "axis_break             = yes\n".
                        "axis_break_style       = 2\n".
                        "<break_style 1>\n".
                        "stroke_color           = black\n".
                        "fill_color             = blue\n".
                        "thickness              = 0.25r\n".
                        "stroke_thickness       = 2\n".
                        "</break>\n".
                        "<break_style 2>\n".
                        "stroke_color           = black\n".
                        "stroke_thickness       = 2\n".
                        "thickness              = 1.5r\n".
                        "</break>\n".
                        "</spacing>\n".
                        "radius                 = 0.66r\n".
                        "thickness              = 100p\n".
                        "fill                   = yes\n".
                        "stroke_color           = dgrey\n".
                        "stroke_thickness       = 1p\n".
                        "show_bands             = yes\n".
                        "fill_bands             = no\n".
                        "band_stroke_thickness  = 1\n".
                        "band_stroke_color      = vlgrey\n".
                        "band_transparency      = 5\n".
                        "show_label             = yes\n".
                        "label_font             = bold\n".
                        "label_radius           = dims(ideogram,radius) + ".$label_radius."r\n".
                        "label_with_tag         = yes\n".
                        "label_size             = $label_size\n".
                        "label_parallel         = yes\n".
                        "label_case             = upper\n".
                        "</ideogram>\n";
    close $output_file;

}

1;
__END__
=head1 NAME
 
 circos_config
 
=head1 SYNOPSIS
 
 use circos_config;
 
=head1 DESCRIPTION
 
 Generates configuration files for circos visualization based on output from compareLineage.
 
=head1 SETTINGS
 
 See circos_config.inf file.
 
=head1 OUTPUT
 
 circos_connections_data.txt
 
 circos_d_variability_data.txt
 
 circos_j_variability_data.txt
 
 circos_v_variability_data.txt
 
 circos_V_mutations_data.txt
 
 circos_V_avg_mutation_data.txt
 
 circos_lineage_populations.txt
 
 circos_lineage_links.txt
 
 circos_lineage_isotype_data.txt
 
 circos_lineage_bands.txt
 
=head1 CAVEATS
 
 
 
=head1 AUTHOR
 
 Christopher Fucile
 University of Rochester Medical Center
 Division of Allergy, Immunology & Rheumatology
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
 
=cut
