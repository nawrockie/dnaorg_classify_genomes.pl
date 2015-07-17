#!/usr/bin/env perl
# EPN, Mon May 18 11:34:28 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# hard-coded-paths:
my $idfetch       = "/netopt/ncbi_tools64/bin/idfetch";
my $esl_fetch_cds = "/panfs/pan1/dnaorg/programs/esl-fetch-cds.pl";
#my $esl_fetch_cds = "/home/nawrocke/notebook/15_0518_dnaorg_virus_compare_script/wd-esl-fetch-cds/esl-fetch-cds.pl";
my $usearch       = "/panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/2015.05.29/usearch";

my $do_uc    = 0; 
my $uc_id    = undef;
my $df_uc_id = 0.9;
my $nuc      = 0;   # number of uclust commands written to uclust script

# The definition of $usage explains the script and usage:
my $usage = "\ndnaorg_classify_genomes.pl\n";
$usage .= "\t<directory created by dnaorg_fetch_dna_wrapper>\n";
$usage .= "\t<list file with all accessions>\n";
$usage .= "\n"; 
$usage .= " This script compares genomes from the same species based\n";
$usage .= " mostly on files containing parsed information from a\n";
$usage .= " 'feature table' file which must already exist. That file is\n";
$usage .= " created by running 'dnaorg_fetch_dna_wrapper.pl -ftable' and\n";
$usage .= " subsequently 'parse-ftable.pl'.\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -t <f>      : fractional length difference threshold for mismatch [default: 0.05]\n";
$usage .= "  -f          : fetch CDS and genomes to per-class fasta files [default: do not]\n";
$usage .= "  -s          : use short names: all CDS seqs will be named as sequence accession\n";
#$usage .= "  -gene       : use 'gene:gene' qualifier value for defining labels\n";
$usage .= "  -product    : add CDS 'product' qualifier value to output sequence files deflines\n";
$usage .= "  -protid     : add CDS 'protein_id' qualifier value to output sequence files deflines\n";
$usage .= "  -codonstart : add CDS 'protein_id' qualifier value to output sequence files deflines\n";
$usage .= "  -uc         : create a shell script to run uclust jobs for all fasta files\n";
$usage .= "  -uc_id <f>  : with -uc, set sequence identity cutoff to <f> [default: $df_uc_id]\n";
$usage .= "  -noexp      : do not print verbose explanation of how classes are defined to stdout\n";
$usage .= "\n";

my ($seconds, $microseconds) = gettimeofday();
my $start_secs    = ($seconds + ($microseconds / 1000000.));
my $executable    = $0;
my $be_verbose    = 1;
my $df_fraclen    = 0.05;
my $fraclen       = undef;
my $do_fetch      = 0; # changed to '1' if -f enabled
my $do_shortnames = 0; # changed to '1' if -s enabled
my $do_gene       = 0; # changed to '1' if -gene enabled
my $do_product    = 0; # changed to '1' if -product enabled
my $do_protid     = 0; # changed to '1' if -protid enabled
my $do_codonstart = 0; # changed to '1' if -codonstart enabled
my $do_noexp      = 0; # changed to '1' if -noexp enabled

&GetOptions( "t=s"        => \$fraclen, 
             "f"          => \$do_fetch,
             "s"          => \$do_shortnames, 
#             "gene"       => \$do_gene,
             "product"    => \$do_product,
             "protid"     => \$do_protid,
             "codonstart" => \$do_codonstart,
             "uc"         => \$do_uc, 
             "uc_id"      => \$uc_id,
             "noexp"      => \$do_noexp) || 
    die "Unknown option";

if(scalar(@ARGV) != 2) { die $usage; }
my ($dir, $listfile) = (@ARGV);

# determine name of output xlist1 and xlist2 files
$dir =~ s/\/*$//; # remove trailing '/' if there is one
my $outdir     = $dir;
my $outdirroot = $outdir;
$outdirroot =~ s/^.+\///;
my $llist_outfile  = "$outdir/$outdirroot.labels.txt";
my $xlist1_outfile = "$outdir/$outdirroot.xlist1.txt";
my $xlist2_outfile = "$outdir/$outdirroot.xlist2.txt";
open(OUTL,  ">" . $llist_outfile)  || die "ERROR unable to open $llist_outfile for writing"; 
open(OUTX1, ">" . $xlist1_outfile) || die "ERROR unable to open $xlist1_outfile for writing"; 
open(OUTX2, ">" . $xlist2_outfile) || die "ERROR unable to open $xlist2_outfile for writing"; 
my @xlist1_output_A    = (); # array of output lines about 'x' characters which indicate CDS that do not map to an obvious ref CDS
my @xlist1_accn_A      = (); # same size as @xlist1_output_A, accessions that pertain to xlist_output_A lines
my @xlist1_label_str_A = (); # same size as @xlist1_output_A, label strings that pertain to xlist_output_A lines
my $nlines_x2 = 0;


# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if(defined $fraclen) { 
  $opts_used_short .= "-t $fraclen";
  $opts_used_long  .= "# option:  setting fractional length threshold to $fraclen [-t]\n";
}
if($do_fetch) { 
  $opts_used_short .= "-f ";
  $opts_used_long  .= "# option:  fetching CDS and genomes to per-class fasta files [-f]\n";
}
if($do_shortnames) { 
  $opts_used_short .= "-s ";
  $opts_used_long  .= "# option:  outputting CDS names as accessions [-s]\n";
}
if($do_gene) { 
  $opts_used_short .= "-gene ";
  $opts_used_long  .= "# option:  using \'gene:gene\' annotation to define classes [-gene]\n";
}
if($do_product) { 
  $opts_used_short .= "-product ";
  $opts_used_long  .= "# option:  adding \'product\' qualifier values to output sequence file deflines [-product]\n";
}
if($do_protid) { 
  $opts_used_short .= "-protid ";
  $opts_used_long  .= "# option:  adding \'protein_id\' qualifier values to output sequence file deflines [-protid]\n";
}
if($do_codonstart) { 
  $opts_used_short .= "-codonstart ";
  $opts_used_long  .= "# option:  adding \'codon_start\' qualifier values to output sequence file deflines [-codonstart]\n";
}
if($do_uc) { 
  $opts_used_short .= "-uc ";
  $opts_used_long  .= "# option:  outputting uclust shell script [-uc]\n";
}
if(defined $uc_id) { 
  $opts_used_short .= "-uc_id ";
  $opts_used_long  .= "# option:  uclust cluster identity set to $uc_id [-uc_id]\n";
}
if(defined $do_noexp) { 
  $opts_used_short .= "-noexp ";
  $opts_used_long  .= "# option:  verbose explanation of class definition not being output to stdout [-noexp]\n";
}

# check for incompatible option values/combinations:
if(defined $fraclen && ($fraclen < 0 || $fraclen > 1)) { 
  die "ERROR with -t <f>, <f> must be a number between 0 and 1."; 
}
if((! $do_uc) && (defined $uc_id)) { 
  die "ERROR -uc_id requires -uc";
}

# set default values if user didn't specify otherwise on the command line
if(! defined $fraclen) { 
  $fraclen = $df_fraclen; 
}
if(! defined $uc_id) { 
  $uc_id = $df_uc_id;
}

###############
# Preliminaries
###############
# check if the $dir exists, and that it contains a .gene.tbl file, and a .length file
if(! -d $dir)      { die "ERROR directory $dir does not exist"; }
if(! -s $listfile) { die "ERROR list file $listfile does not exist, or is empty"; }
my $dir_tail = $dir;
$dir_tail =~ s/^.+\///; # remove all but last dir
my $gene_tbl_file  = $dir . "/" . $dir_tail . ".gene.tbl";
my $cds_tbl_file   = $dir . "/" . $dir_tail . ".CDS.tbl";
my $length_file    = $dir . "/" . $dir_tail . ".length";
my $out_root = $dir . "/" . $dir_tail;
#if(! -s $gene_tbl_file) { die "ERROR $gene_tbl_file does not exist."; }
if(! -s $cds_tbl_file)  { die "ERROR $cds_tbl_file does not exist."; }
if($do_gene) { 
  if(! -s $gene_tbl_file)  { die "ERROR $gene_tbl_file does not exist."; }
}
if(! -s $length_file)   { die "ERROR $length_file does not exist."; }

# output banner
my $script_name = "dnaorg_classify_genomes.pl";
my $script_desc = "Classify genomes based on comparison of their annotations";
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# $script_name: $script_desc\n");
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# command: $executable $opts_used_short $dir $listfile\n");
printf("# date:    %s\n", scalar localtime());
if($opts_used_long ne "") { 
  print $opts_used_long;
}

#########################################################
# Define our explanation of how classes are defined. This
# gets output to several of the output files
#########################################################
my @class_explanation_A = ();
push(@class_explanation_A, sprintf("# ==================================================================================================\n"));
push(@class_explanation_A, sprintf("# Explanation of how classes are defined and labels are determined:\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# Classes are defined by unique labels, such that:\n"));
push(@class_explanation_A, sprintf("# - all accessions in the same class have the same label\n"));
push(@class_explanation_A, sprintf("# - no two accessions in different classes have the same label\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# A label is defined by comparing all CDS in an accession to all CDS in the reference accession.\n"));
push(@class_explanation_A, sprintf("# (The reference accession is the first one listed in this file and is the first accession in the input\n"));
push(@class_explanation_A, sprintf("# <list file> (command line argument 2))\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# An accession's label is (2*<n>+1) characters long, where <n> is the number of CDS that accession has\n"));
push(@class_explanation_A, sprintf("# The first <n> characters are either uppercase letters, lowercase letters or single digit numbers, as explained below.\n"));
push(@class_explanation_A, sprintf("# The <n>+1 character is always '|'.\n"));
push(@class_explanation_A, sprintf("# The final <n> characters (characters <n>+2..(2*<n>+1)) are either '=', '!' or '?', as explained below.\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# For character <i> in 1..<n> in the label, the specific character used is determined by\n"));
push(@class_explanation_A, sprintf("# comparing CDS <i> for the accession to each of the reference CDS to see if the following match or not:\n"));
push(@class_explanation_A, sprintf("#   (1) CDS:product annotation\n"));
push(@class_explanation_A, sprintf("#   (2) length %s\n", ($fraclen > 0) ? "(within $fraclen fraction of reference length)" : "(exact match)"));
push(@class_explanation_A, sprintf("#   (3) strand\n"));
push(@class_explanation_A, sprintf("# If the CDS matches in (1) and (2), a capital letter will be assigned to describe that CDS, e.g. 'A'.\n"));
push(@class_explanation_A, sprintf("# If the CDS matches in (1) but not (2), a lowercase letter will be assigned to describe that CDS, e.g. 'a'\n"));
push(@class_explanation_A, sprintf("# If the CDS matches in (2) but not (1), a number will be assigned to describe that CDS, e.g. '1'\n"));
push(@class_explanation_A, sprintf("# All such CDS are considered 'mappable' to a reference CDS.\n"));
push(@class_explanation_A, sprintf("# The specific character used depends on which reference CDS it matches, 'A', 'a' or '1' is used for a match to\n"));
push(@class_explanation_A, sprintf("# the first CDS, 'B' 'b' or '2' for a match to the second and so on.\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# If, however, the CDS matches neither (1) nor (2) when compared to ALL REFERENCE CDS, an 'x' or 'X' will\n"));
push(@class_explanation_A, sprintf("# be assigned at position <i> of the label, and such CDS are considered 'unmappable' to reference CDS.\n"));
push(@class_explanation_A, sprintf("# Unmappble CDS are listed separately in the following two output files:\n"));
push(@class_explanation_A, sprintf("# - $xlist1_outfile\n"));
push(@class_explanation_A, sprintf("# - $xlist2_outfile\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# The final <n> characters of the label are determined based on the strand of each CDS.\n"));
push(@class_explanation_A, sprintf("# For mappable CDS, if the strand of the CDS matches the strand of the mapped reference CDS, the character\n"));
push(@class_explanation_A, sprintf("# will be a '='. If it does not match it will be a '!'.\n"));
push(@class_explanation_A, sprintf("# For unmappable CDS, the character will always be a '?'.\n"));
push(@class_explanation_A, sprintf("#\n"));
push(@class_explanation_A, sprintf("# Example labels (that may make it easier to understand the label definition rules):\n"));
push(@class_explanation_A, sprintf("# E1.  'ABCD|====':   all 4 CDS match all 4 of the reference CDS in the proper order for all 3 properties: CDS:product, length and strand\n"));
push(@class_explanation_A, sprintf("# E2.  'ABCD|==!=':   same as E1, but CDS 3 is on the opposite strand of reference CDS 3\n"));
push(@class_explanation_A, sprintf("# E3.  'ABcD|====':   same as E1, but CDS 3 has a length that is %s\n", ($fraclen > 0) ? "not within $fraclen fraction of the length of reference CDS 3" : "not identical to the length of reference CDS 3"));
push(@class_explanation_A, sprintf("# E4.  'ABcD|==!=':   same as E3, but CDS 3 is on the opposite strand as reference CDS 3\n"));
push(@class_explanation_A, sprintf("# E5.  'AB3D|====':   same as E1, but CDS 3 has a different CDS:product annotation from the reference CDS\n"));
push(@class_explanation_A, sprintf("# E6.  'AB3D|==!=':   same as E5, but CDS 3 is on the opposite strand as reference CDS 3\n"));
push(@class_explanation_A, sprintf("# E7.  'ABXD|====':   same as E1, but CDS 3 has a different CDS:product annotation AND is %s\n", ($fraclen > 0) ? "not within $fraclen fraction of the length of reference CDS 3" : "not identical to the length of reference CDS 3"));
push(@class_explanation_A, sprintf("# E8.  'ABXD|==!=':   same as E7, but CDS 3 is on the opposite strand as reference CDS 3\n"));
push(@class_explanation_A, sprintf("# E9.  'ABD|===':     only 3 CDS are annotated, and they correspond to reference CDS 1, 2 and 4, in that order, on matching strands\n"));
push(@class_explanation_A, sprintf("# E10. 'BACD|====':   same as E1, but CDS 1 and 2 are in the opposite order in the accession relative to the reference\n"));
push(@class_explanation_A, sprintf("# E11. 'xxx|???':     only 3 CDS are annotated, and none match either length or CDS:product with any of the reference\n"));
push(@class_explanation_A, sprintf("# E12. 'ABCDx|====?': same as E1, but an extra CDS is annotated, after the 4th reference-matching CDS, and it does\n"));
push(@class_explanation_A, sprintf("#                     not match any of the reference CDS in either CDS:product or length\n"));
push(@class_explanation_A, sprintf("# E13. 'AxBCD|=?===': same as E12, but extra CDS is the 2nd CDS, not the 5th\n"));
push(@class_explanation_A, sprintf("# ==================================================================================================\n"));
push(@class_explanation_A, sprintf("#\n"));

#####################
# parse the list file
#####################
my @accn_A = (); # array of accessions
open(IN, $listfile) || die "ERROR unable to open $listfile for reading"; 
my $waccn = 0; # max length of all accessions
while(my $accn = <IN>) { 
  chomp $accn;
  stripVersion(\$accn); # remove version
  push(@accn_A, $accn);
  if(length($accn) > $waccn) { $waccn = length($accn); }
}
close(IN); 

my $head_accn = $accn_A[0];


##################################
# parse the table and length files
##################################
my %gene_tbl_HHA = ();  # Data from .gene.tbl file
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %cds_tbl_HHA = ();   # Data from .cds.tbl file
                        # hash of hashes of arrays, 
                        # 1D: key: accession
                        # 2D: key: column name in gene ftable file
                        # 3D: per-row values for each column
my %totlen_H = (); # key: accession, value length read from length file

parseLength($length_file, \%totlen_H);

parseTable($gene_tbl_file, \%gene_tbl_HHA);
parseTable($cds_tbl_file, \%cds_tbl_HHA);

###########################################################################
# Partition accessions into classes; each class has a unique label string
# This requires the first of two passes through all accessions.
###########################################################################
my $wstrand_str = 0;         # for output formatting, width of max expected length strand string
my $wlabel_str = 0;          # for output formatting, width of max expected length label string
my %label_str2idx_H = ();    # key: label string, class number for this label string 
my %idx2label_str_H = ();    # key: class number,  value: label string
my %ct_label_str_H = ();     # key: label string, value: number of accessions in the class defined by this label string 
my %fa_label_str_H = ();     # key: strand string, value: name of output fasta file for the class defined by this label string 
my %out_label_str_HA = ();   # key: strand string, value: array of output strings for the class defined by this label string 
my @out_fetch_gnm_A = ();    # out_fetch_gnm_A[$c]: for class $c+1, input for idfetch to fetch full genome seqs
my @ct_fetch_gnm_A = ();     # ct_fetch_gnm_A[$c]:  for class $c+1, number of genomes to fetch
my @out_fetch_cds_AA = ();   # out_fetch_cds_AA[$c][$i]: for class $c+1, gene $i+1, input for esl-fetch-cds.pl
my @ct_fetch_cds_AA = ();    # ct_fetch_cds_AA[$c][$i]:  for class $c+1, gene $i+1, number of sequences to fetch 
my %accn2label_str_H = ();   # key: $accn, value: $label_str for that accn
my $label_str = "";          # a label string, e.g. 'ABCD:===='
my $strand_str;              # +/- string for all CDS for an accession: e.g. '+-+': 1st and 3rd CDS are + strand, 2nd is -
my $class_idx = undef;       # class index for current accession
my $nclasses = 0;            # total number of classes
my @ncds_per_class_A = ();   # number of CDS per class
my $max_ncds = 0;            # max number of CDS for any class

# variables related to a reference accession
my $ref_accn      = undef; # changed to <s> with -ref <s>
my $ref_label_str = undef; # label string for reference accn

my $uc_script = $out_root . ".uclust.sh";

# reference information on reference accession, first accession read in ntlist file
my $ref_ncds          = 0;  # number of CDS in reference
my $ref_strand_str    = ""; # strand string for reference 
my @ref_cds_len_A     = (); # [0..$i..$ref_ncds-1]: length of each reference CDS
my @ref_cds_len_tol_A = (); # [0..$i..$ref_ncds-1]: length tolerance, any gene that is within this fraction of the lenght of the ref gene is a match
my @ref_cds_coords_A  = (); # [0..$i..$ref_ncds-1]: CDS coords for reference
my @ref_cds_product_A = (); # CDS:product qualifier data for reference 

my $ncds = 0;              # number of CDS
my $npos = 0;              # number of CDS on positive strand
my $nneg = 0;              # number of CDS on negative strand
my $nunc = 0;              # number of CDS on uncertain strand
my $nbth = 0;              # number of CDS on both strands
my @cds_len_A = ();        # [0..$i..$ncds-1] length of CDS $i
my @cds_coords_A = ();     # [0..$i..$ncds-1] coords of CDS $i
my @cds_product_A = ();    # [0..$i..$ncds-1] CDS:product annotation for CDS $i
my @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
my @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)
my $do_desc = ($do_product || $do_protid || $do_codonstart) ? 1 : 0; # '1' to create a description to add to defline of fetched sequences, '0' not to

my %x_product_len_accn_HHA   = (); # key 1: $product, key 2: $length, value: array of  accessions with this product and length that do not match a ref CDS
my %x_product_len_accn_ct_HH = (); # key 1; $product, key 2: $length, value: number of accessions with this product and length that do not match a ref CDS
my %x_product_accn_ct_H      = (); # key 1; $product, value: number of accessions with this product of ANY length that do not match a ref CDS

# Pass 1 through all accessions, to determine their labels:
my $naccn = scalar(@accn_A);
for(my $a = 0; $a < $naccn; $a++) { 
  my $accn = $accn_A[$a];
  my $nxlist1_this_accn = 0;

  # get ref accession information
  if($a == 0) { # set the reference accession
    $ref_accn = $accn;
    if($do_gene) { 
      if(! exists ($gene_tbl_HHA{$ref_accn})) { die "ERROR no gene information stored for reference accession"; }
    }
    else { 
      if(! exists ($cds_tbl_HHA{$ref_accn})) { die "ERROR no CDS information stored for reference accession"; }
    }
    if($do_gene) { 
      (undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats(\%gene_tbl_HHA, $ref_accn);
      getLengthStatsAndCoordStrings(\%gene_tbl_HHA, $ref_accn, \@ref_cds_len_A, \@ref_cds_coords_A);
      getQualifierValues(\%gene_tbl_HHA, $ref_accn, "gene", \@ref_cds_product_A);
    }
    else { 
      (undef, undef, undef, undef, undef, $ref_strand_str) = getStrandStats(\%cds_tbl_HHA, $ref_accn);
      getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $ref_accn, \@ref_cds_len_A, \@ref_cds_coords_A);
      getQualifierValues(\%cds_tbl_HHA, $ref_accn, "product", \@ref_cds_product_A);
    }
    $ref_ncds = scalar(@ref_cds_len_A);
    # determine length tolerances

    print("#\n");
    printf("# Reference accession: $ref_accn\n");
    print("#\n");
    for(my $rc = 0; $rc < $ref_ncds; $rc++) { 
      $ref_cds_len_tol_A[$rc] = $fraclen * $ref_cds_len_A[$rc];
      printf("# Reference CDS product %d: %-30s reference-length: %6d   length-range-for-match: %8.1f to %8.1f\n", $rc+1, $ref_cds_product_A[$rc], $ref_cds_len_A[$rc], ($ref_cds_len_A[$rc] - $ref_cds_len_tol_A[$rc]), ($ref_cds_len_A[$rc] + $ref_cds_len_tol_A[$rc]));
    }
    $wlabel_str  = (2*$ref_ncds+1) * 2;
    $wstrand_str = $ref_ncds * 2;
    if($wlabel_str  < length("label"))         { $wlabel_str  = length("label"); }
    if($wstrand_str < length("strand-string")) { $wstrand_str = length("strand_string"); }
  }    

  # sanity checks
  if($do_gene) { 
    if($a == 0 && (! exists $gene_tbl_HHA{$accn})) { die "ERROR didn't read any gene table information for first accession in $listfile: $accn\n"; } 
  }
  else { 
    if($a == 0 && (! exists $cds_tbl_HHA{$accn})) { die "ERROR didn't read any CDS table information for first accession in $listfile: $accn\n"; } 
  }
  if(! exists $totlen_H{$accn}) { die "ERROR accession $accn does not exist in the length file $length_file"; }

  # set defaults that will stay if we don't have any CDS information
  $ncds = 0; 
  $npos = 0;
  $nneg = 0;
  $nunc = 0;
  $nbth = 0; 
  $strand_str = "";
  @cds_len_A = ();
  @cds_coords_A = ();
  @cds_product_A = ();    # will remain empty unless $do_product is 1 (-product enabled at cmdline)
  @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
  @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)

  if((   $do_gene  && exists ($gene_tbl_HHA{$accn})) || 
     ((! $do_gene) && exists ($cds_tbl_HHA{$accn}))) { 
    if($do_gene) { 
      ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%gene_tbl_HHA, $accn);
      getLengthStatsAndCoordStrings(\%gene_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
      getQualifierValues(\%gene_tbl_HHA, $accn, "gene", \@cds_product_A);
    }
    else { 
      ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
      getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
      getQualifierValues(\%cds_tbl_HHA, $accn, "product", \@cds_product_A);
    }
    if($do_protid) { 
      if($do_gene) { die "ERROR -protid and -gene are incompatible, currently"; }
      getQualifierValues(\%cds_tbl_HHA, $accn, "protein_id", \@cds_protid_A);
    }
    if($do_codonstart) { 
      if($do_gene) { die "ERROR -codonstartd and -gene are incompatible, currently"; }
      getQualifierValues(\%cds_tbl_HHA, $accn, "codon_start", \@cds_codonstart_A);
    }
  }

  # determine which reference CDS each current accn's CDS corresponds to, if any:
  my $match_product = 0;
  my $match_len     = 0;
  my $match_strand  = 0;

  my $label_sep = "|"; # separate the label string from the label-strand string
  # determine label string
  my @lchar_A        = (); # [0..$c..$ncds-1] array of $lchar       characters for each CDS, concatenated together these will make up part of $label_str
  my @lstrand_char_A = (); # [0..$c..$ncds-1] array of $lstrand_str characters for each CDS, concatenated together these will make up part of $lstrand_str

  # for each CDS in current accession, compare against all reference CDS looking for a match
  for(my $c = 0; $c < $ncds; $c++) { 
    my $cur_strand = substr($strand_str, $c, 1);
    my $ref_strand;
    my $lstrand_char = "";
    my $lchar        = "";
    my $in_order     = 0;
    my $found_match_product = 0;
    # does this match any of the reference genes?
    my $cur_lchar = 'A';
    # for each reference CDS
    for(my $rc = 0; $rc < $ref_ncds; $rc++) { 
      if($rc > 0) { $cur_lchar++; }
      if(! $found_match_product) { 
        $ref_strand    = substr($ref_strand_str, $rc, 1);
        $match_product = ($ref_cds_product_A[$rc] eq $cds_product_A[$c])  ? 1 : 0;
        $match_len     = (abs($ref_cds_len_A[$rc] - $cds_len_A[$c]) <= $ref_cds_len_tol_A[$rc]) ? 1 : 0;
        $match_strand  = $cur_strand eq $ref_strand ? 1 : 0;
        # check for case 
        if($match_product) { 
          $found_match_product = 1; # all ref CDS must be unique product names, so we know we've found our label if product matches
          if($c == $rc) { $in_order = 1; }
          if($match_len) { 
            # Either case 1 or case 2:
            # case 1: 'A:=' same name, same length, same strand (can only match 1 ref CDS, b/c ref CDS names must not contain duplicates)
            # case 2: 'A:!' same name, same length, diff strand (can only match 1 ref CDS, b/c ref CDS names must not contain duplicates)
            $lchar = $cur_lchar;
            $lstrand_char = ($match_strand) ? "=" : "!";
          }
          else { # ($match_product) but (!$match_len)
            # Either case 3 or case 4:
            # case 3: 'a:=' same name, diff length, same strand (can only match 1 ref CDS, b/c ref CDS names must not contain duplicates)
            # case 4: 'a:!' same name, diff length, diff strand (can only match 1 ref CDS, b/c ref CDS names must not contain duplicates)
            $lchar = $cur_lchar;
            $lchar =~ tr/A-Z/a-z/;
            $lstrand_char = ($match_strand) ? "=" : "!";
          } 
        } # end of 'if($match_product)'           
        elsif($match_len) { # (!$match_product) but ($match_len)
          # Either case 5 or case 6:
          # case 5: '1:=' diff name, same length, same strand (could happen twice, take highest indexed gene we find this for)
          # case 6: '1:!' diff name, same length, diff strand (could happen twice, take highest indexed gene we find this for)
          # NOTE: we may want to start by examining $c == $rc instead of rc 0..nrc-1, so we will assign this case to 
          #       the in order gene instead of just the first match in the unlikely case of multiple matches
          $lchar        = $cur_lchar;
          $lstrand_char = ($match_strand) ? "=" : "!";
        }
      }
    }
    if($lchar eq "") { 
       # case 7: x:? diff name, diff length from all ref CDS
      $lchar        = "x";
      $lstrand_char = "?";
      # find nearest length gene
      my $min_diff = $cds_len_A[$c];
      my $min_rc   = -1;
      for(my $rc2 = 0; $rc2 < $ref_ncds; $rc2++) { 
        if(abs($ref_cds_len_A[$rc2] - $cds_len_A[$c]) < $min_diff) { 
          $min_diff = abs($ref_cds_len_A[$rc2] - $cds_len_A[$c]);
          $min_rc   = $rc2;
        }
      }
      # save this to an array, wait until full string has been determined then print them all out, with full string
      push(@xlist1_output_A, sprintf("%4d  %-40s  %6d  %10d  %5s  %-40s  %10d  %10d",
                                $c+1,
                                $cds_product_A[$c], 
                                $cds_len_A[$c], 
                                $min_rc+1, 
                                ($min_rc == $c) ? "yes" : "no",
                                $ref_cds_product_A[$min_rc],
                                $ref_cds_len_A[$min_rc],
                                $min_diff));
      push(@xlist1_accn_A,      $accn);
      $nxlist1_this_accn++;

      my $x_product = $cds_product_A[$c];
      my $x_len     = $cds_len_A[$c];
      if(! exists $x_product_len_accn_HHA{$x_product}) { 
        %{$x_product_len_accn_HHA{$x_product}} = (); 
        %{$x_product_len_accn_ct_HH{$x_product}} = ();
      }
      if(! exists $x_product_len_accn_HHA{$x_product}{$x_len}) { 
        @{$x_product_len_accn_HHA{$x_product}{$x_len}} = ();
      }
      push(@{$x_product_len_accn_HHA{$x_product}{$x_len}}, $accn);
      # printf("xct just incremented $x_product $x_len by 1\n");
      $x_product_len_accn_ct_HH{$x_product}{$x_len}++;
      $x_product_accn_ct_H{$x_product}++;
    }
    #printf("lchar: $lchar\n");
    #printf("lstr:  $lstrand_str\n");

    push(@lchar_A, $lchar);
    push(@lstrand_char_A, $lstrand_char);
  } # end of 'for(my $c = 0; $c < $ncds; $c++)'

  # Now that we've examined all the CDS, we can further differentiate between
  # some of the case 7s. Specifically we're checking if the following 3 conditions
  # hold: 
  # - same number of CDS as reference
  # - all but 1 of the CDS match reference as case 
  #   1, 2, 3 or 4 in the proper order
  # - the single non-matching CDS is a case 7.
  # 
  # If these all hold we change the single case 7 to 
  # a capital 'X' and it's associated strand character 
  # from a '?' to a '=' or '!'.
  # 
  my $exp_uc_char = 'A';
  my $exp_lc_char = 'a';
  my $ncase7 = 0;     # number of case 7s for this accession
  my $case7_idx = -1; # set to highest index of a case7 in accession
  my $in_order = 1; # this will remain 1 until we find that $accn's CDS annotations are not in same order as reference
  # e.g. 'ABxD' would be 'in order', 'ADxB' would not.
  for(my $c = 0; $c < $ncds; $c++) { 
    if($lchar_A[$c] eq 'x') { 
      $ncase7++; 
      $case7_idx = $c;
    }
    elsif($lchar_A[$c] ne $exp_uc_char && $lchar_A[$c] ne $exp_lc_char) {
      $in_order = 0;
    }
    $exp_uc_char++;
    $exp_lc_char++;
  }
  if($ncase7 == 1 && $in_order) { 
    $lchar_A[$case7_idx] = "X";
    $lstrand_char_A[$case7_idx] = (substr($strand_str, $case7_idx, 1) eq substr($ref_strand_str, $case7_idx, 1)) ? '=' : '!';
  }

  $label_str = "";
  for(my $c = 0; $c < $ncds; $c++) { $label_str .= $lchar_A[$c]; }
  $label_str .= $label_sep;
  for(my $c = 0; $c < $ncds; $c++) { $label_str .= $lstrand_char_A[$c]; }

  # printf("accn: $accn $label_str\n");

  if(! exists $label_str2idx_H{$label_str}) { 
    $nclasses++;
    $label_str2idx_H{$label_str} = $nclasses;
    $ct_label_str_H{$label_str} = 0;
    $fa_label_str_H{$label_str} = $dir . "/" . $dir . "." . $nclasses . ".fa";
    @{$out_label_str_HA{$label_str}} = ();
  }
  $ct_label_str_H{$label_str}++;

  $accn2label_str_H{$accn} = $label_str;
  if($a == 0) { $ref_label_str = $label_str; }
  # and save the label_string so we can output the xlist1 below
  for(my $z = 0; $z < $nxlist1_this_accn; $z++) { 
    push(@xlist1_label_str_A, $label_str);
  }
} # end of first pass through all accessions 'for(my $a = 0; $a < $naccn; $a++)'

############################################################
# Set order of classes for output, by number of members, 
# except ref class is always first
############################################################
my %already_used_H = ();
my $oc = 1;
my $nused = 0;

# set ref class as first class
$already_used_H{$ref_label_str} = 1;
$label_str2idx_H{$ref_label_str} = $oc;
$idx2label_str_H{$oc} = $ref_label_str;
$nused++;
# printf("set $ref_label_str with $ct_label_str_H{$ref_label_str} members as class $oc\n");
$oc++;
while($nused < $nclasses) { 
  my $max_ct = 0;
  my $max_label_str = "";
  foreach $label_str (sort keys (%label_str2idx_H)) { 
    if((! exists $already_used_H{$label_str}) && ($ct_label_str_H{$label_str} > $max_ct)) { 
      $max_label_str = $label_str;
      $max_ct = $ct_label_str_H{$label_str};
    }
  }
  $label_str2idx_H{$max_label_str} = $oc;
  $idx2label_str_H{$oc} = $max_label_str;
  $already_used_H{$max_label_str} = 1;
  $nused++;
  # printf("set $max_label_str with $max_ct members as class $oc\n");
  $oc++;
}

#################################################
# Output xlist1, one line per x in any accession
# We constructed 'most' of this output in the first
# pass through accessions, but we didn't output 
# it because we didn't have class indices. Now
# we do so we can output them.
#################################################
printf OUTX1 ("%-5s  %-10s  %3s  %-*s  %4s  %-40s  %6s  %10s  %5s  %-40s  %10s  %10s\n",
              "#idx", "accession", "cls", $wlabel_str, "label", "#cds", "product", "length", "#ref-cds", "match", "ref-product", "ref-length", "lendiff");

my $nlines_x1 = 0;
for(my $l = 0; $l < scalar(@xlist1_output_A); $l++) { 
  my $xlist1_line = $xlist1_output_A[$l];
  my $accn        = $xlist1_accn_A[$l];
  my $label_str   = $xlist1_label_str_A[$l];
  my $class_idx   = $label_str2idx_H{$label_str};

  $nlines_x1++;
  printf OUTX1 ("%-5s  %-10s  %3d  %-*s  %s\n", $nlines_x1, $accn, $class_idx, $wlabel_str, $label_str, $xlist1_line); 
}
# output explanation of column headings and other information:
printf OUTX1 ("#\n");
printf OUTX1 ("# The above table lists all 'unmappable' CDS, one per line.\n");
printf OUTX1 ("# An unmappable CDS is a CDS that does not 'map' to any reference CDS in the reference accession ($ref_accn).\n");
printf OUTX1 ("# A CDS 'maps' to a reference CDS if it matches identically in either the CDS:product annotation\n");
#printf OUTX1 ("# An unmappable CDS is a CDS for which our rules for mapping to reference CDS fail. There is more information\n");
#printf OUTX1 ("# on these rules below in the section 'Explanation of how classes are defined and labels are determined'.\n");
printf OUTX1 ("#\n");
printf OUTX1 ("# A condensed list, with all unmappable CDS with identical CDS:product and length collapsed to a single line\n");
printf OUTX1 ("# has been output to $xlist2_outfile\n");
printf OUTX1 ("#\n");
printf OUTX1 ("# Explanation of column headings in this table:\n");
printf OUTX1 ("#\n");
printf OUTX1 ("# idx:           line index\n");
printf OUTX1 ("# accession:     accession the unmappable CDS occurs in\n");
printf OUTX1 ("# cls:           index of class this accession belongs to\n");
printf OUTX1 ("# label:         class 'label', defines which CDS map to reference, see below for more info\n");
printf OUTX1 ("# #cds:          CDS index in this accession that this line pertains to\n");
printf OUTX1 ("# product:       CDS:product annotation for these unmappable CDS\n");
printf OUTX1 ("# length:        length of the unmappable CDS\n");
printf OUTX1 ("# #ref-cds:      the index of the CDS in the reference that has the closest length to this unmappable CDS\n");
printf OUTX1 ("# match:         'yes' if '#cds' and '#ref-cds' values are equal, else 'no'\n");
printf OUTX1 ("# ref-product:   the CDS:product annotation of the CDS in the reference that has the closest length to this unmappable CDS\n");
printf OUTX1 ("# ref-length:    the length of the CDS in the reference that has the closest length to this unmappable CDS\n");
printf OUTX1 ("# lendiff:       the difference in length between this unmappable CDS and the nearest length reference CDS (abs('ref-length' - 'length')\n");
printf OUTX1 ("#\n");
# explain how classes are defined and label strings are constructed
#foreach my $line (@class_explanation_A) { 
#  print OUTX1 $line;
#}

#################################################
# Collapse list of unmappable CDS and output it
#################################################
my $n_x_product_keys = scalar(keys %x_product_accn_ct_H);
my $nprinted = 0;
printf OUTX2 ("%-5s  %-40s  %10s  %6s  %s\n", 
              "#idx", "CDS:product-annotation", "length", "count", "accessions");
while($nprinted < $n_x_product_keys) { 
  my $max_ct = 0;
  my $max_key = undef;
  foreach my $key (reverse sort keys %x_product_accn_ct_H) { 
    if($x_product_accn_ct_H{$key} > $max_ct) { 
      $max_ct = $x_product_accn_ct_H{$key};
      $max_key = $key;
    }
  }
  foreach my $len (sort {$x_product_len_accn_ct_HH{$max_key}{$a} <=> $x_product_len_accn_ct_HH{$max_key}{$b}} keys %{$x_product_len_accn_ct_HH{$max_key}}) {
    my $accn_str = "";
    my $naccn = scalar(@{$x_product_len_accn_HHA{$max_key}{$len}});
    for(my $z = 0; $z < $naccn; $z++) { 
      $accn_str .= $x_product_len_accn_HHA{$max_key}{$len}[$z];
      if($z < ($naccn-1)) { 
        $accn_str .= ",";
      }
    }
    $nlines_x2++;
    printf OUTX2 ("%-5d  %-40s  %10d  %6d  %s\n", $nlines_x2, $max_key, $len, $x_product_len_accn_ct_HH{$max_key}{$len}, $accn_str);
  }
  printf OUTX2 ("\n");
  $x_product_accn_ct_H{$max_key} = 0; # this won't be max again
  $nprinted++;
}
printf OUTX2 ("#\n");
printf OUTX2 ("# The above table lists all unique CDS:product and length combinations for all\n");
printf OUTX2 ("# $nlines_x1 'unmappable' CDS.\n");
printf OUTX2 ("# Each line reports a different CDS:product and length combination, and includes\n");
printf OUTX2 ("# all accessions that have that combination at the end of the line, each separated by a ','.\n");
printf OUTX2 ("# Identical CDS:product values occur on adjacent lines.\n");
printf OUTX2 ("# A blank line separates each different CDS:product value.\n");
printf OUTX2 ("#\n");
printf OUTX1 ("# An unmappable CDS is a CDS that does not 'map' to any reference CDS in the reference accession ($ref_accn).\n");
printf OUTX1 ("# A CDS 'maps' to a reference CDS if it matches identically in either the CDS:product annotation\n");
#printf OUTX2 ("# An unmappable CDS is a CDS for which our rules for mapping to reference CDS fail. There is more information\n");
#printf OUTX2 ("# on these rules below in the section 'Explanation of how classes are defined and labels are determined'.\n");
printf OUTX2 ("#\n");
printf OUTX2 ("# A full list, with each unmappable CDS on a separate line, which includes more information on each\n");
printf OUTX2 ("# has been output to $xlist1_outfile\n");
printf OUTX2 ("#\n");
printf OUTX2 ("# Explanation of column headings in this table:\n");
printf OUTX2 ("#\n");
printf OUTX2 ("# idx:                     line index in this file\n");
printf OUTX2 ("# CDS:product-annotation:  CDS:product annotation for this unmappable CDS\n");
printf OUTX2 ("# length:                  length of the unmappable CDS\n");
printf OUTX2 ("# count:                   number of accessions that include an unmappable CDS of this CDS:product annotation\n");
printf OUTX2 ("#                          and length\n");
printf OUTX2 ("# accessions:              list of all accessions that include an unmappable CDS of this CDS:product annotation\n");
printf OUTX2 ("#                          and length, each separated by a ','\n");
printf OUTX2 ("#\n");
# explain how classes are defined and label strings are constructed
#foreach my $line (@class_explanation_A) { 
#  print OUTX2 $line;
#}


###################################################################
# Now that we know the class indices for all classes, we can output
# class info for each accession. 
# This requires the second of two passes through all accessions
# In this loop we also store information needed to fetch the
# CDS for each accession.
###################################################################
for(my $a = 0; $a < scalar(@accn_A); $a++) { 
  my $accn = $accn_A[$a];

  # set defaults that will stay if we don't have any CDS information
  $ncds = 0; 
  $npos = 0;
  $nneg = 0;
  $nunc = 0;
  $nbth = 0; 
  $strand_str = "";
  @cds_len_A = ();
  @cds_coords_A = ();
  @cds_product_A = ();    # will remain empty unless $do_product is 1 (-product enabled at cmdline)
  @cds_protid_A = ();     # will remain empty unless $do_protid is 1 (-protid enabled at cmdline)
  @cds_codonstart_A = (); # will remain empty unless $do_codonstart is 1 (-codonstart enabled at cmdline)
  if((   $do_gene  && exists ($gene_tbl_HHA{$accn})) || 
     ((! $do_gene) && exists ($cds_tbl_HHA{$accn}))) { 
    if($do_gene) { 
      ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%gene_tbl_HHA, $accn);
      getLengthStatsAndCoordStrings(\%gene_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
      getQualifierValues(\%gene_tbl_HHA, $accn, "product", \@cds_product_A);
    }
    else { 
      ($ncds, $npos, $nneg, $nunc, $nbth, $strand_str) = getStrandStats(\%cds_tbl_HHA, $accn);
      getLengthStatsAndCoordStrings(\%cds_tbl_HHA, $accn, \@cds_len_A, \@cds_coords_A);
      getQualifierValues(\%cds_tbl_HHA, $accn, "product", \@cds_product_A);
    }
    if($do_protid) { 
      if($do_gene) { die "ERROR -protid and -gene are incompatible, currently"; }
      getQualifierValues(\%cds_tbl_HHA, $accn, "protein_id", \@cds_protid_A);
    }
    if($do_codonstart) { 
      if($do_gene) { die "ERROR -codonstartd and -gene are incompatible, currently"; }
      getQualifierValues(\%cds_tbl_HHA, $accn, "codon_start", \@cds_codonstart_A);
    }
  }
  my $label_str = $accn2label_str_H{$accn};
  $class_idx = $label_str2idx_H{$label_str};
  my $c = $class_idx - 1;
  my $ncds = scalar(@cds_len_A);
  $ncds_per_class_A[$c] = $ncds;
  if($ncds > $max_ncds) { $max_ncds = $ncds; }
  if(! exists $out_fetch_cds_AA[$c]) { @{$out_fetch_cds_AA[$c]} = (); }
  if(! exists $ct_fetch_cds_AA[$c])  { @{$ct_fetch_cds_AA[$c]} = (); }
  
  my $outline = sprintf("%-*s  %3d  %-*s  %-*s  %5d  %7d  ", $waccn, $accn, $class_idx, $wlabel_str, $label_str, $wstrand_str, $strand_str, $ncds, $totlen_H{$accn});
  $out_fetch_gnm_A[$c] .= sprintf("%s\n", $accn);
  $ct_fetch_gnm_A[$c]++;

  for(my $i = 0; $i < $ncds; $i++) { 
    my $desc = "";
    # determine DESCRIPTION: string, esl-fetch-cds.pl will parse this and add it as a description in the defline after fetching the sequence
    if($do_desc) { 
      $desc = "\tDESCRIPTION:";
      if($do_product) { 
        $desc .= "product:";
        if(scalar(@cds_product_A) > 0) { 
          if($i >= scalar(@cds_product_A)) { die "ERROR ran out of products too early for $accn\n"; }
          $desc .= $cds_product_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
        if($do_protid || $do_codonstart) { $desc .= " "; }
      }
      if($do_protid) { 
        $desc .= "protein_id:";
        if(scalar(@cds_protid_A) > 0) { 
          if($i >= scalar(@cds_protid_A)) { die "ERROR ran out of protein_ids too early for $accn\n"; }
          $desc .= $cds_protid_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
        if($do_codonstart) { $desc .= " "; }
      }
      if($do_codonstart) { 
        $desc .= "codon_start:";
        if(scalar(@cds_codonstart_A) > 0) { 
          if($i >= scalar(@cds_codonstart_A)) { die "ERROR ran out of codon_starts too early for $accn\n"; }
          $desc .= $cds_codonstart_A[$i];
        }
        else { 
          $desc .= "none-annotated";
        }
      }
    }
    $outline .= sprintf("  %5d", $cds_len_A[$i]);
    # create line of input for esl-fetch-cds.pl for fetching the genes of this genome
    if($do_shortnames) { 
      $out_fetch_cds_AA[$c][$i] .= "$accn\t$cds_coords_A[$i]";
    }
    else { 
      $out_fetch_cds_AA[$c][$i] .= sprintf("%s:%s%d:%s%d\t$cds_coords_A[$i]", $head_accn, "class", $class_idx, "gene", ($i+1));
    }
    $out_fetch_cds_AA[$c][$i] .= $desc . "\n"; # $desc will be "" unless $do_product is true
    $ct_fetch_cds_AA[$c][$i]++;
  }
  $outline .= "\n";
  push(@{$out_label_str_HA{$label_str}}, $outline);
} # end of second pass over all accessions 'for(my $a = 0; $a < scalar(@accn_A); $a++)'

#######################################################
# Output list of all accessions with label information.
# We do this class by class.
#######################################################
# print header line
printf OUTL ("#%-*s  %3s  %-*s  %-*s  %5s  %7s  ", $waccn-1, "accn", "cls", $wlabel_str, "label", $wstrand_str, "strand-string", "#cds", "tot-len");
for(my $c = 0; $c < $max_ncds; $c++) { 
  printf OUTL ("  %5s", sprintf("%s%d", ($do_gene) ? "gene" : "cds", ($c+1)));
}
printf OUTL ("\n");

# output stats
for(my $c = 1; $c <= $nclasses; $c++) { 
  my $label_str = $idx2label_str_H{$c};
  foreach my $outline (@{$out_label_str_HA{$label_str}}) { 
    print OUTL $outline;
  }
  print OUTL "\n";
  @{$out_label_str_HA{$label_str}} = (); # clear it for consise output
}

# output explanation of column headins

printf OUTL ("#\n");
printf OUTL ("# Explanation of column headings:\n");
printf OUTL ("#\n");
printf OUTL ("# accn:          accession\n");
printf OUTL ("# cls:           index of class this accession belongs to\n");
printf OUTL ("# label:         class 'label', defines which CDS map to reference, see below for more info\n");
printf OUTL ("# strand-string: string indicating strand of all CDS for this accession, in order, '+': positive, '-' negative\n");
printf OUTL ("# #cds:          number of CDS for this accession\n");
printf OUTL ("# tot-len:       total length in nucleotides for this accession\n");
printf OUTL ("# cds<n>:        length of CDS <n> in this accession\n");
printf OUTL ("#\n");
foreach my $line (@class_explanation_A) { 
  print OUTL $line;
}

################################
# Print summary output to stdout
################################
#
printf("\n");
printf("# Number-of-classes: $nclasses\n");
printf("# class  #accn  #genes  strand-string\n");
printf("# -----  -----  ------  -------------\n");
my $tot_ncds = 0;
my $tot_ct = 0;
for(my $c = 0; $c < $nclasses; $c++) { 
  my $label_str = $idx2label_str_H{($c+1)};
  printf("%7d  %5d  %6d  %s\n", ($c+1), $ct_label_str_H{$label_str}, $ncds_per_class_A[$c], $label_str);
  $tot_ncds += $ncds_per_class_A[$c];
  $tot_ct   += $ct_label_str_H{$label_str};
}
printf("# -----  -----  ------  -------------\n");
printf("%7s  %5d  %6d  %s\n", "total",     $tot_ct,             $tot_ncds,             "N/A");
printf("%7s  %5.1f  %6.1f  %s\n", "avg",   $tot_ct / $nclasses, $tot_ncds / $nclasses, "N/A");
printf("\n");

if($do_uc) { 
  open(UC, ">" . $uc_script);
}
if(! $do_noexp) { 
  foreach my $line (@class_explanation_A) { 
    print $line;
  }
}

close(OUTX1);
close(OUTX2);
close(OUTL);
print("#\n");
print("# Per-accession list of classes saved to $llist_outfile ($nclasses classes, $naccn accessions)\n"); 
print("# List of all CDS that are unmappable to reference CDS saved to $xlist1_outfile ($nlines_x1 lines)\n"); 
print("# Collapsed list of unmappable CDS saved to $xlist2_outfile ($nlines_x2 lines)\n");

#################################################################
# Output esl-fetch-cds input, and run esl-fetch-cds.pl for each 
# file to fetch all seqs. We do this class by class.
#################################################################
if($do_fetch) { 
  for(my $c = 0; $c < $nclasses; $c++) { 
    # fetch the full genomes
    my $out_fetch_gnm_file   = $out_root . ".c" . ($c+1) . ".fg.idfetch.in";
    my $tmp_out_fetch_gnm_fa = $out_root . ".c" . ($c+1) . ".fg.fa.tmp";
    my $out_fetch_gnm_fa     = $out_root . ".c" . ($c+1) . ".fg.fa";
    open(OUT, ">" . $out_fetch_gnm_file) || die "ERROR unable to open $out_fetch_gnm_file for writing";
    print OUT $out_fetch_gnm_A[$c];
    close OUT;
    sleep(0.1);
    printf("# Fetching %3d full genome sequences for class %2d ... ", $ct_fetch_gnm_A[$c], $c+1);
    my $cmd = "$idfetch -t 5 -c 1 -G $out_fetch_gnm_file > $tmp_out_fetch_gnm_fa";
    runCommand($cmd, 0);
    # now open up the file and change the sequence names manually
    open(IN, $tmp_out_fetch_gnm_fa)    || die "ERROR unable to open $tmp_out_fetch_gnm_fa for reading";
    open(OUT, ">" . $out_fetch_gnm_fa) || die "ERROR unable to open $out_fetch_gnm_fa for writing";
    while(my $line = <IN>) { 
      if($line =~ m/^>/) { 
        chomp $line;
        if($line =~ /^>.*\|(\S+)\|\S*\s+(.+)$/) { 
          print OUT (">$1 $2\n");
        }
        else { die "ERROR unable to parse defline $line in file $tmp_out_fetch_gnm_fa"; }
      }
      else { print OUT $line; }
    }
    close(IN);
    close(OUT);
    unlink $tmp_out_fetch_gnm_fa;
    printf("done. [$out_fetch_gnm_fa]\n");
    
    # fetch the cds'
    for(my $i = 0; $i < scalar(@{$out_fetch_cds_AA[$c]}); $i++) { 
      my $cg_substr = ".c" . ($c+1) . ".g" . ($i+1);
      my $out_fetch_cds_file = $out_root . $cg_substr . ".esl-fetch-cds.in";
      my $out_fetch_cds_fa   = $out_root . $cg_substr . ".fa";
      my $np_out_fetch_cds_fa    = stripPath($out_fetch_cds_fa);
      my $np_out_centroids_fa    = stripPath($out_root . $cg_substr . ".centroids.fa");
      my $np_out_uclust_sum_file = stripPath($out_root . $cg_substr . ".cluster_summary.txt");
      my $np_out_msa_root        = stripPath($out_root . $cg_substr . ".msa_cluster_");
      open(OUT, ">" . $out_fetch_cds_file) || die "ERROR unable to open $out_fetch_cds_file for writing";
      print OUT $out_fetch_cds_AA[$c][$i];
      close OUT;
      sleep(0.1);
      printf("# Fetching %3d CDS sequences for class %2d gene %2d ... ", $ct_fetch_cds_AA[$c][$i], $c+1, $i+1);
      my $cmd = "";
      if($do_shortnames) { 
        $cmd = "perl $esl_fetch_cds -onlyaccn $out_fetch_cds_file > $out_fetch_cds_fa";
        #printf("$cmd\n");
      }
      else { 
        $cmd = "perl $esl_fetch_cds -nocodon $out_fetch_cds_file > $out_fetch_cds_fa";
      }
      runCommand($cmd, 0);
      
      printf("done. [$out_fetch_cds_fa]\n");
      
      if($do_uc){ 
        print UC ("$usearch -cluster_fast $np_out_fetch_cds_fa -id $uc_id -centroids $np_out_centroids_fa -uc $np_out_uclust_sum_file -msaout $np_out_msa_root\n");
        $nuc++;
      }
    }
  }
} # end of 'if($do_fetch)'

if($do_uc) { 
  close(UC);
  printf("#\n# Shell script for running $nuc usearch commands saved to $uc_script.\n");
}
printf("#[ok]\n");

#############
# SUBROUTINES
#############
# Subroutine: runCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_exp = 2;

  my ($cmd, $be_verbose) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: parseLength()
# Synopsis:   Parses a length file and stores the lengths read
#             into %{$len_HR}.
# Args:       $lenfile: full path to a length file
#             $len_HR:  ref to hash of lengths, key is accession
#
# Returns:    void; fills %{$len_HR}
#
# Dies:       if problem parsing $lenfile

sub parseLength {
  my $sub_name = "parseLength()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($lenfile, $len_HR) = @_;

#HM448898.1	2751

  open(LEN, $lenfile) || die "ERROR unable to open $lenfile for reading";

  while(my $line = <LEN>) { 
    chomp $line;
    my ($accn, $length) = split(/\s+/, $line);
    if($length !~ m/^\d+$/) { die "ERROR couldn't parse length file line: $line\n"; } 

    stripVersion(\$accn);
    $len_HR->{$accn} = $length;
  }
  close(LEN);

  return;
}

# Subroutine: parseTable()
# Synopsis:   Parses a table file and stores the relevant info in it 
#             into $values_HAR.
# Args:       $tblfile:      full path to a table file
#             $values_HHAR:  ref to hash of hash of arrays
#
# Returns:    void; fills @{$values_HHAR}
#
# Dies:       if problem parsing $tblfile

sub parseTable {
  my $sub_name = "parseTable()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($tblfile, $values_HHAR) = @_;

##full-accession	accession	coords	strand	min-coord	gene
#gb|HM448898.1|	HM448898.1	129..476	+	129	AV2

  open(TBL, $tblfile) || die "ERROR unable to open $tblfile for reading";

  # get column header line:
  my $line_ctr = 0;
  my @colnames_A = ();
  my $line = <TBL>;
  my $ncols = undef;
  $line_ctr++;
  if(! defined $line) { die "ERROR did not read any lines from file $tblfile"; }
  chomp $line;
  if($line =~ s/^\#//) { 
    @colnames_A = split(/\t/, $line);
    $ncols = scalar(@colnames_A);
  }
  else { 
    die "ERROR first line of $tblfile did not start with \"#\"";
  }
  if($colnames_A[0] ne "full-accession") { die "ERROR first column name is not full-accession"; }
  if($colnames_A[1] ne "accession")      { die "ERROR second column name is not accession"; }
  if($colnames_A[2] ne "coords")         { die "ERROR third column name is not coords"; }

  # read remaining lines
  while($line = <TBL>) { 
    chomp $line;
    $line_ctr++;
    if($line =~ m/^\#/) { die "ERROR, line $line_ctr of $tblfile begins with \"#\""; }
    my @el_A = split(/\t/, $line);
    if(scalar(@el_A) != $ncols) { 
      die "ERROR, read wrong number of columns in line $line_ctr of file $tblfile";
    }
    my $prv_min_coord = 0;
    # get accession
    my $accn = $el_A[1]; 
    stripVersion(\$accn);
    if(! exists $values_HHAR->{$accn}) { 
      %{$values_HHAR->{$accn}} = (); 
    }

    for(my $i = 0; $i < $ncols; $i++) { 
      my $colname = $colnames_A[$i];
      my $value   = $el_A[$i];
      if($colname eq "min-coord") { 
        if($value < $prv_min_coord) { 
          die "ERROR, minimum coordinates out of order at line $line_ctr and previous line of file $tblfile"; 
        }
        $prv_min_coord = $value; 
        # printf("prv_min_coord: $prv_min_coord\n");
      }

      if(! exists $values_HHAR->{$accn}{$colname}) { 
        @{$values_HHAR->{$accn}{$colname}} = ();
      }
      push(@{$values_HHAR->{$accn}{$colname}}, $el_A[$i]);
      #printf("pushed $accn $colname $el_A[$i]\n");
    }
  }
  close(TBL);
  return;
}

# Subroutine: getStrandStats()
# Synopsis:   Retreive strand stats.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      1D key to print for
#
# Returns:    6 values:
#             $nfeatures:  number of features
#             $npos:       number of genes with all segments on positive strand
#             $nneg:       number of genes with all segmenst on negative strand
#             $nunc:       number of genes with all segments on unknown strand 
#             $nbth:       number of genes with that don't fit above 3 categories
#             $strand_str: strand string, summarizing strand of all genes, in order
#
sub getStrandStats {
  my $sub_name = "getStrandStats()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn) = @_;

  my $nfeatures; # number of genes in this genome
  my $npos = 0;  # number of genes on positive strand 
  my $nneg = 0;  # number of genes on negative strand 
  my $nbth = 0;  # number of genes with >= 1 segment on both strands (usually 0)
  my $nunc = 0;  # number of genes with >= 1 segments that are uncertain (usually 0)
  my $strand_str = "";

  if(! exists $tbl_HHAR->{$accn}{"strand"}) { die("ERROR didn't read strand information for accn: $accn\n"); }

  $nfeatures = scalar(@{$tbl_HHAR->{$accn}{"accession"}});
  if ($nfeatures > 0) { 
    for(my $i = 0; $i < $nfeatures; $i++) { 

      # sanity check
      my $accn2 = $tbl_HHAR->{$accn}{"accession"}[$i];
      stripVersion(\$accn2);
      if($accn ne $accn2) { die "ERROR accession mismatch in gene ftable file ($accn ne $accn2)"; }

      if   ($tbl_HHAR->{$accn}{"strand"}[$i] eq "+") { $npos++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "-") { $nneg++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "!") { $nbth++; }
      elsif($tbl_HHAR->{$accn}{"strand"}[$i] eq "?") { $nunc++; }
      else { die("ERROR unable to parse strand for feature %d for $accn\n", $i+1); }
      $strand_str .= $tbl_HHAR->{$accn}{"strand"}[$i];
    }
  }

  return ($nfeatures, $npos, $nneg, $nunc, $nbth, $strand_str);
}


# Subroutine: getLengthStatsAndCoordStrings()
# Synopsis:   Retreive length stats for an accession
#             the length of all annotated genes.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      accession we're interested in
#             $len_AR:    ref to array to fill with lengths of features in %{$tbl_HAR}
#             $coords_AR: ref to array to fill with coordinates for each gene
# Returns:    void; fills @{$len_AR} and @{$coords_AR}
#
sub getLengthStatsAndCoordStrings { 
  my $sub_name = "getLengthStatsAndCoordStrings()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $len_AR, $coords_AR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { die "ERROR in $sub_name, no data for accession: $accn"; }
  if(! exists $tbl_HHAR->{$accn}{"coords"}) { die "ERROR in $sub_name, no coords data for accession: $accn"; }

  my $ngenes = scalar(@{$tbl_HHAR->{$accn}{"coords"}});

  if ($ngenes > 0) { 
    for(my $i = 0; $i < $ngenes; $i++) { 
      push(@{$len_AR},    lengthFromCoords($tbl_HHAR->{$accn}{"coords"}[$i]));
      push(@{$coords_AR}, addAccnToCoords($tbl_HHAR->{$accn}{"coords"}[$i], $accn));
    }
  }

  return;
}

# Subroutine: getQualifierValues()
# Synopsis:   Retreive values for the qualifier $qualifier in the given %tbl_HHAR
#             and return the values in $values_AR.
#             the length of all annotated genes.
# Args:       $tbl_HHAR:  ref to hash of hash of arrays
#             $accn:      accession we're interested in
#             $qualifier: qualifier we're interested in (e.g. 'Product')
#             $values_AR: ref to array to fill with values of $qualifier
# Returns:    void; fills @{$values_AR}
#
sub getQualifierValues {
  my $sub_name = "getQualifierValues()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($tbl_HHAR, $accn, $qualifier, $values_AR) = @_;

  if(! exists $tbl_HHAR->{$accn}) { die "ERROR in $sub_name, no data for accession: $accn"; }

  if(! exists $tbl_HHAR->{$accn}{$qualifier}) { return; } # no annotation for $qualifier, do not update arrays

  my $nvalues = scalar(@{$tbl_HHAR->{$accn}{$qualifier}});

  if ($nvalues > 0) { 
    for(my $i = 0; $i < $nvalues; $i++) { 
      push(@{$values_AR},  $tbl_HHAR->{$accn}{$qualifier}[$i]);
    }
  }

  return;
}


# Subroutine: lengthFromCoords()
# Synopsis:   Determine the length of a region give its coords in NCBI format.
#
# Args:       $coords:  the coords string
#
# Returns:    length in nucleotides implied by $coords  
#
sub lengthFromCoords { 
  my $sub_name = "lengthFromCoords()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords) = @_;

  my $orig_coords = $coords;
  # Examples:
  # complement(2173412..2176090)
  # complement(join(226623..226774, 226854..229725))

  # remove 'complement('  ')'
  $coords =~ s/^complement\(//;
  $coords =~ s/\)$//;

  # remove 'join('  ')'
  $coords =~ s/^join\(//;
  $coords =~ s/\)$//;

  my @el_A = split(/\s*\,\s*/, $coords);

  my $length = 0;
  foreach my $el (@el_A) { 
    # rare case: remove 'complement(' ')' that still exists:
    $el =~ s/^complement\(//;
    $el =~ s/\)$//;
    $el =~ s/\<//; # remove '<'
    $el =~ s/\>//; # remove '>'
    if($el =~ m/^(\d+)\.\.(\d+)$/) { 
      my ($start, $stop) = ($1, $2);
      $length += abs($start - $stop) + 1;
    }
    else { 
      die "ERROR unable to parse $orig_coords in $sub_name"; 
    }
  }

  # printf("in lengthFromCoords(): orig_coords: $orig_coords returning length: $length\n");
  return $length;
}

# Subroutine: addAccnToCoords()
# Synopsis:   Add accession Determine the length of a region give its coords in NCBI format.
#
# Args:       $coords:  the coords string
#             $accn:    accession to add
# Returns:    The accession to add to the coords string.
#
sub addAccnToCoords { 
  my $sub_name = "addAccnToCoords()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
 
  my ($coords, $accn) = @_;

  my $ret_coords = $coords;
  # deal with simple case of \d+..\d+
  if($ret_coords =~ /^\<?\d+\.\.\>?\d+/) { 
    $ret_coords = $accn . ":" . $ret_coords;
  }
  # replace 'complement(\d' with 'complement($accn:\d+'
  while($ret_coords =~ /complement\(\<?\d+/) { 
    $ret_coords =~ s/complement\((\<?\d+)/complement\($accn:$1/;
  }
  # replace 'join(\d' with 'join($accn:\d+'
  while($ret_coords =~ /join\(\<?\d+/) { 
    $ret_coords =~ s/join\((\<?\d+)/join\($accn:$1/;
  }
  # replace ',\d+' with ',$accn:\d+'
  while($ret_coords =~ /\,\s*\<?\d+/) { 
    $ret_coords =~ s/\,\s*(\<?\d+)/\,$accn:$1/;
  }

  #print("addAccnToCoords(), input $coords, returning $ret_coords\n");
  return $ret_coords;
}

# Subroutine: stripVersion()
# Purpose:    Given a ref to an accession.version string, remove the version.
# Args:       $accver_R: ref to accession version string
# Returns:    Nothing, $$accver_R has version removed
sub stripVersion {
  my $sub_name  = "stripVersion()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}

# Subroutine: stripPath()
# Purpose:    Given a file path, remove the all directories and leave only the file name.
# Args:       $filename: full path to file
# Returns:    only the file name, without any directory structure
sub stripPath {
  my $sub_name  = "stripPath()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($filename) = (@_);

  $filename =~ s/^.+\///;

  return $filename;
}
