#!/usr/bin/perl

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - H E A D E R - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# MANUAL FOR KaKs_ratios.pl

=pod

=head1 NAME

$0 -- Comprehensive pipeline : Inputs Alignment file (Bam format), Reference file used and GFF file and generates 
a folder of all the different genes in separate files with the KaKs ratios.
: Leverages the KaKs_Calculator command line tool.
  
=head1 SYNOPSIS

pipGECOinserttranscriptome.pl [--help] [--manual]

=head1 DESCRIPTION

Accepts all folders from frnakenstein output.
 
=head1 OPTIONS

=over 3

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   DBI
   DBD::mysql
   Getopt::Long
   Pod::Usage

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2015 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage
=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR 
use strict;
use File::Basename;
use DBI;
use Getopt::Long;
use Time::localtime;
use Pod::Usage;
use Time::Piece;
use File::stat;
use DateTime;


# #CREATING LOG FILES
my $std_out = '/home/modupeore17/.LOG/TransDBmypip-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupeore17/.LOG/TransDBmypip-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';
my $jobid = "GECOpip-".`date +%m-%d-%y_%T`;
my $progressnote = "/home/modupeore17/.LOG/progressnote".`date +%m-%d-%y_%T`; chomp $progressnote; $progressnote = $progressnote.'.txt'; 

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";
 
#ARGUMENTS
my($help,$manual,$in1);
GetOptions (	
                                "g|gff=s"       =>      \$gff,
                                "r|ref=s"       =>      \$ref,
                                "b|bam=s"       =>      \$bam,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );


%POSstart = '';%POSend=''; %CHR='';%ORN='';%GENE='';
#converts gff file to a tab-delimited file for genes
$infolder = $ARGV[0];
open(IN,"<",$gff) or die "$gff can't open";
while (<IN>){
  chomp;
  @all = split("\t", $_);
  @newall = split("\;", $all[8]);
  if($all[2] =~ /gene/){
    foreach my $abc (@newall){
      if ($abc =~ /^gene=.*/){
        @bcd = split("\=",$abc,2);
          if (exists $GENE{$bcd[1]}){
            $GENE{$bcd[1]}=$GENE{$bcd[1]}++;
          }
          else {
            $GENE{$bcd[1]}=1;
          }
          $POSstart{$bcd[1]}{$GENE{$bcd[1]}} = $all[3];
          $POSend{$bcd[1]}{$GENE{$bcd[1]}} = $all[4];
          $CHR{$bcd[1]}{$GENE{$bcd[1]}} = $all[0];
          if ($all[6]=~/^\+/) {
            $ORN{$bcd[1]}{$GENE{$bcd[1]}} = "FOR";
          }
          else {
            $ORN{$bcd[1]}{$GENE{$bcd[1]}} = "REV";
          }
      }
    }
  }
}
close(IN);



open(IN,"<",$ARGV[0]);
open(OUT,">", $ARGV[2]);
@filecontent = <IN>;
close (IN);
@position=split("\-",((split("\:",$ARGV[1]))[1]),2);
$header = substr($filecontent[0],1,-1);
$seq = substr($filecontent[1],$position[0],-1);
print length($seq)."\n\n";

	

