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

KaKs_ratios.pl -b <xxx.bam> -r <ref.fa> -g <gff.gff> -o <OUTPUTDIR> -n name [--help] [--manual]

=head1 DESCRIPTION

Accepts bam file, reference file and gff file to calculate the KaKs_ratios implemented by the KaKs_calculator.
 
=head1 OPTIONS

=over 3

=item B<-b, --bam>=FILE

Mapping/Alignment file.  (Required)

=item B<-r, -ref, --reference>=FILE

Reference fasta file used to create the bam (i.e. mapping file).  (Required)

=item B<-g, --gff>=FILE

GFF file to determine the different genes (version GFF3).  (Required)

=item B<-o, -out, --output>=FOLDER

output directory where all the libraries will be stored.  (Required)

=item B<-n, --name>=NAME

Name of the file for output.  (Optional)

=item B<-h, --help>

Displays the usage message.  (Optional) 

=item B<-man, --manual>

Displays full manual.  (Optional) 

=back

=head1 DEPENDENCIES

Requires the following Perl libraries (all standard in most Perl installs).
   Getopt::Long
   Pod::Usage
   File::Basename;
   Time::localtime;
   Time::Piece;
   File::stat;
   DateTime;
   Data::Dumper;
   Parallel::ForkManager;
   Scalar::Util::Numeric qw(isint);

=head1 AUTHOR

Written by Modupe Adetunji, 
Center for Bioinformatics and Computational Biology Core Facility, University of Delaware.

=head1 REPORTING BUGS

Report bugs to amodupe@udel.edu

=head1 COPYRIGHT

Copyright 2016 MOA.  
License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>.  
This is free software: you are free to change and redistribute it.  
There is NO WARRANTY, to the extent permitted by law.  

Please acknowledge author and affiliation in published work arising from this script's usage

=cut

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - U S E R  V A R I A B L E S- - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# CODE FOR KaKs_ratios.pl
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use File::stat;
use DateTime;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use Parallel::ForkManager;
use Scalar::Util::Numeric qw(isint);

# #CREATING LOG FILES
my $std_out = '/home/modupe/.LOG/KaKs-'.`date +%m-%d-%y_%T`; chomp $std_out; $std_out = $std_out.'.log';
my $std_err = '/home/modupe/.LOG/KaKs-'.`date +%m-%d-%y_%T`; chomp $std_err; $std_err = $std_err.'.err';
my $jobid = "KaKs-".`date +%m-%d-%y_%T`;

open(STDOUT, '>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>', "$std_err") or die "Error file doesn't exist";
 
# #ARGUMENTS
my($help,$manual,$gff,$ref,$bam, $out, $name);

GetOptions (	
                                "g|gff=s"       =>      \$gff,
                                "r|ref|reference=s"       =>      \$ref,
                                "b|bam=s"       =>      \$bam,
                                "o|out|output=s"       =>      \$out,
                                "n|name=s"       =>      \$name,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# #VALIDATE ARGS
pod2usage( -verbose  => 2)  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument(s) are not found.\n", -exitval => 2, -verbose => 1)  if (! $gff || ! $bam || ! $ref || ! $out );

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - G L O B A L  V A R I A B L E S- - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#PATHS
my $KAKS = "~/.software/KaKs_Calculator2.0/src/KaKs_Calculator";

#INITIALIZE 
my %POSstart = ''; my %POSend=''; my %CHR=''; my %ORN=''; my %GENE=''; my %SEQ=''; my %FILEPATH = ''; my %HTALL = '';
my $prefix; my $sumofgenes=0;

#Making OUTPUT DIRECTORY[IES].
`mkdir $out`;
if ($name) { $prefix = $name; }
else { $prefix = fileparse($bam, qr/\.[^.]*(\.bam)?$/); }
my $newpath = $out."/".$prefix;
`mkdir $newpath`; chdir $newpath;
#STARTING
# 1. HTSEQ TO GET GENES ACTIVE
`htseq-count -t gene -i gene -s yes -f bam $bam $gff 1>htseq.log 2>htseq.err`;
#HTSEQ();
# 2. INITIALIZING GFF FILE
GFF_FILE();

# 3. PARSE TO REFERENCE  SCRIPT
open(REFPERL,">refperl.pl");
print REFPERL "#!/usr/bin/perl\nuse Scalar::Util::Numeric qw(isint);";
my $refperlcontent = <<"ENDREFPERL";
open(STDOUT, '>>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>>', "$std_err") or die "Error file doesn't exist";
ENDREFPERL
print REFPERL $refperlcontent."\n";
print REFPERL "\$sum=$sumofgenes;\n";
print REFPERL "chdir \"$newpath\";\n";
print REFPERL Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
print REFPERL Data::Dumper->Dump( [ \%POSstart ], [ qw(*POSstart) ] );
print REFPERL Data::Dumper->Dump( [ \%POSend ], [ qw(*POSend) ] );
print REFPERL Data::Dumper->Dump( [ \%ORN ], [ qw(*ORN) ] );
print REFPERL Data::Dumper->Dump( [ \%CHR ], [ qw(*CHR) ] );
print REFPERL "open(REF,\"<$ref\");\n";
$refperlcontent = <<'ENDREFPERL';
my %SEQ='';
$/= ">";
while (<REF>){
  my @pieces = split /\n/;
  my $header = $pieces[0];
  my $seq = ''; my $qua = '';
  foreach my $num (1.. $#pieces){
    $seq .= $pieces[$num];
  }
  $SEQ{$header}=$seq;
}
$/="\n";
close (REF);

foreach my $genename (sort keys %GENE ){
  foreach my $number (1..$GENE{$genename}){
    open(OUTREF,">", $genename."-".$number."_ref.nuc");
    my $count = $POSend{$genename}{$number}-$POSstart{$genename}{$number}+1;
    my $refseq = substr($SEQ{$CHR{$genename}{$number}},$POSstart{$genename}{$number}-1,$count);
    if ($ORN{$genename}{$number} =~ /REV/) {
      $refseq=reverse($refseq);
      $refseq =~ tr/ACGTacgt/TGCAtgca/;
    }
    #checking divisibility
    my $codons = $count/3;
    unless (isint $codons) {
	my $tocodons =  int($codons) + 1;
	my $newcodons = $tocodons*3;
        my $added = $newcodons - $count;
	for (my $in = 0; $in < $added; $in++){
  	    $refseq .= "-";
	}
	print OUTREF $refseq."\n";
    }
    else {
	print OUTREF $refseq."\n";
    }
    close OUTREF;
  }
}
ENDREFPERL
print REFPERL $refperlcontent."\n";
close REFPERL;

# 4. PARSE TO BAM FILE & SAMPLES
open(BAMPERL,">bamperl.pl");
my $bamperlcontent = <<"ENDBAMPERL";
#!/usr/bin/perl
use Parallel::ForkManager;
use Scalar::Util::Numeric qw(isint);
open(STDOUT, '>>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>>', "$std_err") or die "Error file doesn't exist";
my \$std_out = "$std_out";
my \$std_err = "$std_err";
my \$samsort = "samtools sort $bam -o aln.sorted.bam";
my \$samindex = "samtools index aln.sorted.bam";
#`\$samsort`;
#`\$samindex`;
ENDBAMPERL
print BAMPERL $bamperlcontent."\n";
print BAMPERL "\$sum=$sumofgenes;\n";
print BAMPERL "\$sumofgenes=$sumofgenes;\n";
print BAMPERL "\$path=\"$newpath\";\n";
print BAMPERL "\$ref=\"$ref\";\n";
print BAMPERL "chdir \"$newpath\";\n";
print BAMPERL Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
print BAMPERL Data::Dumper->Dump( [ \%POSstart ], [ qw(*POSstart) ] );
print BAMPERL Data::Dumper->Dump( [ \%POSend ], [ qw(*POSend) ] );
print BAMPERL Data::Dumper->Dump( [ \%ORN ], [ qw(*ORN) ] );
print BAMPERL Data::Dumper->Dump( [ \%CHR ], [ qw(*CHR) ] );
$bamperlcontent = <<'ENDBAMPERL';
my $ident = 0;
my $addi = 0;
my $tag = 0;
OVER: {
while ($addi<$sum) {
  $ident++;
  open(SAMPERL,">samperl-$ident.pl");
  print SAMPERL "#!/usr/bin/perl\n";
  print SAMPERL "use File::Basename;\n";
  print SAMPERL "use Scalar::Util::Numeric qw(isint);\n";
  print SAMPERL "use fastq;\n";
  print SAMPERL "chdir \"$path\";\n";
  print SAMPERL "\$sumofgenes=$sumofgenes;\n";
  print SAMPERL "open(STDOUT, '>>', \"$std_out\") or die \"Log file doesn't exist\";\nopen(STDERR, '>>', \"$std_err\") or die \"Error file doesn't exist\";\n";

  $tag=0;

  foreach my $genename (sort keys %GENE ){
    redo OVER if $tag==200;
    ++$tag;
    foreach my $number (1..$GENE{$genename}){
      $addi++;
      my $extractsam ="~/.software/samtools-0.1.19/samtools mpileup -uf $ref aln.sorted.bam -r $CHR{$genename}{$number}:$POSstart{$genename}{$number}-$POSend{$genename}{$number} | bcftools view -cg - | /home/modupe/.software/samtools-0.1.19/bcftools/vcfutils.pl vcf2fq > $genename-$number.pep";
my $samperlcontent = <<"ENDSAMPERL";
`$extractsam`;
FASTQTOFASTA("$genename-$number.pep", $CHR{$genename}{$number});
open(FASTA, "<$genename-$number.fna");
\@filecontent = <FASTA>; close (FASTA);
\$count = $POSend{$genename}{$number}-$POSstart{$genename}{$number}+1;
\$sampleseq = substr(\$filecontent[1],$POSstart{$genename}{$number}-1,-1);
#this is to add dashes to the end of sequences { works for some but not all }
\$add = \$count - length \$sampleseq;
for (my \$in = 0; \$in < \$add; \$in++){
  \$sampleseq .= "-";
}
ENDSAMPERL
print SAMPERL $samperlcontent."\n";

      if ($ORN{$genename}{$number} =~ /REV/) {
        print SAMPERL "\$sampleseq=reverse(\$sampleseq);\n";
        print SAMPERL "\$sampleseq =~ tr/ACGTacgt/TGCAtgca/;\n";
      }
      $outsamname=$genename."-".$number."_sam.nuc";
$samperlcontent = <<"ENDSAMPERL";
open(OUTSAM,">$outsamname");
my \$codons = \$count/3;
unless (isint \$codons) {
    my \$tocodons =  int(\$codons) + 1;
    my \$newcodons = \$tocodons*3;
    my \$added = \$newcodons - \$count;
    for (my \$in = 0; \$in < \$added; \$in++){
	\$sampleseq .= "-";
    }
    print OUTSAM \$sampleseq."\\n";
}
else {
    print OUTSAM \$sampleseq."\\n";
}
close OUTSAM;
ENDSAMPERL
print SAMPERL $samperlcontent."\n";
    }
    delete $GENE{$genename};
  }
}
}
$samperlcontent = <<"ENDSAMPERL";

ENDSAMPERL
print SAMPERL $samperlcontent."\n";
close SAMPERL;
my $max_procs = 15;
my $pm =  new Parallel::ForkManager($max_procs);
for $inc (1..$ident){
  my $pid = $pm->start and next;
  print "samperl-$inc working\n";
  print system ("perl samperl-$inc.pl");
  $pm->finish; # pass an exit code to finish
}
ENDBAMPERL
print BAMPERL $bamperlcontent."\n";
close BAMPERL;

# 5. FOR COMBINE STAGE && KAKS
open(COMBINE,">combine.pl");
my $combinecontent = <<"ENDCOMBINE";
#!/usr/bin/perl
use Parallel::ForkManager;
open(STDOUT, '>>', "$std_out") or die "Log file doesn't exist";
open(STDERR, '>>', "$std_err") or die "Error file doesn't exist";
my \$std_out = "$std_out";
my \$std_err = "$std_err";
ENDCOMBINE
print COMBINE $combinecontent."\n";
print COMBINE "\$sum=$sumofgenes;\n";
print COMBINE "\$path=\"$newpath\";\n";
print COMBINE "\$KAKS=\"$KAKS\";\n";
print COMBINE "chdir \"$newpath\";\n";
print COMBINE Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
$combinecontent = <<'ENDCOMBINE';
my $ident = 0;
my $addi = 0;
my $tag = 0;
OVER: {
while ($addi<$sum) {
  $ident++;
  open(KAKS,">kaks-$ident.pl");
  print KAKS "#!/usr/bin/perl\n";
  print KAKS "use File::Basename;\n";
  print KAKS "chdir \"$path\";\n";
  print KAKS "\$sum=$sum;\n";
  print KAKS "open(STDOUT, '>>', \"$std_out\") or die \"Log file doesn't exist\";\nopen(STDERR, '>>', \"$std_err\") or die \"Error file doesn't exist\";\n"; 
  $tag=0;
  foreach my $genename (sort keys %GENE ){
    redo OVER if $tag==200;
    ++$tag;
    foreach my $number (1..$GENE{$genename}){
      $addi++;
      my $fileloc = $genename."-".$number;
      my $naming = "echo \"$fileloc\" > $fileloc.name";
      my $concat = "cat $fileloc.name ".$fileloc."_ref.nuc ".$fileloc."_sam.nuc > $fileloc.axt";

      my $runkaks = "$KAKS -i $fileloc.axt -o $fileloc.kaks -m MLWL -m GY -m YN -m MYN -m MA";
my $kaksperlcontent = <<"ENDKAKSPERL";
`$naming`;
# $concat;
`$concat`;
`$runkaks`;
ENDKAKSPERL
print KAKS $kaksperlcontent."\n";
    }
    delete $GENE{$genename};
  }
}
}
$kaksperlcontent = <<'ENDKAKSPERL';
AGAIN: {
  my $countkaks = `ls -l *.kaks | wc -l`; chomp $countkaks;
  while ($countkaks <= $sum){
    redo AGAIN if $countkaks < $sum;
    if($countkaks == $sum) {
      `perl clean.pl`;
    }
  }
}
ENDKAKSPERL
print KAKS $kaksperlcontent."\n";
close KAKS;

my $max_procs = 15;
my $pm =  new Parallel::ForkManager($max_procs);
for $inc (1..$ident){
  my $pid = $pm->start and next;
  print "kaks-$inc working\n";
  print system ("perl kaks-$inc.pl");
  $pm->finish; # pass an exit code to finish
  print "All KAKS scripts queued done\n";
}
ENDCOMBINE
print COMBINE $combinecontent."\n";
close COMBINE;



##FOR FASTQTOFASTA SUBROUTINE
my $fastqcontent = <<'ENDFQPERL';
package fastq;
use strict;
use warnings;
use File::Basename;
use base 'Exporter';

our @EXPORT = qw/ FASTQTOFASTA /;

sub FASTQTOFASTA {
  my $fastqout = fileparse($_[0], qr/\.[^.]*(\.gz)?$/);
  open(OUTF, "> $fastqout.fna") or die $!;
  $/ = "\@".$_[1];
  open (FASTQ, $_[0]) or die "$_[0] can't open";
  my @fastqfile = <FASTQ>;
  shift(@fastqfile);
  foreach my $entry (@fastqfile){
    my @pieces = split(/\n/, $entry);
    my $header = ">".$_[1].$pieces[0];
    my $seq = ''; my $qua = '';
    my $check = 0;
    foreach my $num (1.. $#pieces){
      if ($pieces[$num] =~ /^\+$/) {
        $check = 1; next;
      }
      if ($check == 0) {
        $seq .= $pieces[$num];
      }
      elsif ($check == 1) {
        $qua .= $pieces[$num];
      }
    }
    print OUTF "$header\n$seq\n";
  }
  close FASTQ; close OUTF;
  $/ = "\n";
}
1;
ENDFQPERL
open(FQ,">fastq.pm");
print FQ $fastqcontent."\n";
close FQ;

# 6. CLEANUP
open(CLEAN,">cleanup.pl");
print CLEAN "#!/usr/bin/perl\n";
print CLEAN "chdir \"$newpath\";\n";
my $cleancontent = <<'ENDCLEAN';
#All the nucleotide files
`mkdir nucleotides`;
`mv *nuc nucleotides`;
#All the peptide files
`mkdir fastq`;
`mv *pep fastq`;
#All fasta files
`mkdir fasta`;
`mv *fna fasta`;
#All AXT files
`mkdir axtfiles`;
`mv *axt axtfiles`;
#All KaKs files
`mkdir KAKS`;
`mv *kaks KAKS`;
#All Name files
`mkdir NAME`;
`mv *name NAME`;
#All Scripts
`mkdir SCRIPTS`;
`mv *pl *pm SCRIPTS`;
print "All sorted into respective folders\n";

ENDCLEAN
print CLEAN $cleancontent."\n";
close CLEAN;
close STDOUT; close STDERR;

#7. JOB EXECUTION
print "Job execution: refperl.pl\n";
`perl refperl.pl`;
print "Job completed: refperl.pl\nJob execution: bamperl.pl\n";
`perl bamperl.pl`;

my $moveon = 0;
my $verdict = `ps -u | grep samperl`;
while ($verdict != /samperl-\d+\.pl/){
	$verdict = `ps -u | grep samperl`;
	$moveon = 1;
}

if ($moveon == 1){
	$moveon = 2;
	print "Job completed: bamperl.pl\nJob execution: combine.pl\n";
	`perl combine.pl`;
}

$verdict = `ps -u | grep kaks`;
while ($verdict != /kaks-\d+\.pl/) {
	$verdict = `ps -u | grep kaks`;
	$moveon = 3;
}

if ($moveon == 3){
	print "Job completed: combine.pl\nJob execution: cleanup.pl\n";
	`perl cleanup.pl`;
}

print "Job completed: cleanup.pl\n";
print "Finished!!!!\n";


##SUBROUTINES
#converts gff file to a tab-delimited file for genes
sub GFF_FILE {
  open(GFF,"<",$gff) or die "$gff can't open";
  while (<GFF>){
    chomp;
    my @all = split("\t", $_);
    my @newall = split("\;", $all[8]);
    if($all[2] =~ /gene/){
      foreach my $abc (@newall){
        if ($abc =~ /^Name=.*/){
          my @bcd = split("\=",$abc,2);
					$bcd[1] =~ s/ /\_/g;
#          if (exists $HTALL{$bcd[1]}){
            $sumofgenes++;
            if (exists $GENE{$bcd[1]}){
              $GENE{$bcd[1]}=$GENE{$bcd[1]}+1;
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
#          }
        }
      }
    }
  }
  close (GFF);
}

sub HTSEQ {
  open(HTSEQ, "<htseq.log") or die "can't open file\n";
  while (<HTSEQ>){
    chomp;
    my @all = split("\t",$_,2);
    if ($all[1] > 0 && $all[0] =~ /^[0-9a-zA-Z].*$/){
      $HTALL{$all[0]} = $all[1];
    }
  }
}
