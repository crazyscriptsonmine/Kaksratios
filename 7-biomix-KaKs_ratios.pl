#!/usr/bin/perl

# CODE FOR KaKs_ratios.pl
use strict;
use File::Basename;
use Getopt::Long;
use Time::localtime;
use File::stat;
use DateTime;
use Pod::Usage;
use Data::Dumper qw(Dumper);
use threads;
use Thread::Queue;
use Scalar::Util::Numeric qw(isint);

# #ARGUMENTS
my($help,$manual,$gff,$ref,$bam, $out,$filepath, $name);

GetOptions (	
                                "g|gff=s"       =>      \$gff,
                                "r|ref|reference=s"       =>      \$ref,
                                "b|bam=s"       =>      \$bam,
                                "p|path=s"	=>	\$filepath,
				"o|out|output=s"       =>      \$out,
                                "n|name=s"       =>      \$name,
                                "h|help"        =>      \$help,
                                "man|manual"	=>      \$manual );

# #VALIDATE ARGS
pod2usage( -verbose  => 2)  if ($manual);
pod2usage( -verbose => 1 )  if ($help);
pod2usage( -msg  => "ERROR!  Required argument(s) are not found; -g|gff , -r|ref, -b|bam / -p|path, -o|out.\n$0\n", -exitval => 2, -verbose => 1)  if (! $gff && (! $filepath || !$bam) && ! $ref && ! $out );

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
elsif ($bam) { $prefix = fileparse($bam, qr/\.[^.]*(\.bam)?$/); }
else { die "ERROR: option n|name needs to be specified if there isn't a bam file;" }
if ($filepath) { $bam = "$filepath/aln.sorted.bam"; }
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
print REFPERL "#!/usr/bin/perl\nuse Scalar::Util::Numeric qw(isint);\nuse threads;\nuse Thread::Queue;\n";
my $refperlcontent = <<"ENDREFPERL";
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

#- - - - - - - REFERNCE GENOME GENES SEQUENCE - - - - - - - - - - 
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
my @genez = map {$_; } sort keys %GENE ;

#push variables for threads.
push @VAR, [splice @genez, 0, 200] while @genez;
$queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..15; #execute 15 threads
$builder->join; #join threads
foreach (@threads){$_->join;}

sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
  }
  foreach(1..15) { $queue-> enqueue(undef); }
}

sub processor { my $query; while ($query = $queue->dequeue()){ parseinput(@$query); } }

sub parseinput {
  foreach my $genename (@_) {
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
			} else {
				print OUTREF $refseq."\n";
			}
			close OUTREF;
		}
	}
}

ENDREFPERL
print REFPERL $refperlcontent."\n";
close REFPERL;

## - - - - - - - END REFERNCE GENOME GENES SEQUENCE - - - - - - - - - - 


# 4. PARSE TO BAM FILE & SAMPLES
#- - - - - - - BAM FILE GENES SEQUENCE - - - - - - - - - - 
open(BAMPERL,">bamperl.pl");
my $bamperlcontent = <<"ENDBAMPERL";
#!/usr/bin/perl
use threads;
use Thread::Queue;
use Scalar::Util::Numeric qw(isint);
chdir \"$newpath\";

ENDBAMPERL
print BAMPERL $bamperlcontent."\n";
print BAMPERL "\$sum=$sumofgenes;\n";
print BAMPERL "\$sumofgenes=$sumofgenes;\n";
print BAMPERL "\$path=\"$newpath\";\n";
print BAMPERL "\$ref=\"$ref\";\n";
print BAMPERL "\$bam=\"$bam\";\n";
print BAMPERL "chdir \"$newpath\";\n";
print BAMPERL Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
print BAMPERL Data::Dumper->Dump( [ \%POSstart ], [ qw(*POSstart) ] );
print BAMPERL Data::Dumper->Dump( [ \%POSend ], [ qw(*POSend) ] );
print BAMPERL Data::Dumper->Dump( [ \%ORN ], [ qw(*ORN) ] );
print BAMPERL Data::Dumper->Dump( [ \%CHR ], [ qw(*CHR) ] );
$bamperlcontent = <<'ENDBAMPERL';
my ($ident, $addi, $tag) = (0,0,0);
my (@identi, @VAR);

OVER: {
while ($addi<$sum) {
  $ident++;
	push(@identi, "samperl-$ident");
  open(SAMPERL,">samperl-$ident.pl");
  print SAMPERL "#!/usr/bin/perl\n";
  print SAMPERL "use File::Basename;\n";
  print SAMPERL "use Scalar::Util::Numeric qw(isint);\n";
  print SAMPERL "use fastq;\n";
  print SAMPERL "chdir \"$path\";\n";
  print SAMPERL "\$sumofgenes=$sumofgenes;\n";

  $tag=0;

  foreach my $genename (sort keys %GENE ){
    redo OVER if $tag==200;
    ++$tag;
    foreach my $number (1..$GENE{$genename}){
      $addi++;
      my $extractsam ="samtools mpileup -uf $ref $bam -r $CHR{$genename}{$number}:$POSstart{$genename}{$number}-$POSend{$genename}{$number} | bcftools view -cg - | vcfutils.pl vcf2fq > $genename-$number.pep";
			open (CHRS, ">$genename-$number.chrs"); print CHRS "$genename-$number\t$CHR{$genename}{$number}:$POSstart{$genename}{$number}-$POSend{$genename}{$number}\n"; close (CHRS);
my $samperlcontent = <<"ENDSAMPERL";
print "$extractsam";
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

#push variables for threads.
push @VAR, [splice @identi, 0, 1] while @identi;
$queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..15; #execute 15 threads
$builder->join; #join threads
foreach (@threads){$_->join;}

sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
  }
  foreach(1..15) { $queue-> enqueue(undef); }
}

sub processor { my $query; while ($query = $queue->dequeue()){ parseinput(@$query); } }

sub parseinput {
  foreach my $inc (@_) {
		open (TEMPIT, ">temp-$inc.sh");
		print TEMPIT '#!/bin/bash',"\n";
		print TEMPIT '#SBATCH --job-name=',$inc,"\n";
		print TEMPIT '#SBATCH --ntasks=1',"\n";
		print TEMPIT '#SBATCH --mem=160000',"\n";
		print TEMPIT 'echo "',$inc," working\"\n";
	  print TEMPIT "perl $inc.pl\n";
		print TEMPIT "echo;\necho \"done\";\n";
		close (TEMPIT);
		`sbatch temp-$inc.sh`; 
  }
}

ENDBAMPERL
print BAMPERL $bamperlcontent."\n";
close BAMPERL;

##- - - - - - - ENDOF BAM FILE GENES SEQUENCE - - - - - - - - - -

# 5. FOR COMBINE STAGE && KAKS
#- - - - - - - KAKS CALCULATION GENES SEQUENCE - - - - - - - - - -
open(COMBINE,">combine.pl");
my $combinecontent = <<"ENDCOMBINE";
#!/usr/bin/perl
use threads;
use Thread::Queue;
ENDCOMBINE
print COMBINE $combinecontent."\n";
print COMBINE "\$sum=$sumofgenes;\n";
print COMBINE "\$path=\"$newpath\";\n";
print COMBINE "\$KAKS=\"$KAKS\";\n";
print COMBINE "chdir \"$newpath\";\n";
print COMBINE Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
$combinecontent = <<'ENDCOMBINE';
my ($ident, $addi, $tag) = (0,0,0);
my (@identi, @VAR);

OVER: {
while ($addi<$sum) {
  $ident++;
  push(@identi, "kaks-$ident");
	open(KAKS,">kaks-$ident.pl");
  print KAKS "#!/usr/bin/perl\n";
  print KAKS "use File::Basename;\n";
  print KAKS "chdir \"$path\";\n";
  print KAKS "\$sum=$sum;\n";
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
ENDKAKSPERL
print KAKS $kaksperlcontent."\n";
close KAKS;

#push variables for threads.
push @VAR, [splice @identi, 0, 1] while @identi;
$queue = new Thread::Queue();
my $builder=threads->create(\&main); #create thread for each subarray into a thread
push @threads, threads->create(\&processor) for 1..15; #execute 15 threads
$builder->join; #join threads
foreach (@threads){$_->join;}

sub main {
  foreach my $count (0..$#VAR) {
		while(1) {
			if ($queue->pending() < 100) {
				$queue->enqueue($VAR[$count]);
				last;
			}
		}
  }
  foreach(1..15) { $queue-> enqueue(undef); }
}

sub processor { my $query; while ($query = $queue->dequeue()){ parseinput(@$query); } }

sub parseinput {
  foreach my $inc (@_) {
		open (TEMPIT, ">temp-$inc.sh");
		print TEMPIT '#!/bin/bash',"\n";
		print TEMPIT '#SBATCH --job-name=',$inc,"\n";
		print TEMPIT '#SBATCH --ntasks=1',"\n";
		print TEMPIT '#SBATCH --mem=160000',"\n";
		print TEMPIT 'echo "',$inc," working\"\necho\n";
	  print TEMPIT "perl $inc.pl\n";
		print TEMPIT "echo;\necho \"done\";\n";
		close (TEMPIT);
		`sbatch temp-$inc.sh`;
		
  }
}
ENDCOMBINE
print COMBINE $combinecontent."\n";
close COMBINE;
## - - - - - - - ENDOF KAKS CALCULATION GENES SEQUENCE - - - - - - - - - -

#6. FOR FASTQTOFASTA SUBROUTINE
#- - - - - - - FASTA to FASTQ subroutine - - - - - - - - - -
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
##- - - - - - - ENDOF FASTA to FASTQ subroutine - - - - - - - - - -

# 7. EXTRACT KAKS RATIOS TO ONE FILES
open(EXTRACT,">extractkaksratios.pl");
my $extractcontent = <<"ENDEXTRACT";
#!/usr/bin/perl
use threads;
use Thread::Queue;
ENDEXTRACT
print EXTRACT $extractcontent."\n";
print EXTRACT "\$sum=$sumofgenes;\n";
print EXTRACT "chdir \"$newpath\";\n";
print EXTRACT Data::Dumper->Dump( [ \%GENE ], [ qw(*GENE) ] );
$extractcontent = <<'ENDEXTRACT';
my ($method, $seq, $verdict) = (0,0,0);
my (%Name, %METHOD, @array);
my (@filecontent, @chrscontent);

foreach my $genename (sort keys %GENE ){
  foreach my $number (1..$GENE{$genename}){
    my $fileloc = $genename."-".$number;
    open(KAKSIN, "<$fileloc.kaks"); @filecontent = <KAKSIN>; close (KAKSIN);
		open(CHRS, "<$fileloc.chrs"); @chrscontent = <CHRS>; close (CHRS);
		open(ALLOUT, ">$fileloc.finalkaks");
		
		($method, $seq, $verdict) = (0,0,0);
		undef %Name; undef %METHOD; undef @array;
		@header = split("\t", (shift @filecontent)); #header
		foreach (0..$#header) { #header attributes 
      if ($header[$_] =~ /Ka\/Ks/) { $verdict = $_ ; }
      elsif ($header[$_] =~ /Method/) { $method = $_; }
      elsif ($header[$_] =~ /Sequence/) { $seq = $_; }
		} #end foreach 
		unless ($verdict != 0) { next;} #exit if there is no header

		foreach (@filecontent) { #read content
      chomp;
      my @all = split("\t");
      my @name = split(/-([^-]\d*)$/, $all[$seq]);
      #if (length $name[0] > 1) { $Name{1} = "$name[0]\t$name[1]"; }
      unless (exists $METHOD{$all[$method]}) {
        $Name{1} = "$name[0]\t$name[1]";
        $METHOD{$all[$method]} = $all[$verdict];
      }
		} #end foreach
		my @finalheader = map {$_; } sort keys %METHOD;
		print ALLOUT "Name\tCount\tChrlocation\t";
		print ALLOUT join("\t",@finalheader),"\tAVERAGE\tVERDICT\n";
		print ALLOUT "$Name{1}\t", substr(((split("\t", $chrscontent[0]))[1]),0,-1) ,"\t";
		
		foreach (sort keys %METHOD){ 
			print ALLOUT $METHOD{$_},"\t"; 
      if (($METHOD{$_} > 0)  && ($METHOD{$_} != /NA/i) && ($METHOD{$_} != "50")) {
				push @array, $METHOD{$_};  
      }
		}
		if ($#array > 3) {
      my $mean =  &average(\@array);
      print ALLOUT "$mean";
      if (sprintf("%.0f",$mean) > 1) { print ALLOUT "\tpositive"; }
      elsif (sprintf("%.0f",$mean)  < 1) { print ALLOUT "\tnegative"; }
      elsif (sprintf("%.0f",$mean) == 1) { print ALLOUT "\tneutral"; }
      print ALLOUT "\n";
		} else { print ALLOUT "\n"; }
		close (ALLOUT);
	}
}

sub range {
	my ($data) = @_;
	my $max = (sort {$b <=> $a} @$data)[0];
	my $min = (sort {$a <=> $b} @$data)[0];
	return $max, $min;
}
sub median {
	my $median; my ($data) = @_;
	my $length = @$data;
	if ($length > 1) {
	my @median = (sort {$a <=> $b} @$data)[ int($length/2), ceil($length /2) ];
	$median = sprintf("%.3f", (&sum(\@median)/2));
	} else {
		$median = sprintf("%.3f",@$data[0]);
	}
	return $median;
}
sub sum {
	my ($data) = @_;
	if (not @$data) {
		die("Empty arrayn");
	}
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	return $total;
}
sub average {
  my ($data) = @_;
  if (not @$data) {
     die("Empty arrayn");
  }
  my $total = &sum($data);
  my $average = $total / @$data;
  return $average;
}
sub stdev {
	my($data) = @_;
	if(@$data == 1){
		return 0;
	}
	my $average = &average($data);
	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($average-$_) ** 2;
	}
	my $std = ($sqtotal / (@$data-1)) ** 0.5;
	return $std;
}

ENDEXTRACT
print EXTRACT $extractcontent."\n";
close EXTRACT;


# 8. CLEANUP
#- - - - - - - CLEAN UP files - - - - - - - - - -
open(CLEAN,">cleanup.pl");
print CLEAN "#!/usr/bin/perl\n";
print CLEAN "chdir \"$newpath\";\n";
my $cleancontent = <<'ENDCLEAN';
#All the nucleotide files
`mkdir -p nucleotides; mv *nuc nucleotides`;
#All the peptide files
`mkdir -p fastq; mv *pep fastq`;
#All fasta files
`mkdir -p fasta; mv *fna fasta`;
#All AXT files
`mkdir -p axtfiles; mv *axt axtfiles`;
#All KaKs files
`mkdir -p KAKS; mv *.kaks KAKS`;
#All FinalKAKS files
`mkdir -p FINALKAKS; mv *.finalkaks FINALKAKS`;

#All Name files
`mkdir -p NAME; mv *name NAME`;
#All Scripts
`mkdir -p SCRIPTS; mv samperl*pl kaks*pl *pm SCRIPTS`;
#All Chrs
`mkdir -p CHRS; mv *chrs CHRS`;
print "All sorted into respective folders\n";

ENDCLEAN
print CLEAN $cleancontent."\n";
close CLEAN;
##- - - - - - - ENDOF CLEAN UP files - - - - - - - - - -

# 8. JOB EXECUTION
#- - - - - - - JOB EXECUTION - - - - - - - - - -
print "Job execution: 1. refperl.pl\n";
open (TEMPIT, ">temp-refperl.sh");
print TEMPIT '#!/bin/bash',"\n";
print TEMPIT '#SBATCH --job-name=refperl',"\n";
print TEMPIT '#SBATCH --ntasks=1',"\n";
print TEMPIT '#SBATCH --mem=160000',"\n";
print TEMPIT "echo 'refperl working'\n";
print TEMPIT "perl refperl.pl\n";
print TEMPIT "echo;\necho \"done\";\n";
close (TEMPIT);
`sbatch temp-refperl.sh`;
		
print "Job execution: 2. bamperl.pl\n";
open (TEMPIT, ">temp-bamperl.sh");
print TEMPIT '#!/bin/bash',"\n";
print TEMPIT '#SBATCH --job-name=bamperl',"\n";
print TEMPIT '#SBATCH --ntasks=1',"\n";
print TEMPIT '#SBATCH --mem=160000',"\n";
print TEMPIT "echo 'bamperl working'\n";
print TEMPIT "perl bamperl.pl\n";
print TEMPIT "echo;\necho \"done\";\n";
close (TEMPIT);
`sbatch temp-bamperl.sh`;

print "Job execution: 3. combine.pl\n";
open (TEMPIT, ">temp-combine.sh");
print TEMPIT '#!/bin/bash',"\n";
print TEMPIT '#SBATCH --job-name=combine',"\n";
print TEMPIT '#SBATCH --ntasks=1',"\n";
print TEMPIT '#SBATCH --mem=160000',"\n";
print TEMPIT "echo 'combine working'\n";
print TEMPIT "perl combine.pl\n";
print TEMPIT "echo;\necho \"done\";\n";
close (TEMPIT);
print "Run `sbatch temp-combine.sh`\n";

print "Job execution: 4. extractkaksratios.pl and 5. cleanup.sh \n";
open (TEMPIT, ">temp-extract_cleanup.sh");
print TEMPIT '#!/bin/bash',"\n";
print TEMPIT '#SBATCH --job-name=ex-clean',"\n";
print TEMPIT '#SBATCH --ntasks=1',"\n";
print TEMPIT '#SBATCH --mem=160000',"\n";
print TEMPIT "echo 'extractkaksratios working'\n";
print TEMPIT "perl extractkaksratios.pl\n";
print TEMPIT "perl cleanup.pl\n";
print TEMPIT "echo;\necho \"done\";\n";
close (TEMPIT);
print "Run `sbatch temp-extract_cleanup.sh`\n";

print "Finished!!!!\n";
##- - - - - - - ENDOF JOB EXECUTION - - - - - - - - - -

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
   threads;
   Thread::Queue;
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
