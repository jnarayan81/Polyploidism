# $Id$ All messages
# Perl module for Polyploidism polySubs;
# Author: Jitendra Narayan <jnarayan81@gmail.com>
# Copyright (c) 2015 by Jitendra. All rights reserved.
# You may distribute this module under the same terms as Perl itself

##-------------------------------------------------------------------------##
## POD documentation - main docs before the code
##-------------------------------------------------------------------------##

=head1 NAME

EBALib::CommonSubs  - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=cut

=head1 CONTACT

Jitendra <jnarayan81@gmail.com>

=head1 APPENDIX

The rest of the documentation details each of the object methods.

=cut

##-------------------------------------------------------------------------##
## Let the code begin...
##-------------------------------------------------------------------------##

package polySubs;
use warnings::register;
#use Statistics::R;

use Exporter qw(import);
 
our @EXPORT_OK = qw(yesORno deldir);

## Pause the program, activate when press enter
## returns 1 for y
sub yesORno { 
	my ($min, $ran, $gen)=@_;
	my $yflag = 0;
	if ($ran) {
		print "The contigs/scaffolds below $min bp will be ignored, Wanna continew [y/n] : ";
		}
	elsif ($gen) {
		print "Great!! You are using predicted genes locations. Wanna continew [y/n] : ";
		}
	my $input = '';
	while ( $input !~ /y|n/i ) {
		print "\b";
		$input = <STDIN>;
		if ( $input =~ /^y/i ) {	
		print "Cheers!\n";
		}
		elsif ($input =~ /^n/i ) {	
		print "Thank you! I guess you need manual try -h:\n";
		}
		else {print "Ohh wait !! Did you sleep well last night ;) Try again [y/n]:\n";}
		chomp $input;
	}
	if ( $input =~ /^y/i ) { $yflag = 1; }
	return $yflag;
}

##!/usr/bin/perl
#deldir("test"); # or deldir($ARGV[0]) to make it commandline
 
sub deldir {
  my $dirtodel = pop;
  my $sep = '//'; #change this line to "/" on linux.
  opendir(DIR, $dirtodel);
  my @files = readdir(DIR);
  closedir(DIR);
 
  @files = grep { !/^\.{1,2}/ } @files;
  @files = map { $_ = "$dirtodel$sep$_"} @files;
 
  @files = map { (-d $_)?deldir($_):unlink($_) } @files;
 
  rmdir($dirtodel);
}

sub EBAWelcome {
my $VERSION = shift;
print "\nWelcome to help section of Polyploidism tool, version $VERSION\n"; 
printUsage(); ## display documentation
}

# This function prints the script usage mode
sub printUsage {
my ($message) = @_;
if (defined $message) {
    print STDERR "\nERROR: $message\n";
  }

print << "End_Print_Usage";

Usage:

perl Polyploidism.pl -f sampleSeq.fa -o TESTOUT2 -i 90 -p -g -s caenorhabditis

Mandatory parameters: -g (gene based)

  -f <infile> 		Provide the contigs/scaffods file.

  -o <outfile> 		Provide the outfile. 
  
  -g			Genes based estimates

  -i <num> 		Identify percentage for lastz.

  -z <num> 		Max number of ploidy to plot.

  -s <species>		Trained species name for gene prediction


perl Polyploidism.pl -f sampleSeq.fa -o TESTOUT2 -i 90 -p -r -l 100 -c 5 

Mandatory parameters: -r (random based)

  -f <infile> 		Provide the contigs/scaffods file.

  -o <outfile> 		Provide the outfile. 

  -l <length> 		Length of the string.

  -r			For randome sequence split

  -i <num> 		Identify percentage for lastz.

  -z <num> 		Max number of ploidy to plot.

  -c <count> 		Number of string.


End_Print_Usage

}

sub EBAWho {
my $VERSION = shift;
print "\nPlolyplodism, version $VERSION by Jit and Nico\n"; 

}

sub process_contig
  {
  my($contig, $sequence) = @_;
  my($len, $ACGTbases, $ATbases, $GCbases, $nonACGTbases);

  # Remove Contig name prefix?
  $contig =~ s/^.*([Cc]ontig)/$1/ if $SHORTEN_CONTIG_NAMES;
  $len = length($sequence);
  push @CONTIG_LENGTHS, $len;
  $Total_Bases += $len;
  $Max_Bases = $len if $Max_Bases < $len;
  $Min_Bases = $len if $Min_Bases > $len || $Min_Bases < 0;

  $ATbases = ($sequence =~ tr/aAtT/aAtT/);
  $GCbases = ($sequence =~ tr/cCgG/cCgG/);
  $ACGTbases = $ATbases + $GCbases;
  $nonACGTbases = $len - $ACGTbases;
  if ($ACGTbases)
    {
    $GC_per_cent = sprintf "%.1f", 100 * $GCbases / $ACGTbases;
    }
  else
    {
    $GC_per_cent = '-';
    }
  $Total_GC += $GCbases;
  $Total_ACGT += $ACGTbases;
  if ($nonACGTbases)
    {
    my $more_Max_Nons = ($Max_Nons < $MAX_PATTERN_MIN_RPT) ?
      $Max_Nons + 1 : $MAX_PATTERN_MIN_RPT;
    my @Nons = ($sequence =~ /[^acgtACGT]{$more_Max_Nons,}/g);
    foreach (@Nons)
      {
      my $l = length $_;
      $Max_Nons = $l if ($Max_Nons < $l);
      }
    $Total_Non_ACGT_Ends += length $1 if ($sequence =~ /^([^acgtACGT]+)/);
    if (substr($sequence, -1) =~ /[^acgtACGT]+$/)
      {
      my $rs = reverse $sequence;
      $Total_Non_ACGT_Ends += length $1 if ($rs =~ /^([^acgtACGT]+)[acgtACGT]/);
      }
    my $more_Max_Ns = ($Max_Ns < $MAX_PATTERN_MIN_RPT) ?
      $Max_Ns + 1 : $MAX_PATTERN_MIN_RPT;
    my @Ns = ($sequence =~ /[nN]{$more_Max_Ns,}/g);
    foreach (@Ns)
      {
      my $l = length $_;
      $Max_Ns = $l if ($Max_Ns < $l);
      }
    $Total_N_Ends += length $1 if ($sequence =~ /^([nN]+)/);
    if (uc substr($sequence, -1) eq 'N' && uc substr($sequence, 0, 1) ne 'N')
      {
      my $rs = uc reverse $sequence;
      $Total_N_Ends += length $1 if ($rs =~ /^(N+)/);
      }
    }

    my $StringRes= "$contig\t$len\t$GC_per_cent\t$nonACGTbases";

	return ($StringRes);

  } # end process_contig


sub reformatTRF {
my ($dat, $flag, $id)=@_;

if (($dat eq "") || ($flag eq "")) {
	print "Usage: program filename flag\nPlease provide the required parameters and run again\n";
	exit;
}
my @arr = split("dat", $dat); my $file = $arr[0];
my $filename="$file"."dat.parse";
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print  $fh ("Repeat Start\tRepeat End\tPeriod Size\tCopy No.\tAlignment Score\tConsensus\tsequence_id\n");

system("cat ".$file."*txt.html > ".$file."tmp");

local $/ = "\n";
open (DATA, '<', $dat) or die "Can not open TRF data file";
while(my $line=<DATA>) {
my ($rep_start, $rep_end, $period_size, $copy_no, $pattern_size, $percent_match, $percent_indel, $align_score, $a_percent, $c_percent, $g_percent, $t_percent, $entropy, $consensus, $repeat);			
chomp($line);

if ($line =~ /^[0-9]/) {
	my @arr = split(" ", $line);

	$rep_start = $arr[0];
	$rep_end = $arr[1];
	$period_size = $arr[2];
	$copy_no = $arr[3];
	$pattern_size = $arr[4]; 
	$percent_match = $arr[5];
	$percent_indel = $arr[6];
	$align_score = $arr[7];
	$a_percent = $arr[8];
	$c_percent = $arr[9];
	$g_percent = $arr[10];
	$t_percent = $arr[11];
	$entropy = $arr[12];
	$consensus = $arr[13];
	$repeat = $arr[14];
	print $fh ("$rep_start\t$rep_end\t$period_size\t$copy_no\t$align_score\t$consensus\t$id\n");
	}

}
close $fh;
close DATA;
flanking($file, $flag, $id);
system("paste ".$file."dat.parse ".$file."txt.parse > ".$file."final.parse");

=pod
USAGE:

 trfparser datfilename flag_value{0 or 1}
 
 Enter 0 as flag_value if you need to use only dat file, and 1 if you want to extract information from both dat and txt/html file.


AVAILABLE FIELDS:

By default, following information will be included in the .final.parse file:

.dat file: Repeat Start, Repeat End, Period Size, Copy No., Alignment Score, Consensus
.txt.html file: Left Flanking Sequence, Right Flanking Sequence

However, you can always modify the code to display other information, following is a list of all available fields and description:

.dat file

$rep_start	:		Indices of the repeat relative to the start of the sequence
$rep_end 	
$period_size	:		Period size of the repeat
$copy_no	:		Number of copies aligned with the consensus pattern
$pattern_size	:		Size of consensus pattern (may differ slightly from the period size)
$percent_match	:		Percent of matches between adjacent copies overall
$percent_indel	:		Percent of indels between adjacent copies overall
$align_score	:		Alignment score
$a_percent	:		Percent composition for each of the four nucleotides
$c_percent
$g_percent
$t_percent
$entropy	:		Entropy measure based on percent composition
$consensus	:		Consensus sequence
$repeat		:		Repeat sequence

.txt.html file

$start		:		Indices of the repeat relative to the start of the sequence
$end			
$left_start	:		Indices of the Left flanking sequence relative to the start of the sequence
$left_end		
$right_start	:		Indices of the Right flanking sequence relative to the start of the sequence
$right_end		
$left_seq	:		Left flanking sequence
$right_seq	:		Right flanking sequence

=cut

}

sub flanking {
my ($file, $flag)=@_;
my ($line_txt, $count, $start, $end, $left_start, $left_end, $right_start, $right_end, $left_seq, $right_seq);
$count = 0;	

my $fileN="$file"."txt.parse";
open (my $fh2, '>', $fileN) or die "Could not open file $fileN $!";
return if ($flag==0);
print $fh2 ("Left Flanking Sequence\tRight Flanking Sequence\tsequence_id\n");
local $/ = "\n";
open(TXT, '<', $file."tmp") or die "Can not open trf txt/html file";
while($line_txt=<TXT>) {
	chomp($line_txt);
	my @arr = split("[ :-]", $line_txt);

	if ($line_txt =~ "    Indices:") {
		$start = $arr[6]; 
		$end = $arr[8];
		$count = 1;
	}
			
	elsif  ($line_txt =~ "Left flanking sequence:") {
		$left_start = $arr[5]; 
		$left_end = $arr[9];
		chomp($line_txt = <TXT>);
		until ($line_txt eq "") {
			$left_seq .= $line_txt;
			chomp($line_txt = <TXT>);
		}
	}

	elsif  ($line_txt =~ "Right flanking sequence:") {
		$right_start = $arr[5]; 
		$right_end = $arr[9];
		chomp($line_txt = <TXT>);
		until ($line_txt eq "") {
			$right_seq .= $line_txt;
			chomp($line_txt = <TXT>);
		}
	}

	elsif((($line_txt =~ "Found at i:") || (eof)) && ($count == 1)){
		print  $fh2 ("$left_seq\t$right_seq\t$id\n");
    		$start = $end = $left_start = $left_end = $right_start = $right_end = $count = 0;
    		$left_seq = $right_seq = "";
	}
  }
close TXT;
close $fh2;
}


#subroutine to check tri di bases
sub process_di_tri {
  my($contig, $sequence, $extended) = @_;
  my $len = length($sequence);

  # count mono-nucleotides and total
  $x_counters{'total1'} += $len;
  my $newa = ($sequence =~ tr/a/a/);
  $x_counters{'a'} += $newa;
  my $newc = ($sequence =~ tr/c/c/);
  $x_counters{'c'} += $newc;
  my $newg = ($sequence =~ tr/g/g/);
  $x_counters{'g'} += $newg;
  my $newt = ($sequence =~ tr/t/t/);
  $x_counters{'t'} += $newt;
  $x_counters{'n'} += $len - ($newa + $newc + $newg + $newt);

  # count di- and tri-nucleotides
  for ($i = 0; $i < $len - 2; $i++) {
    print STDERR "$contig ... $i\n"
      if (($i % 500000) == 0) && ($i != 0);
    $x_counters{'total2'}++;
    $x_counters{'total3'}++;
    my $tri = substr($sequence, $i, 3);
    my $di = substr($tri, 0, 2);
    if ($di !~ /n/) {
      $x_counters{$di}++;
      }
    else {
      $x_counters{'nn'}++;
      }
    if ($tri !~ /n/) {
      $x_counters{$tri}++;
      }
    else {
      $x_counters{'nnn'}++;
      }
    } # end for ($i = 0; $i < $len - 2; $i++)
  if ($len > 1)		# take care of last di-nucleotide
    {
    $x_counters{'total2'}++;
    my $last2 = substr($sequence, -2, 2);
    if ($last2 !~ /n/) {
      $x_counters{$last2}++;
      }
    else {
      $x_counters{'nn'}++;
      }
    }

  # remove non-ACGT from counts unless requested
  unless ($extended) {
    $x_counters{'total1'} -= $x_counters{'n'};
#    $x_counters{'n'} = 0;
    $x_counters{'total2'} -= $x_counters{'nn'};
    $x_counters{'nn'} = 0;
    $x_counters{'total3'} -= $x_counters{'nnn'};
    $x_counters{'nnn'} = 0;
    }

  #if ($print) {
    #print STDOUT "\n$CONTIG=$contig\n\n" unless ($TOTALS > 1);
    #print_extended(%x_counters) unless ($TOTALS > 1);
    # update counters
print "$x_counters{'total1'} - $x_counters{'total2'} - $x_counters{'total3'} + $x_counters{'n'} - $x_counters{'nn'} -  $x_counters{'a'} - $x_counters{'c'} - $x_counters{'g'} - $x_counters{'t'}";
    #$x_totals{$_} += $x_counters{$_} foreach keys %x_counters;
    #%x_counters = map { ($_, 0) } keys %x_counters;
    #}
} # end do tri caclculation



sub reformatNplot {
my($outLoc,$file, $hash_ref) = @_;
my %ids=%$hash_ref; my @lines;
my $plotFile="$file.final";

open my $outF, '>', $plotFile or die "Could not create: $!\n";;
print $outF "Name\tPloidy\tHits\tCount\tSeq\tGC\tGC_per\tnon_ATGC\tPercentage\n" if (-z "$plotFile");
foreach my $key (keys %ids) {
$/ = "\n";  # read by \n
open FILE, $file;
   while (<FILE>) {
    	chomp $_;
    	my $line = trim ($_);
	next if $. == 1; #Ignore header
	next if ($line =~ /^\s*$/);
    	my @values = split('\t', $_);
    	if ($key eq $values[0]) { 
		push @lines, $line;
		}
	}
   close FILE;
my @sorted_lines = sort { (split "\t", $a)[8] <=> (split "\t", $b)[8] } @lines;
undef @lines;
my $bestPloidy=trim($sorted_lines[-1]);
print $outF "$bestPloidy\n";
	#foreach my $linewa (@sorted_lines) {
	#	my @val = split('\t', $linewa);
	#	print "$keys\t\n";
	#}
}
close $outF;
}


sub reverse_complement_IUPAC {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}


sub storeAlles {
my ($id, $gcor1, $gcor2, $gname, $size, $ghits, $genecnt, $slen, $expectedPoly, $strand, $outDIR)=@_;
open FH1, '<', $ghits or die "Could not open file $ghits $!";
open FH2, '>>', "$outDIR/alles.loc" or die "Could not open file $!";
if (-z "$outDIR/alles.loc") {print FH2 "SeqName\tLength\tGeneCount\tGeneName\tGeneCor1\tGeneCor2\tStrand\tallesCount\tpreDucplication\tAllesCordinates\n";}
$/ = "\n";  # read by \n
open FH1, $ghits;
   while (<FH1>) {
    	chomp $_;
    	my $line = trim ($_);
	#next if $. == 1; #Ignore header
	next if ($line =~ /^\s*$/);


#qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs
#qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs
# xxx qseqid sstrand qlen qstart qend sallseqi sstrand slen sstart send qlen nident pident nident pident
#104726,seq0-seq0,+,1109,0,1109,seq0,+,1681,359,1468,1109/1109,100.0%,1109/1109,100.0%

    	my @values = split('\t', $_);
	#next if (($values[5]-$values[4])/2) <= ($values[10]-$values[9]); # Next if only 1/2 of the length is matched
	next if $values[-1] <= 50; #Next if only 1/2 of the length is matched
	push @allN, $values[4];
	#my $cor=extractAllelsCor($line);
	my $cor="$values[0]:$values[1]:$values[2]:$values[3]:$values[4]:$values[5]";	
	my $ORes = checkCorOverlaps ($gcor1, $gcor2, $values[5], $values[6]);
	if ((!$ORes) and ($id eq $values[0])) { 
		push @corv, "$cor\t";
		}
	elsif (($ORes) and ($id ne $values[0])) {
		push @corv, "$cor\t";
		}
	elsif ((!$ORes) and ($id ne $values[0])) {
		push @corv, "$cor\t";
		}
	elsif (($ORes) and ($id eq $values[0])) { next;}	
	else { print "Hello .. You forgot some conditions\n"; print "$line\n"; exit;}
   }
   close FH1;
#my $decPal=dupsInArray(\@allN);
my $allesN= scalar @corv;
my $preDup="NO";
if ($allesN > $expectedPoly) {$preDup="DupGenes"}
print FH2 "$id\t$slen\t$genecnt\t$gname\t$gcor1\t$gcor2\t$strand\t$allesN\t$preDup\t";
foreach my $v (@corv) {print FH2 $v;} print FH2 "\n";
undef @corv; undef @allN;
close FH2;
#again re-set it as the loop needs it
$/ = "\n>";  # read by \n
}

# Checks if a provided two coordinates overlaps or not it return 1 if overlaps
sub checkCorOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

sub dupsInArray {
my $array_ref = shift;
my %seen; my @array=@$array_ref;
foreach my $string (@array) {
    $string=reverse_complement_IUPAC($string);
    next unless $seen{$string}++;
    return "PalinGenes"; # if duplicates
}
}


sub extractAllelsCor {
my $hits=shift;
my @line = split("\t", $hits);
#$allesCor="$line[6]:$line[9]:$line[10]:$line[7]";

#qseqid staxid qstart qend sallseqi sstart send evalue length frames qcovs
$allesCor="$line[2]:$line[3]:$line[5]:$line[6]";
return $allesCor;
}

#Store column in array

sub storeInHash {
my $file=shift;
my %hash = ();
my @arr;

local $/ = "\n";  # read by \n
open (my $fh, "<", $file) or die "Can't open the file $file: ";

while (my $line =<$fh>)
{
    $line = trim ($line);
    next if $. == 1; #Ignore header
    my @key = split("\t", $line);
    $hash{$key[0]} = 1;
}
close $fh;
return %hash;
}

#Parse the augustus outfile
sub parseAugustus {
my $seqfile = shift; #"preGenes.txt"
my %genCor;
open $FL, $seqfile;
local $/ = "\n";  # read by \n
while (<$FL>) {
    chomp;
    next if /^#/;  # discard comments
    my $seq = $_;
    my @seqLine = split('\t', $seq);
    next if $seqLine[2] ne 'gene';
    my $region="$seqLine[0]\t$seqLine[3]\t$seqLine[4]\t$seqLine[8]\t$seqLine[6]";
    #print "$seqLine[0]\t$seqLine[1]\t$seqLine[2]\t$seqLine[3]\t$seqLine[4]\t$seqLine[5]\t$seqLine[6]\t$seqLine[8]\n";
    $genCor{$region}=$.;
}
return %genCor;
close $FL;
}

#palindrome_r uses recursion. Return 1 if palindrome and 0 if not.
sub palindrome_r {
    my $s = (@_ ? shift : $_);
    if (length $s <= 1) { return 1; }
    elsif (substr($s, 0, 1) ne substr($s, -1, 1)) { return 0; }
    else { return palindrome_r(substr($s, 1, -1)); }
}

sub Mapper {
use Statistics::R;
my ($InFile)=@_;
print "Creating graph for $InFile .. .\t\n";

my $R = Statistics::R->new();
$R->startR;


# 20-06-2016
# works with R 2.15.2 and ggplot 0.9.3.1
# Check ggplot2 help forums or contact Jitendra Narayan jnarayan81@gmail.com if something doesn't run
# because of updated programs/packages
#Store the data to Alien variable

#Alien <- read.csv("TESTOUT2.csv")

$R->run(qq`
	#Set the current working directory
	#setwd("$ENV{PWD}");
	#pdf("$ENV{PWD}/$InFile.PloidyPlot.pdf") 
	Alien <- read.table( "$InFile", sep="\t", header=TRUE)
	names(Alien);
	library("ggplot2")
	# Basic box plot
	p10 <- ggplot(Alien, aes(x = Ploidy, y = GC_per)) + 
        	geom_boxplot(aes(group = cut_width(Ploidy, 0.25), fill = factor (Ploidy), outlier.colour = "red", outlier.shape = 1,notch = TRUE)) + geom_jitter(width = 0.2)
	p10

dev.off()`);

$R->stopR();

}


sub isatty()
{
    my $tty = `/bin/ps -p $$ -o tty --no-headers`;
    $tty =~ s{[\s?]}{}g;
    return $tty;
}


sub KaKs {
use Bio::SeqIO;
my ($inFile, $sequence)=@_;
open (my $ifh, "<", $inFile) or die "Can't open the file $inFile: ";
local $/ = "\n";  # read by \n
while (<$ifh>) {
    chomp;
    next if /^#/;  # discard comments
    my $line = $_;
    my @hitLine = split('\t', $line);

my $in_file = $sequence;
my $start_pos = $hitLine[4]+1; # Needed to to as it start from 1 not zero as in lastz hit
my $end_pos = $hitLine[5];

my $in = Bio::SeqIO->new ( -file => $in_file, -format => 'fasta');
my $out = Bio::SeqIO->new( -file => ">>abba.out", -format => 'fasta');


while (my $seq = $in->next_seq() ) {

    $seq->display_id( $seq->display_id() . "_$start_pos-$end_pos" );
    $out->write_seq( $seq->trunc($start_pos, $end_pos) );
}
}
printLines(50, '*');
#system ('muscle -in abba.out -out abba.out.aln -quiet -clwstrict -seqtype nucleo');
0 == system ('muscle -in abba.out -out abba.out.aln -quiet') or die "Failed to muscle the sequence\n";   
#system ("awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' abba.out.aln");
#system ('perl -pe '/^>/ ? print "\n" : chomp' abba.out.aln > out.fasta')
#unlink ('abba.out') # Unlink after the first check

}

#Trim the lines
sub trim($) {
 my $string = shift;
 $string =~ s/^[\t\s]+//;
 $string =~ s/[\t\s]+$//;
 $string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
 return $string;
}

sub ManualHelp {
print("Unrecognized option(s)!! Please check manual OR Try --help \n\n")
}

sub printLines {
my ($num,$sign)=@_;
my $seplines = ($sign x $num)."\n";
print "\n$seplines";
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
  print LOG $message;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###FUNCTION TO MOVE THE FILES WITH WILD
sub moveFiles {
    my ( $source_ref, $arc_dir ) = @_;
    my @old_files = @$source_ref;
   foreach my $old_file (@old_files)
         {
    #my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    #my $new_file = $arc_dir . $short_file_name;
    move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
   }
}


#$SIG{INT} = \&polySubs::interrupt;
sub interrupt {
	print STDERR "Caught a control c!\n";
	exit;  # or just about anything else you'd want to do
}

1;
