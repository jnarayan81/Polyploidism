#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Temp qw(tempfile);
use File::Copy;
use Cwd qw();
use FindBin qw( $RealBin );
use lib $RealBin;
use polySubs;
use Getopt::Long;
use File::Remove 'remove';
use File::Path qw(make_path remove_tree);
use Fcntl qw( LOCK_EX LOCK_NB );
use Bio::SeqIO;

#USAGE:perl Polyploidism.pl -i sampleSeq.fa -o TESTOUT2 -s 100 -c 5

print <<'WELCOME';   

                ||----P |
                ||     ||  
	 Polyploidism v0.7 [Oct 15, 2017]  - [+]

Citation - automated estimation of polyploidy nature of the genome

License: Creative Commons Licence
Bug-reports and requests to: jitendra.narayanATunamur.be and nicolas.debortoliATunamur.be

Keep an eye at following line: 
#Set up the path of all required tools - stored in the same folder
$ENV{'PATH'} = "/bin:/usr/bin:/usr/bin/env:$ENV{PWD}/augustus.2.5.5/bin: $ENV{PWD}/augustus.2.5.5/config: $ENV{PWD}/lastz-distrib-1.03.73/src: $ENV{PWD}/KaKs_Calculator2.0/bin/Linux: $ENV{PWD}/muscle: $ENV{PWD}/dotter: $ENV{PWD}/trf";


FOR GENE try: perl Polyploidism.pl -f xaa.fa -o TESTOUT3 -s elegans -q 50 -e 4 -g -i 70 -t 7 -z 5
FOR RANDOM: perl Polyploidism.pl -f sampleSeq.fa -o TESTOUT3 -s elegans -q 50 -e 4 -r -i 70 -t 7 -l 500 -c 2 -z 5

-------------------------------------------------------------------
WELCOME

my (
	$infile,
	$outfile, 
	$length,
	$count,
	$random,
	$species,
	$expected,
	$qcov,
	$genes,
	$thread,
	$zip,
	$identity,
	$logfile,
);

# Default option setting for EBA tool
my $VERSION=0.6;
my $verbose=0; 	# Verbose set to 0;
my $blastEval= '5e-5'; # -evalue $blastEval in blast line
my %options = ();
$length=0;
$count=0;
$thread=1; # Default thread to run

GetOptions(
	\%options,
    	'infile|f=s'    	=> \$infile,        	## Infile
    	'outfile|o=s' 		=> \$outfile,           ## Outfile
	'length|l=i' 		=> \$length,		## Size of sub-string
    	'count|c=i' 		=> \$count, 		## Number of times 
	'identity|i=i' 		=> \$identity, 		## identity percentage
	'blastEval|b=i' 	=> \$blastEval, 	## blast evalues
	'qcovl|q=i' 		=> \$qcov, 	## quesry coverage
	'species|s=s' 		=> \$species, 		## Name of the trained species 
	'expected|e=i' 		=> \$expected, 	## expected polyploidy
	'random|r' 		=> \$random, 	## check random 
	'genes|g' 		=> \$genes, 	## check random 
	'thread|t=i' 		=> \$thread, 	## thread 
	'zip|z=i' 		=> \$zip, 	## zip the ploidy to plot
    	'help|?|h!'     	=> sub { polySubs::EBAWelcome($VERSION) },
   	'who|w!'     		=> sub { polySubs::EBAWho($VERSION) },
	'verbose' 		=> \$verbose,
    	'logfile=s' 		=> \$logfile,		## logfile
	
) or die polySubs::EBAWelcome($VERSION);

if ($random) { if ((!$infile) or (!$outfile) or (!$length) or(!$count) or (!$expected) or (!$thread) or (!$zip)) { polySubs::printUsage(); exit; }}
if ($genes) { if ((!$infile) or (!$outfile) or (!$species) or (!$expected) or (!$thread) or (!$zip)) { polySubs::printUsage(); exit; }}

#Check if allready running same script anywhere.
flock(DATA,LOCK_EX|LOCK_NB)
  or  die "This script ($0) is already running.\n Wait the first instances to finish and try again later\n Sorry for inconvenience\n";

#Check if other instances is running
my $cmd = 'ps -eaf | grep -c "\< '."$0".'\>"'; 
chomp(my $instances = `$cmd`);
if($instances > 1) { print "Other instances of your program is running\n Please let them complete first\n"; exit; }

#Check the threads in ur CPU
chomp(my $cpu_count = `grep -c -P '^processor\\s+:' /proc/cpuinfo`);
if ($cpu_count < $thread) {print "Sorry, you have only $cpu_count thread in your PC\n"; exit;}

#Set up the path of all required tools - stored in the same folder
$ENV{'PATH'} = "/bin:/usr/bin:/usr/bin/env:$ENV{PWD}/augustus.2.5.5/bin: $ENV{PWD}/augustus.2.5.5/config: $ENV{PWD}/lastz-distrib-1.03.73/src: $ENV{PWD}/KaKs_Calculator2.0/bin/Linux: $ENV{PWD}/muscle: $ENV{PWD}/dotter: $ENV{PWD}/trf";

if (!$genes and !$random) { print "ERROR: You might forgot to select one of the estimation method: Gene:g or Random:r\n"; exit(0);}

remove_tree( "RES$outfile"); # Remove the outDIR if exist
#Make OUTDIR 
my $path = Cwd::cwd();
my $outdirname ="RES$outfile";
mkdir $outdirname, 0755;
my $outDIR = "$path/$outdirname";

#Remove tmpRES folder
remove_tree( "tmpRES");

#also delete
if ( -e "$outfile/alles.loc" ) { unlink ("$outfile/alles.loc") or die "$outfile/alles.loc: $!"; }

#if ( -e $outDIR/$outfile ) {
#	print "The $outDIR/$outfile already exists, deleting now\n";
#        unlink($outDIR/$outfile) or die "$outDIR/$outfile: $!"
#    }

my $min_length = ( $count - 1 ) * ( 2 * $length - 1 ) + $length;

my $decision=polySubs::yesORno($min_length, $random, $genes); if ($decision != 1) {exit;}

#To check the palindrome
#if ($palindrome) { 
#	print "You turn palindrome flag on, I will check for the palindrome in your infile\n";
#	polySubs::palindrome($infile,$identity);
#	}


#Create a local database for BLAST
system ("makeblastdb -in $infile -parse_seqids -dbtype nucl -out localDATA");

#Replace the unusual character(|) from fasta file header, if any
#system (perl -pi -e 's/\|/_/g' $infile");

#local $/ = "\n>";  # read by FASTA record
#open my $infh,  '<', $infile;
open my $outfh, '>>', "$outDIR/$outfile";
print $outfh "Name\tPloidy\tHits\tCount\tSeq\tGC\tGC_per\tnon_ATGC\tPercentage\tLength\n";
#my $infh  = \*DATA;
#my $outfh = \*STDOUT;
my $seqio = Bio::SeqIO->new(-file => "$infile", '-format' => 'Fasta');

while(my $string = $seqio->next_seq) {
    	my $seq = $string->seq;
	my $id = $string->display_id;
	my $len=length($seq);
	polySubs::printLines(50, '-');
	print "Working on $id ->>\n\n";
	#print "$id\t$min_length\t$len\n";

#If genes based analysis selected
if ($genes) {
	#print "$id-$pos - $substring\n";
	my $tmpfh = new File::Temp( UNLINK => 1 );
	print $tmpfh ">$id\n$seq\n";
	print "Predicting genes region in $id with Augustus:\n";
	0 == system ("augustus --species=$species $tmpfh --outfile=preGenes.txt --AUGUSTUS_CONFIG_PATH=$ENV{PWD}/augustus.2.5.5/config") or die "Failed to augustus the $id sequence\n"; ;
	my %genesCors=polySubs::parseAugustus("preGenes.txt");
	$count= keys %genesCors; #Store the round/count == no of genes in a sequence	
	#if (!!%genesCors) { next; }
	foreach my $key (sort(keys %genesCors)) {
		my @gval = split('\t', $key);
		my $siz=$gval[2]-$gval[1];
		my $seqstr = substr $seq, $gval[1], $siz;
		my $seqstring="";
		#if ($gval[4] eq "-") {
		#Extracted the sequence and reverse complement it, if minus strand
		#$seqstring = polySubs::reverse_complement_IUPAC ($seqstr);
		#} else {$seqstring=$seqstr;}
		$seqstring=$seqstr;
		
        		#print "$id-$pos - $substring\n";
			my $tmpf = new File::Temp( UNLINK => 1 );
			print $tmpf ">$id-$gval[0]\n$seqstring\n";
			#format=general- supress the header .. The name of outfile should be different
			
			my $myMEGAblast="blastn -task megablast -query $tmpf -db localDATA -perc_identity $identity -num_threads $thread -outfmt '6 qseqid qstart qend sseqid sstart send evalue length frames qcovs' -out $gval[2]_sampleSeq.fa_$id.tmp";

			#my $myLASTZ="lastz $tmpf $infile --chain --output=$gval[2]_sampleSeq.fa_$id.tmp --format=general- --progress --ambiguous=iupac --strand=both --identity=$identity"; 

#xxx qseqid sstrand qlen qstart qend sallseqi sstrand slen sstart send qlen nident pident nident pident
#104726,seq0-seq0,+,1109,0,1109,seq0,+,1681,359,1468,1109/1109,100.0%,1109/1109,100.0%
 	
		0 == system ("$myMEGAblast") or die "Failed to megablast the $id-$gval[0] sequence\n";    
		#system ("cp $gval[2]_sampleSeq.fa_$id.tmp $gval[2]_$gval[3]_$id.txt");
		#Store the allelic coordinates
		polySubs::storeAlles($gval[0], $gval[1], $gval[2], $gval[3], $siz, "$gval[2]_sampleSeq.fa_$id.tmp", $count, $len, $expected, $gval[4], $outDIR);
		#KaKs predictions		
		#polySubs::KaKs("$gval[2]_sampleSeq.fa_$id.tmp" , "$tmpf");
		}
}

#If randome split based analysis selected
elsif ($random) {
    	# Need a long enough sequence for multiple substrings with no overlap
    	if ( $min_length > length $seq ) {
        	warn "Line $., Genome too small:  Must be $min_length, not ", length($seq), "\n";
        	next;
    		}

    	# Save all possible positions for substrings in an array.  This enables us
    	# to remove possibilities after each substring to prevent overlap.
    	my @pos = ( 0 .. length($seq) - 1 - ( $length - 1 ) );
	#foreach (@pos) {print "$_\n";}
    	for ( 1 .. $count ) {
        	my $index = int rand @pos;
        	my $pos   = $pos[$index];

        	# Remove from possible positions
        	my $min = $index - ( $length - 1 );
        	$min = 0 if $min < 0;
        	splice @pos, $min, $length + $index - $min;

        	my $substring = substr $seq, $pos, $length;
		my $countN = () = $substring =~ /N|n/g;
		my $substringLen=length $substring;
		next if $countN >= ($substringLen/3); #If more than one third contain NN or nn ignore string

        	#print "$id-$pos - $substring\n";
		my $tmp_fh = new File::Temp( UNLINK => 1 );
		print $tmp_fh ">$id-$pos\n$substring\n";
		
		#format=general- supress the header

		my $myMEGAblast="blastn -task megablast -query $tmp_fh -db localDATA -perc_identity $identity -num_threads $thread -outfmt '6 qseqid qstart qend sseqid sstart send evalue length frames qcovs' -out $index-sampleSeq.fa_$id.tmp";

		#my $myLASTZ="lastz $tmp_fh $infile --chain --output=$index-sampleSeq.fa_$id.tmp --format=general- --progress --ambiguous=iupac --strand=both --identity=$identity";
        0 == system ("$myMEGAblast") or die "Failed to LastZ the $id-$pos sequence\n";
	}
}

else { print "ERROR: What how this did happed !!"; exit(0);}

	#Now handle all lastz files
	my @old_files = glob "*.tmp";
	my $path = Cwd::cwd();
	my $dirname ="tmpRES";
	mkdir $dirname, 0755;
	my $arc_dir = "$path/$dirname";

	foreach my $old_file (@old_files) {
    		#my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    		#my $new_file = $arc_dir . $short_file_name;
    		move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
	}

	#Lets parse all the outfiles
	my @allFiles = glob "$path/$dirname/*.tmp";
	my @allVal;
	foreach my $fileName (@allFiles) { # If user cant to keep the intermediate use @old_files here

		print "$fileName\n";
		open my $resF, '<', "$fileName";
		my $cnt=0;
		local $/ = "\n";  # read by lines
		while (<$resF>) {
			chomp;
#print "$_ -<>-\n";
			next if /^#/; #Ignore the header if available
			my @values = split('\t', $_);
			# Use conditions here to avoid false positive
			# Ignore the match is less than 50% of the size of query
			if ($values[-1] >= $qcov) { $cnt++;}
			#if (($values[10]-$values[9]) >= ($values[5]-$values[4])/2) { $cnt++; }
			#$values[14] =~ s/\s*\d+%$//;
			#if (($values[14]) >= 50.00) { $cnt++; } #If more than 50 percent covered
			#$cnt++; # if !/^\s+?$/;	
		}
	#print "$id\t$fileName\t$cnt\n";
	push @allVal, $cnt;
	}

	# To ignore if not hits reported
	if ( scalar (@allVal) == 0) { print "Interestingly !! No hits for $id, Skipping !!\n"; remove_tree($dirname); next;}
	my $baseStats=polySubs::process_contig($id, $seq);

	my $prev=0;
	my $ct=0;
	
	foreach (sort {$a <=> $b} @allVal) {
    		if($prev != $_) {
        		if($ct) {
            		print $outfh "$id\t$prev\t$ct\t$count\t$baseStats\t";
			my $per= $ct*100/$count;
			print $outfh "$per\t$len\n";
        		}
        	$prev=$_;
        	$ct=0;
    		}
    	$ct++;
	}
	print $outfh "$id\t$prev\t$ct\t$count\t$baseStats\t" if $ct;

my $per= $ct*100/$count;
print $outfh "$per\t$len\n";

polySubs::deldir("$arc_dir");
}
close $outfh;

#Delete the local blast database files
0 == system("rm -rf localDATA*") or die "Can not delete the files\n";

#Plotting the graph
my %allIds=polySubs::storeInHash("$outDIR/$outfile");
#Print the result in friendly format
polySubs::reformatNplot($outDIR, "$outDIR/$outfile", \%allIds);

#Run the Rscript for Boxplot
#polySubs::Mapper("$outDIR/$outfile.final");
refFile($outDIR, "$outDIR/$outfile.final", "$outDIR/$outfile.final2", $zip);
0 == system("Rscript $RealBin/Rscripts/RBoxplot.R $outDIR/plot.final2 $outDIR/boxPlotted.pdf") or die "Failed to Rscript the boxplot file\n";

0 == system("Rscript $RealBin/Rscripts/marginShinyPlot.R $outDIR/plot.final2") or die "Failed to Rscript the boxplot file\n";
print "\nCongratulation plolyploidism estimation accomplished, Check your $outDIR/$outfile file for result !!!\n";

sub refFile {
my($outL, $infile,$outfile, $zip) = @_;
my $plotF="$outL/plot.final2";

open my $OF, '>', $plotF or die "Could not create: $!\n";;
print $OF "Name\tPloidy\tHits\tCount\tSeq\tGC\tGC_per\tnon_ATGC\tPercentage\tLength\n" if (-z "$plotF");

open FILE, $infile;
   while (<FILE>) {
    	chomp $_;
    	my $line = polySubs::trim ($_);
	next if $. == 1; #Ignore header
	next if ($line =~ /^\s*$/);
    	my @values = split('\t', $_);
    	if ($values[1] <= $zip) { 
		my $newVal = $values[1]+0; #added as it one less in old file
		print $OF "$values[0]\t$newVal\t$values[2]\t$values[3]\t$values[4]\t$values[5]\t$values[6]\t$values[7]\t$values[8]\t$values[9]\n";
		}
	}
   close FILE;
   close $OF;
}

__END__

