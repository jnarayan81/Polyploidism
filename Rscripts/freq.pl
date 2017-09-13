use Statistics::R ;
use List::Util qw(sum);

my $R = Statistics::R->new() ;
$R->startR ;

my $fileN=$ARGV[0];
open (my $fh2, $fileN) or die "Could not open file $fileN $!";
my %cHash;
while(<$fh2>) {
	chomp;
	my @arr = split("\t", $_);
	next if $. == 1;
	#if ($arr[7] eq "-") {
	#	my $len=$arr[5]-$arr[4];
	#	push @allPalLen,$len;
	$cHash{$arr[1]}++;
}

foreach my $val (keys %cHash) {
print "$val\t$cHash{$val}\n";
}
