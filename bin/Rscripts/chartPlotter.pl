
use Statistics::R ;
use List::Util qw(sum);

my $R = Statistics::R->new() ;
$R->startR ;

my $fileN=$ARGV[0];
my $gSizeMB=$ARGV[1];
open (my $fh2, $fileN) or die "Could not open file $fileN $!";
my $total;
while(<$fh2>) {
	chomp;
	my @arr = split("\t", $_);
	next if $. == 1;
	#if ($arr[7] eq "-") {
	#	my $len=$arr[5]-$arr[4];
	#	push @allPalLen,$len;
	$cHash{$arr[1]}++;
}
my @keys = keys %cHash;
my $keys=join ',', @keys;
my @values = map {$cHash{$_}} @keys;
my $total=sum(@values);
my $values=join ',', @values;

$R->run(qq`
# Create data for the graph.
x <-  c($total)
labels <-  c($keys)

piepercent<- round(100*x/sum(x), 1)

# Give the chart file a name.
png(file = "covPalindromic.png", width = 880, height = 880, units = "px", pointsize = 12, bg = "white",  res = NA)

# Plot the chart.
pie(x, labels = piepercent, main = "Palindromic Coverage in MB genome",col = rainbow(length(x)))
legend("topright", c($keys), cex = 0.8,
   fill = rainbow(length(x)))

# Save the file.
dev.off()
`);


$R->stopR() ; 

