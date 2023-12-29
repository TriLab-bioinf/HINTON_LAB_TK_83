#!/usr/local/bin/perl
#
my $usage = "$0 -d <datafile_file> -dri data_return_index [-dqi <data query index [0]> ; -i <query_id_file> ; -ii query index [0] | STDIN ; -h help]\n\n";

my %arg = @ARGV;
die $usage if $arg{-h} || !exists($arg{-d}) || !exists($arg{-dri});
my $DQINDEX = $arg{-dqi} || 0;
my $DRINDEX = $arg{-dri}; 
my $QINDEX = $arg{-ii} || 0;

# Load reference ids to query in the input file
open (DATA, "<$arg{-d}") || die "I cannot find file $arg{-f}: $!\n\n";
my (%h, %counter);
while(<DATA>){
	chomp;
	my @x=split /\t/;
	$h{ $x[$DQINDEX] } = $x[$DRINDEX];
	$counter{ $x[$DQINDEX] }++; ## check for duplicated indexes
	# warn "WARNING, index $x[$DQINDEX] is duplicated $counter{$x[$DQINDEX]} times\n" if $counter{ $x[$DQINDEX] } > 1; 
}
close DATA;

my ($QUERY);

# Do the query
if ($arg{-i}){
	open( $QUERY , "$arg{-i}") || die "I cannot find file $arg{-i}: $!\n\n";
} else {
	$QUERY = *STDIN; 
}

while(<$QUERY>){
	chomp;
	my @x = split /\t/;
        if( $h{ $x[$QINDEX] }){
		warn "WARNING, index $x[$QINDEX] is duplicated $counter{$x[$QINDEX]} times\n" if $counter{ $x[$QINDEX] } > 1;
		print "$_\t$h{$x[$QINDEX]}\n";
	}
	else {
		print "$_\t#N/A\n";
	}
	
}

close $QUERY if $arg{-i};

