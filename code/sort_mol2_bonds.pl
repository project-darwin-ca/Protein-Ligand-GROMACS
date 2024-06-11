#!/usr/bin/perl






unless (scalar(@ARGV)==2)
{
      die "Usage: perl sort_mol2_bonds.pl input.mol2 output.mol2\n";
}


my $input = $ARGV[0];
my $output = $ARGV[1];

open(IN, "<$input") || die "Cannot open $input: $!\n";
my @in = <IN>;
close(IN);

# test for header lines that some scripts produce
unless($in[0] =~ /TRIPOS/)
{
    die "Nonstandard header found: $in[0]. Please delete header lines until the TRIPOS molecule definition.\n"; 
}

open(OUT, ">$output") || die "Cannot open $output: $!\n";

# get number of atoms and number of bonds from mol2 file
my @tmp = split(" ", $in[2]);
my $natom = $tmp[0];
my $nbond = $tmp[1];

# check
print "Found $natom atoms in the molecule, with $nbond bonds.\n";

# print out everything prior to bonds section
my $i=0;
while (!($in[$i] =~ /BOND/))
{
    print OUT $in[$i];
    $i++;
}

# print the bond section header line to output
my $bondfmt = $%6d































