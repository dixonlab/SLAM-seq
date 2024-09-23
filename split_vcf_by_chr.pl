#!/usr/bin/perl

use strict;

MAIN : {

    my ($vcf_file) = @ARGV;
    if (not defined $vcf_file) {
	die ("Usage: ./split_vcf_by_chr.pl <vcf file>\n");
    }

    my $hash;
    my @header;
    open(FILE,$vcf_file);
    while (my $line = <FILE>) {
	chomp $line;
	if ($line =~ m/^\#/) {
	    push(@header,$line);
	} else {
	    my ($chr, $loc, @rest) = split(/\t/,$line);
	    push(@{$hash->{$chr}},$line);
	}	
    }
    close(FILE);
    
    foreach my $chr (keys %$hash) {
	open(OUTPUT,">$chr.vcf");
	foreach my $line (@header) {
	    print OUTPUT $line . "\n";
	}
	foreach my $line (@{$hash->{$chr}}) {
	    print OUTPUT $line . "\n";
	}
	close(OUTPUT);
    }
    
}
