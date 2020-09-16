#!/usr/bin/perl -w

my $usage = "\nUsage: $0 [-h] [file]\n\n" .
    "  Take normal ms output and print out sample stats.\n" .
    "  Make sure the program sample_stats (from ms pkg) is installed.\n".
    "\n  Outputs:\n".
    "   1st-5th columns: from sample_stats (pi ss D thetaH H)\n".
    "   Then\n     pi_within_1, pi_within_2, ... pi_sithin_n\n" .
    "   Then\n     pi_between_1-2, pi_between_1-3, ... pi_between_1-n, pi_between_2-3, pi_between_2-4, ...\n" .
    "   pi_within_i is pi_within of i-th population.\n" .
    "   pi_between_i-j is pi_between of i-th and j-th populations.\n\n";

# Copyright 2010, Naoki Takebayashi <ffnt@uaf.edu>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
# 02110-1301 USA

# Version 20100407
#   about 50-fold speed up by using more references

use IO::File;
use File::Temp;

our($opt_h);
use Getopt::Std;
getopts('hr') || die $usage;

die $usage if (defined($opt_h));

my $headerDone = 0;
my %conf = ();
my @outCache = ();
my @pwDiffVect = ();

$sampleStats = "./sample_stats";
my @allMSOut = ();

while(<>) {
    push @allMSOut, $_;
    chomp;
    if($.== 1) { # extract info from ms cmd line
	if (/-I\s*(\d+.+)/) {
	    my @arg = split /\s+/, $1;
	    $conf{'numPop'} = shift @arg;

	    my @samplesInSubpops = splice @arg, 0, $conf{'numPop'};

	    $conf{'samples'} = \@samplesInSubpops;
	}
	next;
    }

    next if ($. < 5); # ignore 1st 4 lines
    next if (/^\s*$/);

    if ($_ !~ /^\/\//) {
	push @outCache, $_;
	next;
    }

    # Extract number of total seg sites
    my $numSeg;
    my ($segLine, $posiLine) = splice(@outCache, 0, 2);  # removing segsites: and position line
    if ($segLine =~ /segsites:\s*(\d+)/) {
	$numSeg = $1;
    }

    my $matRef = MakeMat(\@outCache);

    my @pwWithin = AllPWDiffWithin(\%conf, $numSeg, $matRef);
    my @pwBetween = AllPWDiffBetween(\%conf, $numSeg, $matRef);
    push @pwDiffVect, join("\t", @pwWithin, @pwBetween);
# reset
    @outCache = ();
}

# for the final set
my $numSeg;
my ($segLine, $posiLine) = splice(@outCache, 0, 2);  # removing segsites: and position line
if ($segLine =~ /segsites:\s*(\d+)/) {
    $numSeg = $1;
}

my $matRef = MakeMat(\@outCache);

my @pwWithin = AllPWDiffWithin(\%conf, $numSeg, $matRef);
my @pwBetween = AllPWDiffBetween(\%conf, $numSeg, $matRef);
push @pwDiffVect, join("\t", @pwWithin, @pwBetween);

# run sample_stats
my ($tmpOut, $tmpOutFH);
do {$tmpOut = tmpnam()} until $tmpOutFH = IO::File->new($tmpOut, O_RDWR|O_CREAT|O_EXCL);
END {
    if (defined($tmpOut && -e $tmpOut)) {
	unlink($tmpOut) || die "Couldn't unlink $tmpOut : $!\n";
    }
};

open TMP, ">$tmpOut" || die "can't open $tmpOut\n";
print TMP join("\n", @allMSOut);
close (TMP);

my @ssOut = `$sampleStats < $tmpOut |cut -f 2,4,6,8,10`;
chomp(@ssOut);

if (@ssOut != @pwDiffVect) {
    die "error: result len unmatch\n";
}
foreach my $row (0..$#ssOut) {
    print "$ssOut[$row]\t$pwDiffVect[$row]\n";
}
#print join("\n", @pwDiffVect), "\n";
exit;


sub AllPWDiffWithin {
    my ($confRef, $numSeg, $ms1RunRef) = @_;

    my $nPop = ${$confRef}{'numPop'};
    my @nArr = @{${$confRef}{'samples'}};

    if ($numSeg == 0) {
	return (map {0} (1..$nPop));
    }

    my $endIndex = -1;
    my $beginIndex = 0;
    my @results = ();
    foreach my $popi (0..($nPop-1)) {
	$endIndex += $nArr[$popi];

	# count 1 in each nucleotide site and make an array
	my @cntArr = ColumnSum($ms1RunRef, $beginIndex, $endIndex);

	$beginIndex =  $endIndex + 1;

	my @avePWDiffSite = map {$_ * ($nArr[$popi] - $_) } @cntArr;
	my $avePWdiff = SumElements(\@avePWDiffSite)/($nArr[$popi] * ($nArr[$popi]-1)/2);

	push @results, $avePWdiff;
    }
    return @results;
}

sub AllPWDiffBetween {
    my ($confRef, $numSeg, $ms1RunRef) = @_;

    my $nPop = ${$confRef}{'numPop'};
    my @nArr = @{${$confRef}{'samples'}};

    if ($numSeg == 0) {
	return (map {0} (1..($nPop * ($nPop-1)/2)));
    }

    my $totSampleSize = ${$confRef}{'totSamples'};

    my @beginIndex = (0);
    my @endIndex = ($nArr[0]-1);

    foreach my $popI (1..($nPop-1)) {
	$endIndex[$popI] = $endIndex[$popI-1] + $nArr[$popI];
    }
    foreach my $popI (1..($nPop-1)) {
	$beginIndex[$popI] = $endIndex[$popI-1] + 1;
    }

    # count 1 for each pop
    my @cntMat = ();
    foreach my $popI (0..($nPop-1)) {
	# count 1 in each nucleotide site and make an array
	my @cntArr = ColumnSum($ms1RunRef, $beginIndex[$popI], $endIndex[$popI]);
	push @cntMat, \@cntArr;
    }

    my @results = ();
    foreach my $i (0..($nPop-2)) {
	foreach my $j (($i+1)..($nPop-1)) {
	    my $sum = 0;
	    foreach my $siteI (0..($numSeg-1)) {
		my $c1 = $cntMat[$i]->[$siteI];
		my $c2 = $cntMat[$j]->[$siteI];
		$sum += $c1 * ($nArr[$j] - $c2) + ($nArr[$i] - $c1) * $c2;
	    }
	    $sum = $sum / ($nArr[$i] * $nArr[$j]);
	    push @results, $sum;
	}
    }
    return @results;
}

# $ms1RunRef is a ref to a matrix: Each row corresponds to each sample
# sequence, each column represents a segregating site.  Value is
# either 0 or 1.
#
# It returns a ref to a matrix with n row x s columns, where n is
# number of subpopulation and s is the number of segregating sites.
# The value is number of 1's at the segregating site for the
# sub-population.

sub SiteFreqForEachSubPop {
    my ($confRef, $numSeg, $ms1RunRef) = @_;

    my $nPop = ${$confRef}{'numPop'};
    my @nArr = @{${$confRef}{'samples'}};

    if ($numSeg == 0) {
	my @tmp = ();
	return \@tmp;
    }

    # Create an array of indeces which correspond to the begin/end
    # of each subpop.
    my @beginIndex = (0);
    my @endIndex = ($nArr[0]-1);
    foreach my $popI (1..($nPop-1)) {
	$endIndex[$popI] = $endIndex[$popI-1] + $nArr[$popI];
    }
    foreach my $popI (1..($nPop-1)) {
	$beginIndex[$popI] = $endIndex[$popI-1] + 1;
    }

    # count 1 for each pop
    my @cntMat = ();
    foreach my $popI (0..($nPop-1)) {
	# my @thisPopSeq = @seqArr[$beginIndex..$endIndex];
	# count 1 in each nucleotide site and make an array
	my @cntArr = ColumnSum($ms1RunRef,
			       $beginIndex[$popI], $endIndex[$popI]);
	push @cntMat, \@cntArr;
    }

    return \@cntMat;
}

# Receives an array ref, each element is a char string.  Each element
# of array gets splitted by split //, and a (rugged) matrix is
# created. The reference the matrix is returned.
sub MakeMat {
    my $arrRef = shift;
    my $lastIndex = scalar(@{$arrRef})-1;
    my @result = ();

    return (\@result) if ($lastIndex == -1);

    foreach my $i (0..$lastIndex) {
	my @tmpArr = split //, $arrRef->[$i];
	push @result, \@tmpArr;
    }

    return \@result;
}

# Take a reference to matrix, 0-offset indices of rows to begin and
# end summing.  Each column is summed up and matrix is returned.
# It could be screwy with the rugged array.
sub ColumnSum {
    my ($matRef, $rowBeginIndex, $rowEndIndex) = @_; # give 0-offset index

    my $nCol = scalar(@{$matRef->[$rowBeginIndex]});

    my @result = ();
    foreach my $col (0..($nCol-1)) {
	my $thisCntr=0;
	foreach my $row ($rowBeginIndex..$rowEndIndex) {
	    $thisCntr += $matRef->[$row]->[$col];
	}
	push @result, $thisCntr;
    }
    return @result;
}

# a ref to a vector is taken, and sum of all elements are returned.
sub SumElements {
    my $arrRef = shift;
    my $sum = 0;
    foreach my $i (@$arrRef) {
	$sum += $i;
    }
    return $sum;
}
