#!/usr/bin/env perl
# Charles R. Kissinger, MIT license
#
# Print atom records from Protein Data Bank (PDB) files that contain fields
# that match a pattern.
# 
# Usage: pdbgrep.pl [options] [file ...]
# Options:
#    -a range        - atom number (ordinal)
#    -b range        - temperature factor
#    -c expression   - chain id
#    -e expression   - element type
#    -h              - heterogen atom
#    -i expression   - insertion code
#    -l expression   - alternate location indicator
#    -n expression   - atom name
#    -o range        - occupancy
#    -p              - pass non-ATOM or HETATM records through to the output
#    -q expression   - atomic charge
#    -r range        - residue number
#    -s expression   - segment id
#    -t expression   - residue type
#    -v              - left-shifted atom name (heavy atom or hydrogen)
#    -x range        - x coordinate
#    -y range        - y coordinate
#    -z range        - z coordinate
#
# The uppercase version of each option prints records that DON'T match
# the range or expression.
#
# range:        a single value, or low:high, or low: (no upper bound)
#               or :high (no lower bound) 
# expression:   any regular expression that matches the field completely
#               (leading and trailing spaces in fields are discarded)

use strict;
use warnings;
use Getopt::Std;

our ($opt_a, $opt_A, $opt_b, $opt_B, $opt_c, $opt_C, $opt_e, $opt_E, $opt_h);
our ($opt_H, $opt_i, $opt_I, $opt_l, $opt_L, $opt_n, $opt_N, $opt_o, $opt_O);
our ($opt_p, $opt_P, $opt_q, $opt_Q, $opt_r, $opt_R, $opt_s, $opt_S, $opt_t);
our ($opt_T, $opt_v, $opt_V, $opt_x, $opt_X, $opt_y, $opt_Y, $opt_z, $opt_Z);
my (@arange, @Arange, @brange, @Brange, @orange, @Orange, @rrange, @Rrange);
my (@xrange, @Xrange, @yrange, @Yrange, @zrange, @Zrange);

getopts('a:A:b:B:c:C:e:E:i:I:l:L:n:N:o:O:q:Q:r:R:s:S:t:T:x:X:y:Y:z:Z:hHpvV')
    or exit 1;

if ($opt_a) { @arange = rangeArg($opt_a); }
if ($opt_A) { @Arange = rangeArg($opt_A); }
if ($opt_b) { @brange = rangeArg($opt_b); }
if ($opt_B) { @Brange = rangeArg($opt_B); }
if ($opt_o) { @orange = rangeArg($opt_o); }
if ($opt_O) { @Orange = rangeArg($opt_O); }
if ($opt_r) { @rrange = rangeArg($opt_r); }
if ($opt_R) { @Rrange = rangeArg($opt_R); }
if ($opt_x) { @xrange = rangeArg($opt_x); }
if ($opt_X) { @Xrange = rangeArg($opt_X); }
if ($opt_y) { @yrange = rangeArg($opt_y); }
if ($opt_Y) { @Yrange = rangeArg($opt_Y); }
if ($opt_z) { @zrange = rangeArg($opt_z); }
if ($opt_Z) { @Zrange = rangeArg($opt_Z); }

while (<>) {
	if (/^HETATM|^ATOM  /) {
		if (isMatchingRecord($_)) { print; }
	}
	elsif ($opt_p) { print; }
}

sub isMatchingRecord {
  # atom number
  if ($opt_a && (! inRange(substr($_[0], 6, 5), @arange))) { return 0; }
  if ($opt_A && (inRange(substr($_[0], 6, 5), @Arange))) { return 0; }
  # B-value
  if ($opt_b && (! inRange(substr($_[0], 60, 6), @brange))) { return 0; }
  if ($opt_B && (inRange(substr($_[0], 60, 6), @Brange))) { return 0; }
  # chain id
  if ($opt_c && (substr($_[0], 21, 1) !~ /[$opt_c]/o)) { return 0; }
  if ($opt_C && (substr($_[0], 21, 1) =~ /[$opt_C]/o)) { return 0; }
  # element type
  if ($opt_e && (pdbElementType(@_) !~ /^($opt_e)$/o)) { return 0; }
  if ($opt_E && (PdbElementType(@_) =~ /^($opt_E)$/o)) { return 0; }
  # heterogen atom
  if ($opt_h && ($_[0] !~ /^HETATM/)) { return 0; }
  if ($opt_H && ($_[0] =~ /^HETATM/)) { return 0; }
  # insertion code
  if ($opt_i && (substr($_[0], 26, 1) !~ /[$opt_i]/o)) { return 0; }
  if ($opt_I && (substr($_[0], 26, 1) =~ /[$opt_I]/o)) { return 0; }
  # alternate location indicator
  if ($opt_l && (substr($_[0], 16, 1) !~ /[$opt_l]/o)) { return 0; }
  if ($opt_L && (substr($_[0], 16, 1) =~ /[$opt_L]/o)) { return 0; }
  # atom name
  if ($opt_n && (pdbAtomName(@_) !~ /^($opt_n)$/o)) { return 0; }
  if ($opt_N && (pdbAtomName(@_) =~ /^($opt_N)$/o)) { return 0; }
  # occupancy
  if ($opt_o && (! inRange(substr($_[0], 54, 6), @orange))) { return 0; }
  if ($opt_O && (inRange(substr($_[0], 54, 6), @Orange))) { return 0; }
  # residue number
  if ($opt_r && (! inRange(substr($_[0], 22, 4), @rrange))) { return 0; }
  if ($opt_R && (inRange(substr($_[0], 22, 4), @Rrange))) { return 0; }
  # segment id (or old entry code)
  if ($opt_s && (pdbSegID(@_) !~ /^($opt_s)$/o)) { return 0; }
  if ($opt_S && (pdbSegID(@_) =~ /^($opt_S)$/o)) { return 0; }
  # heavy or hydrogen (left-shifted) atom name
  if ($opt_v && (substr($_[0], 12, 1) eq ' ')) { return 0; }
  if ($opt_V && (substr($_[0], 12, 1) ne ' ')) { return 0; }
  # residue type
  if ($opt_t && (pdbResName(@_) !~ /^($opt_t)$/o)) { return 0; }
  if ($opt_T && (pdbResName(@_) =~ /^($opt_T)$/o)) { return 0; }
  # atomic charge
  if ($opt_q && (pdbAtomicCharge(@_) !~ /^($opt_q)$/o)) { return 0; }
  if ($opt_Q && (PdbAtomicCharge(@_) =~ /^($opt_Q)$/o)) { return 0; }
  # x, y or z coordinate
  if ($opt_x && (! inRange(substr($_[0], 30, 8), @xrange))) { return 0; }
  if ($opt_X && (inRange(substr($_[0], 30, 8), @Xrange))) { return 0; }
  if ($opt_y && (! inRange(substr($_[0], 38, 8), @yrange))) { return 0; }
  if ($opt_Y && (inRange(substr($_[0], 38, 8), @Yrange))) { return 0; }
  if ($opt_z && (! inRange(substr($_[0], 46, 8), @zrange))) { return 0; }
  if ($opt_Z && (inRange(substr($_[0], 46, 8), @Zrange))) { return 0; }
  return 1;
}

sub rangeArg {
	(my $arg) = @_;
    # the following is not a complete test
	if ($arg !~ /^[:.0-9\-]+$/)
		{ die "pdbgrep: invalid range argument-> $arg\n"; }
	my @range = split(/:/, $arg, 2);
	if ($arg =~ /:/) {
		if (length($range[0]) == 0) { $range[0] = -999999999; }
		if (length($range[1]) == 0) { $range[1] = 999999999; }
	}
	else { $range[1] = $range[0]; }
	return @range;
}

sub inRange {
	(my $test, my $low, my $high) = @_;
	return ($test >= $low && $test <= $high);
}

sub pdbAtomName {
	(my $line) = @_;
	my $name = substr($line, 12, 4);
	$name =~ s/^\s+|\s+$//g;
	return $name;
}

sub pdbSegID {
	(my $line) = @_;
	my $name = substr($line, 72, 4);
	$name =~ s/^\s+|\s+$//g;
	return $name;
}

sub pdbElementType {
  (my $line) = @_;
  my $name = substr($line, 76, 2);
  $name =~ s/^\s+|\s+$//g;
  return $name;
}

sub pdbAtomicCharge {
  (my $line) = @_;
  my $name = substr($line, 78, 2);
  $name =~ s/^\s+|\s+$//g;
  return $name;
}

sub pdbResName {
	(my $line) = @_;
	my $name = substr($line, 17, 3);
	$name =~ s/^\s+|\s+$//g;
	return $name;
}

