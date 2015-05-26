#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# Based on a list of states, this updates the basetex.tex file with tmin, tmax, etc.
# Call style:
# ./spectrum_sigma2_texcopy.pl [ensemble] [states...] 

if (@ARGV < 2)
{
        die "Not enough arguments. Expects \"./spectrum_meas_texcopy.pl [ensemble] [states...]\"\n";
}

my $ensemble = $ARGV[0];

my @states = ();
for (my $i = 1; $i < @ARGV; $i++)
{
	push(@states, $ARGV[$i]);
}

my $run = $ensemble;
# This can change on different systems. 
#my $scrfolder = "/projectnb/qcd/Staggered/Scripts/Spectrum/spec_analysis/";
my %params = name2param($run);

my $path = $params{path};
my $cut = $params{cut};
my $L = $params{L};
my $T = $params{T}; my $nt = $T;
my $beta = $params{beta};
my $m_l;
my $m_h;
my $have_l = 0;
my $have_h = 0;

if ($params{flavor} eq '4plus8') {
    $have_l = 1; $have_h = 1;
    $m_h = $params{mh};
        $m_l = $params{ml};
}
elsif ($params{flavor} eq '4plus4') {
    $have_l = 1; $have_h = 1;
    $m_h = $params{mh};
        $m_l = $params{ml};
}
elsif ($params{flavor} eq '8') {
    $have_l = 1; $have_h = 0;
    $m_l = $params{ml};
}
elsif ($params{flavor} eq '4') {
    $have_l = 1; $have_h = 0;
    $m_l = $params{ml};
}
elsif ($params{flavor} eq '12') {
    $have_l = 1; $have_h = 0;
    $m_l = $params{ml};
}

=begin GHOSTCODE
# Clean up tmax.
if (@ARGV == 6)
{
   $tmax = $ARGV[5]; 
}
else
{
   $tmax = $T/2;
}
=end GHOSTCODE
=cut

# Make sure a directory folder exists! 
if (!(-d "$path/$ensemble/spectrum2"))
{
        print "Error: a spectrum directory for $ensemble does not exist. Exiting.\n";
        exit(200);
}

if (!(-d "$path/$ensemble/spectrum2/stoch"))
{
        print "Error: the raw parsed stochastic data for $ensemble do not exist. Please parse the data before running! Exiting.\n";
        exit(200);
}

# Prepare the sed command!
my %sed_commands = ();

# Next, make sure the state makes sense. These are only 'disconnected' measurements.
my @safe_states = ("sc_stoch", "dc_stoch", "dc_stoch_ppp", "dc_stoch_oscil", "sg_stoch", "sg_stoch_ppp");

# All disconnected measurements happen at the same time. Thus, to get the number of measurements,
# frequency of measurements, etc, we just need to look at one of them. Conveniently, we have one
# sane measurement to grab such things from, too!

# Look at the pion!

open(my $filehandle_thing, "<$path/$ensemble/spectrum2/corr/corr.sc_stoch");
my @pion_file = <$filehandle_thing>;
close($filehandle_thing);

# Get three things---the first config #, the second config #, and the last config #.
# This gives us the number of measurements and the spacing... and a way to get a
# consistency check!

# This can be found from the first, N_t+1'th, and last line!
my $first_line = $pion_file[0];
my $second_line = $pion_file[$nt];
my $last_line = $pion_file[@pion_file-1];

my @tmp_split = split(' ', $first_line); my $first_cfg = $tmp_split[0];
@tmp_split = split(' ', $second_line); my $second_cfg = $tmp_split[0];
@tmp_split = split(' ', $last_line); my $last_cfg = $tmp_split[0];

my $step_size = $second_cfg - $first_cfg;
my $expected_meas = ($last_cfg - $first_cfg)/$step_size + 1;
my $actual_meas = (@pion_file/$nt);
my $missing_meas = $expected_meas - $actual_meas;

# Verified good.
printf("First cfg: %d. Second cfg: %d. Last cfg: %d. Step size: %d. Expected Measurements: %d. Actual measurements: %d. Missing measurements: %d\n", $first_cfg, $second_cfg, $last_cfg, $step_size, $expected_meas, $actual_meas, $missing_meas);


