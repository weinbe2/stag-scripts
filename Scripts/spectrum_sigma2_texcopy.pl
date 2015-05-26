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
#printf("First cfg: %d. Second cfg: %d. Last cfg: %d. Step size: %d. Expected Measurements: %d. Actual measurements: %d. Missing measurements: %d\n", $first_cfg, $second_cfg, $last_cfg, $step_size, $expected_meas, $actual_meas, $missing_meas);

$sed_commands{'ACT_MEAS'} = int($actual_meas);
$sed_commands{'MISS_MEAS'} = int($missing_meas);
$sed_commands{'FREQ_MEAS'} = int($step_size);
$sed_commands{'THERM'} = int($first_cfg);
$sed_commands{'LAST_MEAS'} = int($last_cfg);

# Cool! Next, we need to go through all the inputed states and grab their t_min, t_max, and blocksize. (Item 2, 4, 5).


foreach my $indiv_state (@states)
{
	# First, make sure it's an approved state.
	if (!($indiv_state ~~ @safe_states))
	{
		die "$indiv_state is not an approved state. Exiting.\n";
	}
	
	# Next, make sure a fitparams exists. 
	if (!(-f "$path/$ensemble"."/spectrum2/fitparams/fitparam.$indiv_state"))
	{
		die "A fitparams file for state $indiv_state does not exist. Exiting.\n";
	}
	
	# We're safe!
	open(my $f_fitparam, "<$path/$ensemble"."/spectrum2/fitparams/fitparam.$indiv_state");
	my @info_line = <$f_fitparam>;
	close($f_fitparam);
	
	my @info_splitter = split(' ', $info_line[0]);
	
	my $tmin = $info_splitter[1];
	my $tmax = $info_splitter[3];
	my $blocksize = $info_splitter[4];
	
	# Verified!
	#printf("State $indiv_state has tmin $tmin, tmax $tmax, and blocksize $blocksize\n");
	$sed_commands{uc($indiv_state)."_TMIN"} = $tmin;
	$sed_commands{uc($indiv_state)."_TMAX"} = $tmax;
	$sed_commands{uc($indiv_state)."_MEASPERBLOCK"} = $blocksize;
}

# Verified!
#print %sed_commands;

# Build a sed command!

my @split_direc = split('/', $path);
my $direc = $split_direc[1];

my $command = "sed -e 's/:ENSEMBLE:/$direc\\/$ensemble/g' "; # < basetex_disc.tex > ../$path/$ensemble/spectrum2/tex/f4plus8l24t48b40m005m080_disc.tex

foreach my $a_key (keys %sed_commands)
{
	$command = $command."-e 's/:".$a_key.":/".$sed_commands{$a_key}."/g' ";
}

$command = $command." < /projectnb/qcd/Staggered/All/spectrum_disconnected_24/tex/basetex_disc.tex > $path/$ensemble/spectrum2/tex/$ensemble"."_disc.tex";

#print $command."\n";

`$command`;

=begin

my $flag = 0;
foreach my $indiv_state (@states)
{
	if ($indiv_state eq $state)
	{
		$flag = 1;
	}
}

if ($flag == 0)
{
	die "State $state is not a predefined state.\n";
}

# A few states get special options. Prepare for this.
my $extra_args = " ";

=end
=cut

=begin 
if ($state eq "sg_stoch") # If we're looking at the sigma...
{
	# We expect one additional argument, the tmin_pref to use for the analytic subtraction of connected.
	
	if (@ARGV != 7)
	{
		die "The state $state expects an extra argument, [analytic subtraction tmin-pref].\n";
	}
	
	$extra_args = $extra_args."$ARGV[6] ";
}

=end
=cut

=begin
# Now that we have all of this together, we can print out a fit interval file.

open(my $filehandle, ">$path/$ensemble"."/spectrum2/fitparams/fitparam.$state");

my $output_string = "$tmin_min $tmin_pref $tmin_max $tmax $blocksize $extra_args";
print $filehandle $output_string;
close($filehandle);
=end
=cut

