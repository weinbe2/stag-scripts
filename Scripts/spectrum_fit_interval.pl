#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# This creates a fitparams file in a given ensemble.
# Call style:
# ./set_fit_interval.pl [ensemble] [state] {[tmin-min] [tmin-pref] [tmin-max] [tmax] [blocksize]}
# This implies the fitter will fit to a range of tmins, tmin-min to tmin-max, to a fixed tmax (default half-Nt)
# If a function calls for a fixed tmin, it'll use tmin-pref.
# If only passed with two arguments, it returns the existing file.

# Update 2015-02-23: be able to only update some fields.

# Check and see if there are any - flags. If so, process, if not, go to default behavior.
my $check_flags = 0;

foreach my $arg (@ARGV)
{
   if (substr($arg, 0, 1) eq "-")
   {
      $check_flags = 1;
   }
}

my $ensemble = $ARGV[0];
my $state = $ARGV[1];


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
if (!(-d "$path/$ensemble"))
{
        print "Error: a directory for $ensemble does not exist. Exiting.\n";
        exit(200);
}



# Make necessary fodlers!
if (!(-d "$path/$ensemble"."/spectrum2"))
{
        my $command = "mkdir $path/$ensemble"."/spectrum2";
        `$command`;
}

# This is where we put fit intervals!
if (!(-d "$path/$ensemble"."/spectrum2/fitparams"))
{
        my $command = "mkdir $path/$ensemble"."/spectrum2/fitparams";
        `$command`;
}

my @states = ();
# Load safe states from file.
open(my $states_handle, "<./EvanSpectrum/states.txt");
@states = <$states_handle>;
close($states_handle);
for (my $i=0;$i<@states;$i++)
{
	# I don't know why I have to do this.
	my @tmp_arr = split(' ', $states[$i]);
	$states[$i] = $tmp_arr[0];
	my $tmp2 = length($states[$i]);
}

# 2015-02-09 Do a check for a wf_state
my $check_state = "";
my @states_split = split("_", $state);
my $split_count = @states_split;
if ($split_count > 2 && substr($states_split[$split_count-2], 0, 2) eq "wt") # It's a wf state.
{
   splice(@states_split, $split_count-1, 1);
   splice(@states_split, $split_count-2, 1);
   $check_state = join("_", @states_split);
   #print "wf state\n";
}
else
{
   $check_state = $state;
}

my $flag = 0;
foreach my $indiv_state (@states)
{
	if ($indiv_state eq $check_state)
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

# Now that we have all of this together, we can print out a fit interval file.

if ($check_flags == 0) # Old behavior.
{

if (@ARGV != 7 && @ARGV != 2)
{
        die "Not enough arguments. Expects \"./set_fit_interval.pl [ensemble] [state] [tmin-min] [tmin-pref] [tmin-max] [tmax] [blocksize]\"\n";
}


if (@ARGV == 7)
{
my $tmin_min = $ARGV[2];
my $tmin_pref = $ARGV[3];
my $tmin_max = $ARGV[4];
my $tmax = $ARGV[5];
my $blocksize = $ARGV[6];


open(my $filehandle, ">$path/$ensemble"."/spectrum2/fitparams/fitparam.$state");

my $output_string = "$tmin_min $tmin_pref $tmin_max $tmax $blocksize $extra_args";
print $filehandle $output_string;
close($filehandle);
}
else
{
   if (-f "$path/$ensemble"."/spectrum2/fitparams/fitparam.$state")
   {
      open (my $filehandle, "<$path/$ensemble"."/spectrum2/fitparams/fitparam.$state");
      my @temp_contents = <$filehandle>;
      print "./spectrum_fit_interval.pl $ensemble $state $temp_contents[0]\n";
      close ($filehandle);
   }
   else
   {
     die "A fitparams file for $ensemble state $state does not exist!\n";
   }

}
}
else # New flag checking behavior.
{

my $tmin_min = -1;
my $tmin_pref = -1;
my $tmin_max = -1;
my $tmax = -1;
my $blocksize = -1;

# See if a fitparam file exists.
if (-f "$path/$ensemble"."/spectrum2/fitparams/fitparam.$state")
{
   open (my $filehandle, "<$path/$ensemble"."/spectrum2/fitparams/fitparam.$state");
   my @temp_contents = <$filehandle>;
   close ($filehandle);
   my @splitvals = split(' ', $temp_contents[0]);
   $tmin_min = $splitvals[0];
   $tmin_pref = $splitvals[1];
   $tmin_max = $splitvals[2];
   $tmax = $splitvals[3];
   $blocksize = $splitvals[4];
}


for (my $i = 2; $i < @ARGV; $i++)
{
	if (substr($ARGV[$i], 0, 1) ne "-")
	{
		die "Expected a flag, got another argument.\n";
	}
	
	substr($ARGV[$i], 0, 1) = "";
	
	switch ($ARGV[$i])
	{
		case "tmin-min" {
			$tmin_min = $ARGV[$i+1];
			$i++;
		}
		case "tmin" {
			$tmin_pref = $ARGV[$i+1];
			$i++;
		}
		case "tmin-max" {
			$tmin_max = $ARGV[$i+1];
			$i++;
		}
		case "tmax" {
			$tmax = $ARGV[$i+1];
			$i++;
		}
		case "blocksize" {
			$blocksize = $ARGV[$i+1];
			$i++;
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

# Make sure nothing is -1.
if ($tmin_min == -1 || $tmin_pref == -1 || $tmin_max == -1 || $tmax == -1 || $blocksize == -1)
{
   die "Not all parameters are set!\n";
}

open(my $filehandle, ">$path/$ensemble"."/spectrum2/fitparams/fitparam.$state");

my $output_string = "$tmin_min $tmin_pref $tmin_max $tmax $blocksize ";
print $filehandle $output_string;
close($filehandle);


}



