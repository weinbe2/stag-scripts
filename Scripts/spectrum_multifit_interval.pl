#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# This creates a multifitparams file in a given ensemble.
# Call style:
# ./set_multifit_interval.pl [ensemble] [state] {tmin tmax blocksize}
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

# Next, make sure the state makes sense.
my @states = ("ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "sc_stoch", "dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "sg_111", "sg_211", "rvt", "rpv", "pll", "plc", "pcc", "all", "acc", "rll", "rcc");

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

if (@ARGV != 5 && @ARGV != 2)
{
        die "Not enough arguments. Expects \"./set_multifit_interval.pl [ensemble] [state] [tmin] [tmax] [blocksize]\"\n";
}


if (@ARGV == 5)
{
my $tmin_pref = $ARGV[2];
my $tmax = $ARGV[3];
my $blocksize = $ARGV[4];


open(my $filehandle, ">$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state");

my $output_string = "$tmin_pref $tmax $blocksize ";
print $filehandle $output_string;
close($filehandle);
}
else
{
   if (-f "$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state")
   {
      open (my $filehandle, "<$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state");
      my @temp_contents = <$filehandle>;
      print "./spectrum_multifit_interval.pl $ensemble $state $temp_contents[0]\n";
      close ($filehandle);
   }
   else
   {
     die "A multifitparams file for $ensemble state $state does not exist!\n";
   }

}
}
else # New flag checking behavior.
{

my $tmin_pref = -1;
my $tmax = -1;
my $blocksize = -1;

# See if a multifitparam file exists.
if (-f "$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state")
{
   open (my $filehandle, "<$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state");
   my @temp_contents = <$filehandle>;
   close ($filehandle);
   my @splitvals = split(' ', $temp_contents[0]);
   $tmin_pref = $splitvals[0];
   $tmax = $splitvals[1];
   $blocksize = $splitvals[2];
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
		case "tmin" {
			$tmin_pref = $ARGV[$i+1];
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
if ($tmin_pref == -1 || $tmax == -1 || $blocksize == -1)
{
   die "Not all parameters are set!\n";
}

open(my $filehandle, ">$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$state");

my $output_string = "$tmin_pref $tmax $blocksize ";
print $filehandle $output_string;
close($filehandle);


}



