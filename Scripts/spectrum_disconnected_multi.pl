#! /usr/bin/perl 
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_disconnected_multi.pl -do {build|analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] {-blocksize [numbers for blocksizes]} {-noerrors}

if (@ARGV < 4) {
   die "Inappropriate number of arguments: expected ./spectrum_disconnected_multi.pl -do {analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] {-blocksize [numbers for blocksizes]}\n";
}

# Handle the list of arguments.

my @ensemble_list = ();
my @state_list = ();
my @todo_list = ();
my @blocksize_list = ();

my @safe_states = ("dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "dc_stoch_three", "sg_stoch_three");

my $build = 0;
my $analyze = 0;
my $plots = 0;
my $blocksize = 0;
my $noerrors = 0;

my $perl_command = ""; # for plotting
my $matlab_command = "cd('../Scripts/EvanSpectrum/matlab');"; # for analysis.

for (my $i = 0; $i < @ARGV; $i++)
{
	if (substr($ARGV[$i], 0, 1) ne "-")
	{
		die "Expected a flag, got another argument.\n";
	}
	
	substr($ARGV[$i], 0, 1) = "";
	
	switch ($ARGV[$i])
	{
		case "do" {
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				if ($ARGV[$i+1] eq "plots")
				{
					$plots = 1;
				}
				if ($ARGV[$i+1] eq "analyze")
				{
					$analyze = 1;
				}
				if ($ARGV[$i+1] eq "build")
				{
					$build = 1;
				}
				
				$i++;
			}
		}
		case "ensembles" {
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
			
				# Prepare pre-parsed info.
				
				my %params = name2param($ARGV[$i+1]);
				$params{ensemble} = $ARGV[$i+1];
				
				if ($params{flavor} eq '4plus8') {
					$params{have_l} = 4;
					$params{have_h} = 8;
				}
				elsif ($params{flavor} eq '4plus4') {
					$params{have_l} = 4;
					$params{have_h} = 4;
				}
				elsif ($params{flavor} eq '8') {
					$params{have_l} = 8;
					$params{have_h} = 0;
					$params{mh} = 0;
				}
				elsif ($params{flavor} eq '4') {
					$params{have_l} = 4;
					$params{have_h} = 0;
					$params{mh} = 0;
				}
				elsif ($params{flavor} eq '12') {
					$params{have_l} = 12;
					$params{have_h} = 0;
					$params{mh} = 0;
				}
			
				push (@ensemble_list, \%params);
				$i++;
			}
		}
		case "states" {
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				if (!($ARGV[$i+1] ~~ @safe_states)) # Remember, number code at end.
				{
					die "State $ARGV[$i+1] is not an approved state.\n";
				}
				push (@state_list, $ARGV[$i+1]);
				$i++;
			}
		}
		case "blocksize" {
			$blocksize = 1;
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				if ($ARGV[$i+1] eq "~")
				{
					push(@blocksize_list, 1);
				}
				else
				{
					push(@blocksize_list, $ARGV[$i+1]);
				}
				$i++;
			}
		}
		case "noerrors" {
			$noerrors = 1;
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

# Do a quick check on the blocksize list.
if ($blocksize == 1 && (0+@blocksize_list) != (0+@ensemble_list))
{
	die "Blocksize list and ensemble list are not the same length!\n";
}

# Build up commands! Both the plots and analyze commands, since it's cheap.

# Update this for disconnected.

my $ensemble_counter = 0; # For blocksize lookup.

foreach my $params_ref (@ensemble_list)
{
	my %params = %{$params_ref};
	my $direc = $params{ensemble};
	my $path = $params{path};
	
	if ($build == 1)
	{
		if ($blocksize == 1)
		{
			$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, $blocksize_list[$ensemble_counter]);"

		}
		else
		{
			if (-f "$path/$direc/spectrum2/fitparams/fitparam.dc_stoch")
			{
				my @split_lines = split(' ',`cat ./$path/$direc/spectrum2/fitparams/fitparam.dc_stoch`);
				$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, $split_lines[4]);"
			}
			else
			{
				$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, 1);"
			}
		}
		$perl_command = $perl_command."./EvanSpectrum/scripts/run_disconnected.pl $path/$direc;";
	}
	
	$ensemble_counter++;
	
	foreach my $state (@state_list)
	{
		if ($analyze == 1)
		{
			if ($noerrors == 0)
			{
				$matlab_command = $matlab_command." study_disconnected(\'../../$path/$direc\', \'$state\', 6);";
			}
			else # Ignore error analysis
			{
				$matlab_command = $matlab_command." study_disconnected(\'../../$path/$direc\', \'$state\', 6, 1);";
			}
		}
		
		if ($plots == 1)
		{
			if ($state eq "dc_stoch_three")
			{
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc dc_stoch_oscil;";
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc dc_stoch_deriv";
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc dc_stoch_wconst;";
			}
			elsif ($state eq "sg_stoch_three")
			{
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc sg_stoch;";
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc sg_stoch_deriv";
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc sg_stoch_wconst;";
			}
			else
			{
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc $state;";
			}
		}
	}
}

print $matlab_command;

# Only run matlab if we need to.
if ($matlab_command ne "cd('../Scripts/EvanSpectrum/matlab');")
{
	my $matlab = "matlab -softwareopengl -nosplash -nodisplay -r \"$matlab_command quit;\"";
	system($matlab);
}


`$perl_command`;