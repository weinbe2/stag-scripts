#! /usr/bin/perl 

# 2015-04-10 Extended to include calling disconnected functions.
# 2015-04-17 Added the diagonal only flag
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_connected_multi.pl -do {build|analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states]

if (@ARGV < 4) {
   die "Inappropriate number of arguments: expected ./spectrum_connected_multi.pl -do {build|analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] {-noerrors}\n";
}

# Handle the list of arguments.

my @ensemble_list = ();
my @state_list = ();
my @todo_list = ();

my @safe_states = ("sc_stoch", "ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "rvt", "rpv", "dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "dc_stoch_three", "sg_stoch_three", "dc_stoch_zcen", "sg_stoch_zcen");

my @old_disconnected_filter = ("dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "dc_stoch_three", "sg_stoch_three");
my @new_disconnected_filter = ("dc_stoch_zcen", "sg_stoch_zcen");

my $build = 0;
my $analyze = 0;
my $plots = 0;
my $noerrors = 0;
my $diagonly = 0;

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
		case "noerrors" {
			$noerrors = 1;
		}
		case "diagonly" {
			$diagonly = 1;
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

# Build up commands! Both the plots and analyze commands, since it's cheap.

foreach my $params_ref (@ensemble_list)
{
	my %params = %{$params_ref};
	my $direc = $params{ensemble};
	my $path = $params{path};
	
	if ($build == 1)
	{
		if (-f "$path/$direc/spectrum2/fitparams/fitparam.dc_stoch")
		{
			my @split_lines = split(' ',`cat ./$path/$direc/spectrum2/fitparams/fitparam.dc_stoch`);
			$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, $split_lines[4]);"
		}
		else
		{
			# If we have no blocksize to refer to, just assume 1.
			$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, 1);"
		}
		
		$perl_command = $perl_command."./EvanSpectrum/scripts/run_disconnected.pl $path/$direc;";
	}
	
	foreach my $state (@state_list)
	{
		print $state."\n";
		# Is this an old disconnected measurement?
		if ($state ~~ @old_disconnected_filter)
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
		elsif ($state ~~ @new_disconnected_filter) # or a new one?
		{
			if ($noerrors == 0)
			{
				$matlab_command = $matlab_command." study_disconnected_zcen(\'../../$path/$direc\', \'$state\');";
			}
			else # Ignore error analysis
			{
				$matlab_command = $matlab_command." study_disconnected_zcen(\'../../$path/$direc\', \'$state\', 1);";
			}
		}
		else # otherwise, it's a connected business as usual.
		{
			if ($noerrors == 0 && $diagonly == 0)
			{
				$matlab_command = $matlab_command." study_connected_new(\'../../$path/$direc\', \'$state\');";
			}
			else # Ignore error analysis
			{
				my $flag = $noerrors + 2*$diagonly;
				$matlab_command = $matlab_command." study_connected_new(\'../../$path/$direc\', \'$state\', $flag);";
			}
			
			if ($plots == 1)
			{
				$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc $state;";
			}
		}
	}
}

# And run it!
	print $matlab_command."\n";
	print $perl_command."\n";

if ($build == 1 || $analyze == 1)
{
	my $matlab = "matlab -softwareopengl -nosplash -nodisplay -r \"$matlab_command quit;\"";

	system($matlab);
}

if ($build == 1 || $plots == 1)
{

	`$perl_command`;
}

