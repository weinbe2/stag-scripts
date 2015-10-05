#! /usr/bin/perl -w

# 2015-04-10 Extended to include calling disconnected functions.
# 2015-04-17 Added the diagonal only flag
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_connected_multi.pl -do {build|analyze|plots|tex|multitex} -ensembles {f4plus8l24t48b30m005m020} -states [list of states]

if (@ARGV < 4) {
   die "Inappropriate number of arguments: expected ./spectrum_connected_multi.pl -do {build|analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] {-noerrors}\n";
}

# Handle the list of arguments.

my @ensemble_list = ();
my @state_list = ();
my @todo_list = ();

# Load safe states from file.
open(my $states_handle, "<./EvanSpectrum/states.txt");
my @safe_states = <$states_handle>;
close($states_handle);
#my @safe_states = ("sc_stoch", "ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "rvt", "rpv", "dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "dc_stoch_three", "sg_stoch_three", "dc_stoch_zcen", "sg_stoch_zcen");

my @old_disconnected_filter = ("dc_stoch", "dc_stoch_oscil", "dc_stoch_ppp", "sg_stoch", "dc_stoch_three", "sg_stoch_three");
my @new_disconnected_filter = ("dc_stoch_zcen", "sg_stoch_zcen");

# Update this per machine.
my $main_dir = "/projectnb/qcd/Staggered/";

my $build = 0;
my $analyze = 0;
my $plots = 0;
my $tex = 0;
my $multitex = 0;
my $noerrors = 0;
my $diagonly = 0;

# Lots of things trigger a disconnected build.
# Only do it once!
my $did_build_disc = 0;

my $perl_command = ""; # for plotting
my $matlab_command = ""; # for analysis.

$| = 1;

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
				if ($ARGV[$i+1] eq "tex")
				{
					$tex = 1;
				}
				if ($ARGV[$i+1] eq "multitex")
				{
					$multitex = 1;
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
	
	print "$direc $path\n"; 
	
	# Reset the disc counter.
	$did_build_disc = 0;
	
	foreach my $state (@state_list)
	{
		print $state."\n";
		# Is this an old disconnected measurement?
		if ($state ~~ @old_disconnected_filter)
		{
			if ($build == 1 && $did_build_disc == 0)
			{
				# We should add an option for connected states to build just the "sum." file.
				# Disconnected should build finite difference files instead. 
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
				
				$did_build_disc = 1;
			}
		
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
			if ($build == 1 && $did_build_disc == 0)
			{
				# We should add an option for connected states to build just the "sum." file.
				# Disconnected should build finite difference files instead. 
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
				
				$did_build_disc = 1;
			}
		
			if ($analyze == 1)
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
		}
		else # otherwise, it's a connected business as usual.
		{
			# Build sc_stoch if needed.
			if ($state eq 'sc_stoch' && $build == 1 && $did_build_disc == 0)
			{
				# We should add an option for connected states to build just the "sum." file.
				# Disconnected should build finite difference files instead. 
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
				
				$did_build_disc = 1;
			}
			elsif ($build == 1)
			{
				if (-f "$path/$direc/spectrum2/fitparams/fitparam.$state")
				{
					my @split_lines = split(' ',`cat ./$path/$direc/spectrum2/fitparams/fitparam.$state`);
					$matlab_command = $matlab_command." build_correlator(\'../../$path/$direc\', \'$state\', $split_lines[4]);"
				}
				else
				{
					# If we have no blocksize to refer to, just assume 1.
					$matlab_command = $matlab_command." build_correlator_scalar(\'../../$path/$direc\', 6, 1);"
				}
			}
		
			if ($analyze == 1)
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
			}
			
			if ($plots == 1)
			{
				# If a multi-exp fit exists...
				if (-f "$path/$direc/spectrum2/multifits/multifits.$state")
				{
					$perl_command = $perl_command." ./EvanSpectrum/scripts/run_multiconn.pl -do plots -ensemble $path/$direc -state $state;";
				}
				else # business as usual.
				{
					$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc $state;";
				}
			}
		}
	}
	
}

# Prepare a texcopy if we need to!
my $ensemble_state_str = "";
if ($tex == 1 || $multitex == 1)
{
	$ensemble_state_str = "-ensembles ";
	foreach my $params_ref (@ensemble_list)
	{
		my %params = %{$params_ref};
		my $direc = $params{ensemble};
		$ensemble_state_str = $ensemble_state_str.$direc." ";
	}
	
	$ensemble_state_str = $ensemble_state_str."-states ";
	
	foreach my $state (@state_list)
	{
		$ensemble_state_str = $ensemble_state_str.$state." ";
	}
}

if ($tex == 1)
{
	$perl_command = $perl_command." ./spectrum_meas_texcopy.pl ".$ensemble_state_str.";";
}

if ($multitex == 1)
{
	$perl_command = $perl_command." ./spectrum_meas_multitexcopy.pl ".$ensemble_state_str.";";
}


# And run it!
	print $perl_command."\n";

if ($matlab_command ne "")
{
	$matlab_command = "cd('../Scripts/EvanSpectrum/matlab');".$matlab_command."\n";
	
	my $matlab = "matlab -softwareopengl -nosplash -nodisplay -r \"$matlab_command quit;\"";

	system($matlab);
}

if ($build == 1 || $plots == 1 || $tex == 1 || $multitex == 1)
{

	`$perl_command`;
}

