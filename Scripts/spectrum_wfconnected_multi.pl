#! /usr/bin/perl 
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_connected_multi.pl -do {analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] -wf [list of wt[time]]

if (@ARGV < 4) {
   die "Inappropriate number of arguments: expected ./spectrum_connected_multi.pl -do {analyze|plots} -ensembles {f4plus8l24t48b30m005m020} -states [list of states] -wf [list of wf times] {-noerrors}\n";
}

# Handle the list of arguments.

my @ensemble_list = ();
my @state_list = ();
my @wf_list = ();
my @todo_list = ();

my @safe_states = ("sc_stoch", "ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "rvt", "rpv");

my $analyze = 0;
my $plots = 0;
my $noerrors = 0;
my $wftime = 0;

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
                case "wf" {
			$wftime = 1;
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
                        {
                                push (@wf_list, $ARGV[$i+1]);
                                $i++;
                        }
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

if ($wftime == 0)
{
	push (@wf_list, "wt0_00");
}

# Build up commands! Both the plots and analyze commands, since it's cheap.

foreach my $params_ref (@ensemble_list)
{
	my %params = %{$params_ref};
	my $direc = $params{ensemble};
	my $path = $params{path};
	
	foreach my $state (@state_list)
	{
		foreach my $wf_entry (@wf_list)
		{
			my $wf_string = "";
			if ($wf_entry ne "wt0_00")
			{
				$wf_string = "_".$wf_entry;
			}

			if ($noerrors == 0)
			{
				$matlab_command = $matlab_command." study_connected(\'../../$path/$direc\', \'$state"."$wf_string\');";
			}
			else # Ignore error analysis
			{
				$matlab_command = $matlab_command." study_connected(\'../../$path/$direc\', \'$state"."$wf_string\', 1);";
			}
		
			$perl_command = $perl_command." ./EvanSpectrum/scripts/run_connected.pl $path/$direc $state"."$wf_string;";
		}
	}
}

# And run it!

if ($analyze == 1)
{
	my $matlab = "matlab -softwareopengl -nosplash -nodisplay -r \"$matlab_command quit;\"";
	system($matlab);
}

if ($plots == 1)
{
	`$perl_command`;
}

