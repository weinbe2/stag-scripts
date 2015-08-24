#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# Based on a list of states, this updates the basetex_multi.tex file with tmin, tmax, etc.
# Call style:
# ./spectrum_meas_texcopy.pl -ensembles [list of ensembles] {-ensemblelist [relative path to ensemble list]} -states [list of state]

if (@ARGV < 4)
{
        die "Not enough arguments. Expects \"./spectrum_meas_texcopy.pl -ensembles [list of ensembles] {-ensemblelist [relative path to ensemble list]} -states [list of state]\"\n";
}

# Handle the list of arguments.
my @ensemble_list = ();
my @state_list = ();

# A variable to hold file handles in.
my $file_handle;

# A list of safe states. These are only 'connected' measurements.
my @safe_states = ("ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "rvt", "rpv", "sc_stoch");

# Loop over all arguments.
for (my $i = 0; $i < @ARGV; $i++)
{
	if (substr($ARGV[$i], 0, 1) ne "-")
	{
		die "Expected a flag, got another argument.\n";
	}
	
	substr($ARGV[$i], 0, 1) = "";
	
	switch ($ARGV[$i])
	{
		case "ensembles" {
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				# Save parsed information on all of the ensembles.
				
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
		case "ensemblelist" {
			if (-e $ARGV[$i+1])
			{
				open($file_handle, "<".$ARGV[$i+1]);
				my @ensemble_file = <$file_handle>;
				close($file_handle);
				
				foreach my $indiv_ensemble (@ensemble_file)
				{
					chomp($indiv_ensemble);
					#print $indiv_ensemble."\n";
					
					my %params = name2param($indiv_ensemble);
					$params{ensemble} = $indiv_ensemble;
					
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
				}
			}
			else
			{
				die "Ensemble list $ARGV[$i+1] does not exist!\n";
			}
			
			$i++;
		}
		case "states" {
			# Loop over all subsequent entries until we hit the end or another flag.
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				if (!($ARGV[$i+1] ~~ @safe_states)) # Remember, number code at end.
				{
					die "State $ARGV[$i+1] is not an approved state type. Exiting.\n";
				}
				push (@state_list, $ARGV[$i+1]);
				$i++;
			}
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}


# Loop over all ensembles!
foreach my $params_ref (@ensemble_list)
{
	my %params = %{$params_ref};
	my $ensemble = $params{ensemble};
	my $path = $params{path};
	my $nt = $params{T};
	
	# Make sure a directory folder exists! 
	if (!(-d "$path/$ensemble/spectrum2"))
	{
		print "Error: a spectrum directory for $ensemble does not exist. Skipping to next state.\n";
		next;
	}

	if (!(-d "$path/$ensemble/spectrum2/corr"))
	{
		print "Error: the raw parsed correlators for $ensemble do not exist. Please parse the data before running! Skipping to next state.\n";
		next; 
	}

	# Prepare the sed command!
	my %sed_commands = ();

	# All connected measurements happen at the same time. Thus, to get the number of measurements,
	# frequency of measurements, etc, we just need to look at one of them.

	# Look at the pion!

	open($file_handle, "<$path/$ensemble/spectrum2/corr/corr.ps");
	my @pion_file = <$file_handle>;
	close($file_handle);

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

	$sed_commands{'ACT_MEAS'} = $actual_meas;
	$sed_commands{'MISS_MEAS'} = $missing_meas;
	$sed_commands{'FREQ_MEAS'} = $step_size;
	$sed_commands{'THERM'} = $first_cfg;
	$sed_commands{'LAST_MEAS'} = $last_cfg;
	
	# We also need to grab the stochastic ones, if they exist.
	if (-e "$path/$ensemble/spectrum2/corr/corr.sc_stoch")
	{
		open($file_handle, "<$path/$ensemble/spectrum2/corr/corr.sc_stoch");
		@pion_file = <$file_handle>;
		close($file_handle);

		# Get three things---the first config #, the second config #, and the last config #.
		# This gives us the number of measurements and the spacing... and a way to get a
		# consistency check!

		# This can be found from the first, N_t+1'th, and last line!
		$first_line = $pion_file[0];
		$second_line = $pion_file[$nt];
		$last_line = $pion_file[@pion_file-1];

		@tmp_split = split(' ', $first_line); $first_cfg = $tmp_split[0];
		@tmp_split = split(' ', $second_line); $second_cfg = $tmp_split[0];
		@tmp_split = split(' ', $last_line); $last_cfg = $tmp_split[0];

		$step_size = $second_cfg - $first_cfg;
		$expected_meas = ($last_cfg - $first_cfg)/$step_size + 1;
		$actual_meas = (@pion_file/$nt);
		$missing_meas = $expected_meas - $actual_meas;

		# Verified good.
		#printf("First cfg: %d. Second cfg: %d. Last cfg: %d. Step size: %d. Expected Measurements: %d. Actual measurements: %d. Missing measurements: %d\n", $first_cfg, $second_cfg, $last_cfg, $step_size, $expected_meas, $actual_meas, $missing_meas);

		$sed_commands{'ACT_MEAS_STOCH'} = int($actual_meas);
		$sed_commands{'MISS_MEAS_STOCH'} = int($missing_meas);
		$sed_commands{'FREQ_MEAS_STOCH'} = int($step_size);
		$sed_commands{'LAST_MEAS_STOCH'} = int($last_cfg);
	}
	
	# Assume things like sc_stoch, nu, rij, and the taste states don't exist, and only
	# fix if they do!
	$sed_commands{'STOCH_YES_BEGIN'} = "\\\\begin{comment}";
	$sed_commands{'STOCH_YES_END'} = "\\\\end{comment}";
	$sed_commands{'NU_YES_BEGIN'} = "\\\\begin{comment}";
	$sed_commands{'NU_YES_END'} = "\\\\end{comment}";
	$sed_commands{'RIJ_YES_BEGIN'} = "\\\\begin{comment}";
	$sed_commands{'RIJ_YES_END'} = "\\\\end{comment}";
	$sed_commands{'TASTE_YES_BEGIN'} = "\\\\begin{comment}";
	$sed_commands{'TASTE_YES_END'} = "\\\\end{comment}";
	$sed_commands{'STOCH_COMM'} = "%";
	$sed_commands{'NU_COMM'} = "%";
	$sed_commands{'RIJ_COMM'} = "%";
	$sed_commands{'TASTE_COMM'} = "%";
	
	# Cool! Next, we need to go through all the inputed states and grab their t_min, t_max, and blocksize. (Item 2, 4, 5).

	my $milc_flag = 0;
	
	# Prepare to set notify strings.
	
	# Save pion mass from each source.
	my $cps = 0;
	my $eps = 0;
	my $cps2 = 0;
	my $eps2 = 0;
	
	# Save a0 mass from each source.
	my $csc = 0;
	my $esc = 0;
	my $cscs = 0;
	my $escs = 0;

	foreach my $indiv_state (@state_list)
	{
		# First, make sure a multifitparams exists for this state. Otherwise, print
		# a warning, and skip on to the next state!
		if (-f "$path/$ensemble"."/spectrum2/multifitparams/multifitparam.$indiv_state" && -f "$path/$ensemble/spectrum2/corr/corr.$indiv_state" && -f "$path/$ensemble/spectrum2/multicentral/multicentral.$indiv_state")
		{
			# Next, check the milc flag.
			if ($indiv_state eq "rvt" || $indiv_state eq "rpv")
			{
				$milc_flag = 1;
			}

			# Check for the new rho flag.
			if ($indiv_state eq "rij")
			{
				$sed_commands{'RIJ_YES_BEGIN'} = "";
				$sed_commands{'RIJ_YES_END'} = "";
				$sed_commands{'RIJ_COMM'} = "";
			}
			
			if ($indiv_state eq "sc_stoch")
			{
				$sed_commands{'STOCH_YES_BEGIN'} = "";
				$sed_commands{'STOCH_YES_END'} = "";
				$sed_commands{'STOCH_COMM'} = "";
			}
			
			if ($indiv_state eq "nu")
			{
				$sed_commands{'NU_YES_BEGIN'} = "";
				$sed_commands{'NU_YES_END'} = "";
				$sed_commands{'NU_COMM'} = "";
			}
			
			if ($indiv_state eq "ij" || $indiv_state eq "i5")
			{
				$sed_commands{'TASTE_YES_BEGIN'} = "";
				$sed_commands{'TASTE_YES_END'} = "";
				$sed_commands{'TASTE_COMM'} = "";
			}
			
			# Next, make sure a fitparams exists. 
			
			
			# We're safe!
			open($file_handle, "<$path/$ensemble/spectrum2/multifitparams/multifitparam.$indiv_state");
			my @info_line = <$file_handle>;
			close($file_handle);
			
			my @info_splitter = split(' ', $info_line[0]);
			
			my $tmin = $info_splitter[0];
			my $tmax = $info_splitter[1];
			my $blocksize = $info_splitter[2];
			
			open($file_handle, "<$path/$ensemble/spectrum2/multicentral/multicentral.$indiv_state");
			@info_line = <$file_handle>;
			close($file_handle);
			
			@info_splitter = split(' ', $info_line[0]);
			my @info_splitter2 = split(' ', $info_line[1]);
			
			my $ndir = $info_splitter[0];
			my $cdir = $info_splitter[1];
			my $edir = $info_splitter[2];
			
			my $nosc = $info_splitter2[0];
			my $cosc = $info_splitter2[1];
			my $eosc = $info_splitter2[2];
			
			if ($indiv_state eq "ps")
			{
				$cps = $cdir;
				$eps = $edir;
			}
			if ($indiv_state eq "ps2")
			{
				$cps2 = $cdir;
				$eps2 = $edir;
			}
			if ($indiv_state eq "sc")
			{
				$csc = $cosc;
				$esc = $eosc;
			}
			if ($indiv_state eq "sc_stoch")
			{
				$cscs = $cdir;
				$escs = $edir;
			}
			
			# Verified!
			#printf("State $indiv_state has tmin $tmin, tmax $tmax, and blocksize $blocksize\n");
			$sed_commands{uc($indiv_state)."_TMIN"} = $tmin;
			$sed_commands{uc($indiv_state)."_TMAX"} = $tmax;
			$sed_commands{uc($indiv_state)."_MEASPERBLOCK"} = $blocksize;
			$sed_commands{uc($indiv_state)."_NDIR"} = $ndir;
			$sed_commands{uc($indiv_state)."_NOSC"} = $nosc;
			$sed_commands{uc($indiv_state)."_DIR_MASS"} = parenFormWiki($cdir,$edir);
			if ($cosc > 0)
			{
				$sed_commands{uc($indiv_state)."_OSC_MASS"} = parenFormWiki($cosc,$eosc);
			}
		}
		else
		{
			print "A fitparams file for state $indiv_state in ensemble $ensemble does not exist. Skipping to the next state.\n";
		}
	}
	
	# Prepare the notify string.
	my $notify_string = "";
	
	# Compare pion masses.
	my $ps_dev = 0;
	if ($cps > $cps2)
	{
		$ps_dev = ceil(($cps-$cps2)/($eps+$eps2));
	}
	else
	{
		$ps_dev = ceil(($cps2-$cps)/($eps+$eps2));
	}
	
	if ($ps_dev > 2)
	{
		$notify_string = $notify_string."The pion masses agree at only \\alert{$ps_dev sigma}, which is large.";
	}
	
	# Compare a0 masses.
	my $a0_dev = 0;
	if ($cscs != 0)
	{
		if ($csc > $cscs)
		{
			$a0_dev = ceil(($csc-$cscs)/($esc+$escs));
		}
		else
		{
			$a0_dev = ceil(($cscs-$csc)/($escs+$esc));
		}
		if ($a0_dev > 2)
		{
			$notify_string = $notify_string."The \$a_0\$ masses agree at only \\alert{$a0_dev sigma}, which is large.";
		}
	}
	else
	{
		$notify_string = $notify_string."There are no stochastic \$a_0\$ measurements."
	}
	
	if ($notify_string eq "")
	{
		$notify_string = "None.";
	}
	
	$sed_commands{"NOTIFY_STRING"} = $notify_string;
	
	# The sed_commands hashtable contains all of the sed replacements to perform.
	# Let's build the executable command!
	
	my @split_direc = split('/', $path);
	my $direc = $split_direc[1];

	my $command = "sed -e 's/:ENSEMBLE:/$direc\\/$ensemble/g' "; # < basetex_disc.tex > ../$path/$ensemble/spectrum2/tex/f4plus8l24t48b40m005m080_disc.tex

	foreach my $a_key (keys %sed_commands)
	{
		$command = $command."-e 's/:".$a_key.":/".$sed_commands{$a_key}."/g' ";
	}

	if ($milc_flag == 0)
	{
		$command = $command." < /projectnb/qcd/Staggered/Scripts/EvanSpectrum/tex/basetex_multi.tex > $path/$ensemble/spectrum2/tex/$ensemble".".tex";
	}
	elsif ($milc_flag == 1) # milc_flag set.
	{
		$command = $command." < /projectnb/qcd/Staggered/All/spectrum_connected_24/tex/basetex_lsd.tex > $path/$ensemble/spectrum2/tex/$ensemble".".tex";
	}


	`$command`;
}



sub parenFormWiki
{
   my ($the_value, $the_error) = @_;
   
   my $neg_flag = 1;
   if ($the_value < 0)
   {
		$neg_flag = -1;
		$the_value = -$the_value;
   }
   
   if (abs($the_value) < 1e-10)
   {
	 return "0";
   }

   my $err_size = floor(log($the_error)/log(10));

   my $value_size = floor(log(abs($the_value))/log(10));

   $the_error *= 10**(-$err_size+1);

        if (floor($the_error + 0.5) == 100) # Special case.
        {
                $the_error/=10;
                $err_size--;
        }

   $the_value *= 10**(-$err_size+1);


   $the_value = floor($the_value + 0.5);
   $the_error = floor($the_error + 0.5);

   #print $the_value." ".$the_error."\n";

   # If the value is O(10^1) to O(10^-2), don't exponential notation it.
   
   my $negative_string = "";
   if ($neg_flag == -1)
   {
	$negative_string = "-";
	}
   
   if ($value_size <= 1 && $value_size >= -2)
   {
		switch ($value_size)
		{
			case 1 {
				return $negative_string.substr($the_value, 0, 2).".".substr($the_value, 2)."(".$the_error.")";
			}
			case 0 {
				return $negative_string.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.")";
			}
			case -1 {
				return $negative_string."0.".$the_value."(".$the_error.")";
			}
			case -2 {
				return $negative_string."0.0".$the_value."(".$the_error.")";
			}
		}
   }
   else
   {
		return $negative_string.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.") x 10 ^$value_size^";
	}
   


        #return 0;

}
