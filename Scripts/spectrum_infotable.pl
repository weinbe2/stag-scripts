#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# Create a table of data! Output to either wiki format or to a tex file. 
# Optionally, you can specify if errors should be split or not.
# Call style:
# ./spectrum_infotable.pl [wiki|tex] [-ensembles] {list of ensembles} [-states] {list of states}

if (@ARGV < 5)
{
        die "Not enough arguments. Expects \"./spectrum_infotable.pl [wiki|tex] [-ensembles] {list of ensembles} [-states] {list of states} {-errors}\"\n";
}

# Figure out wiki or text.

my $output_flag = "";
if ($ARGV[0] eq "wiki")
{
	$output_flag = "wiki";
}
elsif ($ARGV[0] eq "tex")
{
	$output_flag = "tex";
}
else
{
	die "Invalid output format. Expected a first argument of \"wiki\" or \"states\".\n";
}

# Handle the rest of the arguments.

my @ensemble_list = ();
my @state_list = ();
my $error_form = 1; # 1 = parenthasis, 2 = separate.

# Extra things: "wflow", "nconfig", "nmeas_conn", "nmeas_disc"
my @safe_states = ("sc_stoch", "dc_stoch", "dc_stoch_ppp", "dc_stoch_oscil", "sg_stoch", "sg_stoch_ppp", "ps", "ps2", "sc", "ij", "i5", "ri5", "ris", "rij", "r0", "nu", "de", "rvt", "rpv", "wflow", "nconfig", "nmeas_conn", "nmeas_disc");

my %safe_state_names_wiki = (
	sc_stoch1 => "aM ~a0~",
	sc_stoch2 => "aM ~&pi; sc~",
	dc_stoch1 => "aM ~&sigma; D~",
	dc_stoch_ppp1 => "aM ~&sigma; D~",
	dc_stoch_oscil1 => "aM ~&sigma; D~",
	dc_stoch_oscil2 => "aM ~&pi; sc D~",
	sg_stoch1 => "aM ~&sigma; D-C~",
	sg_stoch2 => "aM ~&pi; D-C~",
	ps1 => "aM ~&pi; GB~",
	ps21 => "aM ~&pi; GB~",
	ps22 => "aF ~&pi;~",
	sc1 => "aM ~&pi; sc~",
	sc2 => "aM ~a0~",
	ij1 => "aM ~&pi; ij~",
	i51 => "aM ~&pi; i5~",
	i52 => "aM ~a0 i5~",
	ri51 => "aM ~&rho; i5~",
	ri52 => "aM ~axial~",
	ris1 => "aM ~&rho; is~",
	ris2 => "aM ~tensor~",
	rij1 => "aM ~&rho; ij~",
	rij2 => "aM ~tensor~",
	r01 => "aM ~&rho; 0~",
	r02 => "aM ~axial~",
	nu1 => "aM ~nucl+~",
	nu2 => "aM ~nucl-~",
	de1 => "aM ~delta+~",
	de2 => "aM ~delta-~",
	rvt1 => "aM ~&rho; i~",
	rvt2 => "aM ~tensor~",
	rpv1 => "aM ~&rho; i0",
	rpv2 => "aM ~axial~",
	wflow => "a ^-1^ sqrt(8t ~0~ )",
);

my $command = "";

for (my $i = 1; $i < @ARGV; $i++)
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
				if (!(substr($ARGV[$i+1],0,-1) ~~ @safe_states) && ($ARGV[$i+1] ne "wflow")) # Remember, number code at end.
				{
					die "State $ARGV[$i+1] is not an approved state.\n";
				}
				push (@state_list, $ARGV[$i+1]);
				$i++;
			}
		}
		case "errors" {
			$error_form = 2; # Have errors separate!
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

foreach my $state (@state_list)
{
#	print $state."\n";
}

foreach my $ensemble (@ensemble_list)
{
#	print $ensemble."\n";
}

print "Note: For the time being, only \"wiki\" output is supported.\n";


# Okay! Now we go ensemble by ensemble, state by state, and build a table.

my $meas_size = @state_list;
print "|\\8=. |\\2=. Connected |\\2=. Disconnected |\\".($meas_size*$error_form)."=. |\n";
print "|_. &beta; |_. L |_. T |_. N ~fl~ |_. N ~fh~ |_. am ~l~ |_. am ~h~ |_. N ~config~ |_. N ~meas~ |_. MDTU/meas |_. N ~meas~ |_. MDTU/meas |"; 

foreach my $state (@state_list)
{
	print "_. ".$safe_state_names_wiki{$state}." |";
	if ($error_form == 2)
	{
		print "_. ".$safe_state_names_wiki{$state}." err |";
	}
}

print "\n";

foreach my $item (@ensemble_list)
{
	my %params = %{$item};
	
	# Print the intro info.

	print "| $params{beta} | $params{L} | $params{T} | $params{have_l} | $params{have_h} | $params{ml} | $params{mh} | ";
	
	# Get nconfig:
	
	{
		my $traj_count = 0;
		my $cut = $params{cut};
		if (-f "$params{path}/$params{ensemble}/info.txt")
		{
			# Grab the TrajectoryNo.
			$command = "grep \"trajectories\" $params{path}/$params{ensemble}/info.txt";
			#print $command."\n";
			
			my @trajinfo = `$command`;
			
			if (@trajinfo > 0)
			{
			   my @trajsplit = split("->", $trajinfo[0]);
			   if (@trajsplit > 0)
			   {
			   	$traj_count = $trajsplit[1];
			   }
			}
		
		}
		
		if (($traj_count - $cut) >= 0)
		{
			print " ".(($traj_count - $cut)/10)." | ";
		}
		else
		{
			print " --- | ";
		}
	}
	
	# Get nmeas_conn
	{
		# Let's just get creative here, shall we?
		if (-d "$params{path}/$params{ensemble}/spectrum2/corr")
		{
			$command = "cat $params{path}/$params{ensemble}/spectrum2/corr/corr.ps | wc -l";
			my $nlines = `$command`;
			my $nmeas = $nlines/$params{T};
			
			# get the MDTU between things.
			$command = "head -1 $params{path}/$params{ensemble}/spectrum2/corr/corr.ps | tail -1";
			my @split_start = split(' ', `$command`);
			
			$command = "head -".(1+$params{T})." $params{path}/$params{ensemble}/spectrum2/corr/corr.ps | tail -1";
			my @split_end = split(' ', `$command`);
			
			my $MDTU_block = $split_end[0]-$split_start[0];
			
			if (-f "$params{path}/$params{ensemble}/spectrum2/fitparams/fitparam.ps")
			{
				my $command = "cat $params{path}/$params{ensemble}/spectrum2/fitparams/fitparam.ps";
				my @splitting = split(' ', `$command`);
				print " ".$nmeas."/".$splitting[4]." | ".$MDTU_block." MDTU | ";
			}
			else
			{
				print " ".$nmeas." | --- | ";
			}
		}
		else
		{
			print " --- | --- | ";
		}
	}
	
	# get nmeas_disc
	{
		if (-d "$params{path}/$params{ensemble}/spectrum2/stoch")
		{
			$command = "cat $params{path}/$params{ensemble}/spectrum2/stoch/SPECTRUM.dat | wc -l";
			my $nlines = `$command`;
			my $nmeas = $nlines/$params{T};
			
			# get the MDTU between things.
			$command = "head -1 $params{path}/$params{ensemble}/spectrum2/stoch/SPECTRUM.dat | tail -1";
			my @split_start = split(' ', `$command`);
			
			$command = "head -".(1+$params{T})." $params{path}/$params{ensemble}/spectrum2/stoch/SPECTRUM.dat | tail -1";
			my @split_end = split(' ', `$command`);
			
			my $MDTU_block = $split_end[0]-$split_start[0];
			
			if (-f "$params{path}/$params{ensemble}/spectrum2/fitparams/fitparam.dc_stoch")
			{
				my $command = "cat $params{path}/$params{ensemble}/spectrum2/fitparams/fitparam.dc_stoch";
				my @splitting = split(' ', `$command`);
				print " ".$nmeas."/".$splitting[4]." | ".$MDTU_block." MDTU | ";
			}
			else
			{
				print " ".$nmeas." | --- | ";
			}
		}
		else
		{
			print " --- | --- | ";
		}
	}
	
	# Make sure the central folder exists!
	if (!(-d "$params{path}/$params{ensemble}/spectrum2/central"))
	{
		print "Error: a spectrum or central directory for $params{ensemble} does not exist.\n";
	}
	
	# Go through all states.
	foreach my $state (@state_list)
	{
		switch ($state)
		{
			case "wflow" {
				if (-f "$params{path}/$params{ensemble}/plots/wflow.txt")
				{
					open (my $temphandle, "<$params{path}/$params{ensemble}/plots/wflow.txt");
					my @in_file = <$temphandle>;
					close($temphandle);
					
					my @split_value = split(' ', $in_file[0]);
					
					if ($error_form == 1)
					{
						my $mass_val = parenFormWiki($split_value[2], $split_value[3]);
						print $mass_val." | ";
					}
					elsif ($error_form == 2)
					{
						print $split_value[2]." | ".$split_value[3]." | ";
					}
				}
				else
				{
					print " --- | " x $error_form; # print $error_form times.
				}
			}
			else { # Regular spectrum state.
				my $state_code = substr($state, 0, -1);
				my $state_num = substr($state, -1);
				
				if (-f "$params{path}/$params{ensemble}/spectrum2/central/central_new.".$state_code)
				{
					open (my $temphandle, "<$params{path}/$params{ensemble}/spectrum2/central/central_new.".$state_code);
					my @in_file = <$temphandle>;
					close($temphandle);
					
					my @split_value = split(' ', $in_file[$state_num-1]);
					
					if ($error_form == 1)
					{
						my $mass_val = parenFormWiki($split_value[0], $split_value[1]);
						print $mass_val." | ";
					}
					elsif ($error_form == 2)
					{
						print $split_value[0]." | ".$split_value[1]." | ";
					}
				}
				else
				{
					print " --- | " x $error_form;
				}
			}
		}
	}
	
	print "\n";
}








sub parenForm
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

   if ($neg_flag == -1)
   {
		return '-'.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.")\\mbox{e".$value_size.'} ';
	}
	else
	{
		return ''.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.")\\mbox{e".$value_size.'} ';
	}
   


        #return 0;

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

=begin GHOSTCODE


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


