#! /usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# Copy data into a git folder.
# call: ./spectrum_put_git.pl -dir [relative path to Data folder] -ensembles [list of ensembles] -data {meas_conn|meas_disc|pion|pion_i5|pion_ij|pion_sc|fpi|rho|a0|axial}

if (@ARGV < 4)
{
        die "Not enough arguments. Expects \"./spectrum_put_git.pl -dir [relative path to Data folder] -ensembles [list of ensembles] -data {meas_conn|meas_disc|pion|pion_i5|pion_ij|pion_sc|fpi|rho|a0|axial}\"\n";
}

# Handle the list of arguments.
my @ensemble_list = ();
my @data_list = ();

my $data_directory = "";

# A variable to hold file handles in.
my $file_handle;

# Data list. These all correspond to subfolders of the data directory that I manage.

my @safe_data = ("meas_conn", "meas_disc", "pion", "pion_i5", "pion_ij", "pion_sc", "fpi", "rho", "a0", "axial");

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
		case "dir" {
			
			# Check and make sure the Data directory exists.
			if (!(-d $ARGV[$i+1]))
			{
				die "Data directory does not exist.\n";
			}
			
			# If so, store it, and move along!
			$data_directory = $ARGV[$i+1];
			$i++;
			
		}
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
		case "data" {
			# Loop over all subsequent entries until we hit the end or another flag.
			while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 1) ne "-")
			{
				if (!($ARGV[$i+1] ~~ @safe_data)) # Remember, number code at end.
				{
					die "Data $ARGV[$i+1] is not an approved data type.\n";
				}
				push (@data_list, $ARGV[$i+1]);
				$i++;
			}
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

# Make sure the user has recently pulled from the git.
print "Have you recently performed a git pull (y/n): ";
my $pullcheck = <STDIN> ;
chomp ($pullcheck);

# Crash unless the return starts with 'y'.
if (substr(lc($pullcheck), 0, 1) ne 'y')
{
	die "Perform a git pull and try again!\n";
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
			die "Error: a spectrum directory for $ensemble does not exist. Exiting.\n";
	}

	if (!(-d "$path/$ensemble/spectrum2/corr"))
	{
			die "Error: the raw parsed correlators for $ensemble do not exist. Please parse the data before running! Exiting.\n";
	}
	
	# Loop over all possible data measurements.
	foreach my $data_pick (@data_list)
	{
		switch ($data_pick)
		{
			case "meas_conn" {
				# We need the number of measurements and the bin size.
				my $num_meas = 0;
				my $bin_size = 0;
				
				# Borrow this trick from "spectrum_meas_texcopy.pl"
				
				if (-f "$path/$ensemble/spectrum2/corr/corr.ps")
				{
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

					$num_meas = (@pion_file/$nt);
					
					# Next, get the bin size. 
					# Make sure a fitparams exists. 
					if (-f "$path/$ensemble/spectrum2/fitparams/fitparam.ps")
					{
						# Load the fitparams file.
						open($file_handle, "<$path/$ensemble/spectrum2/fitparams/fitparam.ps");
						my @info_line = <$file_handle>;
						close($file_handle);
						
						# The binsize is the last (5th) entry.
						my @info_splitter = split(' ', $info_line[0]);
						$bin_size = $info_splitter[4];
						
						# Excellent! We have all the information, now spit it out to file.
						open($file_handle, ">$data_directory/meas_conn/meas_conn_$ensemble");
						print $file_handle "$ensemble $num_meas $bin_size\n";
						close($file_handle);
					}
					else
					{
						print "No way to find the preferred binsize. Run spectrum_fit_inverval.pl for $ensemble and state ps.\n";
					}
				}
				else
				{
					print "No way to count the number of connected measurements exists! Run spectrum_meas_parse.pl for $ensemble.\n";
				}
			}
			case "meas_disc" {
				# We need the number of measurements and the bin size.
				my $num_meas = 0;
				my $bin_size = 0;
				
				# Borrow this trick from "spectrum_meas_texcopy.pl"
				
				if (-f "$path/$ensemble/spectrum2/corr/corr.sc_stoch")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/corr/corr.sc_stoch");
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

					$num_meas = (@pion_file/$nt);
					
					# Next, get the bin size. 
					# Make sure a fitparams exists. 
					if (-f "$path/$ensemble/spectrum2/fitparams/fitparam.dc_stoch")
					{
						# Load the fitparams file.
						open($file_handle, "<$path/$ensemble/spectrum2/fitparams/fitparam.dc_stoch");
						my @info_line = <$file_handle>;
						close($file_handle);
						
						# The binsize is the last (5th) entry.
						my @info_splitter = split(' ', $info_line[0]);
						$bin_size = $info_splitter[4];
						
						# Excellent! We have all the information, now spit it out to file.
						open($file_handle, ">$data_directory/meas_disc/meas_disc_$ensemble");
						print $file_handle "$ensemble $num_meas $bin_size\n";
						close($file_handle);
					}
					else
					{
						print "No way to find the preferred binsize. Run spectrum_fit_inverval.pl for $ensemble and state dc_stoch.\n";
					}
				}
				else
				{
					print "No way to count the number of disconnected measurements exists! Run spectrum_sigma2_parse.pl for $ensemble, and then spectrum_connected_multi.pl -do build.\n";
				}
				
			}
			case "pion" {
				# This and subsequent entries are based on the following shell script:
				#direc=$1
				#state=$2
				#stateline=$3
				#
				#while read line; do
				#	
				#	echo -n $line" " > /usr3/graduate/weinbe2/fourpluseight/Data/$direc/${direc}_${line};
				#	cat /projectnb/qcd/Staggered/FourPlusEight/$line/spectrum2/central/central_new.${state} | sed -n "${stateline}p" >> /usr3/graduate/weinbe2/fourpluseight/Data/$direc/${direc}_${line}
				#	echo "$line"
				#done < $4
				
				if (-f "$path/$ensemble/spectrum2/central/central_new.ps")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.ps");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[0]);
					open($file_handle, ">$data_directory/pion/pion_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the pion exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states ps.\n";
				}
			}			
			case "pion_i5" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.i5")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.i5");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[0]);
					open($file_handle, ">$data_directory/pion_i5/pion_i5_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the pion_i5 exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states i5.\n";
				}
			}
			case "pion_ij" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.ij")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.ij");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[0]);
					open($file_handle, ">$data_directory/pion_ij/pion_ij_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the pion_i5 exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states ij.\n";
				}
			}
			case "pion_sc" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.sc")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.sc");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[0]);
					open($file_handle, ">$data_directory/pion_sc/pion_sc_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the pion_sc exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states sc.\n";
				}
			}
			case "fpi" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.ps2")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.ps2");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[1]);
					open($file_handle, ">$data_directory/fpi/fpi_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of fpi exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states ps2.\n";
				}
			}
			case "rho" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.ri5")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.ri5");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[0]);
					open($file_handle, ">$data_directory/rho/rho_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the rho exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states ri5.\n";
				}
			}
			case "a0" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.sc")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.sc");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[1]);
					open($file_handle, ">$data_directory/a0/a0_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the a0 exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states sc.\n";
				}
			}
			case "axial" {
				if (-f "$path/$ensemble/spectrum2/central/central_new.ri5")
				{
					open($file_handle, "<$path/$ensemble/spectrum2/central/central_new.ri5");
					my @central_file = <$file_handle>;
					close($file_handle);
					
					my @tmp_split = split(' ', $central_file[1]);
					open($file_handle, ">$data_directory/axial/axial_$ensemble");
					print $file_handle "$ensemble $tmp_split[0] $tmp_split[1]\n";
					close($file_handle);
				}
				else
				{
					# Don't crash, just point it out.
					die "No measurement of the axial exists! Run ./spectrum_connected_multi.pl -do analyze plots -ensembles $ensemble -states ri5.\n";
				}
			}
		}
		
	}

}

