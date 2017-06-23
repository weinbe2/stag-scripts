#! /usr/bin/perl -w
use warnings;
#use strict;
use POSIX;
use Switch;

# call: ./spectrum_parse.pl --path ../FourPlusEight/f4plus8l24t48b30m005m020 --parse {conn|disc|decay} --start 20 --stop 400 --sep 20
# Folks editing this code---Ctrl+F for "ADD" to find where you should add new measurements.

if (@ARGV < 8) {
   #ADD
   die "Inappropriate number of arguments: expected ./spectrum_parse.pl --path [path] --parse {conn|disc|decay} --start 20 {--stop 400} --sep 20 {--delete-bad}\n";
}

# Handle the list of arguments.

my $pathname; 

# MEASUREMENT TYPES. ADD EXTENSIONS HERE.
my $parse_conn = 0;
my $parse_disc = 0;
my $parse_decay = 0;

#END MEASUREMENT TYPES.

my $start = 0;
my $sep = 0;
my $stop = -1;
my $delete_bad = 0;

for (my $i = 0; $i < @ARGV; $i++)
{
  if (substr($ARGV[$i], 0, 2) ne "--")
  {
    die "Expected a flag, got another argument.\n";
  }
  
  substr($ARGV[$i], 0, 2) = "";
  
  switch ($ARGV[$i])
  {
    case "path" {
      $pathname = $ARGV[$i+1]; 
      
      # Make sure path exists.
      if (!(-d $pathname))
      {
        die "Ensemble directory ".$pathname." does not exist.\n";
      }
      $i++;
    }
    case "parse" {
      # ADD
      while (($i+1) < @ARGV && substr($ARGV[$i+1], 0, 2) ne "--")
      {
        if ($ARGV[$i+1] eq "conn")
        {
          $parse_conn = 1;
        }
        elsif ($ARGV[$i+1] eq "disc")
        {
          $parse_disc = 1;
        }
        elsif ($ARGV[$i+1] eq "decay")
        {
          $parse_decay = 1;
        }
        else
        {
          die "Unsupported measurement type ".$ARGV[$i+1]."\n";
        }
        $i++
      }
    }
    case "start" {
      $start = $ARGV[$i+1];
      $i++;
    }
    case "sep" {
      $sep = $ARGV[$i+1];
      $i++;
    }
    case "stop" {
      $stop = $ARGV[$i+1];
      $i++;
    }
  }
}

if ($stop < $start && $stop != -1)
{
  die "Stop configuration cannot be less than start configuration.\n"
}

#####################################################################################################
# Learn some things about the ensemble name. There's probably a smart regex for this, but whatever. #
#####################################################################################################

my $ns = 0;
my $nt = 0;
my $beta = 0;
my @flavors = ();
my @mass = ();
my $stream = -1;
my $is_stream = 0;

if (substr($pathname, -1) eq "/")
{
  chop($pathname); # remove last character.
}

# Get last incidence of / in string.
$tmp = rindex($pathname, "/");
my $ensemble = substr($pathname, $tmp+1);

# Do an inefficient run through the string.
if (substr($ensemble, 0, 1) ne "f")
{
  die "Improperly formatted ensemble name: 'f' comes first.\n";
}

# Learn about flavor structure. 
if (substr($ensemble, 0, 7) eq "f4plus8")
{
  $flavors[0] = 4;
  $flavors[1] = 8;
  substr($ensemble, 0, 7) = "";
}
elsif (substr($ensemble, 0, 7) eq "f4plus4")
{
  $flavors[0] = 4;
  $flavors[1] = 4;
  substr($ensemble, 0, 7) = "";
}
elsif (substr($ensemble, 0, 7) eq "f8plus4")
{
  $flavors[0] = 8;
  $flavors[1] = 4;
  substr($ensemble, 0, 7) = "";
}
elsif (substr($ensemble, 0, 3) eq "f12")
{
  $flavors[0] = 12;
  substr($ensemble, 0, 3) = "";
}
elsif (substr($ensemble, 0, 2) eq "f8")
{
  $flavors[0] = 8;
  substr($ensemble, 0, 2) = "";
}
elsif (substr($ensemble, 0, 2) eq "f4")
{
  $flavors[0] = 4;
  substr($ensemble, 0, 2) = "";
}
else
{
  die "Improper number of flavors.\n";
}

# Learn about ns.
if (substr($ensemble, 0, 1) ne "l")
{
  die "Improperly formatted ensemble name: 'l' comes after 'f'.\n";
}

$tmp = index($ensemble, "t");
if ($tmp > 4)
{
  die "Improperly formatted ensemble name: I doubt N_s is in the thousands.\n";
}

$ns = substr($ensemble, 1, $tmp-1);
substr($ensemble, 0, $tmp) = "";

# Learn about nt.
if (substr($ensemble, 0, 1) ne "t")
{
  die "Improperly formatted ensemble name: 't' comes after 'l'.\n";
}

$tmp = index($ensemble, "b");
if ($tmp > 4)
{
  die "Improperly formatted ensemble name: I doubt N_t is in the thousands.\n";
}
if ($tmp == -1)
{
  die "Improperly formatted ensemble name: Beta should come after n_t.\n";
}

$nt = substr($ensemble, 1, $tmp-1);
substr($ensemble, 0, $tmp) = "";

# Learn about beta.

$tmp = index($ensemble, "m");
if ($tmp > 4)
{
  die "Improperly formatted ensemble name: I imagine beta isn't longer than three digits.\n";
}
if ($tmp == -1)
{
  die "Improperly formatted ensemble name: Masses should come after beta.\n";
}

$beta = substr($ensemble, 1, $tmp-1)/10.0;
substr($ensemble, 0, $tmp) = "";

# Learn about masses. 

$tmp = index($ensemble, "m");

if ($tmp == -1)
{
  die "Improperly formatted ensemble name: No masses supplied.\n";
}

while ($tmp != -1)
{
  substr($ensemble, 0, 1) = "";
  
  # see if there's another mass afterwards.
  $tmp = index($ensemble, "m");
  if ($tmp != -1)
  {
    my $tmp_mass = substr($ensemble, 0, $tmp);
    push(@mass, $tmp_mass/(10.0 ** length($tmp_mass)));
    substr($ensemble, 0, $tmp) = "";
  }
  else
  {
    # see if there's a stream afterwards.
    $tmp = index($ensemble, "s");
    if ($tmp != -1)
    {
      $is_stream = 1;
      my $tmp_mass = substr($ensemble, 0, $tmp);
      push(@mass, $tmp_mass/(10.0 ** length($tmp_mass)));
      substr($ensemble, 0, $tmp) = "";
      $tmp = -1;
    }
    else
    {
      # See if it's just a numeric. If so, it's the end.
      if ($ensemble =~ /^[0-9]+$/)
      {
        push(@mass, $ensemble/(10.0 ** length($ensemble)));
      }
      else
      {
        die "Improperly formatted ensemble name: the only valid flag after 'm' is 'm' or 's'.\n";
      }
    }
  } 
}

if ($is_stream == 1)
{
  substr($ensemble, 0, 1) = "";
  if ($ensemble =~ /^[0-9]+$/)
  {
    $stream = $ensemble;
  }
  else
  {
    die "Improperly formatted ensemble name: There should not be any flags after stream.\n";
  }
}

foreach my $item (@flavors)
{
  print "Flavor: $item.\n";
}

print "$ns x $ns x $ns x $nt Volume.\n";
print "Beta = $beta.\n";

foreach my $item (@mass)
{
  print "Mass: $item.\n";
}

if ($is_stream == 1)
{
  print "Stream: $stream.\n";
}


die "End of file.\n"; 

#########################################
# Make sure some directories exist. ADD #
#########################################



### CONTINUE FROM HERE ###
my $direc = $ARGV[0];
my $start = $ARGV[1];
my $sep = $ARGV[2];
my $del_flag = 0;

# Old code to delete runs that don't match the spacing.

#if (@ARGV == 5)
#{
#	$del_flag = $ARGV[4];
#}


my $run = $direc;
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

if ($start == 0)
{
   $start = $cut + $sep;
}

# First, make sure the 'meas' folder exists. If not, error out.

if (!(-d "$path/$direc"."/meas"))
{
	print "Error: a meas directory does not exist. Exiting.\n";
	exit(200);
}

# Next, let's make a 'corr' folder. 

if (!(-d "$path/$direc"."/spectrum2"))
{
        my $command = "mkdir $path/$direc"."/spectrum2";
        `$command`;
}


if (!(-d "$path/$direc"."/spectrum2/corr"))
{
	my $command = "mkdir $path/$direc"."/spectrum2/corr";
	`$command`;
}

if (!(-d "$path/$direc"."/spectrum2/sum"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/sum";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/uwerr"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/uwerr";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/effmass"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/effmass";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/fits"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/fits";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/tex"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/tex";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/central"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/central";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/fits/plots"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/fits/plots";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/effmass/plots"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/effmass/plots";
        `$command`;
}

if (!(-d "$path/$direc"."/spectrum2/sum/plots"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/sum/plots";
        `$command`;
}




# Push out the block size and p-val cut to the data folder.
# This is so we have it on record.
#{
#	open(my $spectrum_info, ">"."$path/$direc"."/spectrum.txt");
#	print $spectrum_info "$blocksize $pval_cut";
#	close($spectrum_info);
#}


#my $tag = $ARGV[0];
#chomp($tag);

#print $tag."\n";

my $a_path = "$path/$direc"."/meas/run*out";

print $a_path."\n";

my @allfiles = glob( $a_path );

my @filecontents = ();

# There are several different correlators.

my @m_ps = ();
my @m_sc = ();
my @m_i5 = ();
my @m_ij = ();
my @m_r0 = ();
my @m_ris = ();
my @m_rij = ();
my @m_ri5 = ();
my @m_ps2 = ();
my @m_rvt = (); # MILC only
my @m_rpv = (); # MILC only
my @b_nu = ();
my @b_de = ();
my @m_times = ();

my $milc_flag = 0;

foreach my $file (@allfiles)
{
	#print $file."\n";
	my @splitpieces = split('/', $file);
	#print $splitpieces[1]."\n";

	my @argh = split('\.', $splitpieces[4]);

	my $n_config = 0;
	my $machine = "";

	# Check for the scc or vulcan.
	if ($argh[1] eq "scc" || $argh[1] eq "vulcan")
	{
		my @argh2 = split('-', $argh[0]);
		print $argh2[1]."\n";
		$n_config = $argh2[1];
		$machine = $argh[1];   
	}
	else
	{
		my @argh2 = split('_', $argh[0]); # $argh2[1] is the number.
		my @argh3 = split('-', $argh2[1]);
		#print $argh3[0]."\n";
		$n_config = $argh3[0];
		$machine = $argh[2];
	}

	my $orig_n_config = $n_config;
	# Assume a separation of 10.
	$n_config = ($n_config-$start)/$sep;

	if (($orig_n_config - $start) % $sep != 0 || $n_config < 0)
        {
           print "NOPE: $orig_n_config\n";
           next; 
	}


	open(my $filething, "<".$file);
	my @parseit = <$filething>;
	close($filething);

	my @m_time_lines = ();

	# There are four chunks to look for.
	my @rwall_lines = ();
	my @pion_lines = ();
	my @rho_lines = ();
	my @rho_milc_lines = (); # MILC
	my @baryon_lines = ();

	for (my $i = 0; $i < @parseit; $i++)
	{
		chomp($parseit[$i]);
		if ($parseit[$i] eq "SINKS: POINT_KAON_5 WALL_KAON_5")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@rwall_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: PION_PS PION_SC PION_i5 PION_ij")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@pion_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: RHO_0 RHO_is RHO_ij RHO_i5")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@rho_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: PION_PS PION_SC RHO_VT RHO_PV")
		{
			$milc_flag = 1;
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@rho_milc_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: NUCLEON DELTA")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@baryon_lines,$parseit[$i+$j+1]);
			}
		}
		else
		{
			my @splitted_vals = split(' ', $parseit[$i]);

			if (@splitted_vals > 1)
			{
				if ($splitted_vals[0] eq "Job")
				{
					push(@m_time_lines, "# time-start ".$parseit[1]);
				}
				elsif ($splitted_vals[1] eq "time-start")
				{
					#print $line."\n";
					push(@m_time_lines, $parseit[$i]);
				}
				elsif ($splitted_vals[1] eq "time-finish")
				{
					#print $line."\n";
					push(@m_time_lines, $parseit[$i]);
				}
			}
		}

	}
	
	if (!(exists $pion_lines[0]) || !((exists $rho_lines[0]) || (exists $rho_milc_lines[0])) || !(exists $baryon_lines[0]) || !(exists $rwall_lines[0])) # || !(@m_time_lines == 2))
	{
		print "File $orig_n_config doesn't exist.\n";
		if ($del_flag != 0)
		{
			my $tmp_command = "rm -rf $file";
			print $tmp_command."\n";
			`$tmp_command`;
		}
	}
	else
	{
		# Next, process each item.
		my @indiv_m_ps = ();
		my @indiv_m_sc = ();
		my @indiv_m_i5 = ();
		my @indiv_m_ij = ();
		my @indiv_m_r0 = ();
		my @indiv_m_ris = ();
		my @indiv_m_rij = ();
		my @indiv_m_ri5 = ();
		my @indiv_m_rvt = (); # MILC
		my @indiv_m_rpv = (); # MILC
		my @indiv_m_ps2 = ();
		my @indiv_b_nu = ();
		my @indiv_b_de = ();
		
		my $local_good_flag = 1;
		
		for (my $i = 0; $i < $nt; $i++)
		{
			# Chomp the line.
			chomp($rwall_lines[$i]);
			chomp($pion_lines[$i]);
			if (exists $rho_lines[0])
			{
				chomp($rho_lines[$i]);
			}
			else
			{
				chomp($rho_milc_lines[$i]);
			}
			chomp($baryon_lines[$i]);
			
			# First, make sure there's the right number of columns for each element.
			my @temparray = split(' ', $rwall_lines[$i]);
			my $tempsize = @temparray;
			if ($tempsize != 5 && $local_good_flag == 1)
			{
				print "NOPE: $orig_n_config Rwall pion\n";
				$local_good_flag = 0;
				next;
			}
			
			@temparray = split(' ', $pion_lines[$i]);
			$tempsize = @temparray;
			if ($tempsize != 9 && $local_good_flag == 1)
			{
				print "NOPE: $orig_n_config Wall pion\n";
				$local_good_flag = 0;
				next;
			}
			
			if (exists $rho_lines[0])
			{
				@temparray = split(' ', $rho_lines[$i]);
				$tempsize = @temparray;
				if ($tempsize != 9 && $local_good_flag == 1)
				{
					print "NOPE: $orig_n_config Wall rho\n";
					$local_good_flag = 0;
					next;
				}
			}
			else
			{
				@temparray = split(' ', $rho_milc_lines[$i]);
				$tempsize = @temparray;
				if ($tempsize != 9 && $local_good_flag == 1)
				{
					print "NOPE: $orig_n_config MILC rho\n";
					$local_good_flag = 0;
					next;
				}
			}
			
			@temparray = split(' ', $baryon_lines[$i]);
			$tempsize = @temparray;
			if ($tempsize != 5 && $local_good_flag == 1)
			{
				print "NOPE: $orig_n_config Wall Baryon\n";
				$local_good_flag = 0;
				next;
			}
			
			@temparray = split(' ', $rwall_lines[$i]);
			push(@indiv_m_ps2, "$i $temparray[1]\n");
			
			
			@temparray = split(' ', $pion_lines[$i]);
			push(@indiv_m_ps, "$i $temparray[1]\n");
			push(@indiv_m_sc, "$i $temparray[3]\n");
			push(@indiv_m_i5, "$i $temparray[5]\n");
			push(@indiv_m_ij, "$i $temparray[7]\n");
			
			if (exists $rho_lines[0])
			{
				@temparray = split(' ', $rho_lines[$i]);
				push(@indiv_m_r0, "$i $temparray[1]\n");
				push(@indiv_m_ris, "$i $temparray[3]\n");
				push(@indiv_m_rij, "$i $temparray[5]\n");
				push(@indiv_m_ri5, "$i $temparray[7]\n");
			}
			else # exists $rho_milc_lines
			{
				@temparray = split(' ', $rho_milc_lines[$i]);
				push(@indiv_m_rvt, "$i $temparray[5]\n");
				push(@indiv_m_rpv, "$i $temparray[7]\n");
			}
			
			
			@temparray = split(' ', $baryon_lines[$i]);
			push(@indiv_b_nu, "$i $temparray[1]\n");
			push(@indiv_b_de, "$i $temparray[3]\n");
		}
		
		if ($local_good_flag == 1)
		{
			
			$m_ps2[$n_config] = \@indiv_m_ps2;
			$m_ps[$n_config] =  \@indiv_m_ps;
			$m_sc[$n_config]= \@indiv_m_sc;
			$m_i5[$n_config]= \@indiv_m_i5;
			$m_ij[$n_config]= \@indiv_m_ij;
			$m_r0[$n_config]= \@indiv_m_r0;
			$m_ris[$n_config]= \@indiv_m_ris;
			$m_rij[$n_config]= \@indiv_m_rij;
			$m_ri5[$n_config]= \@indiv_m_ri5;
			$m_rpv[$n_config]= \@indiv_m_rpv;
			$m_rvt[$n_config]= \@indiv_m_rvt;
			$b_nu[$n_config]= \@indiv_b_nu;
			$b_de[$n_config]= \@indiv_b_de;

			if (@m_time_lines == 2)
			{
				my @begin_split = split(' ', $m_time_lines[0]);
				my @end_split = split(' ', $m_time_lines[1]);

				my @begin_break = split(':', $begin_split[5]);
				my @end_break = split(':', $end_split[5]);

				# Get a difference of seconds.
				my $begin_sec = $begin_break[2];
				my $begin_min = $begin_break[1];
				my $begin_hrs = $begin_break[0];
				my $end_sec = $end_break[2];
				my $end_min = $end_break[1];
				my $end_hrs = $end_break[0];

				my $diff_sec = $end_sec - $begin_sec;
				my $diff_min = $end_min - $begin_min;
				my $diff_hrs = $end_hrs - $begin_hrs;

				if ($diff_sec < 0)
				{
					$diff_min -= 1;
					$diff_sec += 60;
				}

				if ($diff_min < 0)
				{
					$diff_hrs -= 1;
					$diff_min += 60;
				}

				if ($diff_hrs < 0)
				{
					$diff_hrs += 24;
				}

				my $diff_time = 60*60*$diff_hrs + 60*$diff_min + $diff_sec;

				$m_times[$n_config] = $diff_time." ".$machine;
				#print $begin_split[5]." ".$end_split[5]."\n";

				#print $diff_time."\n";
			}
			else
			{
				$m_times[$n_config] = "0 ".$machine;
			}
			
		}
	}

}


open(my $tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ps2");
for (my $i = 0; $i < @m_ps2; $i++)
{
	if (exists $m_ps2[$i])
	{
		foreach my $line (@{$m_ps2[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
	else
	{
		print "NOPE! ".($i*$sep+$start)."\n";
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ps");
for (my $i = 0; $i < @m_ps; $i++)
{
	if (exists $m_ps[$i])
	{
		foreach my $line (@{$m_ps[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.sc");
for (my $i = 0; $i < @m_sc; $i++)
{
	if (exists $m_sc[$i])
	{
		foreach my $line (@{$m_sc[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.i5");
for (my $i = 0; $i < @m_i5; $i++)
{
	if (exists $m_i5[$i])
	{
		foreach my $line (@{$m_i5[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ij");
for (my $i = 0; $i < @m_ij; $i++)
{
	if (exists $m_ij[$i])
	{
		foreach my $line (@{$m_ij[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

if ($milc_flag == 0)
{
	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.r0");
	for (my $i = 0; $i < @m_r0; $i++)
	{
		if (exists $m_r0[$i])
		{
			foreach my $line (@{$m_r0[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);

	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ris");
	for (my $i = 0; $i < @m_ris; $i++)
	{
		if (exists $m_ris[$i])
		{
			foreach my $line (@{$m_ris[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);

	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rij");
	for (my $i = 0; $i < @m_rij; $i++)
	{
		if (exists $m_rij[$i])
		{
			foreach my $line (@{$m_rij[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);

	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ri5");
	for (my $i = 0; $i < @m_ri5; $i++)
	{
		if (exists $m_ri5[$i])
		{
			foreach my $line (@{$m_ri5[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);
}
else # $milc_flag == 1
{
	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rvt");
	for (my $i = 0; $i < @m_rvt; $i++)
	{
		if (exists $m_rvt[$i])
		{
			foreach my $line (@{$m_rvt[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);

	open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rpv");
	for (my $i = 0; $i < @m_rpv; $i++)
	{
		if (exists $m_rpv[$i])
		{
			foreach my $line (@{$m_rpv[$i]})
			{
				print $tmpfile ($i*$sep+$start)." ".$line;
			}
		}
	}
	close($tmpfile);
}

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.nu");
for (my $i = 0; $i < @b_nu; $i++)
{
	if (exists $b_nu[$i])
	{
		foreach my $line (@{$b_nu[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.de");
for (my $i = 0; $i < @b_de; $i++)
{
	if (exists $b_de[$i])
	{
		foreach my $line (@{$b_de[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

#open($tmpfile, ">$path/$direc"."/corr/times.dat");
#for (my $i = 0; $i < @m_times; $i++)
#{
#	if (exists $m_times[$i])
#	{
#		my $time_thing = $i*$sep+$start;
#		print $tmpfile "$time_thing $m_times[$i]\n";
#	}
#}
#close($tmpfile);
