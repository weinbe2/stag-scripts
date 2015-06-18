#! /usr/bin/perl 
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_meas_parse.pl f4plus8l24t48b30m005m020 20 20

if (@ARGV != 3) {
   die "Inappropriate number of arguments: expected connected_check.pl [directory] [start] [space]\n";
}

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

my $a_path = "$path/$direc"."/decay/decay*";

print $a_path."\n";

my @allfiles = glob( $a_path );

my @filecontents = ();

# There are several different correlators.

my @p_ll = ();
my @p_lc = ();
my @p_cc = ();
my @r_ll = ();
my @r_ll_x = ();
my @r_ll_y = ();
my @r_ll_z = ();
my @r_cc = ();
my @r_cc_x = ();
my @r_cc_y = ();
my @r_cc_z = ();
my @a_ll = ();
my @a_ll_x = ();
my @a_ll_y = ();
my @a_ll_z = ();
my @a_cc = ();
my @a_cc_x = ();
my @a_cc_y = ();
my @a_cc_z = ();

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

	# There are three chunks to look for.
	my @pion_lines = ();
	my @ll_lines = ();
	my @cc_lines = ();

	for (my $i = 0; $i < @parseit; $i++)
	{
		chomp($parseit[$i]);
		if ($parseit[$i] eq "SINKS: PION_LL PION_LC PION_CC")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@pion_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: RHO_LL RHO_LL_X RHO_LL_Y RHO_LL_Z AXIAL_LL AXIAL_LL_X AXIAL_LL_Y AXIAL_LL_Z")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@ll_lines,$parseit[$i+$j+1]);
			}
		}
		elsif ($parseit[$i] eq "SINKS: RHO_CC RHO_CC_X RHO_CC_Y RHO_CC_Z AXIAL_CC AXIAL_CC_X AXIAL_CC_Y AXIAL_CC_Z")
		{
			for (my $j = 0; $j < $nt; $j++)
			{
				push(@cc_lines,$parseit[$i+$j+1]);
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
	
	if (!(exists $pion_lines[0]) || !(exists $ll_lines[0]) || !(exists $cc_lines[0])) # || !(@m_time_lines == 2))
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
		my @indiv_p_ll = ();
      my @indiv_p_lc = ();
      my @indiv_p_cc = ();
      my @indiv_r_ll = ();
      my @indiv_r_ll_x = ();
      my @indiv_r_ll_y = ();
      my @indiv_r_ll_z = ();
      my @indiv_r_cc = ();
      my @indiv_r_cc_x = ();
      my @indiv_r_cc_y = ();
      my @indiv_r_cc_z = ();
      my @indiv_a_ll = ();
      my @indiv_a_ll_x = ();
      my @indiv_a_ll_y = ();
      my @indiv_a_ll_z = ();
      my @indiv_a_cc = ();
      my @indiv_a_cc_x = ();
      my @indiv_a_cc_y = ();
      my @indiv_a_cc_z = ();
		
		for (my $i = 0; $i < $nt; $i++)
		{
			chomp($pion_lines[$i]);
			my @temparray = split(' ', $pion_lines[$i]);
			push(@indiv_p_ll, "$i $temparray[1]\n");
			push(@indiv_p_lc, "$i $temparray[3]\n");
			push(@indiv_p_cc, "$i $temparray[5]\n");
			
			chomp($ll_lines[$i]);
			@temparray = split(' ', $ll_lines[$i]);
			push(@indiv_r_ll, "$i $temparray[1]\n");
			push(@indiv_r_ll_x, "$i $temparray[3]\n");
			push(@indiv_r_ll_y, "$i $temparray[5]\n");
			push(@indiv_r_ll_z, "$i $temparray[7]\n");
			push(@indiv_a_ll, "$i $temparray[9]\n");
			push(@indiv_a_ll_x, "$i $temparray[11]\n");
			push(@indiv_a_ll_y, "$i $temparray[13]\n");
			push(@indiv_a_ll_z, "$i $temparray[15]\n");
			
			chomp($cc_lines[$i]);
			@temparray = split(' ', $cc_lines[$i]);
			push(@indiv_r_cc, "$i $temparray[1]\n");
			push(@indiv_r_cc_x, "$i $temparray[3]\n");
			push(@indiv_r_cc_y, "$i $temparray[5]\n");
			push(@indiv_r_cc_z, "$i $temparray[7]\n");
			push(@indiv_a_cc, "$i $temparray[9]\n");
			push(@indiv_a_cc_x, "$i $temparray[11]\n");
			push(@indiv_a_cc_y, "$i $temparray[13]\n");
			push(@indiv_a_cc_z, "$i $temparray[15]\n");
			
		}
		
		$p_ll[$n_config] = \@indiv_p_ll;
      $p_lc[$n_config] = \@indiv_p_lc;
      $p_cc[$n_config] = \@indiv_p_cc;
      $r_ll[$n_config] = \@indiv_r_ll;
      $r_ll_x[$n_config] = \@indiv_r_ll_x;
      $r_ll_y[$n_config] = \@indiv_r_ll_y;
      $r_ll_z[$n_config] = \@indiv_r_ll_z;
      $r_cc[$n_config] = \@indiv_r_cc;
      $r_cc_x[$n_config] = \@indiv_r_cc_x;
      $r_cc_y[$n_config] = \@indiv_r_cc_y;
      $r_cc_z[$n_config] = \@indiv_r_cc_z;
      $a_ll[$n_config] = \@indiv_a_ll;
      $a_ll_x[$n_config] = \@indiv_a_ll_x;
      $a_ll_y[$n_config] = \@indiv_a_ll_y;
      $a_ll_z[$n_config] = \@indiv_a_ll_z;
      $a_cc[$n_config] = \@indiv_a_cc;
      $a_cc_x[$n_config] = \@indiv_a_cc_x;
      $a_cc_y[$n_config] = \@indiv_a_cc_y;
      $a_cc_z[$n_config] = \@indiv_a_cc_z;
	

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


open(my $tmpfile, ">$path/$direc"."/spectrum2/corr/corr.pll");
for (my $i = 0; $i < @m_p_ll; $i++)
{
	if (exists $m_p_ll[$i])
	{
		foreach my $line (@{$m_p_ll[$i]})
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

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.plc");
for (my $i = 0; $i < @m_p_lc; $i++)
{
	if (exists $m_p_lc[$i])
	{
		foreach my $line (@{$m_p_lc[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.pcc");
for (my $i = 0; $i < @m_p_cc; $i++)
{
	if (exists $m_p_cc[$i])
	{
		foreach my $line (@{$m_p_cc[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rll");
for (my $i = 0; $i < @m_r_ll; $i++)
{
	if (exists $m_r_ll[$i])
	{
		foreach my $line (@{$m_r_ll[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rllx");
for (my $i = 0; $i < @m_r_ll_x; $i++)
{
	if (exists $m_r_ll_x[$i])
	{
		foreach my $line (@{$m_r_ll_x[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rlly");
for (my $i = 0; $i < @m_r_ll_y; $i++)
{
	if (exists $m_r_ll_y[$i])
	{
		foreach my $line (@{$m_r_ll_y[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rllz");
for (my $i = 0; $i < @m_r_ll_z; $i++)
{
	if (exists $m_r_ll_z[$i])
	{
		foreach my $line (@{$m_r_ll_z[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rcc");
for (my $i = 0; $i < @m_r_cc; $i++)
{
	if (exists $m_r_cc[$i])
	{
		foreach my $line (@{$m_r_cc[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rccx");
for (my $i = 0; $i < @m_r_cc_x; $i++)
{
	if (exists $m_r_cc_x[$i])
	{
		foreach my $line (@{$m_r_cc_x[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rccy");
for (my $i = 0; $i < @m_r_cc_y; $i++)
{
	if (exists $m_r_cc_y[$i])
	{
		foreach my $line (@{$m_r_cc_y[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.rccz");
for (my $i = 0; $i < @m_r_cc_z; $i++)
{
	if (exists $m_r_cc_z[$i])
	{
		foreach my $line (@{$m_r_cc_z[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.all");
for (my $i = 0; $i < @m_a_ll; $i++)
{
	if (exists $m_a_ll[$i])
	{
		foreach my $line (@{$m_a_ll[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.allx");
for (my $i = 0; $i < @m_a_ll_x; $i++)
{
	if (exists $m_a_ll_x[$i])
	{
		foreach my $line (@{$m_a_ll_x[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.ally");
for (my $i = 0; $i < @m_a_ll_y; $i++)
{
	if (exists $m_a_ll_y[$i])
	{
		foreach my $line (@{$m_a_ll_y[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.allz");
for (my $i = 0; $i < @m_a_ll_z; $i++)
{
	if (exists $m_a_ll_z[$i])
	{
		foreach my $line (@{$m_a_ll_z[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.acc");
for (my $i = 0; $i < @m_a_cc; $i++)
{
	if (exists $m_a_cc[$i])
	{
		foreach my $line (@{$m_a_cc[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.accx");
for (my $i = 0; $i < @m_a_cc_x; $i++)
{
	if (exists $m_a_cc_x[$i])
	{
		foreach my $line (@{$m_a_cc_x[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.accy");
for (my $i = 0; $i < @m_a_cc_y; $i++)
{
	if (exists $m_a_cc_y[$i])
	{
		foreach my $line (@{$m_a_cc_y[$i]})
		{
			print $tmpfile ($i*$sep+$start)." ".$line;
		}
	}
}
close($tmpfile);

open($tmpfile, ">$path/$direc"."/spectrum2/corr/corr.accz");
for (my $i = 0; $i < @m_a_cc_z; $i++)
{
	if (exists $m_a_cc_z[$i])
	{
		foreach my $line (@{$m_a_cc_z[$i]})
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
