#! /usr/bin/perl -w

# Parse the new multiexponential fits to produce plots and tables, depending on the flags.

use POSIX;
use strict;
use Switch;
use BSMrun;

# Define where various parameters are saved in save files.
use constant {
	TMIN => 0,
	TMAX => 1,
	FITTYPE => 2,
	DIRAMP1 => 3,
	DIRMASS1 => 5,
	OSCAMP1 => 15,
	OSCMASS1 => 17,
	PVAL => 28,
	AUX1 => 41,
	TMIN_OLD => 0,
	TMAX_OLD => 1,
	DIRAMP1_OLD => 2,
	DIRMASS1_OLD => 4,
	OSCAMP1_OLD => 10,
	OSCMASS1_OLD => 12,
	PVAL_OLD => 19,
	AUX1_OLD => 21,
};


# Call as:
# ./run_multiconn.pl -do {plots tables} -ensemble [ensemble] -state [state]

# Let's grab all the data we can!

if (@ARGV < 4)
{
	die "Call as ./run_multiconn.pl -do {plots tables} -ensemble [ensemble] -state [state]\n";
}

my $ensemble = "";
my $state = "";

my $plots = 0;
my $tables = 0;

my %ensemble_params = ();

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
				if ($ARGV[$i+1] eq "tables")
				{
					$tables = 1;
				}
				$i++;
			}
		}
		case "ensemble" {
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
		
			%ensemble_params = %params;
			$i++;
		}
		case "state" {
			$state = $ARGV[$i+1];
			$i++;
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

# Prepare to run some shell commands later.
my $command = "";
my $path = "/projectnb/qcd/Staggered/Scripts";

# Prepare some state names.
my $dir_label = "";
my $osc_label = "";
my $is_fpi = 0;
switch ($state)
{
	case "ps" {
		$dir_label = "M_\\\\pi";
	}
	case "ps2" {
		$dir_label = "M_\\\\pi";
		$osc_label = "F_\\\\pi";
		$is_fpi = 1;
	}
	case "pll" {
		$dir_label = "M_\\\\pi";
		$osc_label = "F_\\\\pi";
		$is_fpi = 1;
	}
	case "ri5" {
		$dir_label = "M_\\\\rho";
		$osc_label = "M_{a_1}";
	}
	case "sc" {
		$dir_label = "M_{\\\\pi_{sc}}";
		$osc_label = "M_{a_0}";
	}
	case "sc_stoch" {
		$dir_label = "M_{a_0}";
		$osc_label = "M_{\\\\pi_{sc}}";
	}
	case "nu" {
		$dir_label = "M_{N+}";
		$osc_label = "M_{N-}";
	}
	else: {
		$dir_label = "M_{dir}";
		$osc_label = "M_{osc}";
	}
}

printf("Dir lab. %s\n", $dir_label);
printf("Osc lab. %s\n", $osc_label);

# First, check and see if the multifit file exists.
my $direc = $ensemble_params{ensemble};
my $relpath = $ensemble_params{path};
my @multifit_file;
if (-f "$path/$relpath/$direc/spectrum2/multifits/multifits.$state")
{
	open (my $fit_handle, "<$path/$relpath/$direc/spectrum2/multifits/multifits.$state");
	@multifit_file = <$fit_handle>;
	close($fit_handle);
}
else
{
	die "No multifit file exists!\n";
}

# Good! Now let's learn about it. We'll split it up, too.
my @fit_info = ();
my $num_dir = 0; # What's the maximum number of direct states?
my $num_osc = 0; # What's the maximum number of oscillating states?
my @num_change = (); # Find out where the number of states changes.
my $tmin_min = 1000;
my $tmin_max = 0;
my $tmax = 0;
my $prev_ln_dir = 0;
my $prev_ln_osc = 0;
my $num_fits = $#multifit_file;
for (my $i = 1; $i <= $#multifit_file; $i++)
{
	my @splitted = split(' ', $multifit_file[$i]);
	push(@fit_info, \@splitted);
	
	if ($splitted[TMIN] < $tmin_min)
	{
		$tmin_min = $splitted[TMIN];
	}
	
	if ($splitted[TMIN] > $tmin_max)
	{
		$tmin_max = $splitted[TMIN];
	}
	
	my $tmax = $splitted[TMAX]; # This should always be the same.
	
	my $ln_dir = ($splitted[FITTYPE])%4;
	my $ln_osc = ($splitted[FITTYPE]-$ln_dir)/4;
	
	if ($ln_dir > $num_dir)
	{
		$num_dir = $ln_dir;
	}
	
	if ($ln_osc > $num_osc)
	{
		$num_osc = $ln_osc;
	}
	
	if ($ln_dir != $prev_ln_dir || $ln_osc != $prev_ln_osc)
	{
		if ($i > 1) # Obviously there's a change at the first one.
		{
			push(@num_change, $splitted[TMIN]+0.5);
		}
	}
	
	$prev_ln_dir = $ln_dir;
	$prev_ln_osc = $ln_osc;
	
}

# Good! Print this out.

printf("Tmin_min = %d\n", $tmin_min);
printf("Tmin_max = %d\n", $tmin_max);
printf("Num dir. = %d\n", $num_dir);
printf("Num osc. = %d\n", $num_osc);
printf("Num chg. at ");
foreach my $tchg (@num_change)
{
	printf("%.1f ", $tchg);
}
printf("\n");

# Alright! Learn something about the mean and error of ground states.
# This helps us pick y-intervals.

# Get the avg ground state mass, avg ground state error.
my $mean_dir = 0;
my $err_dir = 0;
my $mean_osc = 0;
my $err_osc = 0;

foreach my $elem (@fit_info)
{
	my @vals = @{$elem};
	if ($num_dir > 0)
	{
		$mean_dir += $vals[DIRMASS1]/$num_fits;
		$err_dir += $vals[DIRMASS1+1]/$num_fits;
	}
	if ($num_osc > 0)
	{
		$mean_osc += $vals[OSCMASS1]/$num_fits;
		$err_osc += $vals[OSCMASS1+1]/$num_fits;
	}
	if ($is_fpi == 1) # We have fpi, reuse oscil.
	{
		$mean_osc += $vals[AUX1]/$num_fits;
		$err_osc += $vals[AUX1+1]/$num_fits;
	}
}

printf("Avg dir. = %.6f\n", $mean_dir);
printf("Err dir. = %.6f\n", $err_dir);
printf("Avg osc. = %.6f\n", $mean_osc);
printf("Err osc. = %.6f\n", $err_osc);

# Bound the error to be < 0.01. This keeps the y-range on plots
# from being too large.
if ($err_dir > 0.01)
{
	$err_dir = 0.01;
}
if ($err_osc > 0.01)
{
	$err_osc = 0.01;
}

# Before building plots, if it exists, grab the single-state fits.
my $singlefit_flag = 0;
my @singlefit_file;
if (-f "$path/$relpath/$direc/spectrum2/fits/fit_new.$state")
{
	open (my $fit_handle, "<$path/$relpath/$direc/spectrum2/fits/fit_new.$state");
	@singlefit_file = <$fit_handle>;
	$singlefit_flag = 1;
	close($fit_handle);
}

my @singlefit_info = ();
if ($singlefit_flag == 1)
{
	for (my $i = 0; $i <= $#singlefit_file; $i++)
	{
		my @splitted = split(' ', $singlefit_file[$i]);
		push(@singlefit_info, \@splitted);
	}
}

# Save a local copy of data we need to file.
if ($num_dir > 0)
{
	for (my $i = 0; $i < $num_dir; $i++)
	{
		open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/multifit_dir$i");
		foreach my $a_row (@fit_info)
		{
			my @the_row = @{$a_row};
			if ($the_row[DIRMASS1+4*$i] > 1e-6)
			{
				print $outfile3 $the_row[TMIN]." ".$the_row[DIRMASS1+4*$i]." ".$the_row[DIRMASS1+4*$i+1]."\n";
			}
		}
		close($outfile3);
	}
	
	if ($singlefit_flag == 1)
	{
		open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/singlefit_dir");
		foreach my $a_row (@singlefit_info)
		{
			my @the_row = @{$a_row};
			if ($the_row[DIRMASS1_OLD] > 1e-6 && $the_row[TMIN_OLD] < $num_change[0])
			{
				print $outfile3 $the_row[TMIN_OLD]." ".$the_row[DIRMASS1_OLD]." ".$the_row[DIRMASS1_OLD+1]."\n";
			}
		}
		close($outfile3);
	}
}

if ($num_osc > 0)
{
	for (my $i = 0; $i < $num_osc; $i++)
	{
		open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/multifit_osc$i");
		foreach my $a_row (@fit_info)
		{
			my @the_row = @{$a_row};
			if ($the_row[OSCMASS1+4*$i] > 1e-6)
			{
				print $outfile3 $the_row[TMIN]." ".$the_row[OSCMASS1+4*$i]." ".$the_row[OSCMASS1+4*$i+1]."\n";
			}
		}
		close($outfile3);
	}
	
	if ($singlefit_flag == 1)
	{
		open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/singlefit_osc");
		foreach my $a_row (@singlefit_info)
		{
			my @the_row = @{$a_row};
			if ($the_row[DIRMASS1_OLD] > 1e-6 && $the_row[TMIN_OLD] < $num_change[0])
			{
				print $outfile3 $the_row[TMIN_OLD]." ".$the_row[OSCMASS1_OLD]." ".$the_row[OSCMASS1_OLD+1]."\n";
			}
		}
		close($outfile3);
	}
}

if ($is_fpi == 1)
{

	open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/multifit_osc0");
	foreach my $a_row (@fit_info)
	{
		my @the_row = @{$a_row};
		if ($the_row[AUX1] > 1e-6)
		{
			print $outfile3 $the_row[TMIN]." ".$the_row[AUX1]." ".$the_row[AUX1+1]."\n";
		}
	}
	close($outfile3);

	
	if ($singlefit_flag == 1)
	{
		open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/singlefit_osc");
		foreach my $a_row (@singlefit_info)
		{
			my @the_row = @{$a_row};
			if ($the_row[AUX1_OLD] > 1e-6 && $the_row[TMIN_OLD] < $num_change[0])
			{
				print $outfile3 $the_row[TMIN_OLD]." ".$the_row[AUX1_OLD]." ".$the_row[AUX1_OLD+1]."\n";
			}
		}
		close($outfile3);
	}
}

# Save pvalues!
open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/pvalue");
foreach my $a_row (@fit_info)
{
	my @the_row = @{$a_row};
	print $outfile3 $the_row[TMIN]." ".$the_row[PVAL]."\n";
}
close($outfile3);

if ($singlefit_flag == 1)
{
	open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/singlepvalue");
	foreach my $a_row (@singlefit_info)
	{
		my @the_row = @{$a_row};
		if ($the_row[TMIN_OLD] < $num_change[0])
		{
			print $outfile3 $the_row[TMIN_OLD]." ".$the_row[PVAL_OLD]."\n";
		}
	}
	close($outfile3);
}
# I guess it's time to build some plots!

open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/plots.plt");

# Define colors. Even are for multiexp, odd are for single if they exist.
my $ptsize = 1.5;
print $outfile3 "set style line 11 lt 1 lc rgb '#000000' pt 10 ps $ptsize\n";
print $outfile3 "set style line 12 lt 1 lc rgb '#B8860B' pt 8 ps $ptsize\n";
print $outfile3 "set style line 13 lt 1 lc rgb '#CC0000' pt 6 ps $ptsize\n";
print $outfile3 "set style line 14 lt 1 lc rgb '#0000CC' pt 4 ps $ptsize\n";
print $outfile3 "set style line 15 lt 1 lc rgb '#009900' pt 12 ps $ptsize\n";
print $outfile3 "set style line 1 lt 1 lc rgb '#555555' pt 11 ps $ptsize\n";
print $outfile3 "set style line 2 lt 1 lc rgb '#D0A015' pt 9 ps $ptsize\n";
print $outfile3 "set style line 3 lt 1 lc rgb '#885555' pt 7 ps $ptsize\n";
print $outfile3 "set style line 4 lt 1 lc rgb '#5555FF' pt 5 ps $ptsize\n";
print $outfile3 "set style line 5 lt 1 lc rgb '#558855' pt 13 ps $ptsize\n";
print $outfile3 "set style line 21 lt 1 lc rgb '#000000' lw 3 pt 10 ps $ptsize\n";

# First, prepare the plot of all states.

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/allstates.tex\"\n";
print $outfile3 "set bar 3\n"; # make top and bottom of errorbars bigger.
print $outfile3 "set clip points\n"; # don't plot points on boundaries.
# Set labels.
print $outfile3 "set xlabel \"\$t_{min}\$\"\n";
print $outfile3 "set ylabel \"\$aM\$\"\n";

# Set range. For the all plot, just go from 0 to 1.
print $outfile3 "set xrange [".($tmin_min-1).":".($tmin_max+1)."]\n";
print $outfile3 "set yrange [0:1]\n";
print $outfile3 "set key top right Left reverse\n";

my $sm = 10;

# Plot some arrows.
my $iter = 1;
foreach my $tchg (@num_change)
{
	print $outfile3 "set arrow $iter from ".($tchg).",0 to ".($tchg).",1 nohead ls 21\n";
	$iter++;
}

# Prepare to plot!

print $outfile3 "plot ";

if ($num_dir > 0)
{
	for (my $i = 0; $i < $num_dir; $i++)
	{
		my $local_dir_label = $dir_label;
		# Excited states.
		if ($i > 0)
		{
			$local_dir_label = $local_dir_label."^{".("*" x $i)."}";
		}
		
		if ($i > 0)
		{
			print $outfile3 ", ";
		}
		
		print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/multifit_dir$i\" using 1:2:3 with yerrorbars ls ".($i+1)." title \"\$ ".$local_dir_label." \$\"";
	
	}
	
	#if ($singlefit_flag == 1)
	#{
	#	print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_dir\" using (\$1+0.2):2:3 with yerrorbars ls ".(1+10)." notitle";
	#}
}

if ($num_osc > 0)
{
	for (my $i = 0; $i < $num_osc; $i++)
	{
		my $local_osc_label = $osc_label;
		# Excited states.
		if ($i > 0)
		{
			$local_osc_label = $local_osc_label."^{".("*" x $i)."}";
		}
		
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/multifit_osc$i\" using 1:2:3 with yerrorbars ls ".($i+3+1)." title \"\$ ".$local_osc_label." \$\"";
	
	}
	
	#if ($singlefit_flag == 1)
	#{
	#	print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_osc\" using (\$1+0.2):2:3 with yerrorbars ls ".(3+1+10)." notitle";
	#}
}

if ($is_fpi == 1)
{
	my $local_osc_label = $osc_label;
	
	print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/multifit_osc0\" using 1:2:3 with yerrorbars ls ".(3+1)." title \"\$ ".$local_osc_label." \$\"";
	
	#if ($singlefit_flag == 1)
	#{
	#	print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_osc\" using (\$1+0.2):2:3 with yerrorbars ls ".(3+1+10)." notitle";
	#}
}

print $outfile3 "\n";
print $outfile3 "set output\n";

# Reset arrows.
$iter = 1;
foreach my $tchg (@num_change)
{
	print $outfile3 "unset arrow $iter\n";
	$iter++;
}

# Next, create zoomed in plots.

# Direct channel.

if ($num_dir > 0)
{
	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/dirstate.tex\"\n";
	print $outfile3 "set bar 3\n"; # make top and bottom of errorbars bigger.
	print $outfile3 "set clip points\n"; # don't plot points on boundaries.
	# Set labels.
	print $outfile3 "set xlabel \"\$t_{min}\$\"\n";
	print $outfile3 "set ylabel \"\$a".$dir_label."\$\"\n";

	# Set range. For the all plot, just go from 0 to 1.
	print $outfile3 "set xrange [".($tmin_min-1).":".($tmin_max+1)."]\n";
	print $outfile3 "set yrange [".($mean_dir-$sm*$err_dir).":".($mean_dir+$sm*$err_dir)."]\n";
	print $outfile3 "set key top right Left reverse\n";

	# Plot some arrows.
	$iter = 1;
	foreach my $tchg (@num_change)
	{
		print $outfile3 "set arrow $iter from ".($tchg).",".($mean_dir-$sm*$err_dir)." to ".($tchg).",".($mean_dir+$sm*$err_dir)." nohead ls 21\n";
		$iter++;
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/multifit_dir0\" using 1:2:3 with yerrorbars ls ".(1)." title \"\$ ".$dir_label." \$\"";
		
		
	if ($singlefit_flag == 1)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_dir\" using (\$1+0.2):2:3 with yerrorbars ls ".(1+10)." notitle";
	}

	print $outfile3 "\n";
	print $outfile3 "set output\n";
}

if ($num_osc > 0 || $is_fpi == 1)
{
	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/oscstate.tex\"\n";
	print $outfile3 "set bar 3\n"; # make top and bottom of errorbars bigger.
	print $outfile3 "set clip points\n"; # don't plot points on boundaries.
	# Set labels.
	print $outfile3 "set xlabel \"\$t_{min}\$\"\n";
	print $outfile3 "set ylabel \"\$a".$osc_label."\$\"\n";

	# Set range. For the all plot, just go from 0 to 1.
	print $outfile3 "set xrange [".($tmin_min-1).":".($tmin_max+1)."]\n";
	print $outfile3 "set yrange [".($mean_osc-$sm*$err_osc).":".($mean_osc+$sm*$err_osc)."]\n";
	print $outfile3 "set key top right Left reverse\n";

	# Plot some arrows.
	$iter = 1;
	foreach my $tchg (@num_change)
	{
		print $outfile3 "set arrow $iter from ".($tchg).",".($mean_osc-$sm*$err_osc)." to ".($tchg).",".($mean_osc+$sm*$err_osc)." nohead ls 21\n";
		$iter++;
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/multifit_osc0\" using 1:2:3 with yerrorbars ls ".(1+3)." title \"\$ ".$osc_label." \$\"";
		
		
	if ($singlefit_flag == 1)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_osc\" using (\$1+0.2):2:3 with yerrorbars ls ".(1+3+10)." notitle";
	}

	print $outfile3 "\n";
	print $outfile3 "set output\n";
}

# Last, create a p-value plot.

if ($num_dir > 0)
{
	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/pvalue.tex\"\n";
	print $outfile3 "set bar 3\n"; # make top and bottom of errorbars bigger.
	print $outfile3 "unset clip points\n"; # don't plot points on boundaries.
	# Set labels.
	print $outfile3 "set xlabel \"\$t_{min}\$\"\n";
	print $outfile3 "set ylabel \"p value\"\n";

	# Set range. For the all plot, just go from 0 to 1.
	print $outfile3 "set xrange [".($tmin_min-1).":".($tmin_max+1)."]\n";
	print $outfile3 "set yrange [0:1]\n";
	print $outfile3 "set key top left Left reverse\n";

	# Plot some arrows.
	$iter = 1;
	foreach my $tchg (@num_change)
	{
		print $outfile3 "set arrow $iter from ".($tchg).",0 to ".($tchg).",1 nohead ls 21\n";
		$iter++;
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/pvalue\" using 1:2 ls ".(1)." notitle";
		
		
	if ($singlefit_flag == 1)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlepvalue\" using (\$1+0.2):2 ls ".(1+10)." notitle";
	}
	
	print $outfile3 ", 0.05 ls 15 title \"5 percent\", 0.01 ls 5 title \"1 percent\"";

	print $outfile3 "\n";
	print $outfile3 "set output\n";
}


`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/allstates-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/allstates.tex`;
`convert ./allstates.pdf ./allstates.png`;

$command = "mv allstates.pdf $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state".".pdf";
`$command`;
$command = "mv allstates.png $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state".".png";
`$command`;

if ($num_dir > 0)
{
	`epstopdf ./EvanSpectrum/scripts/tmp_space/dirstate-inc.eps`;
	`pdflatex ./EvanSpectrum/scripts/tmp_space/dirstate.tex`;
	`convert ./dirstate.pdf ./dirstate.png`;
	
	$command = "mv dirstate.pdf $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-dir.pdf";
	`$command`;
	$command = "mv dirstate.png $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-dir.png";
	`$command`;
}

if ($num_osc > 0 || $is_fpi == 1)
{
	`epstopdf ./EvanSpectrum/scripts/tmp_space/oscstate-inc.eps`;
	`pdflatex ./EvanSpectrum/scripts/tmp_space/oscstate.tex`;
	`convert ./oscstate.pdf ./oscstate.png`;
	
	$command = "mv oscstate.pdf $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-osc.pdf";
	`$command`;
	$command = "mv oscstate.png $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-osc.png";
	`$command`;
}
`epstopdf ./EvanSpectrum/scripts/tmp_space/pvalue-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/pvalue.tex`;
`convert ./pvalue.pdf ./pvalue.png`;
$command = "mv pvalue.pdf $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-pv.pdf";
`$command`;
$command = "mv pvalue.png $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-pv.png";
`$command`;



`rm *.aux`;
`rm *.log`;

close($outfile3);

return;

my $ensemble = $ARGV[0];
my $state = $ARGV[1];
my @ens_splitted = split('t', $ensemble);
my @ens_splitted2 = split('b', $ens_splitted[@ens_splitted-1]);
my $nt = $ens_splitted2[0];
chomp($nt);
my $tmin_min = 0;
my $tmin_pref = 0;
my $tmin_max = 0;
my $tmax_fit = 0;

if (@ARGV == 2) # look for a fitparams file.
{
	if (-f "$path/$ensemble/spectrum2/fitparams/fitparam.$state")
	{
		open(my $fit_handle, "<$path/$ensemble/spectrum2/fitparams/fitparam.$state");
		my @fit_lines = <$fit_handle>;
		# Grab +/- 1 around the preferred value.
		my @splitted_values = split(' ', $fit_lines[0]);
		$tmin_pref = $splitted_values[1];
		$tmin_min = $tmin_pref-1;
		$tmin_max = $tmin_pref+1;
		$tmax_fit = $splitted_values[3];
		close($fit_handle);
	}
	else
	{
		die "No t ranges were specified and a fitparams file doesn't exist!\n";
	}
}
else
{
	$tmin_min = $ARGV[2];
	$tmin_pref = $ARGV[3];
	$tmin_max = $ARGV[4];
	$tmax_fit = $ARGV[5];
}

my $tmax_nontriv = 0;
if ($tmax_fit != $nt/2) # If t_max on the fit isn't the center of the lattice, we put an extra dashed line in.
{
	$tmax_nontriv = 1; 
}



# First, build a table of fit results. This also builds us up an array of masses and whatnot.
# Load the right file.

open(my $infile, "<$path/$ensemble/spectrum2/fits/fit_new.".$state);
my @inlines = <$infile>;
close($infile);

# Split and iterate through each line, building some tex!

my @textlines = ();

my @tmins = ();
my @mass_1s = ();
my @mass_1_errs = ();
my @mass_2s = (); # Gets hijacked for f_pi if need be.
my @mass_2_errs = ();
my @pvals = ();
my @amp_1s = ();
my @amp_2s = ();

my $fpi_exists = 0;
my $oscil_exists = 0;

foreach my $line (@inlines)
{
	my @splitted = split(' ', $line);
	
	my $tmin = $splitted[0];
	my $tmax = $splitted[1];
	
	my $amp1 = $splitted[2];
	my $amp1err = $splitted[3];
	my $mass1 = $splitted[4];
	my $mass1err = $splitted[5];
	my $amp2 = $splitted[6];
	my $amp2err = $splitted[7];
	my $mass2 = $splitted[8];
	my $mass2err = $splitted[9];
	
	my $amp3 = $splitted[10];
	my $amp3err = $splitted[11];
	my $mass3 = $splitted[12];
	my $mass3err = $splitted[13];
	my $amp4 = $splitted[14];
	my $amp4err = $splitted[15];
	my $mass5 = $splitted[16];
	my $mass5err = $splitted[17];
	
	my $chisqdof = $splitted[18];
	my $pval = $splitted[19];
	my $condnum = $splitted[20];
	
	my $fpi = 0; # If we have it!
	my $fpierr = 0; # If we have it!
	if (@splitted == 23)
	{
		$fpi = $splitted[21];
		$fpierr = $splitted[22];
	}

	# If the error is zero, the fit failed, either in the central value or the error.
	if (abs($mass1err) < 1e-10)
	{
		next;
	}

	
	push(@tmins, $tmin);
	push(@mass_1s, $mass1);
	push(@mass_1_errs, $mass1err);
	push(@amp_1s, $amp1);
	
	# If there's an fpi, hijack it into mass2!

	if (@splitted == 23)
	{
		push(@mass_2s, $fpi);
		push(@mass_2_errs, $fpierr);
		$fpi_exists = 1;
	}
	else
	{
		push(@mass_2s, $mass3);
		push(@mass_2_errs, $mass3err);
		push(@amp_2s, $amp3);
	}
	
	if ($fpi_exists == 0 && abs($amp3) > 1e-10)
	{
		$oscil_exists = 1;
	}
	
	push(@pvals, $pval);
	

	# Get parenthetical forms of the errors.
	
	my $amp1par = "";
	my $mass1par = "";
	my $amp2par = ""; # Gets hijacked into fpi.
	my $mass2par = "";
	
	if ($fpi_exists == 0)
	{
		$amp1par = parenForm($amp1, $amp1err);
		$mass1par = parenForm($mass1, $mass1err);
		if ($oscil_exists == 1)
		{
			$amp2par = parenForm($amp3, $amp3err);
			$mass2par = parenForm($mass3, $mass3err);
		}
	}
	else
	{
		$amp1par = parenForm($amp1, $amp1err);
		$mass1par = parenForm($mass1, $mass1err);
		$amp2par = parenForm($fpi, $fpierr);
	}
	
	my $a_line = "";
	
	if ($tmin == $tmin_pref) # Bold our preferred value.
	{
               if ($oscil_exists == 1) # cosh+oscil.
                {
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & \$$mass2par\$ & %.2f  \\\\", $tmin, $pval);
                }
                elsif ($fpi_exists == 1) # also need to print f_pi.
                {
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & %.2f  \\\\", $tmin, $pval);
                }
                else # Just a cosh.
                {
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & %.2f  \\\\", $tmin, $pval);
                }

	}
	else
	{
		if ($oscil_exists == 1) # cosh+oscil.
		{
			$a_line = sprintf("%i & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & \$$mass2par\$ & %.2f  \\\\", $tmin, $pval);
		}
		elsif ($fpi_exists == 1) # also need to print f_pi.
		{
			$a_line = sprintf("%i & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & %.2f  \\\\", $tmin, $pval);
		}
		else # Just a cosh.
		{
			$a_line = sprintf("%i & \$$amp1par\$ & \$$mass1par\$ & %.2f  \\\\", $tmin, $pval);
		}
	}
	
	push(@textlines, $a_line);

}


open (my $outfile, ">$path/$ensemble/spectrum2/tex/table_new.".$ARGV[1]);
foreach my $line (@textlines)
{
	print $outfile $line."\n\\hline\n";
}
close($outfile);


# Next, get central values.

my $mass1_mean = 0;
my $mass2_mean = 0; # fpi goes here, too, for ps2 correlator.
my $mass1_weight = 0;
my $mass2_weight = 0;
my $mass1_err = 0;
my $mass2_err = 0;
my $mass1_max = 0;
my $mass1_min = 1000;
my $mass2_max = 0;
my $mass2_min = 1000;
my $mass1_center = 0;
my $mass2_center = 0;
my $mass1_center_err = 0;
my $mass2_center_err = 0;
my $mass1_sys_err = 0;
my $mass2_sys_err = 0;
my $amp1_center = 0;
my $amp2_center = 0;
my $mass_pval = 0;
my $count = 0;

for (my $i = 0; $i < @tmins; $i++)
{
	# Make sure we don't accidentally get weird rounding on the t value.
	$tmins[$i] = floor($tmins[$i]+0.5);
	if ($tmins[$i] >= $tmin_min && $tmins[$i] <= $tmin_max)
	{
		$mass1_mean += $mass_1s[$i]/($mass_1_errs[$i]**2);
		$mass1_weight += 1.0/($mass_1_errs[$i]**2);
		$mass1_err += $mass_1_errs[$i];
		$count++;
		if ($mass1_max < $mass_1s[$i])
		{
			$mass1_max = $mass_1s[$i];
		}
		if ($mass1_min > $mass_1s[$i])
		{
			$mass1_min = $mass_1s[$i];
		}
		
		if ($tmins[$i] == $tmin_pref)
		{
			$mass1_center = $mass_1s[$i];
			$mass1_center_err = $mass_1_errs[$i];
			$mass_pval = $pvals[$i];
			$amp1_center = $amp_1s[$i];
		}
		
		if ($fpi_exists == 1 || $oscil_exists == 1)
		{
			$mass2_mean += $mass_2s[$i]/($mass_2_errs[$i]**2);
			$mass2_weight += 1.0/($mass_2_errs[$i]**2);
			$mass2_err += $mass_2_errs[$i];
			if ($mass2_max < $mass_2s[$i])
			{
				$mass2_max = $mass_2s[$i];
			}
			if ($mass2_min > $mass_2s[$i])
			{
				$mass2_min = $mass_2s[$i];
			}
			
			if ($tmins[$i] == $tmin_pref)
			{
				$mass2_center = $mass_2s[$i];
				$mass2_center_err = $mass_2_errs[$i];
				if ($fpi_exists == 0)
				{
					$amp2_center = $amp_2s[$i];
				}
			}
		}
	}
}

$mass1_mean /= $mass1_weight;
$mass1_err /= $count;

$mass1_max = abs($mass1_max - $mass1_mean);
$mass1_min = abs($mass1_mean - $mass1_min);
$mass1_sys_err = ($mass1_max > $mass1_min) ? $mass1_max : $mass1_min;

if ($fpi_exists == 1 || $oscil_exists == 1)
{
	$mass2_mean /= $mass2_weight;
	$mass2_err /= $count;

	$mass2_max = abs($mass2_max - $mass2_mean);
	$mass2_min = abs($mass2_mean - $mass2_min);
	$mass2_sys_err = ($mass2_max > $mass2_min) ? $mass2_max : $mass2_min;
}

# Print some errors!
print "Avg Mass1 value: ".(parenFormTwo($mass1_mean, $mass1_err, $mass1_sys_err))."\n";
print "Cen Mass1 value: ".(parenForm($mass1_center, $mass1_center_err))."\n";

if ($fpi_exists == 1 || $oscil_exists == 1)
{
	print "Avg Mass2 value: ".(parenFormTwo($mass2_mean, $mass2_err, $mass2_sys_err))."\n";
	print "Cen Mass2 value: ".(parenForm($mass2_center, $mass2_center_err))."\n";
}

$mass_pval = sprintf("%.2f", $mass_pval);

open (my $an_outfile, ">$path/$ensemble/spectrum2/tex/central_new.".$ARGV[1]);
if ($fpi_exists == 1 || $oscil_exists == 1)
{
	my $output_str = $tmin_min." to ".$tmin_max." & ".(parenFormTwo($mass1_mean, $mass1_err, $mass1_sys_err))." & ".(parenFormTwo($mass2_mean, $mass2_err, $mass2_sys_err))." & ".($tmin_pref)." & ".(parenForm($mass1_center, $mass1_center_err))." & ".(parenForm($mass2_center, $mass2_center_err))." & $mass_pval\n";
	print $an_outfile $output_str;
}
else
{
	print $an_outfile "$tmin_min to $tmin_max & ".(parenFormTwo($mass1_mean, $mass1_err, $mass1_sys_err))." & ".($tmin_pref)." & ".(parenForm($mass1_center, $mass1_center_err))." & $mass_pval\n";
}
close($an_outfile);

# Spit out some values the big table script can use.

open($an_outfile, ">$path/$ensemble/spectrum2/central/central_new.".$ARGV[1]);
if ($fpi_exists == 1 || $oscil_exists == 1)
{
	print $an_outfile "$mass1_center $mass1_center_err\n";
	print $an_outfile "$mass2_center $mass2_center_err\n";
}
else
{
	print $an_outfile "$mass1_center $mass1_center_err\n";
}
close($an_outfile);

# Now create an effective mass plot. This is new!

# First, grab the effective mass files and put them locally.

$command = "cp $path/$ensemble/spectrum2/effmass/effmass.$state"."* ./EvanSpectrum/scripts/tmp_space";
#print $command."\n";
`$command`;



open(my $outfile3, ">EvanSpectrum/scripts/tmp_space/plots.plt");
print $outfile3 "set style line 1 lt 1 lc rgb '#0000FF' pt 6 ps 1.5\n";
print $outfile3 "set style line 2 lt 1 lc rgb '#D70A52' pt 8 ps 1.5\n";
print $outfile3 "set style line 3 lt 1 lc rgb '#228B22' pt 4 ps 1.5\n";
print $outfile3 "set style line 4 lt 2 lc rgb 'black' lw 3 ps 1.5\n";
print $outfile3 "set style line 5 lt 1 lc rgb '#00008B' lw 3 pt 6 ps 1.5\n";
print $outfile3 "set style line 6 lt 1 lc rgb 'black' lw 4 ps 1.5\n";
print $outfile3 "set style line 7 lt 1 lc rgb '#228B22' lw 2 pt 6 ps 1.5\n";
print $outfile3 "set style line 8 lt 1 lc rgb '#DAA520' lw 2 pt 8 ps 1.5\n";
print $outfile3 "set style line 9 lt 1 lc rgb '#77EB77' lw 2 pt 6 ps 1.5\n";
print $outfile3 "set style line 10 lt 1 lc rgb '#FCCF70' lw 2 pt 8 ps 1.5\n";


# Effective Mass

# We have a new ymin/ymax---it's the preferred central value +/- 5 standard errors. We make one plot for each state.

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Time-dependent Mass\"\n";
print $outfile3 "set xrange [-0.4:".($nt/2+0.4)."]\n";
print $outfile3 "set yrange [0:1]\n";
#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",0 to ".($tmin_pref-0.5).",1 nohead ls 4\n";
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass1.tex\"\n";

my $sm = 20;

# Oscillating state, if it exists.
if ($oscil_exists == 1)
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmax_fit+0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	}
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.$state"."2\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Effective Mass\"\n";

	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmax_fit+0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	}
	print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
	print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Oscil Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.$state"."1\" using 1:2:3 with yerrorbars ls 2 title \"Oscil Effective Mass\"\n";
}
else # Non-oscil
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmax_fit+0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	}
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Effective Mass\"\n";
}
print $outfile3 "set output\n";
#print $outfile3 "set terminal x11\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass1-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass1.tex`;
`convert ./effmass1.pdf ./effmass1.png`;
if ($oscil_exists == 1)
{
	`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass2-inc.eps`;
	`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass2.tex`;
	`convert ./effmass2.pdf ./effmass2.png`;
}
#`cd ..`;

$command = "mv effmass1.pdf $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."1.pdf";
`$command`;
$command = "mv effmass1.png $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."1.png";
`$command`;

if ($oscil_exists == 1)
{
	$command = "mv effmass2.pdf $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."2.pdf";
	`$command`;
	$command = "mv effmass2.png $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."2.png";
	`$command`;
}

`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;

#Now create a fit mass plot. This is new!


$command = "cp $path/$ensemble/spectrum2/fits/fit.$state"."* ./EvanSpectrum/scripts/tmp_space";
#print $command."\n";
`$command`;



open($outfile3, ">EvanSpectrum/scripts/tmp_space/plots.plt");
print $outfile3 "set style line 1 lt 1 lc rgb '#0000FF' pt 6 ps 1.5\n";
print $outfile3 "set style line 2 lt 1 lc rgb '#D70A52' pt 8 ps 1.5\n";
print $outfile3 "set style line 3 lt 1 lc rgb '#228B22' pt 4 ps 1.5\n";
print $outfile3 "set style line 4 lt 2 lc rgb 'black' lw 3 ps 1.5\n";
print $outfile3 "set style line 5 lt 1 lc rgb '#00008B' lw 3 pt 6 ps 1.5\n";
print $outfile3 "set style line 6 lt 1 lc rgb 'black' lw 4 ps 1.5\n";
print $outfile3 "set style line 7 lt 1 lc rgb '#228B22' lw 2 pt 6 ps 1.5\n";
print $outfile3 "set style line 8 lt 1 lc rgb '#DAA520' lw 2 pt 8 ps 1.5\n";
print $outfile3 "set style line 9 lt 1 lc rgb '#77EB77' lw 2 pt 6 ps 1.5\n";
print $outfile3 "set style line 10 lt 1 lc rgb '#FCCF70' lw 2 pt 8 ps 1.5\n";


# Effective Mass

# We have a new ymin/ymax---it's the preferred central value +/- 5 standard errors. We make one plot for each state.

print $outfile3 "set xlabel \"\$t_{min}\$\"\n";
print $outfile3 "set ylabel \"Non-Linear Fit Mass\"\n";
print $outfile3 "set xrange [-0.4:".($nt/2+0.4)."]\n";
print $outfile3 "set yrange [0:1]\n";
#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",0 to ".($tmin_pref-0.5).",1 nohead ls 4\n";
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass1.tex\"\n";

# Oscillating state, if it exists.
if ($oscil_exists == 1)
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmax_fit+0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	}
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Preferred Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Non-Linear Fit Mass\"\n";

	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmax_fit+0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	}
	print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
	print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Preferred Oscil Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."2\" using 1:2:3 with yerrorbars ls 2 title \"Oscil Non-Linear Fit Mass\"\n";
}
elsif ($fpi_exists)
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
        print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
        if ($tmax_nontriv == 1)
		{
			print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmax_fit+0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
		}
		print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Preferred Cosh Non-linear Fit Mass\", \\\n";
        print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Non-Linear Fit Mass\"\n";

	print $outfile3 "set ylabel \"Non-Linear Fit \$f_{\\\\pi}\$\"\n";
        print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
        print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
		if ($tmax_nontriv == 1)
		{
			print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmax_fit+0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
		}
        print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
        print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Preferred \$f_{\\\\pi}\$ Non-linear Fit Value\", \\\n";
        print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."2\" using 1:2:3 with yerrorbars ls 2 title \"\$f_{\\\\pi}\$ Non-Linear Fit Value\"\n";

}
else # Non-oscil
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	if ($tmax_nontriv == 1)
	{
		print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmax_fit+0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	}
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Preferred Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Non-Linear Fit Mass\"\n";
}
print $outfile3 "set output\n";
#print $outfile3 "set terminal x11\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass1-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass1.tex`;
`convert ./effmass1.pdf ./effmass1.png`;
if ($oscil_exists == 1 || $fpi_exists == 1)
{
	`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass2-inc.eps`;
	`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass2.tex`;
	`convert ./effmass2.pdf ./effmass2.png`;
}
#`cd ..`;

$command = "mv effmass1.pdf $path/$ensemble/spectrum2/fits/plots/fit-$state"."1.pdf";
`$command`;
$command = "mv effmass1.png $path/$ensemble/spectrum2/fits/plots/fit-$state"."1.png";
`$command`;

if ($oscil_exists == 1 || $fpi_exists == 1)
{
	$command = "mv effmass2.pdf $path/$ensemble/spectrum2/fits/plots/fit-$state"."2.pdf";
	`$command`;
	$command = "mv effmass2.png $path/$ensemble/spectrum2/fits/plots/fit-$state"."2.png";
	`$command`;
}

`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;


# Correlator plot... ehhh...
# Now just a log plot!

my $tmp_state = $state;
if ($state eq "dc_stoch_oscil")
{
	$state = "dc_stoch";
}

# First, load the sum file.
open (my $corr_handle, "<$path/$ensemble/spectrum2/sum/sum.$state");
my @corr_contents = <$corr_handle>;
close($corr_handle);

$state = $tmp_state;

# Split things into tmin/tmax/whatnot.
my @corr_t = ();
my @corr_val = ();
my @corr_err = ();
my $max_corr_val = 0;
my $min_corr_val = 1e100;

foreach my $corr_line (@corr_contents)
{
	my @temp_splitter = split(' ', $corr_line);
	push(@corr_t, $temp_splitter[0]);
	push(@corr_val, $temp_splitter[1]);
	push(@corr_err, $temp_splitter[2]);
	
	if (abs($temp_splitter[1])+$temp_splitter[2] > $max_corr_val)
	{
		$max_corr_val = abs($temp_splitter[1])+$temp_splitter[2];
	}
	
	if (abs(abs($temp_splitter[1])-$temp_splitter[2]) < $min_corr_val && abs($temp_splitter[1]) > 1e-9)
	{
		$min_corr_val = abs(abs($temp_splitter[1])-$temp_splitter[2]);
	}
	
}



# Write the correlator out to file. Maybe this is redundant now!

open(my $outcorr_handle, ">./EvanSpectrum/scripts/tmp_space/asinhcorr.dat");
for (my $i = 0; $i < @corr_t; $i++)
{
	if ($state eq "sc_stoch") # Because it's all negative
	{
		print $outcorr_handle "$corr_t[$i] ".(-1*$corr_val[$i])." $corr_err[$i]\n";
	}
	else
	{
		print $outcorr_handle "$corr_t[$i] $corr_val[$i] $corr_err[$i]\n";
	}
}
close($outcorr_handle);

# Make a plot of the correlator!

open($outfile3, ">EvanSpectrum/scripts/tmp_space/plots.plt");
print $outfile3 "set style line 1 lc rgb 'blue' pt 6 ps 1.5\n";
print $outfile3 "set style line 2 lt 1 lc rgb '#D70A52' lw 3 pt 8 ps 1.5\n";
print $outfile3 "set style line 3 lt 1 lc rgb '#228B22' pt 4 ps 1.5\n";
print $outfile3 "set style line 4 lt 2 lc rgb 'black' lw 3 ps 1.5\n";
print $outfile3 "set style line 5 lt 1 lc rgb '#00008B' lw 3 pt 6 ps 1.5\n";
print $outfile3 "set style line 6 lt 1 lc rgb 'black' lw 4 ps 1.5\n";

my $maxt = ($nt/2+0.4);
if ($state eq "nu" || $state eq "de")
{
	$maxt = $nt;
}

# Alrighty. This sucks. We need to plot the asinh form of the fit curve correlator.

# Clean up our cheating with fpi.
if ($fpi_exists)
{
	$mass2_center = 0;
	$mass2_center_err = 0;
}

my $curve_line = "";

if ($state eq "nu" || $state eq "de")
{
	# Nucleons are special.
	# Obnoxious, too.
	$curve_line = $amp1_center."*(exp(-".$mass1_center."*x)-cos(3.1415926535*x)*exp(-".$mass1_center."*(".$nt."-x))) - ".$amp2_center."*cos(3.1415926535*x)*(exp(-".$mass2_center."*x)-cos(3.1415926535*x)*exp(-".$mass2_center."*(".$nt."-x)))";
}
elsif ($state eq "sc_stoch")
{
	$curve_line = "-".$amp1_center."*cosh(".$mass1_center."*(".$nt."/2-x)) - ".$amp2_center."*cos(3.1415926535*x)*cosh(".$mass2_center."*(".$nt."/2-x))";
}
else
{
	$curve_line = $amp1_center."*cosh(".$mass1_center."*(".$nt."/2-x)) + ".$amp2_center."*cos(3.1415926535*x)*cosh(".$mass2_center."*(".$nt."/2-x))";
}

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Log Plot of Correlator\"\n";
print $outfile3 "set xrange [-0.4:".$maxt."]\n";
print $outfile3 "set log y\n";
print $outfile3 "set format y \"%2.0e\"\n";
#print $outfile3 "set format y \"%2.0tx10^{%L}\"\n";
print $outfile3 "set yrange [".(0.9*$min_corr_val).":".(1.1*$max_corr_val)."]\n";
#print $outfile3 "set yrange [-".($mass1_mean*30).":".($mass1_mean*30)."]\n";
print $outfile3 "set samples 1000\n";
print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".(0.9*$min_corr_val)." to ".($tmin_pref-0.5).",".(1.1*$max_corr_val)." nohead ls 4\n";
if ($state eq "nu" || $state eq "de")
{
	print $outfile3 "set arrow 2 from ".($nt-($tmin_pref-0.5)).",".(0.9*$min_corr_val)." to ".($nt-($tmin_pref-0.5)).",".(1.1*$max_corr_val)." nohead ls 4\n";
}
if ($state eq "sg_111" || $state eq "sg_211" || $state eq "dc_stoch")
{
	print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".(0.9*$min_corr_val)." to ".($tmax_fit+0.5).",".(1.1*$max_corr_val)." nohead ls 4\n";
}
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/asinhcorr.tex\"\n";
#if ($state eq "sc_stoch")
#{
#	print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/asinhcorr.dat\" using 1:(-\$2):3 with yerrorbars ls 1 notitle, \\\n";
#}
#else
#{
	print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/asinhcorr.dat\" using 1:2:3 with yerrorbars ls 1 notitle, \\\n";
#}
print $outfile3 "   $curve_line ls 2 notitle\n";
print $outfile3 "set output\n";
#print $outfile3 "set terminal wxt\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/asinhcorr-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/asinhcorr.tex`;
`convert ./asinhcorr.pdf ./asinhcorr.png`;

#`cd ..`;

$command = "mv asinhcorr.pdf $path/$ensemble/spectrum2/sum/plots/sum-$state.pdf";
`$command`;
$command = "mv asinhcorr.png $path/$ensemble/spectrum2/sum/plots/sum-$state.png";
`$command`;


`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;


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

sub parenFormTwo
{
   my ($the_value, $the_error_1, $the_error_2) = @_;
   
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
   
   my $the_error = 0;
   my $other_error = 0;
   my $error_no = 0;
   
   if ($the_error_1 > $the_error_2)
   {
		$the_error = $the_error_1;
		$other_error = $the_error_2;
		$error_no = 1;
	}
	else
	{
		$the_error = $the_error_2;
		$other_error = $the_error_1;
		$error_no = 2;
	}

   my $err_size = floor(log($the_error)/log(10));
   my $other_err_size = floor(log($other_error)/log(10));

   my $value_size = floor(log(abs($the_value))/log(10));

   $the_error *= 10**(-$err_size+1);
   $other_error *= 10**(-$err_size+1);

        if (floor($the_error + 0.5) == 100) # Special case.
        {
                $the_error/=10;
				$other_error/=10;
                $err_size--;
        }
		
		if (floor($other_error + 0.5) == 100)
		{
			$other_error /= 10;
			$other_err_size--;
		}
	
	

   $the_value *= 10**(-$err_size+1);


   $the_value = floor($the_value + 0.5);
   $the_error = floor($the_error + 0.5);
   $other_error = floor($other_error + 0.5);
   
   if ($other_error < 10)
   {
		$other_error = "0".$other_error;
	}
	

   #print $the_value." ".$the_error."\n";
   
   if ($error_no == 1)
   {
		$the_error_1 = $the_error;
		$the_error_2 = $other_error;
   }
   else
   {
		$the_error_2 = $the_error;
		$the_error_1 = $other_error;
	}

   if ($neg_flag == -1)
   {
		return '-'.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error_1.")(".$the_error_2.")\\mbox{e".$value_size.'} ';
	}
	else
	{
		return ''.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error_1.")(".$the_error_2.")\\mbox{e".$value_size.'} ';
	}
   


        #return 0;

}

sub parenFormFake
{
   my ($the_value) = @_;
   
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
   
   my $the_error = $the_value/100;

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
		return '-'.substr($the_value, 0, 1).".".substr($the_value, 1)."\\mbox{e".$value_size.'} ';
	}
	else
	{
		return ''.substr($the_value, 0, 1).".".substr($the_value, 1)."\\mbox{e".$value_size.'} ';
	}
   


        #return 0;

}

