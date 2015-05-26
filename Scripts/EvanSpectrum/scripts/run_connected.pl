#! /usr/bin/perl -w

use POSIX;
use strict;
use Switch;

# Call as:
# ./run_connected.pl [ensemble] [state] {tmin-min} {tmin-pref} {tmin-max} {tmax}
# Last params are if there is no fitparam file/if we want to ignore it.

# Let's grab all the data we can!

if (@ARGV != 6 && @ARGV != 2)
{
	die "Need a [ensemble] [state] {tmin-min} {tmin-pref} {tmin-max} {tmax}.\n";
}

my $path = "/projectnb/qcd/Staggered/Scripts";
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

my $temp_state = $state;
if ($state eq "dc_stoch_oscil" || $state eq "dc_stoch_deriv" || $state eq "dc_stoch_wconst")
{
	$state = "dc_stoch";
}

if ($state eq "sg_stoch" || $state eq "sg_stoch_deriv" || $state eq "sg_stoch_wconst")
{
	$state = "sg_stoch";
}

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

$state = $temp_state;

my $tmax_nontriv = 0;
if ($tmax_fit != $nt/2) # If t_max on the fit isn't the center of the lattice, we put an extra dashed line in.
{
	$tmax_nontriv = 1; 
}

# Prepare to run some shell commands later.
my $command = "";

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
my @const_1s = ();
my @const_1_errs = ();

my $fpi_exists = 0;
my $oscil_exists = 0;
my $const_exists = 0;

my $mass1_center = 0;
my $mass2_center = 0;
my $mass1_center_err = 0;
my $mass2_center_err = 0;
my $amp1_center = 0;
my $amp2_center = 0;
my $const_center = 0;
my $const_center_err = 0;

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
	# 2015-01-12: Update, we may just be doing a no errors run.
	#if (abs($mass1err) < 1e-10)
	#{
	#	next;
	#}

	
	push(@tmins, $tmin);
	push(@mass_1s, $mass1);
	push(@mass_1_errs, $mass1err);
	push(@amp_1s, $amp1);
	push(@const_1s, $amp2);
	push(@const_1_errs, $amp2err); 
	
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
	
	if ($state eq "dc_stoch_wconst" || $state eq "sg_stoch_wconst")
	{
		$const_exists = 1;
	}
	
	push(@pvals, $pval);
	

	# Get parenthetical forms of the errors.
	
	my $amp1par = "";
	my $mass1par = "";
	my $amp2par = ""; # Gets hijacked into fpi.
	my $mass2par = "";
	my $constpar = "";
	
	if ($fpi_exists == 0)
	{
		$amp1par = parenForm($amp1, $amp1err);
		$mass1par = parenForm($mass1, $mass1err);
		if ($oscil_exists == 1)
		{
			$amp2par = parenForm($amp3, $amp3err);
			$mass2par = parenForm($mass3, $mass3err);
		}
		if ($const_exists == 1)
		{
			$constpar = parenForm($amp2, $amp2err);
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
				if ($const_exists == 1)
				{
					$a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & \$$mass2par\$ & \$$constpar\$ & %.2f  \\\\", $tmin, $pval);
					$const_center = $amp2;
					$const_center_err = $amp2err;
					$amp1_center = $amp1;
						$mass1_center = $mass1;
						$mass1_center_err = $mass1err;
						$amp2_center = $amp3;
						$mass2_center = $mass3;
						$mass2_center_err = $mass3err;
				}
				elsif ($oscil_exists == 1) # cosh+oscil.
                {
					
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & \$$mass2par\$ & %.2f  \\\\", $tmin, $pval);
						$amp1_center = $amp1;
						$mass1_center = $mass1;
						$mass1_center_err = $mass1err;
						$amp2_center = $amp3;
						$mass2_center = $mass3;
						$mass2_center_err = $mass3err;
                }
                elsif ($fpi_exists == 1) # also need to print f_pi.
                {
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & %.2f  \\\\", $tmin, $pval);
						$amp1_center = $amp1;
						$mass1_center = $mass1;
						$mass1_center_err = $mass1err;
						$mass2_center = $fpi;
						$mass2_center_err = $fpierr;
                }
                else # Just a cosh.
                {
                        $a_line = sprintf("\\textbf{%i} & \$$amp1par\$ & \$$mass1par\$ & %.2f  \\\\", $tmin, $pval);
						$amp1_center = $amp1;
						$mass1_center = $mass1;
						$mass1_center_err = $mass1err;
                }

	}
	else
	{
		if ($const_exists == 1)
		{
			$a_line = sprintf("%i & \$$amp1par\$ & \$$mass1par\$ & \$$amp2par\$ & \$$mass2par\$ & \$$constpar\$ & %.2f  \\\\", $tmin, $pval);
		}
		elsif ($oscil_exists == 1) # cosh+oscil.
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

=begin GHOSTCODE
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
#print "Avg Mass1 value: ".(parenFormTwo($mass1_mean, $mass1_err, $mass1_sys_err))."\n";
print "Cen Mass1 value: ".(parenForm($mass1_center, $mass1_center_err))."\n";

if ($fpi_exists == 1 || $oscil_exists == 1)
{
	#print "Avg Mass2 value: ".(parenFormTwo($mass2_mean, $mass2_err, $mass2_sys_err))."\n";
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
=end GHOSTCODE
=cut


# Spit out some values the big table script can use.

open(my $an_outfile, ">$path/$ensemble/spectrum2/central/central_new.".$ARGV[1]);
if ($const_exists == 1)
{
	print $an_outfile "$mass1_center $mass1_center_err\n";
	print $an_outfile "$mass2_center $mass2_center_err\n";
	print $an_outfile "$const_center $const_center_err\n";
}
elsif ($fpi_exists == 1 || $oscil_exists == 1)
{
	print $an_outfile "$mass1_center $mass1_center_err\n";
	print $an_outfile "$mass2_center $mass2_center_err\n";
}
else
{
	print $an_outfile "$mass1_center $mass1_center_err\n";
}
close($an_outfile);

# Before making plots, let's make sure the errors aren't zero! If they are, just set them such that plots go from
# 0 to 2*the value of the mass.

my $noerr_flag = 0;

if ($mass1_center_err == 0)
{
	$mass1_center_err = abs($mass1_center)/20.0;
	$noerr_flag = 1;
}

if ($fpi_exists == 1 || $oscil_exists == 1)
{
	if ($mass2_center_err == 0)
	{
		$mass2_center_err = abs($mass2_center)/20.0;
	}
}



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

if ($noerr_flag == 1) # turn off errorbar colors
{
	print $outfile3 "set style line 9 lt 1 lc rgb '#FFFFFF' lw 2 pt 6 ps 1.5\n";
	print $outfile3 "set style line 10 lt 1 lc rgb '#FFFFFF' lw 2 pt 8 ps 1.5\n";

}


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

if ($noerr_flag == 1) # turn off errorbar colors
{
	print $outfile3 "set style line 9 lt 1 lc rgb '#FFFFFF' lw 2 pt 6 ps 1.5\n";
	print $outfile3 "set style line 10 lt 1 lc rgb '#FFFFFF' lw 2 pt 8 ps 1.5\n";

}

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
if ($state eq "dc_stoch_oscil" || $state eq "dc_stoch_deriv" || $state eq "dc_stoch_wconst")
{
	$state = "dc_stoch";
}
if ($state eq "sg_stoch" || $state eq "sg_stoch_deriv" || $state eq "sg_stoch_wconst")
{
	$state = "sg_stoch";
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
elsif ($state eq "dc_stoch_wconst" || $state eq "sg_stoch_wconst")
{
	$curve_line = $amp1_center."*cosh(".$mass1_center."*(".$nt."/2-x)) + ".$amp2_center."*cos(3.1415926535*x)*cosh(".$mass2_center."*(".$nt."/2-x)) + ".$const_center;
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


=begin GHOSTCODE
# We use Claudio's arcsinh method!

# First, load the sum file.
open (my $corr_handle, "<../ensembles/$path"."$ensemble/sum/sum.$state");
my @corr_contents = <$corr_handle>;
close($corr_handle);

# Next, we need to figure out the asinh scale factor.
# This requires us to find the maximum absolute value of the correlator,
# add the error, then add some free room.
# Split things into tmin/tmax/whatnot.
my @corr_t = ();
my @corr_val = ();
my @corr_err = ();
my $max_corr_val = 0;

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
}

# We want to solve, somewhat arbitrarily, 4 = asinh(sc*$max_corr_val*1.1).
# The 1.1 gives a bit of wiggle room.
my $sc_factor = sinh($mass1_mean*30)/($max_corr_val*1.1);

my @corr_rescale = ();
my @corr_rescale_uperr = ();
my @corr_rescale_downerr = ();

for (my $i = 0; $i < @corr_t; $i++)
{
	my $asinh_val = my_asinh($corr_val[$i]*$sc_factor);
	my $asinh_up = my_asinh(($corr_val[$i]+$corr_err[$i])*$sc_factor);
	my $asinh_down = my_asinh(($corr_val[$i]-$corr_err[$i])*$sc_factor);
	
	push(@corr_rescale, $asinh_val);
	push(@corr_rescale_uperr, $asinh_up - $asinh_val);
	push(@corr_rescale_downerr, -$asinh_down + $asinh_val);
}

# Write the correlator out to file.

open(my $outcorr_handle, ">./EvanSpectrum/scripts/tmp_space/asinhcorr.dat");
for (my $i = 0; $i < @corr_t; $i++)
{
	print $outcorr_handle "$corr_t[$i] $corr_rescale[$i] $corr_rescale_downerr[$i] $corr_rescale_uperr[$i]\n";
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
	$curve_line = "asinh(".$sc_factor."*(".$amp1_center."*(exp(-".$mass1_center."*x)-cos(3.1415926535*x)*exp(-".$mass1_center."*(".$nt."-x))) - ".$amp2_center."*cos(3.1415926535*x)*(exp(-".$mass2_center."*x)-cos(3.1415926535*x)*exp(-".$mass2_center."*(".$nt."-x)))))";
}
else
{
	$curve_line = "asinh(".$sc_factor."*(".$amp1_center."*cosh(".$mass1_center."*(".$nt."/2-x)) + ".$amp2_center."*cos(3.1415926535*x)*cosh(".$mass2_center."*(".$nt."/2-x))))";
}

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Weighted Arcsinh of Correlator\"\n";
print $outfile3 "set xrange [-0.4:".$maxt."]\n";
print $outfile3 "set yrange [-".($mass1_mean*30).":".($mass1_mean*30)."]\n";
print $outfile3 "set samples 1000\n";
print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".(-$mass1_mean*30)." to ".($tmin_pref-0.5).",".($mass1_mean*30)." nohead ls 4\n";
if ($state eq "nu" || $state eq "de")
{
	print $outfile3 "set arrow 2 from ".($nt-($tmin_pref-0.5)).",".(-$mass1_mean*30)." to ".($nt-($tmin_pref-0.5)).",".($mass1_mean*30)." nohead ls 4\n";
}
if ($state eq "sg_111" || $state eq "sg_211")
{
	print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".(-$mass1_mean*30)." to ".($tmax_fit+0.5).",".($mass1_mean*30)." nohead ls 4\n";
}
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/asinhcorr.tex\"\n";
print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/asinhcorr.dat\" using 1:2:(\$2-\$3):(\$2+\$4) with yerrorbars ls 1 notitle, \\\n";
print $outfile3 "   $curve_line ls 2 notitle\n";
print $outfile3 "set output\n";
print $outfile3 "set terminal wxt\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/asinhcorr-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/asinhcorr.tex`;
`convert ./asinhcorr.pdf ./asinhcorr.png`;

#`cd ..`;

$command = "mv asinhcorr.pdf ../ensembles/$path"."$ensemble/sum/plots/sum-$state.pdf";
`$command`;
$command = "mv asinhcorr.png ../ensembles/$path"."$ensemble/sum/plots/sum-$state.png";
`$command`;


`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;

sub my_asinh
{
	if (abs($_[0]+sqrt(1.0+$_[0]**2)) > 1e-8)
	{
		return log($_[0]+sqrt(1.0+$_[0]**2));
	}
	else
	{
		return -1000000;
	}
}

=end GHOSTCODE
=cut

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
   
   # Check and see if the error is zero. If it is, we set a dummy error to round it to 2 digits.
   if (abs($the_error) < 1e-20)
   {
	 return parenFormFake($the_value);
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

