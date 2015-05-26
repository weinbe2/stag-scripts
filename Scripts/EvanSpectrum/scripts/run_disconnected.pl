#! /usr/bin/perl -w

use POSIX;
use strict;
use Switch;

# Call as:
# ./run_connected.pl [ensemble] {tmin-min} {tmin-pref} {tmin-max} {tmax}
# Last params are if there is no fitparam file/if we want to ignore it.

# Let's grab all the data we can!

if (@ARGV != 5 && @ARGV != 1)
{
	die "Need a [ensemble] {tmin-min} {tmin-pref} {tmin-max} {tmax}.\n";
}

my $path = "/projectnb/qcd/Staggered/Scripts";
my $ensemble = $ARGV[0];
my $state = "sg_211"; #$ARGV[1];
my @ens_splitted = split('t', $ensemble);
my @ens_splitted2 = split('b', $ens_splitted[@ens_splitted-1]);
my $nt = $ens_splitted2[0];
chomp($nt);
my $tmin_min = 0;
my $tmin_pref = 0;
my $tmin_max = 0;
my $tmax_fit = 0;

if (@ARGV == 1) # look for a fitparams file.
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
	$tmin_min = $ARGV[1];
	$tmin_pref = $ARGV[2];
	$tmin_max = $ARGV[3];
	$tmax_fit = $ARGV[4];
}


# Prepare to run some shell commands later.
my $command = "";

# First, do some things with the vev.

# We're being a little lazy here---have this create the vev folder.

if (!(-d "$path/$ensemble/spectrum2/vev/plots"))
{
	$command = "mkdir $path/$ensemble/spectrum2/vev/plots";
	`$command`;
}

# Load the right file. Do both at once.
open(my $pbpfile, "<$path/$ensemble/spectrum2/vev/history.vev1"); #pbp
my @pbp_lines = <$pbpfile>;
close($pbpfile);

open(my $vevfile, "<$path/$ensemble/spectrum2/vev/history.vev2"); #vev
my @vev_lines = <$vevfile>;
close($vevfile);

# Split and interate through each line, building some tex.
my @texlines = ();
my $pbp_center = 0;
my $pbp_err = 0;
my $vev_center = 0;
my $vev_err = 0;
my $max_cfg = 0;
my $sm = 20;

my $centertexline = "";
for (my $i = 0; $i < @pbp_lines; $i++)
{
	my @splitted_pbp = split(' ', $pbp_lines[$i]);
	my @splitted_vev = split(' ', $vev_lines[$i]);
	
	my $vevpar = parenForm($splitted_pbp[1], $splitted_pbp[2]);
	my $pbppar = parenForm($splitted_vev[1], $splitted_vev[2]);
	
	# Build a line of tex.
	my $a_line = "\$$splitted_pbp[0]\$ & \$$pbppar\$ & \$$vevpar\$ \\\\";
	
	push(@texlines, $a_line);

	if (($i+1) == @pbp_lines)
	{
		$centertexline = "\$$splitted_pbp[0]\$ & \$$pbppar\$ & \$$vevpar\$ \n";
		$pbp_center = $splitted_pbp[1];
		$pbp_err = $splitted_pbp[2];
		$vev_center = $splitted_vev[1];
		$vev_err = $splitted_vev[2];
		$max_cfg = $splitted_pbp[0];
	}
	
}

open(my $vevoutfile, ">$path/$ensemble/spectrum2/tex/table_new.vev");
foreach my $line (@texlines)
{
	print $vevoutfile, $line."\n\\hline\n";
}
close($vevoutfile);

# Print some nice full stats overview out.

open (my $an_outfile, ">$path/$ensemble/spectrum2/tex/central_new.vev");
print $an_outfile $centertexline;
close($an_outfile);

# And something nice the bib table script can use.

open($an_outfile, ">$path/$ensemble/spectrum2/central/central_new.vev");
print $an_outfile "$pbp_center $pbp_err\n";
print $an_outfile "$vev_center $vev_err\n";
close($an_outfile);

# Plot the vev!

# First, grab it put them locally.
$command = "cp $path/$ensemble/spectrum2/vev/history.vev"."* ./EvanSpectrum/scripts/tmp_space";
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


print $outfile3 "set xlabel \"Number of Blocks\"\n";
print $outfile3 "set ylabel \"Value of \$\\\\langle \\\\bar{\\\\psi} \\\\psi \\\\rangle(t)\$\"\n";
print $outfile3 "set xrange [".(4.0/sqrt(2.0)).":".($max_cfg*sqrt(2.0))."]\n";
print $outfile3 "set log x\n";
print $outfile3 "set yrange [0:1]\n";
#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",0 to ".($tmin_pref-0.5).",1 nohead ls 4\n";
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass1.tex\"\n";

print $outfile3 "set yrange [".($pbp_center-$sm*$pbp_err).":".($pbp_center+$sm*$pbp_err)."]\n";
print $outfile3 "plot ".($pbp_center+$pbp_err)." ls 9 notitle, ".($pbp_center-$pbp_err)." ls 9 notitle, $pbp_center ls 7 title \"Preferred \$\\\\langle \\\\bar{\\\\psi} \\\\psi \\\\rangle(t)\$\", \\\n";
print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/history.vev"."1\" using 1:2:3 with yerrorbars ls 1 title \"History of \$\\\\langle \\\\bar{\\\\psi} \\\\psi \\\\rangle(t)\$\"\n";

print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
print $outfile3 "set yrange [".($vev_center-$sm*$vev_err).":".($vev_center+$sm*$vev_err)."]\n";
print $outfile3 "plot ".($vev_center+$vev_err)." ls 9 notitle, ".($vev_center-$vev_err)." ls 9 notitle, $vev_center ls 7 title \"Preferred VEV\", \\\n";
print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/history.vev"."2\" using 1:2:3 with yerrorbars ls 1 title \"History of VEV\"\n";


print $outfile3 "set output\n";
#print $outfile3 "set terminal x11\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass1-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass1.tex`;
`convert ./effmass1.pdf ./effmass1.png`;

`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass2-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass2.tex`;
`convert ./effmass2.pdf ./effmass2.png`;


$command = "mv effmass1.pdf $path/$ensemble/spectrum2/vev/plots/history-vev1.pdf";
`$command`;
$command = "mv effmass1.png $path/$ensemble/spectrum2/vev/plots/history-vev1.png";
`$command`;

$command = "mv effmass2.pdf $path/$ensemble/spectrum2/vev/plots/history-vev2.pdf";
`$command`;
$command = "mv effmass2.png $path/$ensemble/spectrum2/vev/plots/history-vev2.png";
`$command`;


`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;



# Now create an effective mass plot. This is new! Disconnected, specifically.

# First, grab the effective mass files and put them locally.
# This also grabs the new kuti and kmi effmasses.
$command = "cp $path/$ensemble/spectrum2/effmass/effmass.dc_stoch"."* ./EvanSpectrum/scripts/tmp_space";
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

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Time-dependent Mass\"\n";
print $outfile3 "set xrange [-0.4:".($nt/2+0.4)."]\n";
print $outfile3 "set yrange [0:1]\n";
#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",0 to ".($tmin_pref-0.5).",1 nohead ls 4\n";
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass1.tex\"\n";


=begin GHOSTCODE
# Oscillating state, if it exists.
if ($oscil_exists == 1)
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.$state"."2\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Effective Mass\"\n";

	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
	print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Oscil Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.$state"."1\" using 1:2:3 with yerrorbars ls 2 title \"Oscil Effective Mass\"\n";
}
else # Non-oscil
{
=end GHOSTCODE
=cut
	print $outfile3 "set yrange [0:1]\n";
	#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
	#print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/effmass.dc_stoch"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Effective Mass\"\n";
	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass_spec1.tex\"\n";
# Kuti and KMI effmasses.
print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/effmass.dc_stoch_kuti"."1\" using 1:2:3 with yerrorbars ls 1 title \"Kuti Effective Mass\", \\\n";
print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/effmass.dc_stoch_kmi"."1\" using 1:2:3 with yerrorbars ls 2 title \"KMI Effective Mass\"\n";

print $outfile3 "set output\n";
#print $outfile3 "set terminal x11\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass1-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass1.tex`;
`convert ./effmass1.pdf ./effmass1.png`;

`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass_spec1-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass_spec1.tex`;
`convert ./effmass_spec1.pdf ./effmass_spec1.png`;


=begin GHOSTCODE
if ($oscil_exists == 1)
{
	`epstopdf ./EvanSpectrum/scripts/tmp_space/effmass2-inc.eps`;
	`pdflatex ./EvanSpectrum/scripts/tmp_space/effmass2.tex`;
	`convert ./effmass2.pdf ./effmass2.png`;
}
=end GHOSTCODE
=cut
#`cd ..`;

$command = "mv effmass1.pdf $path/$ensemble/spectrum2/effmass/plots/effmass-dc_stoch"."1.pdf";
`$command`;
$command = "mv effmass1.png $path/$ensemble/spectrum2/effmass/plots/effmass-dc_stoch"."1.png";
`$command`;

$command = "mv effmass_spec1.pdf $path/$ensemble/spectrum2/effmass/plots/effmass-dc_stoch_spec.pdf";
`$command`;
$command = "mv effmass_spec1.png $path/$ensemble/spectrum2/effmass/plots/effmass-dc_stoch_spec.png";
`$command`;

=begin GHOSTCODE
if ($oscil_exists == 1)
{
	$command = "mv effmass2.pdf $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."2.pdf";
	`$command`;
	$command = "mv effmass2.png $path/$ensemble/spectrum2/effmass/plots/effmass-$state"."2.png";
	`$command`;
}
=end GHOSTCODE
=cut

`rm *.aux`;
`rm *.log`;
`rm ./EvanSpectrum/scripts/tmp_space/*`;

=begin GHOSTCODE
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
	print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Preferred Cosh Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Non-Linear Fit Mass\"\n";

	print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
	print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
	print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Preferred Oscil Non-linear Fit Mass\", \\\n";
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."2\" using 1:2:3 with yerrorbars ls 2 title \"Oscil Non-Linear Fit Mass\"\n";
}
elsif ($fpi_exists)
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
        print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
        print $outfile3 "plot ".($mass1_center+$mass1_center_err)." ls 9 notitle, ".($mass1_center-$mass1_center_err)." ls 9 notitle, $mass1_center ls 7 title \"Preferred Cosh Non-linear Fit Mass\", \\\n";
        print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."1\" using 1:2:3 with yerrorbars ls 1 title \"Cosh Non-Linear Fit Mass\"\n";

	print $outfile3 "set ylabel \"Non-Linear Fit \$f_{\\\\pi}\$\"\n";
        print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/effmass2.tex\"\n";
        print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass2_center-$sm*$mass2_center_err)." to ".($tmin_pref-0.5).",".($mass2_center+$sm*$mass2_center_err)." nohead ls 4\n";
        print $outfile3 "set yrange [".($mass2_center-$sm*$mass2_center_err).":".($mass2_center+$sm*$mass2_center_err)."]\n";
        print $outfile3 "plot ".($mass2_center+$mass2_center_err)." ls 10 notitle, ".($mass2_center-$mass2_center_err)." ls 10 notitle, $mass2_center ls 8 title \"Preferred \$f_{\\\\pi}\$ Non-linear Fit Value\", \\\n";
        print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/fit.$state"."2\" using 1:2:3 with yerrorbars ls 2 title \"\$f_{\\\\pi}\$ Non-Linear Fit Value\"\n";

}
else # Non-oscil
{
	print $outfile3 "set yrange [".($mass1_center-$sm*$mass1_center_err).":".($mass1_center+$sm*$mass1_center_err)."]\n";
	print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".($mass1_center-$sm*$mass1_center_err)." to ".($tmin_pref-0.5).",".($mass1_center+$sm*$mass1_center_err)." nohead ls 4\n";
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
=end GHOSTCODE
=cut

# Correlator plot... ehhh...
# Now just a log plot!

# First, load the sum file.
open (my $corr_handle, "<$path/$ensemble/spectrum2/sum/sum.sc_stoch");
my @corr_contents = <$corr_handle>;
close($corr_handle);


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
	
	if (abs(abs($temp_splitter[1])-$temp_splitter[2]) < $min_corr_val)
	{
		$min_corr_val = abs(abs($temp_splitter[1])-$temp_splitter[2]);
	}
	
}



# Write the correlator out to file. Maybe this is redundant now!

open(my $outcorr_handle, ">./EvanSpectrum/scripts/tmp_space/asinhcorr.dat");
for (my $i = 0; $i < @corr_t; $i++)
{
		print $outcorr_handle "$corr_t[$i] ".(-1*$corr_val[$i])." $corr_err[$i]\n";
}
close($outcorr_handle);

# Grab the disconnected correlator, too!

# First, load the sum file.
open ($corr_handle, "<$path/$ensemble/spectrum2/sum/sum.dc_stoch");
@corr_contents = <$corr_handle>;
close($corr_handle);


# Split things into tmin/tmax/whatnot.
@corr_t = ();
@corr_val = ();
@corr_err = ();
#$max_corr_val = 0;
#$min_corr_val = 1e100;

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

        if (abs(abs($temp_splitter[1])-$temp_splitter[2]) < $min_corr_val)
        {
                $min_corr_val = abs(abs($temp_splitter[1])-$temp_splitter[2]);
        }

}



# Write the correlator out to file. Maybe this is redundant now!

open($outcorr_handle, ">./EvanSpectrum/scripts/tmp_space/asinhcorr2.dat");
for (my $i = 0; $i < @corr_t; $i++)
{
                print $outcorr_handle "$corr_t[$i] $corr_val[$i] $corr_err[$i]\n";
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

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Log Plot of Correlator\"\n";
print $outfile3 "set xrange [-0.4:".$maxt."]\n";
print $outfile3 "set log y\n";
print $outfile3 "set format y \"%2.0e\"\n";
#print $outfile3 "set format y \"%2.0tx10^{%L}\"\n";
print $outfile3 "set yrange [".(0.9*$min_corr_val).":".(1.1*$max_corr_val)."]\n";
#print $outfile3 "set yrange [-".($mass1_mean*30).":".($mass1_mean*30)."]\n";
print $outfile3 "set samples 1000\n";
#print $outfile3 "set arrow 1 from ".($tmin_pref-0.5).",".(0.9*$min_corr_val)." to ".($tmin_pref-0.5).",".(1.1*$max_corr_val)." nohead ls 4\n";
#%if ($state eq "nu" || $state eq "de")
#{
#	print $outfile3 "set arrow 2 from ".($nt-($tmin_pref-0.5)).",".(0.9*$min_corr_val)." to ".($nt-($tmin_pref-0.5)).",".(1.1*$max_corr_val)." nohead ls 4\n";
#}
#if ($state eq "sg_111" || $state eq "sg_211")
#{
#	print $outfile3 "set arrow 2 from ".($tmax_fit+0.5).",".(0.9*$min_corr_val)." to ".($tmax_fit+0.5).",".(1.1*$max_corr_val)." nohead ls 4\n";
#}
#print $outfile3 "set arrow 2 from ".($minimum_incl-0.5).",0.8 to ".($minimum_incl-0.5+4).",0.8 ls 6\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/asinhcorr.tex\"\n";
#if ($state eq "sc_stoch")
#{
#	print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/asinhcorr.dat\" using 1:(-\$2):3 with yerrorbars ls 1 notitle, \\\n";
#}
#else
#{
	print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/asinhcorr.dat\" using 1:2:3 with yerrorbars ls 1 title \"Connected\", \\\n";
	print $outfile3 "   \"./EvanSpectrum/scripts/tmp_space/asinhcorr2.dat\" using 1:2:3 with yerrorbars ls 2 title \"Disconnected\"\n";
#}
print $outfile3 "set output\n";
#print $outfile3 "set terminal wxt\n";

close($outfile3);

`gnuplot -e "load \\"EvanSpectrum/scripts/tmp_space/plots.plt\\""`; 
#`cd ./EvanSpectrum/scripts/tmp_space`;
`epstopdf ./EvanSpectrum/scripts/tmp_space/asinhcorr-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/asinhcorr.tex`;
`convert ./asinhcorr.pdf ./asinhcorr.png`;

#`cd ..`;

$command = "mv asinhcorr.pdf $path/$ensemble/spectrum2/sum/plots/sum-scdc.pdf";
`$command`;
$command = "mv asinhcorr.png $path/$ensemble/spectrum2/sum/plots/sum-scdc.png";
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

