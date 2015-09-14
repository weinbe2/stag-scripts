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
	
	$tmax = $splitted[TMAX]; # This should always be the same.
	
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

if (@num_change == 0)
{
	push(@num_change, 0);
}

# Good! Print this out.

printf("Tmin_min = %d\n", $tmin_min);
printf("Tmin_max = %d\n", $tmin_max);
printf("Tmax = %d\n", $tmax);
printf("Num dir. = %d\n", $num_dir);
printf("Num osc. = %d\n", $num_osc);
printf("Num chg. at ");
foreach my $tchg (@num_change)
{
	printf("%.1f ", $tchg);
}
printf("\n");

# If there's a multifitinterval file, grab it and learn something!
my $tmin_pref = 0;
my $tmax_pref = $tmax;
my $bin_size = 0;
if (-f "$path/$relpath/$direc/spectrum2/multifitparams/multifitparam.$state")
{
	open(my $fit_handle, "<$path/$relpath/$direc/spectrum2/multifitparams/multifitparam.$state");
	my @fit_lines = <$fit_handle>;
	my @splitted_values = split(' ', $fit_lines[0]);
	$tmin_pref = $splitted_values[0];
	$bin_size = $splitted_values[2];
	close($fit_handle);
	printf("Prf tmin %d\n", $tmin_pref);
}
else
{
	$tmin_pref = 5;
	$bin_size = 1;
	printf("Prf tmin %d\n", 5);
}

# Push updated tmax anyway.
open(my $fit_handle2, ">$path/$relpath/$direc/spectrum2/multifitparams/multifitparam.$state");
print $fit_handle2 "$tmin_pref ".int($tmax+0.5)." $bin_size"; # some defaults.
close($fit_handle2);


# Alright! Learn something about the mean and error of ground states.
# This helps us pick y-intervals.

# Get the avg ground state mass, avg ground state error.
my $mean_dir = 0;
my $err_dir = 0;
my $mean_osc = 0;
my $err_osc = 0;

# Get preferred ones too, if we can.
my $pref_dir = 0;
my $pref_dir_err = 0;
my $pref_osc = 0;
my $pref_osc_err = 0;
my $pref_num_dir = 0;
my $pref_num_osc = 0;

# Save masses, coefficients.
my @mdc_list = ();
my @mde_list = ();
my @moc_list = ();
my @moe_list = ();

my @cdc_list = ();
my @cde_list = ();
my @coc_list = ();
my @coe_list = ();

foreach my $elem (@fit_info)
{
	my @vals = @{$elem};
	if ($num_dir > 0)
	{
		$mean_dir += $vals[DIRMASS1]/$num_fits;
		$err_dir += $vals[DIRMASS1+1]/$num_fits;
		if ($vals[TMIN] == $tmin_pref)
		{
			$pref_dir = $vals[DIRMASS1];
			$pref_dir_err = $vals[DIRMASS1+1];
			$pref_num_dir = $vals[FITTYPE]%4;
			
			for (my $j = 0; $j < $pref_num_dir; $j++)
			{
				push(@cdc_list, $vals[DIRAMP1+4*$j]);
				push(@cde_list, $vals[DIRAMP1+4*$j+1]);
				push(@mdc_list, $vals[DIRMASS1+4*$j]);
				push(@mde_list, $vals[DIRMASS1+4*$j+1]);
			}
		}
	}
	if ($num_osc > 0)
	{
		$mean_osc += $vals[OSCMASS1]/$num_fits;
		$err_osc += $vals[OSCMASS1+1]/$num_fits;
		if ($vals[TMIN] == $tmin_pref)
		{
			$pref_osc = $vals[OSCMASS1];
			$pref_osc_err = $vals[OSCMASS1+1];
			$pref_num_osc = ($vals[FITTYPE]-$pref_num_dir)/4;
			
			for (my $j = 0; $j < $pref_num_osc; $j++)
			{
				push(@coc_list, $vals[OSCAMP1+4*$j]);
				push(@coe_list, $vals[OSCAMP1+4*$j+1]);
				push(@moc_list, $vals[OSCMASS1+4*$j]);
				push(@moe_list, $vals[OSCMASS1+4*$j+1]);
			}
		}
	}
	if ($is_fpi == 1) # We have fpi, reuse oscil.
	{
		$mean_osc += $vals[AUX1]/$num_fits;
		$err_osc += $vals[AUX1+1]/$num_fits;
		if ($vals[TMIN] == $tmin_pref)
		{
			$pref_osc = $vals[AUX1];
			$pref_osc_err = $vals[AUX1+1];
			
			push(@moc_list, $vals[AUX1]);
			push(@moe_list, $vals[AUX1+1]);
		}
	}
}

printf("Avg dir. = %.6f\n", $mean_dir);
printf("Err dir. = %.6f\n", $err_dir);
printf("Avg osc. = %.6f\n", $mean_osc);
printf("Err osc. = %.6f\n", $err_osc);

# Save the preferred info to a multicentral file.
# Save #dir, #osc, and center+errs of each.
open (my $an_outfile, ">$path/$relpath/$direc/spectrum2/multicentral/multicentral.".$state);
print $an_outfile "$pref_num_dir $pref_dir $pref_dir_err\n";
print $an_outfile "$pref_num_osc $pref_osc $pref_osc_err\n";
close($an_outfile);

# Also, save all preferred coefficients to a tex table.

my @table_lines = ();
if ($num_dir > 0 && $num_osc == 0 && $is_fpi == 0)
{
	(my $dirl = $dir_label) =~ s/\\\\/\\/;
	push(@table_lines, "\\begin{tabular}{|c|c|}");
	push(@table_lines, "\$A\$ & \$".$dirl."\$ \\\\");
}
elsif ($num_osc > 0)
{
	(my $dirl = $dir_label) =~ s/\\\\/\\/;
	(my $oscl = $osc_label) =~ s/\\\\/\\/;
	push(@table_lines, "\\begin{tabular}{|c|c|c|c|}");
	push(@table_lines, "\$A\$ & \$".$dirl."\$ & \$B\$ & \$".$oscl."\$\\\\");
}
elsif ($is_fpi == 1)
{
	(my $dirl = $dir_label) =~ s/\\\\/\\/;
	(my $oscl = $osc_label) =~ s/\\\\/\\/;
	push(@table_lines, "\\begin{tabular}{|c|c|c|}");
	push(@table_lines, "\$A\$ & \$".$dirl."\$ & \$".$oscl."\$\\\\");
}

for (my $j = 0; $j < 3; $j++)
{
	my $local_string = "";
	# There's both a direct and oscillating state.
	if ((defined $mdc_list[$j]) && (defined $moc_list[$j]))
	{
		if ($is_fpi == 1)
		{
			$local_string = sprintf("%s & %s & %s \\\\", parenForm($cdc_list[$j], $cde_list[$j]), parenForm($mdc_list[$j], $mde_list[$j]), parenForm($moc_list[$j], $moe_list[$j]));
		}
		else
		{
			$local_string = sprintf("%s & %s & %s & %s \\\\", parenForm($cdc_list[$j], $cde_list[$j]), parenForm($mdc_list[$j], $mde_list[$j]), parenForm($coc_list[$j], $coe_list[$j]), parenForm($moc_list[$j], $moe_list[$j]));
		}
	}
	elsif (defined $mdc_list[$j]) # there's just a direct state.
	{
		if ($num_osc > 0) # there just isn't an osc. state here.
		{
			$local_string = sprintf("%s & %s & --- & --- \\\\", parenForm($cdc_list[$j], $cde_list[$j]), parenForm($mdc_list[$j], $mde_list[$j]));
		}
		elsif ($is_fpi == 1)
		{
			$local_string = sprintf("%s & %s & --- \\\\", parenForm($cdc_list[$j], $cde_list[$j]), parenForm($mdc_list[$j], $mde_list[$j]));
		}
		else # it's the wall source pion
		{
			$local_string = sprintf("%s & %s \\\\", parenForm($cdc_list[$j], $cde_list[$j]), parenForm($mdc_list[$j], $mde_list[$j]));
		}
	}
	elsif (defined $moc_list[$j]) # there's just an oscil state.
	{
		# there's no weird conditionals here.
		$local_string = sprintf("--- & --- & %s & %s \\\\", parenForm($coc_list[$j], $coe_list[$j]), parenForm($moc_list[$j], $moe_list[$j]));
	}
	else # Ain't got nuttin!
	{
		if ($num_osc > 0) # there just isn't an osc. state here.
		{
			$local_string = "--- & --- & --- & --- \\\\";
		}
		elsif ($is_fpi == 1)
		{
			$local_string = "--- & --- & --- \\\\";
		}
		else # it's the wall source pion
		{
			$local_string = "--- & --- \\\\";
		}
	}
		
	push(@table_lines, $local_string);
}

# Write it out to file!
open (my $outfile5, ">$path/$relpath/$direc/spectrum2/tex/multitable.$state");
foreach my $line (@table_lines)
{
	print $outfile5 $line."\n\\hline\n";
}
print $outfile5 "\\end{tabular}\n";
close($outfile5);

##############################
# Prepare the fit curve now! #
##############################

my $curve_line = "";
for (my $j = 0; $j < $pref_num_dir; $j++)
{
	if ($state eq "nu" || $state eq "de")
	{
		$curve_line = $curve_line."+".$cdc_list[$j]."*(exp(".$mdc_list[$j]."*(".($ensemble_params{T}/2)."-x))-cos(3.1415926535*x)*exp(-".$mdc_list[$j]."*(".($ensemble_params{T}/2)."-x)))/2.0";
	}
	else
	{
		$curve_line = $curve_line."+".$cdc_list[$j]."*cosh(".$mdc_list[$j]."*(".($ensemble_params{T}/2)."-x))";
	}
}

for (my $j = 0; $j < $pref_num_osc; $j++)
{
	if ($state eq "nu" || $state eq "de")
	{
		$curve_line = $curve_line."+".$coc_list[$j]."*cos(3.1415926535*x)*(exp(".$moc_list[$j]."*(".($ensemble_params{T}/2)."-x))-cos(3.1415926535*x)*exp(-".$moc_list[$j]."*(".($ensemble_params{T}/2)."-x)))/2.0";
	}
	else
	{
		$curve_line = $curve_line."+".$coc_list[$j]."*cos(3.1415926535*x)*cosh(".$moc_list[$j]."*(".($ensemble_params{T}/2)."-x))";
	}
}

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
open(my $outfile4, ">EvanSpectrum/scripts/tmp_space/pvalue");
foreach my $a_row (@fit_info)
{
	my @the_row = @{$a_row};
	print $outfile4 $the_row[TMIN]." ".$the_row[PVAL]."\n";
}
close($outfile4);

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
############################################################
# Copy the correlator and learn something about its range.
############################################################

# First, load the sum file.
open (my $corr_handle, "<$path/$relpath/$direc/spectrum2/sum/sum.$state");
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
	
	if (abs(abs($temp_splitter[1])-$temp_splitter[2]) < $min_corr_val && abs($temp_splitter[1]) > 1e-9)
	{
		$min_corr_val = abs(abs($temp_splitter[1])-$temp_splitter[2]);
	}
	
}



# Write the correlator out to file. Maybe this is redundant now!

open(my $outcorr_handle, ">./EvanSpectrum/scripts/tmp_space/corr.dat");
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

############################################
# I guess it's time to build some plots!
############################################

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
print $outfile3 "set style line 22 lt 1 lc rgb '#CC0066' lw 3 pt 10 ps $ptsize\n";
print $outfile3 "set style line 31 lt 2 lc rgb '#555555' pt 11 ps $ptsize\n";

##########################################
# First, prepare the plot of all states.
##########################################

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
	print $outfile3 "set arrow $iter from ".($tchg).",0 to ".($tchg).",1 nohead ls 31\n";
	$iter++;
}

# If tmin_pref is set, window the preferred value.
if ($tmin_pref != 0)
{
	print $outfile3 "set arrow 10 from ".($tmin_pref-0.3).",0 to ".($tmin_pref-0.3).",1 nohead ls 22\n";
	print $outfile3 "set arrow 11 from ".($tmin_pref+0.3).",0 to ".($tmin_pref+0.3).",1 nohead ls 22\n";
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

#################################
# Next, create zoomed in plots.
#################################

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
		print $outfile3 "set arrow $iter from ".($tchg).",".($mean_dir-$sm*$err_dir)." to ".($tchg).",".($mean_dir+$sm*$err_dir)." nohead ls 31\n";
		$iter++;
	}
	
	if ($tmin_pref != 0)
	{
		print $outfile3 "set arrow 10 from ".($tmin_pref-0.3).",".($mean_dir-$sm*$err_dir)." to ".($tmin_pref-0.3).",".($mean_dir+$sm*$err_dir)." nohead ls 22\n";
		print $outfile3 "set arrow 11 from ".($tmin_pref+0.3).",".($mean_dir-$sm*$err_dir)." to ".($tmin_pref+0.3).",".($mean_dir+$sm*$err_dir)." nohead ls 22\n";
		
		print $outfile3 "set arrow 12 from ".($tmin_min-1).",".($pref_dir+$pref_dir_err)." to ".($tmin_max+1).",".($pref_dir+$pref_dir_err)." nohead ls 22\n";
		print $outfile3 "set arrow 13 from ".($tmin_min-1).",".($pref_dir-$pref_dir_err)." to ".($tmin_max+1).",".($pref_dir-$pref_dir_err)." nohead ls 22\n";
		
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/multifit_dir0\" using 1:2:3 with yerrorbars ls ".(1)." title \"\$ ".$dir_label." \$\"";
		
		
	if ($singlefit_flag == 1 && $num_change[0] != 0)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_dir\" using (\$1+0.2):2:3 with yerrorbars ls ".(2+10)." notitle";
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
		print $outfile3 "set arrow $iter from ".($tchg).",".($mean_osc-$sm*$err_osc)." to ".($tchg).",".($mean_osc+$sm*$err_osc)." nohead ls 31\n";
		$iter++;
	}
	
	if ($tmin_pref != 0)
	{
		print $outfile3 "set arrow 10 from ".($tmin_pref-0.3).",".($mean_osc-$sm*$err_osc)." to ".($tmin_pref-0.3).",".($mean_osc+$sm*$err_osc)." nohead ls 22\n";
		print $outfile3 "set arrow 11 from ".($tmin_pref+0.3).",".($mean_osc-$sm*$err_osc)." to ".($tmin_pref+0.3).",".($mean_osc+$sm*$err_osc)." nohead ls 22\n";
		
		print $outfile3 "set arrow 12 from ".($tmin_min-1).",".($pref_osc+$pref_osc_err)." to ".($tmin_max+1).",".($pref_osc+$pref_osc_err)." nohead ls 22\n";
		print $outfile3 "set arrow 13 from ".($tmin_min-1).",".($pref_osc-$pref_osc_err)." to ".($tmin_max+1).",".($pref_osc-$pref_osc_err)." nohead ls 22\n";
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/multifit_osc0\" using 1:2:3 with yerrorbars ls ".(1+3)." title \"\$ ".$osc_label." \$\"";
		
		
	if ($singlefit_flag == 1 && $num_change[0] != 0)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlefit_osc\" using (\$1+0.2):2:3 with yerrorbars ls ".(2+3+10)." notitle";
	}

	print $outfile3 "\n";
	print $outfile3 "set output\n";
}

################################
# Last, create a p-value plot.
################################

if ($tmin_pref != 0)
{
	print $outfile3 "unset arrow 12\n";
	print $outfile3 "unset arrow 13\n";
}

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
		print $outfile3 "set arrow $iter from ".($tchg).",0 to ".($tchg).",1 nohead ls 31\n";
		$iter++;
	}
	
	if ($tmin_pref != 0)
	{
		print $outfile3 "set arrow 10 from ".($tmin_pref-0.3).",0 to ".($tmin_pref-0.3).",1 nohead ls 22\n";
		print $outfile3 "set arrow 11 from ".($tmin_pref+0.3).",0 to ".($tmin_pref+0.3).",1 nohead ls 22\n";
	}

	# Prepare to plot!

	print $outfile3 "plot ";
	
	print $outfile3 "\"./EvanSpectrum/scripts/tmp_space/pvalue\" using 1:2 ls ".(1)." notitle";
		
		
	if ($singlefit_flag == 1 && $num_change[0] != 0)
	{
		print $outfile3 ", \"./EvanSpectrum/scripts/tmp_space/singlepvalue\" using (\$1+0.2):2 ls ".(1+10)." notitle";
	}
	
	print $outfile3 ", 0.05 ls 15 title \"5 percent\"";

	print $outfile3 "\n";
	print $outfile3 "set output\n";
}

#####################################################
# And now, the most fun possible, correlator plots!
#####################################################

# Reset arrows.
$iter = 1;
foreach my $tchg (@num_change)
{
	print $outfile3 "unset arrow $iter\n";
	$iter++;
}

print $outfile3 "set xlabel \"\$t\$\"\n";
print $outfile3 "set ylabel \"Log Plot of Correlator\"\n";
print $outfile3 "set xrange [-0.4:".($ensemble_params{T}/2+0.4)."]\n";
print $outfile3 "set log y\n";
print $outfile3 "set format y \"%2.0e\"\n";
print $outfile3 "set yrange [".(0.9*$min_corr_val).":".(1.1*$max_corr_val)."]\n";
print $outfile3 "set samples 1000\n";
print $outfile3 "set arrow 10 from ".($tmin_pref-0.5).",".(0.9*$min_corr_val)." to ".($tmin_pref-0.5).",".(1.1*$max_corr_val)." nohead ls 22\n";
print $outfile3 "set arrow 11 from ".($tmax_pref+0.5).",".(0.9*$min_corr_val)." to ".($tmax_pref+0.5).",".(1.1*$max_corr_val)." nohead ls 22\n";

print $outfile3 "set terminal epslatex color standalone\n";
print $outfile3 "set output \"./EvanSpectrum/scripts/tmp_space/corr.tex\"\n";
print $outfile3 "plot \"./EvanSpectrum/scripts/tmp_space/corr.dat\" using 1:2:3 with yerrorbars ls 1 notitle, \\\n";
print $outfile3 "   $curve_line ls 2 notitle\n";
print $outfile3 "set output\n";


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

`epstopdf ./EvanSpectrum/scripts/tmp_space/corr-inc.eps`;
`pdflatex ./EvanSpectrum/scripts/tmp_space/corr.tex`;
`convert ./corr.pdf ./corr.png`;
$command = "mv corr.pdf $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-corr.pdf";
`$command`;
$command = "mv corr.png $path/$relpath/$direc/spectrum2/multifits/plots/multifits-$state"."-corr.png";
`$command`;

`rm *.aux`;
`rm *.log`;

close($outfile3);


sub parenForm
{
   my ($the_value, $the_error) = @_;
   
   my $neg_flag = 1;
   if ($the_value < 0)
   {
		$neg_flag = -1;
		$the_value = -$the_value;
   }
   
   if (abs($the_value) < 1e-20)
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
