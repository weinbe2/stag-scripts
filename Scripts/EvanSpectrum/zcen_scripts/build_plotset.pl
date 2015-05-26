#! /usr/bin/perl
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;


# Build the set of plots for the new analysis automatically.

# Make sure there's a sufficient number of arguments.

my $num_args = @ARGV;
if ($num_args < 2)
{
   die "Invalid number of arguments!";
}

my %params = ();
my $amp = 5e-3;
my $fittype = "sg";
my $tmin = 10;

# See if we can grab something from the infofile.
my $command = "";

my $got_ens = 0;
my $amp_changed = 0;
my $type_changed = 0;
my $tmin_changed = 0;

# Prepare for tex outut
my $mass_describe = "";
my $flavor_describe = "";

for (my $i = 0; $i < @ARGV; $i++)
{
	if (substr($ARGV[$i], 0, 1) ne "-")
	{
		die "Expected a flag, got another argument.\n";
	}
	
	substr($ARGV[$i], 0, 1) = "";
	
	switch ($ARGV[$i])
	{
		case "amp" {
			my $temp_amp = $ARGV[$i+1];
			if ($temp_amp != $amp)
			{
				$amp_changed = 1;
				$amp = $temp_amp;
			}
			$i++;
		}
		case "type" {
			my $temp_type = $ARGV[$i+1];
			if ($temp_type ne $fittype && ($temp_type eq "sg" || $temp_type eq "dc"))
			{
				$type_changed = 1;
				$fittype = $temp_type;
			}
			$i++;
		}
		case "tmin" {
			my $temp_tmin = $ARGV[$i+1];
			if ($temp_tmin != $tmin)
                        {
                                $tmin_changed = 1;
                                $tmin = $temp_tmin;
                        }
                        $i++;
                }
		case "ensemble" {
			# Prepare pre-parsed info.
			
			$got_ens = 1;		
			%params = name2param($ARGV[$i+1]);
			$params{ensemble} = $ARGV[$i+1];
				
			if ($params{flavor} eq '4plus8') {
				$params{have_l} = 4;
				$params{have_h} = 8;
				$mass_describe = "\$m_\\\\ell = $params{ml}, m_h = $params{mh}\$";
				$flavor_describe = "4+8";
			}
			elsif ($params{flavor} eq '4plus4') {
				$params{have_l} = 4;
				$params{have_h} = 4;
				$mass_describe = "\$m_\\\\ell = $params{ml}, m_h = $params{mh}\$";
                                $flavor_describe = "4+4";
			}
			elsif ($params{flavor} eq '8') {
				$params{have_l} = 8;
				$params{have_h} = 0;
				$params{mh} = 0;
				$mass_describe = "\$m_q = $params{ml}\$";
                                $flavor_describe = "8";

			}
			elsif ($params{flavor} eq '4') {
				$params{have_l} = 4;
				$params{have_h} = 0;
				$params{mh} = 0;
				$mass_describe = "\$m_q = $params{ml}\$";
                                $flavor_describe = "4";

			}
			elsif ($params{flavor} eq '12') {
				$params{have_l} = 12;
				$params{have_h} = 0;
				$params{mh} = 0;
				$mass_describe = "\$m_q = $params{ml}\$";
                                $flavor_describe = "12";

			}

			# Look for an entry in the infofile.
			if (-f "infofile.dat") {
			$command = "grep $params{ensemble} \"infofile.dat\" | tail -n 1";
			my $infofile_info = `$command`;
			if (length($infofile_info) > 0)
			{
				my @infos = split(' ', $infofile_info);
				if ($amp_changed == 0) { $amp = $infos[1]; }
				if ($type_changed == 0) { $fittype = $infos[2]; }
				if ($tmin_changed == 0) { $tmin = $infos[3]; }
			}
			}
			$i++;
		}
		else: {
			die "Unknown flag.\n";
		}
	}

}

if ($got_ens == 0)
{
	die "Didn't get an ensemble!\n";
}

if ($amp_changed == 1 || $type_changed == 1 || $tmin_changed == 1)
{
	# Spit out this info to the tracker file. Sloppily.
	open(my $info_file, ">>./infofile.dat");
	print $info_file "$params{ensemble} $amp $fittype $tmin\n";
	close($info_file);
}

my $basepath = "../../".$params{path}."/".$params{ensemble}."/spectrum2";

# What files do we need?
# uwerr/uwerr.vev
# sum/sum.dc_stoch
# sum/sum.sg_stoch
# fits/fit_new.dc_stoch_zcen
# fits/fit_new.sg_stoch_zcen
# effmass/effmass.dc_stoch_kuti1
# effmass/effmass.sg_stoch_kuti1

# Copy all these files locally.


$command = "cp $basepath/uwerr/uwerr.vev ./uwerr.vev";
`$command`;

$command = "cp $basepath/sum/sum.dc_stoch ./sum.dc_stoch";
`$command`;
$command = "cp $basepath/sum/sum.sg_stoch ./sum.sg_stoch";
`$command`;
$command = "./fix_orig.sh $basepath/fits/fit_new.dc_stoch_zcen > ./dc_stoch_zcen.dat";
`$command`;
$command = "./fix_orig.sh $basepath/fits/fit_new.sg_stoch_zcen > ./sg_stoch_zcen.dat";
`$command`;
$command = "cp $basepath/effmass/effmass.dc_stoch_kuti1 ./effmass.dc_stoch_kuti1";
`$command`;
$command = "cp $basepath/effmass/effmass.sg_stoch_kuti1 ./effmass.sg_stoch_kuti1";
`$command`;
$command = "cp $basepath/fitparams/fitparam.dc_stoch ./fitparam.dc_stoch";
`$command`;

# Learn the vacuum!
open (my $vev_file, "<./uwerr.vev");
my @vev_info = <$vev_file>;
close($vev_file);

# Get the vev, the vacuum, and the autocorr.
my @vev1 = split(' ', $vev_info[0]);
my $vev_center = parenForm($vev1[1], $vev1[2]);
my $vev_autocorr = floor($vev1[4]+0.5);
my @vev2 = split(' ', $vev_info[1]);
my $vev_info = parenForm($vev2[1], $vev2[2]);

# Get an updated autocorr.
open (my $acorr_file, "<./fitparam.dc_stoch");
my @acorr_info = <$acorr_file>;
close($acorr_file);

my @acorr = split(' ', $acorr_info[0]);
$vev_autocorr = $acorr[4];

# Get tau in MDTU

$command = "sed '1q;d' $basepath/stoch/SPECTRUM.dat | cut -d\" \" -f1";
my $config_first = `$command`;

$command = "sed '".(1+$params{T})."q;d' $basepath/stoch/SPECTRUM.dat | cut -d\" \" -f1";
my $config_second = `$command`;

my $MDTU_per_meas = $config_second - $config_first;


# Get N_meas, N.
$command = "wc -l $basepath/stoch/SPECTRUM.dat";
my $N_meas_str = `$command`;
my @N_meas_split = split(' ', $N_meas_str);
my $N_meas = $N_meas_split[0]/$params{T};
my $N = floor($N_meas/$vev_autocorr);

# Get M_pi, M_rho.
$command = "cut -d\" \" -f2,3 ../../../Data/pion/pion_$params{ensemble}";
my @pion_info = split(' ', `$command`);
my $pions = parenFormWiki($pion_info[0], $pion_info[1]);

$command = "cut -d\" \" -f2,3 ../../../Data/rho/rho_$params{ensemble}";
my @rho_info = split(' ', `$command`);
my $rhos = parenFormWiki($rho_info[0], $rho_info[1]);

# Get more!
$command = "cut -d\" \" -f2,3 ../../../Data/fpi/fpi_$params{ensemble}";
@pion_info = split(' ', `$command`);
my $fpis = parenFormWiki($pion_info[0], $pion_info[1]);

$command = "cut -d\" \" -f2,3 ../../../Data/a0/a0_$params{ensemble}";
@rho_info = split(' ', `$command`);
my $a0s = parenFormWiki($rho_info[0], $rho_info[1]);

$command = "cut -d\" \" -f2,3 ../../../Data/axial/axial_$params{ensemble}";
@pion_info = split(' ', `$command`);
my $axials = parenFormWiki($pion_info[0], $pion_info[1]);

$command = "cut -d\" \" -f2,3 ../../../Data/nucleon/nucleon_$params{ensemble}";
@rho_info = split(' ', `$command`);
my $nucleons = parenFormWiki($rho_info[0], $rho_info[1]);

# Prepare a fit curve in the data!
# There's a special level of hell for this scripting.

my $corrname = ($fittype eq "sg") ? ($params{have_l}/4)."D-C" : ($params{have_l}/4)."D";
my $corrtype = $fittype;
$command = "awk ' { if (int(\$1+0.5) == ".($params{T}/2).") { print \$2 } }' sum.$fittype"."_stoch";
my $corrcen = `$command`; chomp($corrcen);
$command = "awk ' { if (int(\$1+0.5) == $tmin) { print \$14 } }' $fittype"."_stoch_zcen.dat";
my $corrconst = `$command`; chomp($corrconst);
$command = "awk ' { if (int(\$1+0.5) == $tmin) { print \$6 } }' $fittype"."_stoch_zcen.dat";
my $corrmin = abs(`$command`)*0.5;
my $corrmax = 1e1;
my $corrtmin = $tmin;
$command = "./get_vals.sh $fittype"."_stoch_zcen.dat $tmin";
my $corrfunc = `$command`; chomp($corrfunc); $corrfunc = $corrfunc."-($corrconst)";


# Create a custom sigma_plots.plt

$command = "sed -e 's/%NF%/".($params{have_l}/4)."/g' -e 's/%VEV%/\$".($vev_info)."\$/g' -e 's/%AMP%/$amp/g' -e 's/%CORRNAME%/$corrname/g' -e 's/%CORRTYPE%/$corrtype/g' -e 's/%CORRCEN%/$corrcen/g' -e 's/%CORRCONST%/$corrconst/g' -e 's/%CORRMIN%/$corrmin/g' -e 's/%CORRMAX%/$corrmax/g' -e 's/%CORRTMIN%/$corrtmin/g' -e 's/%CORRFUNC%/$corrfunc/g' -e 's/%NTD2%/".(floor($params{T}/2+0.5))."/g' sigma_plot.plt > sigma_plots_run.plt";
`$command`;

# Create a custom text_base.txt




$command = "sed -e 's/%NF%/".($params{have_l}/4)."/g' -e 's/%NMEAS%/$N_meas/g' -e 's/%N%/$N/g' -e 's/%TAU%/$vev_autocorr/g' -e 's/%TAUMDTU%/".($MDTU_per_meas*$vev_autocorr)."/g' -e 's/%4NF%/$flavor_describe/g' -e 's/%BETA%/$params{beta}/g' -e 's/%L%/$params{L}/g' -e 's/%MASSDESC%/$mass_describe/g' -e 's/%PION%/$pions/g' -e 's/%RHO%/$rhos/g' -e 's/%FPI%/$fpis/g' -e 's/%A0%/$a0s/g' -e 's/%AXIAL%/$axials/g' -e 's/%NUCLEON%/$nucleons/g' -e 's/%CORRNAME%/$corrname/g' -e 's/%CORRTMIN%/$corrtmin/g' tex_base.tex > tex_copy.tex";
`$command`;

$command = "sed -e 's/%NF%/".($params{have_l}/4)."/g' -e 's/%NMEAS%/$N_meas/g' -e 's/%N%/$N/g' -e 's/%TAU%/$vev_autocorr/g' -e 's/%TAUMDTU%/".($MDTU_per_meas*$vev_autocorr)."/g' -e 's/%4NF%/$flavor_describe/g' -e 's/%BETA%/$params{beta}/g' -e 's/%L%/$params{L}/g' -e 's/%MASSDESC%/$mass_describe/g' -e 's/%PION%/$pions/g' -e 's/%RHO%/$rhos/g' -e 's/%FPI%/$fpis/g' -e 's/%A0%/$a0s/g' -e 's/%AXIAL%/$axials/g' -e 's/%NUCLEON%/$nucleons/g' -e 's/%ENSEMBLE%/$params{ensemble}/g' -e 's/%CORRNAME%/$corrname/g' -e 's/%CORRTMIN%/$corrtmin/g' tex_base_input.tex > tex_copy_input.tex";
`$command`;


# Execute it!
$command = "gnuplot -e \"load './sigma_plots_run.plt'\"";
`$command`;

$command = "epstopdf m0pp_zcen_amp-inc.eps; pdflatex m0pp_zcen_amp.tex";
`$command`;

$command = "epstopdf m0pp_zcen_amp_dc-inc.eps; pdflatex m0pp_zcen_amp_dc.tex";
`$command`;

$command = "epstopdf m0pp_zcen_cmp-inc.eps; pdflatex m0pp_zcen_cmp.tex";
`$command`;

$command = "epstopdf m0pp_zcen_eff-inc.eps; pdflatex m0pp_zcen_eff.tex";
`$command`;

$command = "epstopdf m0pp_zcen_amp_orig-inc.eps; pdflatex m0pp_zcen_amp_orig.tex";
`$command`;

$command = "epstopdf m0pp_zcen_fitlines-inc.eps; pdflatex m0pp_zcen_fitlines.tex";
`$command`;

# Build a big tex file.

$command = "pdflatex tex_copy.tex; pdflatex tex_copy.tex; pdflatex tex_copy.tex;";
`$command`;


# Clean up at end.

#$command = "rm ./tex_copy.aux; rm ./tex_copy.out; rm ./tex_copy.log; rm ./tex_copy.tex;";
$command = "rm ./tex_copy.aux; rm ./tex_copy.out; rm ./tex_copy.log;";
`$command`;

$command = "rm ./uwerr.vev";
`$command`;
$command = "rm ./sum.dc_stoch";
`$command`;
$command = "rm ./sum.sg_stoch";
`$command`;
$command = "rm ./dc_stoch_zcen.dat";
`$command`;
$command = "rm ./sg_stoch_zcen.dat";
`$command`;
$command = "rm ./effmass.dc_stoch_kuti1";
`$command`;
$command = "rm ./effmass.sg_stoch_kuti1";
`$command`;

$command = "rm -f m0pp*.eps; rm -f m0pp*.tex; rm -f m0pp*.log; rm -f m0pp*.aux; rm -f m0pp*-inc.pdf;";
`$command`;

$command = "rm sigma_plots_run.plt";
`$command`; 

# Move everything into the package folder.
$command = "mv m0pp*.pdf package/plots; mv tex_copy.pdf package; mv tex_copy.tex package; mv tex_copy_input.tex package;";
`$command`;

# Copy things
$command = "mkdir -p /projectnb/qcd/Staggered/All/spectrum_disconnected_lsd/zcen_file/$params{ensemble}";
`$command`;
$command = "cp -r package/* /projectnb/qcd/Staggered/All/spectrum_disconnected_lsd/zcen_file/$params{ensemble}/";
`$command`;

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
                return '-'.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.")\\\\\\\\mbox{e".$value_size.'} ';
        }
        else
        {
                return ''.substr($the_value, 0, 1).".".substr($the_value, 1)."(".$the_error.")\\\\\\\\mbox{e".$value_size.'} ';
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


