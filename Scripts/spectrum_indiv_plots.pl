#! /usr/bin/perl 
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_indiv_plots.pl f4plus8l24t48b30m005m

if ($#ARGV < 0) {
   die "Inappropriate number of arguments: expected spectrum.pl [directory]\n";
}

my $direc = $ARGV[0];
my $run = $direc;
# This can change on different systems. 
my $scrfolder = "/projectnb/qcd/Staggered/Scripts/Spectrum/spec_analysis/";

my $direc_star = namestar($run);
my %params = name2param($direc_star);

my $path = $params{path};
my $cut = $params{cut};
my $L = $params{L};
my $T = $params{T};
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

# Find all folders matching the star match pattern.

my @match_directories = <./$direc_star>;

# Learn what we can from each directory,
# building a hash for a combined analysis.

my %combined_data;

foreach my $indiv_path (@match_directories)
{
	# Get info about the directory.
	my %indiv_params = name2param(substr($indiv_path, 2));
	
	my $indiv_L = $params{L};
	my $indiv_T = $params{T};
	my $indiv_beta = $params{beta};
	my $indiv_m_l;
	my $indiv_m_h;
	my $indiv_have_l = 0;
	my $indiv_have_h = 0;

	if ($params{flavor} eq '4plus8') {
		$indiv_have_l = 1; $indiv_have_h = 1;
		$indiv_m_h = $indiv_params{mh};
		$indiv_m_l = $indiv_params{ml};
	}
	elsif ($params{flavor} eq '4plus4') {
		$indiv_have_l = 1; $indiv_have_h = 1;
		$indiv_m_h = $indiv_params{mh};
		$indiv_m_l = $indiv_params{ml};
	}
	elsif ($params{flavor} eq '8') {
		$indiv_have_l = 1; $indiv_have_h = 0;
		$indiv_m_l = $indiv_params{ml};
	}
	elsif ($params{flavor} eq '4') {
		$indiv_have_l = 1; $indiv_have_h = 0;
		$indiv_m_l = $indiv_params{ml};
	}
	
	# Make sure a hadrons folder exists. If it doesn't, continue.
	
	if (!(-d $path."/".$indiv_path."/hadrons"))
	{
		next;
	}
	
	print "$indiv_path has a hadrons folder!\n";
	
	# Get all files in the hadrons folder.
	my $temp_path = $indiv_path."/hadrons/*";
	my @spec_files = <./$temp_path>;
	
	foreach my $hadron_file (@spec_files)
	{
		my @splitter = split("/", $hadron_file);
		my @splitter2 = split(/\./, $splitter[@splitter-1]);
		my $spectype = $splitter2[0];
		
		chomp($spectype);
		
		open(my $filehandle, "<$hadron_file");
		my @filecontents = <$filehandle>;
		close($filehandle);
		if (@filecontents == 0)
		{
			next;
		}
		
		chomp($filecontents[0]);
		
		if (exists $combined_data{$spectype})
		{
			if ($indiv_have_h == 1)
			{
				push(@{$combined_data{$spectype}}, "$indiv_m_h $filecontents[0]\n");
			}
			else
			{
				push(@{$combined_data{$spectype}}, "$indiv_m_l $filecontents[0]\n");
			}
		}
		else
		{
			if ($indiv_have_h == 1)
			{
				push(@{$combined_data{$spectype}}, "$indiv_m_h $filecontents[0]\n");
			}
			else
			{
				push(@{$combined_data{$spectype}}, "$indiv_m_l $filecontents[0]\n");
			}
		}
	}
	
}

# Make a hash of output file names.
my %outputs_hash;

# Write the combined file.
foreach my $key (keys %combined_data)
{
	open(my $outfile, ">$path"."/hadrons_data/$direc"."_$key.dat");
	$outputs_hash{$key} = "$path"."/hadrons_data/$direc"."_$key.dat";
	chomp($key);
	print $key."\n";
	foreach my $line (@{$combined_data{$key}})
	{
		print $outfile $line;
	}
	close($outfile);
}

# Next, we can make plots!

my $M_script = "cd('../Scripts');";
my $label;

if ($have_h == 1)
{
	$label = "m_h";
}
else
{
	$label = "m_l";
}

foreach my $key (keys %outputs_hash)
{
	my $outputs_eps = "$path"."/hadrons_plots/$direc"."_$key.eps";
	$M_script = $M_script."myplot3_spec('$outputs_hash{$key}', '$outputs_eps', 1, 2, 3, 'id', '$label', 'm_{".$key."}', 'm_{".$key."}');";
}

# Also, add plots of all the pions and rhos together!
my $combined_dat = "$path"."/hadrons_data/$direc";
my $combined_plot = "$path"."/hadrons_plots/$direc";

$M_script = $M_script."myplot_many_spec('$combined_plot"."_pions.eps', 1, '$combined_dat"."_ps1.dat', 2, 3, '$combined_dat"."_sc1.dat', 2, 3, '$combined_dat"."_ij1.dat', 2, 3, '$combined_dat"."_i51.dat', 2, 3, 'id', '$label', 'm_{pion}', {'m_{ps1}', 'm_{sc1}', 'm_{ij1}', 'm_{i51}'});";
$M_script = $M_script."myplot_many_spec('$combined_plot"."_rhos.eps', 1, '$combined_dat"."_r01.dat', 2, 3, '$combined_dat"."_ris1.dat', 2, 3, '$combined_dat"."_rij1.dat', 2, 3, '$combined_dat"."_ri51.dat', 2, 3, 'id', '$label', 'm_{rho}', {'m_{r01}', 'm_{ris1}', 'm_{rij1}', 'm_{ri51}'});";

#print $M_script;

my $matlab = "matlab -softwareopengl -nosplash -nodisplay -r \"$M_script quit;\"";
system($matlab);

# convert all eps to pdf (in plot directory and remove eps)
# needed to preserve bounding box
my @epsfiles =  glob "hadrons_plots/*eps";
foreach (@epsfiles) {
    my $h = `head -1 $_`;
    if ($h ne '%!PS-Adobe-3.0') {
	my $f = `cat $_`;
	my $fix = "%!PS-Adobe-3.0\n%%LanguageLevel: 2\n%%Pages: (atend)\n%%BoundingBox: (atend)\n%%EndComments\n\n";
	open(my $ofix, '>', $_) or die "Could not fix file";
	print $ofix $fix;
	print $ofix $f;
	close($ofix);
    }
    print $_."\n";
    system("eps2goodpdf.pl $_");
    system("rm -f $_"); 
}
system("chmod 664 hadrons_plots/*pdf");

# Next, make the PDF.

my $localcommand = "cat $scrfolder"."../../plot5_spec.tex";
my $tex = `$localcommand`;
$tex =~ s/:DIREC:/$direc/g;

open (my $Fpdf, ">$path"."/hadrons_pdf/$direc.tex") or die "Could not open file: $!"; chmod 0664, $Fpdf;
print $Fpdf "\\documentclass[letterpaper,10pt]{article}\n";
print $Fpdf "\\usepackage[pdftex]{color, graphicx, xcolor}\n";
print $Fpdf "\\usepackage[colorlinks=true, linkcolor=blue, filecolor=blue, urlcolor=blue, citecolor=blue, pdftex=true, plainpages=false]{hyperref}\n";
print $Fpdf "\\setlength{\\textheight}{24cm}\n\\setlength{\\textwidth}{17cm}\n\\setlength{\\hoffset}{-2.5cm}\\setlength{\\voffset}{-2.5cm}\n";
print $Fpdf "\\begin{document}\n";
print $Fpdf CreateInfo($direc, 2);
print $Fpdf "\\clearpage\n";
print $Fpdf $tex."\n";
print $Fpdf "\\end{document}";
close $Fpdf;

exit(0);

my $pval_cut = 0.0;

# Get all output files in the specified folder.
my @f_list;
my @outputstr = ("$path/$direc"."/wflow*"."_*.sh*.out", "$path/$direc"."/wflow*"."*.bungee.bu.edu.out");

foreach (@outputstr) {
    @f_list = glob $_ unless $#f_list >= 0;
}
if ($#f_list <0) {
    print "Error: no files found. Exiting.\n";
    exit(9);
}

my $found = 0;
my %file_list;
foreach my $file (@f_list) {
    #print $file."\n";
    $_ = `grep "configuration:" $file`;
    if (/_(\d*)_scidac/) {
        if (exists $file_list{$1}) {
            print "Trajectory numbers not unique! Multiple offsets $1. Remove conflicting output files. Exiting\n"; exit(4);
        }
        else {
            $file_list{$1} =  $file;
            $found++;
        }
    }
}

if ($found == 0) {
    print "Error no configurations found!\n";
    exit(66);
}
my %M;
foreach my $i_file (keys(%file_list)) {
   my $indiv_file = $file_list{$i_file};
   #print "Parsing $indiv_file\n";
   $M{$i_file} = ParseFile($indiv_file, $i_file);
}
# sort hash of hashes by trajecoty number
my @config = sort {$a <=>$b} (keys(%M));


# Make a folder if it does not exist.
if(! (-d "$path/$run"."/data")) {
   system("mkdir -m 775 $path/$run"."/data");
}



# open output files
open (my $Fwflow, ">$path/$run"."/data/wflow.dat")  or die "Could not open file: $!"; chmod 0664, $Fwflow;

# write measurement data
foreach (@config) {
    if ($_ >= $cut) {
       print $Fwflow $M{$_}."\n";
    }
} 
close $Fwflow;


#$cut = floor($Ntraj/2) if ($Ntraj - $cut < $cut);


system("tar -cjf  $path/$direc"."/plots/$direc.tar.bz2 $path/$direc"."/README $path/$direc"."/plots/*pdf $path/$direc"."/plots/*txt; chmod 664 $path/$direc"."/plots/*.bz2 $path/$direc"."/plots/*.pdf $path/$direc"."/plots/*.txt");
print "Done2.\n";



##############################################################################################
##############################################################################################

##############################################################################################

sub getMass {
	my ($basedir, $direc, $nt) = @_;
	
	# Build a list of extensions to look for.
	my @states = ("ps1", "sc1", "sc2", "i51", "ij1", "r01", "ris1", "rij1", "ri51", "nu1", "de1");
	
	my $mpi = 0.0;
	my $dmpi = 0.0;
	
	my %spec_data = ();
	# Prepare a hash of outputs.
	foreach my $item (@states)
	{
	   $spec_data{ $item } = ();
	}
	
	my $info_path = $basedir."spectrum/spectrum.";
	
	my $t_min = 0; 
	my $t_max = $nt/2;
	
	switch ($nt)
	{
		case "32" {$t_min = 5; $t_max = 13;}
		case "40" {$t_min = 6; $t_max = 17;}
		case "48" {$t_min = 7; $t_max = 21;}
		case "64" {$t_min = 9; $t_max = 29;}
		else {$t_min = $nt/4-3; $t_max = $nt/2-3;}
	}
	
	foreach my $state (@states)
	{
		#print $info_path.$state."\n";
		open(my $file, "<".$info_path.$state);
		my @file_contents = <$file>;
		close($file);

		# This part reproduces Anna's get_masses16 script,
		# but with better control over the output.

		my $avg = 0.0;
		my $err = 0.0;
		my $cnt = 0;
		foreach my $line (@file_contents)
		{
			my @infos = split(' ', $line);
			# if 5 < t_min < 13 and we don't have nan's and err < val/3 and
			# our pval isn't too low.
			if (@infos > 3 && $infos[0] > $t_min && $infos[0] < $t_max && $infos[2] ne "nan" && $infos[2] ne "-nan" && $infos[2] < $infos[1]/3.0 && $infos[4] > $pval_cut)
			{
				$avg += $infos[1];
				$err += $infos[2]**2;
				$cnt++;
			}
		}
		if ($cnt > 0)
		{
			$avg /= $cnt;
			$err = sqrt($err/$cnt);
			
			if ($state eq "ps1")
			{
				$mpi = $avg;
				$dmpi = $err;
			}
			
			push(@{$spec_data{$state}}, "$avg $err");
		}
		else
		{
			delete $spec_data{$state};
		}
	}
	
	foreach my $key (keys %spec_data)
	{
		open(my $fileout, ">$basedir"."hadrons/$key".".txt");
		if (@{$spec_data{$key}} > 0)
		{
			print "$key $spec_data{$key}[0]\n";
			print $fileout $spec_data{$key}[0];
		}
		else
		{
			print $fileout "0.0 0.0";
		}
		close($fileout);
	}
	
	return ($mpi, $dmpi);
}

sub getFpi {
	my ($basedir, $direc, $nt, $mass) = @_;
	
	# Build a list of extensions to look for.
	my @states = ("fp", "pp1");
	
	my $mpi = 0.0;
	my $dmpi = 0.0;
	
	my %spec_data = ();
	# Prepare a hash of outputs.
	foreach my $item (@states)
	{
	   $spec_data{ $item } = ();
	}
	
	my $info_path = $basedir."spectrum/spectrum.";
	
	my $t_min = 0; 
	my $t_max = $nt/2;
	
	switch ($nt)
	{
		case "32" {$t_min = 9; $t_max = 13;}
		case "40" {$t_min = 13; $t_max = 17;}
		case "48" {$t_min = 17; $t_max = 21;}
		case "64" {$t_min = 25; $t_max = 29;}
		else {$t_min = $nt/2-7; $t_max = $nt/2-3;}
	}
	
	foreach my $state (@states)
	{
		#print $info_path.$state."\n";
		open(my $file, "<".$info_path.$state);
		my @file_contents = <$file>;
		close($file);

		# This part reproduces Anna's get_masses16 script,
		# but with better control over the output.

		my $avg = 0.0;
		my $err = 0.0;
		my $avgsq = 0.0;
		my $cnt = 0;
		foreach my $line (@file_contents)
		{
			my @infos = split(' ', $line);
			# if 5 < t_min < 13 and we don't have nan's and err < val/3 and
			# our pval isn't too low.
			if ($state eq "fp")
			{
				if (@infos > 3 && $infos[0] > $t_min && $infos[0] < $t_max && $infos[2] ne "nan" && $infos[2] ne "-nan" && $infos[2] < $infos[1]/3.0 && $infos[11] > $pval_cut)
				{
					$avg += $infos[1];
					$err += $infos[2]**2;
					$avgsq += $infos[1]**2;
					$cnt++;
				}
			}
			elsif ($state eq "pp1")
			{
				if (@infos > 3 && $infos[0] > $t_min && $infos[0] < $t_max && $infos[4] ne "nan" && $infos[4] ne "-nan" && $infos[4] < $infos[3]/3.0 && $infos[5] > $pval_cut)
				{
					$avg += $infos[3];
					$err += $infos[4]**2;
					$cnt++;
				}
			}
		}
		
		if ($cnt > 0)
		{
			if ($state eq "fp")
			{
				$avg /= $cnt;
				$err = sqrt(($avgsq/$cnt-$avg**2)+$err/$cnt);
				push(@{$spec_data{$state}}, "$avg $err");
			}
			elsif ($state eq "pp1")
			{
				my $c = $mass*$mass*(($nt/2)**3);
				$avg = $c*$avg/$cnt;
				$err = $c*sqrt($err/$cnt);
				my $newavg = $avg**(0.2);
				$err = ($avg**0.2)/$avg*$err/5.0;
				push(@{$spec_data{$state}}, "$newavg $err");
			}
		}
	}
	
	foreach my $key (keys %spec_data)
	{
		open(my $fileout, ">$basedir"."hadrons/$key".".txt");
		if (@{$spec_data{$key}} > 0)
		{
			print "$key $spec_data{$key}[0]\n";
			print $fileout $spec_data{$key}[0];
		}
		else
		{
			print $fileout "0.0 0.0";
		}
		close($fileout);
	}
}

sub ParseFile{
   my ($indiv_file, $config_no) =  @_;     
   # only print Wilson flow time and symmE
   $_ = `grep Step $indiv_file | awk \'{ print "$config_no:"\$2":" \$7 }\' `;
   s/:/ /g;
   chomp($_);
   return $_;
}
