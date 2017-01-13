#! /usr/bin/perl -w 
use warnings;
use strict;
use POSIX;
use BSMrun;
use Switch;

# call: ./spectrum_sigma2_parse.pl f4plus8l24t48b30m005m020 20 20

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


# First, make sure the 'sigma2' folder exists. If not, error out.

if (!(-d "$path/$direc"."/sigma2"))
{
	print "Error: a sigma2 directory does not exist. Exiting.\n";
	exit(200);
}

# Next, make sure the 'spectrum2' folder exists. If not, error out.

if (!(-d "$path/$direc"."/spectrum2"))
{
	print "Error: a spectrum2 directory does not exist. Run the connected analysis first. Exiting.\n";
	exit(200);
}

# Next, let's make a 'stoch' folder to put what we parsed.

if (!(-d "$path/$direc"."/spectrum2/stoch"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/stoch";
        `$command`;
}

# And a 'vev' folder!

if (!(-d "$path/$direc"."/spectrum2/vev"))
{
        my $command = "mkdir $path/$direc"."/spectrum2/vev";
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


my $a_path = "$path/$direc"."/sigma2/glu*";

print $a_path."\n";

my @allfiles = glob( $a_path );

my @filecontents = ();

my @m_pbp = ();
my @m_conn = ();
my @m_time = ();

foreach my $file (@allfiles)
{
   #print $file."\n";
   my @splitpieces = split('/', $file);
   #print $splitpieces[1]."\n";

   my @argh = split('\.', $splitpieces[4]);
   
   my $n_config = 0;
   my $machine = "";

   # Check for the scc.
   if ($argh[1] eq "scc" || $argh[1] eq "vulcan")
   {
      my @argh2 = split('-', $argh[0]);
      #print $argh2[1]."\n";
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
      print "Out of order: $orig_n_config\n";
      next; 
   }

   open(my $filething, "<".$file);
   my @parseit = <$filething>;
   close($filething);

   my @m_pbp_lines = ();
   my @m_conn_lines = ();
   my @m_time_lines = ();

   my @lines = ();
   my $flag_pbp = 0;
   my $flag_conn = 0;

   foreach my $line (@parseit)
   {
      chomp($line);
      if ($line eq "START_PBPPART" || $line eq "BEGIN_PBPPART")
      {
         #print "Got!\n";
         $flag_pbp = 1;
      }
      elsif ($line eq "START_SPECTRUM" || $line eq "BEGIN_SPECTRUM")
      {
         $flag_conn = 1;
      }
      elsif ($flag_pbp == 1)
      {
         if ($line eq "END_PBPPART")
         {
            $flag_pbp = 2;
         }
         else
         {
            push(@m_pbp_lines, $line);
         }
      }
      elsif ($flag_conn == 1)
      {
         if ($line eq "END_SPECTRUM")
         {
            $flag_conn = 2;
         }
         else
         {
            push(@m_conn_lines, $line);
         }
      }
      else
      {
         my @splitted_vals = split(' ', $line);
         
         if (@splitted_vals > 1)
         {
            if ($splitted_vals[0] eq "Job")
            {
               push(@m_time_lines, "# time-start ".$parseit[1]);
            }
            elsif ($splitted_vals[1] eq "time-start")
            {
               #print $line."\n";
               push(@m_time_lines, $line);
            }
            elsif ($splitted_vals[1] eq "time-finish")
            {
               #print $line."\n";
               push(@m_time_lines, $line);
            }
         }
      }

   }

   if (!($flag_pbp == 2 && $flag_conn == 2)) # && @m_time_lines == 2))
   {
      print "Incomplete file: $orig_n_config\n";
   }
   else
   {
      $m_pbp[$n_config] = \@m_pbp_lines;
      $m_conn[$n_config] = \@m_conn_lines;
      
      if (@m_time_lines == 2)
      {
      my @begin_split = split(' ', $m_time_lines[0]);
      my @end_split = split(' ', $m_time_lines[1]);

      my @begin_break = split(':', $begin_split[5]);
      my @end_break = split(':', $end_split[5]);

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

      $m_time[$n_config] = $diff_time." ".$machine;
      }
      else
      {
         $m_time[$n_config] = "0 ".$machine;
      }
      #print $begin_split[5]." ".$end_split[5]."\n";

      #print $diff_time."\n";
   }

}


open(my $pbpfile, ">$path/$direc"."/spectrum2/stoch/PBPPART.dat");
open(my $connfile, ">$path/$direc"."/spectrum2/stoch/SPECTRUM.dat");
open(my $timefile, ">$path/$direc"."/spectrum2/stoch/times.dat");

for (my $i = 0; $i < @m_pbp; $i++)
{
   if (exists $m_pbp[$i])
   {
      
	  my $good_flag = 1;
	  foreach my $line (@{$m_pbp[$i]})
	  {
	    my @tmp = split('\t', $line);
		my $tmp_size = @tmp;
		if ($tmp_size != 4 && $good_flag == 1)
		{
		  $good_flag = 0;
		  print "Nope! ".($i*$sep+$start)." disc pbppart\n";
		}
	  }
	  
	  foreach my $line (@{$m_conn[$i]})
	  {
	    my @tmp = split('\t', $line);
		my $tmp_size = @tmp;
		if ($tmp_size != 2 && $good_flag == 1)
		{
		  $good_flag = 0;
		  print "Nope! ".($i*$sep+$start)." disc spectrum\n";
		}
	  }
	  
	  if ($good_flag == 1)
	  {
	  
		foreach my $line (@{$m_pbp[$i]})
		{
		 print $pbpfile ($i*$sep+$start)." ".$line."\n";
		}

		foreach my $line (@{$m_conn[$i]})
		{
		 print $connfile ($i*$sep+$start)." ".$line."\n";
		}

		my $time_thing = $i*$sep+$start;
		print $timefile "$time_thing $m_time[$i]\n";
		}
   }
   else
   {
      print "NOPE! ".($i*$sep+$start)."\n";
   }
}

close($pbpfile); close($connfile); close($timefile);

