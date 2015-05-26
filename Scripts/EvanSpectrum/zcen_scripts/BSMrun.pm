package BSMrun;
use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

our $VERSION = '0.1';
use base 'Exporter';
our @EXPORT = qw(name2param namestar CreateInfo CreateComments);

sub name2param{
#   Sets hash with values: flavor, L, T, beta, ml, mh, and path
#   may return empty values for L, T. beta or masses
    my ($name) = @_;
    my $check = -1;
    my %params;

    if ($name =~ m/k(.+)l(\d*)t(\d*)b(\d*)/) {
	$check = 0;
        $params{group} = 8; # LatKMI
        $params{flavor} =  $1; 
        $params{L} = $2;
        $params{T} = $3;
        $params{beta} = looks_like_number($4) ? sprintf("%3.1f",$4/10) : '';
        if ($params{flavor} eq '8') {
	    $check += 1;
            $params{path} = '../Eight';
            if ($name =~ m/m(\d*)/) {
               $params{ml} = "0.".$1;
               $check += 2;
            }
	}
    }
    elsif ($name =~ m/f(.+)l(\d*)t(\d*)b(\d*)/) {
	$check = 0;
        $params{flavor} =  $1; 
        $params{L} = $2;
        $params{T} = $3;
        $params{beta} = looks_like_number($4) ? sprintf("%3.1f",$4/10) : '';

        if ($params{flavor} eq '4plus8') {
	    $check += 1;
	    #print "T=$params{T} == $params{L} =L\n";
            if (looks_like_number($params{T}) && ($params{T} == $params{L}/2)) {
                print "FiniteT\n";
		$params{path} = '../FiniteT'; 
	    }
	    else {
		$params{path} = '../FourPlusEight';
	    }
            $params{group} = 1; # 4+8
            if ($name =~ m/m(\d*)m(\d*)/) {
               $params{ml} = "0.".$1;
               $params{mh} = "0.".$2;
               $check += 2;
           }
        }
        elsif ($params{flavor} eq '4plus4') {
	    $check += 1;
            $params{path} = '../FourPlusFour';
            $params{group} = 1; # 4+8
            if ($name =~ m/m(\d*)m(\d*)/) {
               $params{ml} = "0.".$1;
               $params{mh} = "0.".$2;
               $check += 2;
            }
        }
        elsif ($params{flavor} eq '12') {
	    $check += 1;
            $params{path} = '../Twelve';
            $params{group} = 2; # 4+8 + Boulder
            if ($name =~ m/m(\d*)/) {
               $params{ml} = "0.".$1;
               $check += 2;
            }
        }
        elsif ($params{flavor} eq '8') {
	    $check += 1;
            $params{path} = '../Eight';
            $params{group} = 4; # LSD
            if ($name =~ m/m(\d*)/) {
               $params{ml} = "0.".$1;
               $check += 2;
            }
        }
        elsif ($params{flavor} eq '4') {
           $check += 1;
           $params{path} = '../Four';
            $params{group} = 5; # 4+8 + LSD
            if ($name =~ m/m(\d*)/) {
               $params{ml} = "0.".$1;
               $check += 2;
            }
        }
    }
    if (-e $params{path}."/".$name."/cut") {
	my $cut = `cat $params{path}/$name/cut`;
        chomp($cut);
        $params{cut} = $cut;
    }
    else {
        $params{"cut"} = 200;
    }
#    print  "$params{'flavor'}   $params{'L'} \n";
    return %params unless $check < 3;
    if ($check == -1) {
	print "[name2param] Error: cannot decompose name.\n";
    }
    elsif ($check == 0) {
	print "[name2param] Error: cannot process flavors.\n";
    }
    elsif ($check == 1) {
	print "[name2param] Error: cannot determine masses.\n";
    }
    exit(90+$check);
}

sub namestar{
    # insert '*' for omitted parameters
    my ($star) =  @_;
    $star =~ s/lt/l*t/;
    $star =~ s/tb/t*b/;
    $star =~ s/bm/b*m/;
    $star =~ s/mm/m*m/;
    $star =~ s/m$/m*/;
    return $star;
}

sub CreateInfo {
    my ($name, $ana) = @_;
    my %params = name2param($name);
    my $str =  <<END;
\\begin{figure}[p]
\\vspace{-5mm} {\\Large{\\textbf{$params{flavor} flavor at \$\\beta = $params{beta}\$\\\\[3mm]}\\pdfbookmark[1]{Overview}{Overview}}}
\\begin{minipage}{\\textwidth}
\\begin{tabular}{ll} 
 \$\\bullet\$& Staggered fermions with nHyp smeared gauge links (\$\\alpha = (0.4,0.5,0.5)\$ [FUEL convention])\\\\ 
 \$\\bullet\$& adjoint plaquette action with \$\\beta = $params{beta},\\, \\beta_a=-\\beta/4\$\\\\ 
 \$\\bullet\$& \$$params{L}^3 \\times $params{T}\$\\\\ 
END
    if ($params{flavor} =~ m/4plus8/) {
	$str = $str.<<END;
 \$\\bullet\$& \$m_l = $params{ml}\$: mass of the four light flavors \\\\ 
 \$\\bullet\$& \$m_h\$: mass of the eight heavy flavors \\\\ 
END
    }
    else {
	$str = $str." \$\\bullet\$& \$m_l\$: mass of the light flavors \\\\ \n";
    }
    $str = $str." \$\\bullet\$& Simulations performed with \\texttt{qhmc} [FUEL] \\\\ \n";
    if ($ana == 1) {
	$str = $str."\$\\bullet\$& Data analyzed with UWerr.m (cf.~ Ulli Wolff, hep-lat/0306017) on \\texttt{scc}\\\\ \n";
    }
    elsif ($ana == 2) {
        $str = $str."\$\\bullet\$& Data analyzed with Anna's analysis scripts on \\texttt{scc}\\\\ \n";
    }
    $str = $str.<<END;
\\end{tabular}
\\end{minipage}
\\end{figure}
END
    return $str;
}


sub CreateComments{
    my ($str) =  @_;
    my @comments = split('#', $str);
    $str = '';
    foreach my $f (@comments) {
       my $run = $2 if ($f =~ m+/(.*)/(.*)/+);
       my ($cut, $trajectories, $acceptance, $machine, $init);
       my $hmc = '';
       my @items;
       open (my $info, "<$f") or die "cannot open $!";
       while (<$info>) {
           chomp($_);
	   if (/trajectory/) {
	       @items = split('->', $_);
               $hmc = $hmc.$items[1]."#";
               $hmc =~ s+ml/Ml+\$m_l/M_l\$+;
               $hmc =~ s+mh/Mh+\$m_h/M_h\$+;
           }
           if (/cut/) {
     	       @items = split('->', $_);
               $cut = $items[1];
	   }
           if (/trajectories/) {
     	       @items = split('->', $_);
               $trajectories = $items[1];
	   }
           if (/acceptance/) {
     	       @items = split('->', $_);
               $acceptance = $items[1];
               $acceptance =~ s/%/\\%/;
	   }
           if (/machine/) {
     	       @items = split('->', $_);
               $machine = $items[1];
	   }
           if (/InitialConfig/) {
     	       @items = split('->', $_);
               $init = $items[1];
	   }

       }
       close($info);
       $str = $str.$run." & ".$acceptance." & ".$trajectories." & ".$cut." & ".$machine."\\\\ & &\\multicolumn{3}{l}{".$init."}\\\\ \n";
       @items = split('#',$hmc);
       foreach (@items) {
          $str = $str." & &\\multicolumn{3}{l}{".$_."}\\\\ \n";
       }
    }
    my $str1 = <<END;
\\begin{figure}[p]
\\vspace{-5mm} {\\Large{\\textbf{Comments}\\\\[3mm]} \\pdfbookmark[1]{Comments}{Comments}}
 \\begin{minipage}{\\textwidth}
 \\begin{center}
 \\begin{tabular}{ccccl}
 \\hline
run     & acc. &  traj. & cut & machine  \\\\ \\hline
END
    my $str2 = <<END;
\\hline
\\end{tabular}
\\end{center}
\\end{minipage}
\\end{figure}
END
    return $str1.$str.$str2;
}

1;
