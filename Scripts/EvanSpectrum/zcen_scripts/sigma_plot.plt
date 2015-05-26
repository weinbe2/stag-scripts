set style line 1 lt 1 lc rgb '#000000' pt 11 #ps 1.5
set style line 2 lt 1 lc rgb '#B8860B' pt 9 #ps 1.5
set style line 3 lt 1 lc rgb '#FF0000' pt 7 #ps 1.5
set style line 4 lt 1 lc rgb '#0000FF' pt 5 #ps 1.5
set style line 5 lt 1 lc rgb '#008800' pt 13 #ps 1.5
set style line 11 lt 1 lc rgb '#555555' pt 10 #ps 1.5
set style line 12 lt 1 lc rgb '#D0A015' pt 8 #ps 1.5
set style line 13 lt 1 lc rgb '#FF5555' pt 6 #ps 1.5
set style line 14 lt 1 lc rgb '#5555FF' pt 4 #ps 1.5
set style line 15 lt 1 lc rgb '#558855' pt 12 #ps 1.5

set key bottom left Left reverse

set terminal epslatex color standalone

set output "m0pp_zcen_amp.tex"
set yrange [-%AMP%:%AMP%]
set ylabel "Amplitude"
set title "Vacuum: %VEV%"
plot "./sg_stoch_zcen.dat" using ($1-0.1):6:7 with yerrorbars ls 3 title "%NF%D-C $c_{0^{++}}$", "./sg_stoch_zcen.dat" using ($1+0.1):(-$14):15 with yerrorbars ls 4 title "%NF%D-C $-v$"
set output


set output "m0pp_zcen_amp_dc.tex"
set yrange [-%AMP%:%AMP%]
set ylabel "Amplitude"
set title "Vacuum: %VEV%"
plot "./dc_stoch_zcen.dat" using ($1-0.1):6:7 with yerrorbars ls 3 title "%NF%D $c_{0^{++}}$", "./dc_stoch_zcen.dat" using ($1+0.1):(-$14):15 with yerrorbars ls 4 title "%NF%D $-v$"
set output

unset title

set yrange [0.0001:0.75]
set key bottom left Left reverse
set terminal epslatex color standalone
set output "m0pp_zcen_cmp.tex"
set xlabel "$t_{min}$"
set ylabel "$aM_{0^{++}}$"
plot "./sg_stoch_zcen.dat" using ($1-0.1):8:9 with yerrorbars ls 2 title "%NF%D-C $M_{0^{++}}$", "./dc_stoch_zcen.dat" using ($1+0.1):8:9 with yerrorbars ls 5 title "%NF%D $M_{0^{++}}$
set output

set output "m0pp_zcen_amp_orig.tex"
set yrange [0:1]
set key top right Left reverse
set ylabel "Amplitude at $t=0$"
plot "./sg_stoch_zcen.dat" using ($1-0.1):16:17 with yerrorbars ls 2 title "%NF%D-C $c_{0^{++}}$cosh$(m_{0^{++}}T/2)$", "./dc_stoch_zcen.dat" using ($1+0.1):16:17 with yerrorbars ls 5 title "%NF%2D $c_{0^{++}}$cosh$(m_{0^{++}}T/2)$"
set output


set xrange [0.0001:16.5]
set yrange [0.0001:0.75]
set key bottom left Left reverse
set terminal epslatex color standalone
set output "m0pp_zcen_eff.tex"
set xlabel "$t_{min}$"
set ylabel "$aM_{0^{++}}^{eff}$"
plot "./effmass.sg_stoch_kuti1" using ($1-1-0.1):2:3 with yerrorbars ls 2 title "%NF%D-C $M_{0^{++}}$", "./effmass.dc_stoch_kuti1" using ($1-1+0.1):2:3 with yerrorbars ls 5 title "%NF%D $M_{0^{++}}$
set output


set log y


set output "m0pp_zcen_fitlines.tex"
set format y "$10^{%T}$"
set ylabel "$%CORRNAME%$ Correlator Plot"
set xlabel "t"
set xrange [0.5:%NTD2%.5]
set yrange [%CORRMIN%:%CORRMAX%]
plot "./sum.%CORRTYPE%_stoch" using 1:($2-(%CORRCEN%)-(%CORRCONST%)):3 with yerrorbars ls 1 title "Correlator", %CORRFUNC% ls 2 title "Fit Curve, $t_{min} = %CORRTMIN%$"
set output



