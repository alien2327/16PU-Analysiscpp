reset

set terminal pngcairo size 1600, 1200
set output 'moment.png'

set multiplot
set datafile separator ","

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.1
set tmargin screen 0.45
set xlabel "turn"
set ylabel "dipole moment[mm]"
set format x
set format y
plot 'dipole_15.csv' using 0:1 with line title 'x #15', '' using 0:2 with line title 'y #15'

set lmargin screen 0.1
set rmargin screen 0.9
set bmargin screen 0.55
set tmargin screen 0.9
set xlabel "turn"
set ylabel "dipole moment[mm]"
set format x
set format y
plot 'dipole_13.csv' using 0:1 with line title 'x #13', '' using 0:2 with line title 'y #13'

unset multiplot
set terminal qt
set output