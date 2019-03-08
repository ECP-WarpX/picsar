set term pngcairo size 800,800 font "Helvetica,24"

set xlabel "x"
set ylabel "y"
set zlabel "z"

set cblabel "optical depth"
set cbrange [0:12]

set view equal xyz
set xyplane 0

set xrange [-6:6]
set yrange [-6:6]
set zrange [-6:6]

unset key

do for [t=0:100] {
  outfile = sprintf('out_%03.0f.png',t)
  tit = sprintf('time : %04.3f',t*0.05)
  set title tit
  set output outfile
  splot 'out_phot1_'.t.'.dat' u 1:2:3:13 w p ps 0.8 pt 7 lc palette
}
