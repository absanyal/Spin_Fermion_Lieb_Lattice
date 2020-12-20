#set terminal pdfcairo size 12, 6
#set output "path_B_x_lR_x.pdf"

set xtics \
('{/Symbol G}' 0, \
'X' 200, \
'M' 400, \
'{/Symbol G}' 600)

set key off

set xr [1:10]
set yr [-3:-2]

set arrow from 200, graph 0 to 200, graph 1 nohead dt 2
set arrow from 400, graph 0 to 400, graph 1 nohead dt 2

plot \
'band_path.dat' u 1:4 w l lt rgb "black" lw 3, \
'band_path.dat' u 1:5 w l lt rgb "black" lw 3, \
'band_path.dat' u 1:6 w l lt rgb "black" lw 3, \
'band_path.dat' u 1:7 w l lt rgb "black" lw 3, \
'band_path.dat' u 1:8 w l lt rgb "black" lw 3, \
'band_path.dat' u 1:9 w l lt rgb "black" lw 3

unset arrow
unset xr
unset yr
unset xtics
#unset terminal
#unset output
unset key