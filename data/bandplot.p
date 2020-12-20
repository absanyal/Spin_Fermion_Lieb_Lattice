#set terminal pdf
#set output "band.pdf"
#set pm3d
#set palette define (0 "red", 1 "blue")
#set pm3d corners2color c1

unset key

set xr [0.5:1.5]
set yr [-1.0:1.0]

plot 'band_ky_pi' u 1:2 w l, 'band_ky_pi' u 1:3 w l, 'band_ky_pi' u 1:4 w l, 'band_ky_pi' u 1:5 w l, 'band_ky_pi' u 1:6 w l, 'band_ky_pi' u 1:7 w l

#unset output
unset xr
unset yr
