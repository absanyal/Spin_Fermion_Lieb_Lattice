set terminal pdf
set output "dos.pdf"

unset key

plot 'dos.dat' u 1:2 w l

unset output

