set terminal pdfcairo size 12, 6
set output "dos_k_fixed_B_x_lR_x.pdf"
set key

set xrange[-5:5]

set xlabel "{/Symbol w}"
set ylabel "{/Symbol r}"

plot 'dos_k.dat' u 1:($2+$3+$4+$5+$6+$7) w l title "DoS"

unset terminal
unset output
unset key
unset xrange