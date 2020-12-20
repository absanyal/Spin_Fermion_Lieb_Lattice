set terminal pdf
set output "Sq_test.pdf"

unset key

#set pm3d map
set palette define (0 "red", 1 "blue")
set pm3d corners2color c1

set view 50, 100

splot "Sq.txt" u 1:2:3 w pm3d

unset output
unset view
unset pm3d
