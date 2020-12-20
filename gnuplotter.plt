set terminal pdf
set output "sq_6x6_AFM.pdf"

unset key

#set palette defined (0 "white", "1" red)
set pm3d

splot "Sq.txt" u 1:2:3 w pm3d

unset output
unset pm3d