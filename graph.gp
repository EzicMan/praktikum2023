set terminal png size 2048,1024
set output "graphxx.png"
set title "sigma_{xx} angle = 45"
plot [0:20]"graphxx.txt" u 1:2 with l lw 3 t "Numerical", "graphxxanalit.txt" u 1:2 with l lw 3 t "Analytical"

set terminal png size 2048,1024
set output "graphxy.png"
set title "sigma_{xy} angle = 45"
plot [0:20][-6:1]"graphxy.txt" u 1:2 with l lw 3 t "Numerical", "graphxyanalit.txt" u 1:2 with l lw 3 t "Analytical"

set terminal png size 2048,1024
set output "graphyy.png"
set title "sigma_{yy} angle = 45"
plot [0:20][-5:15]"graphyy.txt" u 1:2 with l lw 3 t "Numerical", "graphyyanalit.txt" u 1:2 with l lw 3 t "Analytical"
