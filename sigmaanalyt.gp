set terminal png size 2048,1024
set output "sigmaxxanalyt.png"
set title "sigma_{xx} angle = 0"
plot [0:20][0:4] "sigmaxx.txt" u 1:2 with l lw 3 t "Numerical", "analytxx.txt" u 1:2 with l lw 3 t "Analytical"

set terminal png size 2048,1024
set output "sigmaxyanalyt.png"
set title "sigma_{xy} angle = 0"
plot [0:20][-0.2:0.2] "sigmaxy.txt" u 1:2 with l lw 3 t "Numerical", "analytxy.txt" u 1:2 with l lw 3 t "Analytical"

set terminal png size 2048,1024
set output "sigmayyanalyt.png"
set title "sigma_{yy} angle = 0"
plot [0:20][9:35] "sigmayy.txt" u 1:2 with l lw 3 t "Numerical", "analytyy.txt" u 1:2 with l lw 3 t "Analytical"
