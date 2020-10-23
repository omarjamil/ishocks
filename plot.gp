reset
set pm3d map
set terminal postscript eps enhanced color
set output "BLF_diff.eps"
#set key bottom right
set cbrange[0:0.3]
set xlabel 'BLF 1'
set ylabel 'BLF 2'
splot "gamma_1_3.dat" u 1:2:($4-$3) title 'BLF(approx.) - BLF(exact)'
set term wxt

reset
set pm3d map
set terminal postscript eps enhanced color
set output "Efficiency.eps"
set cbrange[0:0.3]
#set key bottom right
set xlabel 'BLF 1'
set ylabel 'BLF 2'
splot "gamma_1_3.dat" u 1:2:5 title 'Efficieny'
set term wxt

reset
set pm3d map
set terminal postscript eps enhanced color
set output "Efficiency_diff.eps"
#set key bottom right
set cbrange[0:0.14]
set xlabel 'BLF 1'
set ylabel 'BLF 2'
splot "gamma_1_3.dat" u 1:2:($5-$6) title 'Efficieny(exact) - Efficiency(approx.)'
set term wxt

reset
set pm3d map
set terminal postscript eps enhanced color
set output "Efficiency_1075_113.eps"
#set key bottom right
set cbrange[0:0.0018]
set xlabel 'BLF 1'
set ylabel 'BLF 2'
splot "gamma_1075_113.dat" u 1:2:5 title 'Efficieny'
set term wxt