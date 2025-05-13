set terminal pngcairo
set output 'histogram_2.png'
set title 'Поток машин'
set xlabel 'x_i'
set ylabel 'n_i/n'
set style fill solid
set yrange[0:*]
set boxwidth 0.75
plot '-' using 1:2 with boxes title 'Частота'
0.375 0.72
1.125 0.23
1.875 0.047
2.625 0.003
e
