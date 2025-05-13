set terminal pngcairo
set output 'histogram_3.png'
set title 'Время ожидания'
set xlabel 'x_i'
set ylabel 'n_i/n'
set style fill solid
set yrange[0:*]
set boxwidth 3
plot '-' using 1:2 with boxes title 'Частота'
1.5 0.51506
4.5 0.283133
7.5 0.0933735
10.5 0.0451807
13.5 0.063253
e
