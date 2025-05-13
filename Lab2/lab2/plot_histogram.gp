set terminal pngcairo
set output 'histogram.png'
set title 'Гистограмма частот'
set xlabel 'x_i'
set ylabel 'n_i/n'
set style fill solid
set boxwidth 380.575
plot '-' using 1:2 with boxes title 'Частота'
7814.88 0.003
8195.46 0.011
8576.03 0.066
8956.61 0.142
9337.18 0.214
9717.76 0.255
10098.3 0.189
10478.9 0.087
10859.5 0.029
11240.1 0.004
e
