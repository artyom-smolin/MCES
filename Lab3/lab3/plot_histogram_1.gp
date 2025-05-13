set terminal pngcairo
set output 'histogram_1.png'
set title 'Время обслуживания'
set xlabel 'x_i'
set ylabel 'n_i/n'
set style fill solid
set yrange[0:*]
set boxwidth 0.0198975
plot '-' using 1:2 with boxes title 'Частота'
0.9107 0.106317
0.930597 0.103236
0.950495 0.0909091
0.970392 0.118644
0.99029 0.0939908
1.01019 0.100154
1.03008 0.11094
1.04998 0.0847458
1.06988 0.0909091
1.08978 0.100154
e
