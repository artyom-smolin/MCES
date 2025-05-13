set terminal pngcairo
set output 'histogram.png'
set title 'Гистограмма частот'
set xlabel 'x_i'
set ylabel 'n_i/n'
set style fill solid
set boxwidth 0.845
plot '-' using 1:2 with boxes title 'Частота'
2.5525 0.0789474
3.3975 0.157895
4.2425 0.210526
5.0875 0.289474
5.9325 0.157895
6.7775 0.105263
e
