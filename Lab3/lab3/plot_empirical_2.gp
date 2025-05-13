set terminal pngcairo
set output 'empirical_function_2.png'
set title 'Поток машин (накопительная функция)'
set xlabel 'x_i'
set ylabel 'F*'
plot '-' with lines title 'F*'
0 0.72
0.75 0.95
1.5 0.997
2.25 1
e
