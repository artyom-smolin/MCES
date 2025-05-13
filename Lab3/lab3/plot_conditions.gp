set terminal pngcairo
set output 'conditions.png'
set title 'Количество сломанных машин'
set xlabel 'count'
set ylabel 'p_i'
set style fill solid
set xrange[-1:12]
set yrange[0:0.525]
set boxwidth 1
plot '-' using 1:2 with boxes title 'Вероятность'
0.5 0.475
1.5 0.221
2.5 0.16
3.5 0.062
4.5 0.028
5.5 0.02
6.5 0.014
7.5 0.013
8.5 0.007
9.5 0
10.5 0
e
