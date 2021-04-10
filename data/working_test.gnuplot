set terminal pngcairo size 800,600 font 'Noto Serif'
set output 'test.png'
set xrange [0:10]
set yrange [0:10]
set xtics out nomirror
set ytics out nomirror
plot '-' with lines title 'data'
    1 7
    5 1
    9 4
    4 5
    3 3
    e
