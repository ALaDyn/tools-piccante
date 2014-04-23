#!/bin/bash

echo "set term pngcairo size 1000,1000">gnuplot.scr
echo "set size ratio -1">>gnuplot.scr

for f in TEST/EMfield_*.bin
	do ./reader $f
	echo "set output '${f}_Ex.png'">>gnuplot.scr
	echo "plot '${f}.txt' u 1:2:4 w image">>gnuplot.scr 

        echo "set output '${f}_Ey.png'">>gnuplot.scr
        echo "plot '${f}.txt' u 1:2:5 w image">>gnuplot.scr 

        echo "set output '${f}_Ez.png'">>gnuplot.scr
        echo "plot '${f}.txt' u 1:2:6 w image">>gnuplot.scr 

        echo "set output '${f}_Bx.png'">>gnuplot.scr
        echo "plot '${f}.txt' u 1:2:7 w image">>gnuplot.scr 

        echo "set output '${f}_By.png'">>gnuplot.scr
        echo "plot '${f}.txt' u 1:2:8 w image">>gnuplot.scr 

        echo "set output '${f}_Bz.png'">>gnuplot.scr
        echo "plot '${f}.txt' u 1:2:9 w image">>gnuplot.scr 
done

gnuplot gnuplot.scr
