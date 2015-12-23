set term png
set output "visual.png"


set logscale x 2
set logscale y 

p "res" u 1:2 w l, "res" u 1:3 w l, "res" u 1:4 w l
