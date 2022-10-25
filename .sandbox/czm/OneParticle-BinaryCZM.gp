i=1
my_terminal = 'aqua'
result_file="OneParticle-BinaryCZM.dat"
#my_terminal = 'wxt'
unset multiplot

set xrange [*:*]


threshold(x) = 1e-04

set term my_terminal i
set multiplot
set size 0.8,0.4		
set origin 0.1,0.1
set yrange [-1e-04:3e-4]
plot result_file u 1:2 w l t 'x position ball',  threshold(x)
#set origin 0.1,0.5
#plot result_file u 1:3 w l t 'y position ball'
set origin 0.1,0.5
set yrange [-1e-04:1e-03]
plot result_file u 1:4 w l t 'x velocity ball'
#set origin 0.5,0.5
#plot result_file u 1:5 w l t 'y velocity ball'
unset multiplot

delta_c = 1e-04
threshold(x) = delta_c

i=i+1
set term my_terminal i
set multiplot
set size 0.8,0.25		
set origin 0.1,0.1
set yrange [-1e-4: 2e-4]
plot result_file u 1:7 w l t 'y_N', threshold(x)
set origin 0.1,0.3
set yrange [-2e-4: 1e-4]
plot result_file u 1:9 w l t '$lambda_N$'

set origin 0.1,0.5
set yrange [-0.1:1.1]
plot result_file u 1:11 w l t 'beta'

set origin 0.1,0.7
set yrange [-1:2.0]
plot result_file u 1:12 w l t 'applied force'




unset multiplot


i=i+1
set term my_terminal i
set multiplot
set size 0.4,0.4		
set origin 0.1,0.1
set yrange [-0.1:1.1]
plot result_file u 1:11 w l t 'beta'

set origin 0.5,0.1
plot result_file u 7:11 w l t 'beta vs y_N'
		
set origin 0.1,0.5
#set yrange [0:1.1]
set yrange [-1:2.0]
plot result_file u 1:12 w l t 'applied force'



set origin 0.0,0.0
set size 1.0,1.0	

unset multiplot
