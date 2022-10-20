i=1
my_terminal = 'aqua'
result_file="BouncingBallTS-BinaryCZM.dat"
my_terminal = 'wxt'
unset multiplot

set term my_terminal i
set multiplot
set size 0.4,0.4		
set origin 0.1,0.1
plot result_file u 1:2 w l t 'x position ball'
set origin 0.1,0.5
plot result_file u 1:3 w l t 'y position ball'
set origin 0.5,0.1
plot result_file u 1:4 w l t 'x velocity ball'
set origin 0.5,0.5
plot result_file u 1:5 w l t 'y velocity ball'
unset multiplot


i=i+1
set term my_terminal i
set multiplot
set size 0.4,0.4		
set origin 0.1,0.1
plot result_file u 1:7 w l t 'y_N'
set origin 0.1,0.5
plot result_file u 1:8 w l t 'y_T'
set origin 0.5,0.1
plot result_file u 1:9 w l t 'lambda_N'
set origin 0.5,0.5
plot result_file u 1:10 w l t 'lambda_T'


unset multiplot


i=i+1
set term my_terminal i
set multiplot
set size 0.4,0.4		
set origin 0.1,0.1
set yrange [-0.1:1.1]
plot result_file u 1:11 w l t 'beta'
		
set origin 0.1,0.5
#set yrange [0:1.1]
set yrange [-12:3]
plot result_file u 1:12 w l t 'applied force'


set origin 0.0,0.0
set size 1.0,1.0	

unset multiplot
