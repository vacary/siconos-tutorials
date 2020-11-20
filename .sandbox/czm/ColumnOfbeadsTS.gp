i=1
my_terminal = 'aqua'

result_file="ColumnOfbeadsTS-BinaryCZM.ref"
result_file="ColumnOfbeadsTS-BinaryCZM.dat"


my_terminal = 'wxt'
set term my_terminal i
plot result_file u 1:2 w l t 'position beads #1'
i=i+1
set term my_terminal i
plot result_file u 1:3 w l t 'velocity beads #1'
i=i+1
set term my_terminal i
plot result_file u 1:4 w l t 'position beads #2'
i=i+1
set term my_terminal i
plot result_file u 1:5 w l t 'velocity beads #2'
i=i+1
set term my_terminal i

set multiplot
set size 0.4,0.4		
set origin 0.1,0.1
plot result_file u 1:6 w l t 'lambda inter floor/bead'
set origin 0.1,0.5
plot result_file u 1:7 w l t 'y inter floor/bead'
set origin 0.5,0.1
plot result_file u 1:8 w l t 'beta inter floor/bead'

unset multiplot
i=i+1
set term my_terminal i
set multiplot
set size 0.4,0.4		
set origin 0.1,0.1	
plot result_file u 1:9 w l t 'lambda inter[0]'
set origin 0.1,0.5
plot result_file u 1:10 w l t 'y inter[0]'
set origin 0.5,0.1
plot result_file u 1:11 w l t 'beta inter[0]'
set origin 0.5,0.5
plot result_file u 9:10 w l t 'lambda vs y inter[0]'
unset multiplot


set origin 0.0,0.0
set size 1.0,1.0	


