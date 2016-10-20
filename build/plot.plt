Jeff = 1.1064
Tc = 2 / log(1.0+sqrt(2.0))
TcK = Tc * Jeff
set parametric
set xrange [0:3]
set trange [0:1]
plot "local.out"  u 1:2 w lp, "local.out"  u 1:(10*$3/$1) w lp,\
     "localK.out" u 1:2 w lp, "localK.out" u 1:(10*$3/$1) w lp,\
     "wolff.out"  u 1:2 w lp, "wolff.out"  u 1:(10*$3/$1) w lp,\
     "wolffK.out" u 1:2 w lp, "wolffK.out" u 1:(10*$3/$1) w lp,\
     Tc, t lw 4 lc 0,\
     TcK, t lw 4 lc 0


