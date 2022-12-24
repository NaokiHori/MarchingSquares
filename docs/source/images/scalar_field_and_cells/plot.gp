reset

mrgn = 0.125

lx = 2.
ly = 1.

nx = 8
ny = 4

dx_ = lx / nx
dy_ = ly / ny

set terminal epslatex standalone color size lx+2.*mrgn,ly+2.*mrgn font ',12'
set output 'result.tex'

unset border

set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

unset xlabel
unset ylabel

set xrange [-mrgn:lx+mrgn]
set yrange [-mrgn:ly+mrgn]

unset xtics
unset ytics

array xf[nx+1]
array yf[ny+1]
array xc[nx]
array yc[ny]

do for [i = 1 : nx+1 : 1] {
  xf[i] = 1.*(i-1)*dx_
}
do for [j = 1 : ny+1 : 1] {
  yf[j] = 1.*(j-1)*dy_
}
do for [i = 1 : nx : 1] {
  xc[i] = 0.5*(xf[i]+xf[i+1])
}
do for [j = 1 : ny : 1] {
  yc[j] = 0.5*(yf[j]+yf[j+1])
}

# xf and yf grids
set style line 1 lc rgb '#000000' lw 3
set style arrow 1 nohead ls 1
do for [i = 1 : nx+1 : 1]{
  set arrow from xf[i], yf[1] to xf[i], yf[ny+1] as 1
}
do for [j = 1 : ny+1 : 1]{
  set arrow from xf[1], yf[j] to xf[nx+1], yf[j] as 1
}
# hatch
ddx_ = lx / nx / 4
do for [i = 1 : 4*nx : 1] {
  x_ = 1.*i*ddx_
  ln = 0.0625
  set arrow from x_   , yf[   1] to x_-ln, yf[   1]-ln as 1
  set arrow from x_-ln, yf[ny+1] to x_   , yf[ny+1]+ln as 1
}

# cell center
min_(a, b) = a < b ? a : b
do for [i = 1 : nx : 1]{
  set object circle at xc[i], yf[   1] size first 0.0625*min_(dx_, dy_) fs solid fc rgb '#000000' lc rgb '#000000'
  set object circle at xc[i], yf[ny+1] size first 0.0625*min_(dx_, dy_) fs solid fc rgb '#000000' lc rgb '#000000'
  do for [j = 1 : ny : 1]{
    set object circle at xc[i], yc[j] size first 0.0625*min_(dx_, dy_) fs solid fc rgb '#000000' lc rgb '#000000'
  }
}

# cell in the current context
set style line 2 lc rgb '#000000' lw 1 dt 3
set style arrow 2 nohead ls 2
do for [i = 1 : nx : 1]{
  set arrow from xc[i], yf[1] to xc[i], yf[ny+1] as 2
}
do for [j = 1 : ny : 1]{
  set arrow from xc[1], yc[j] to xc[nx], yc[j  ] as 2
}

plot \
  NaN notitle

