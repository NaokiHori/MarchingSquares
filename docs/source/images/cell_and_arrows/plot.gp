reset

ln = 2.0
mrgn = 0.05
ln_ = ln-2.*mrgn

xmin =    0.-mrgn
xmax = 4.*ln+mrgn
ymin =    0.-mrgn
ymax = 4.*ln+mrgn

lx = xmax-xmin
ly = ymax-ymin

set terminal epslatex standalone color size lx,ly font ',17.28'
set output 'result.tex'

unset border

set lmargin 0
set rmargin 0
set bmargin 0
set tmargin 0

unset xlabel
unset ylabel

set xrange [xmin:xmax]
set yrange [ymin:ymax]

unset xtics
unset ytics

do for [case = 0 : 15 : 1] {

  col = case % 4
  row = case / 4

  ox = col * ln + mrgn
  oy = row * ln + mrgn

  v0 = case % 2
  v1 = (case / 2) % 2
  v2 = ((case / 2) / 2) % 2
  v3 = (((case / 2) / 2) / 2) % 2
  set label sprintf('$v: 0b%d%d%d%d$', v3, v2, v1, v0) center at (col + 0.5) * ln, (row + 0.60) * ln front
  if(v0 == v1){
    e0 = 0
  }else{
    e0 = 1
  }
  if(v1 == v2){
    e1 = 0
  }else{
    e1 = 1
  }
  if(v2 == v3){
    e2 = 0
  }else{
    e2 = 1
  }
  if(v3 == v0){
    e3 = 0
  }else{
    e3 = 1
  }
  set label sprintf('$e: 0b%d%d%d%d$', e3, e2, e1, e0) center at (col + 0.5) * ln, (row + 0.40) * ln front

  set object rectangle from ox, oy to ox+ln_, oy+ln_ fc rgb '#FFFFFF' fs empty border rgb '#000000' lw 10

  set style line 1 lc rgb '#FF0000' lw 10
  set style line 2 lc rgb '#0000FF' lw 10
  set style arrow 1 head size 0.5, 15 filled front ls 1
  set style arrow 2 head size 0.5, 15 filled front ls 2
  e0x = ox+0.5*ln_
  e0y = oy+0.0*ln_
  e1x = ox+1.0*ln_
  e1y = oy+0.5*ln_
  e2x = ox+0.5*ln_
  e2y = oy+1.0*ln_
  e3x = ox+0.0*ln_
  e3y = oy+0.5*ln_
  if(case == 1){
    set arrow from e3x, e3y to e0x, e0y as 1 back
  }
  if(case == 2){
    set arrow from e0x, e0y to e1x, e1y as 1 back
  }
  if(case == 3){
    set arrow from e3x, e3y to e1x, e1y as 1 back
  }
  if(case == 4){
    set arrow from e1x, e1y to e2x, e2y as 1 back
  }
  if(case == 5){
    set arrow from e3x, e3y to e0x, e0y as 1 back
    set arrow from e1x, e1y to e2x, e2y as 1 back
    set arrow from e3x, e3y to e2x, e2y as 2 back
    set arrow from e1x, e1y to e0x, e0y as 2 back
  }
  if(case == 6){
    set arrow from e0x, e0y to e2x, e2y as 1 back
  }
  if(case == 7){
    set arrow from e3x, e3y to e2x, e2y as 1 back
  }
  if(case == 8){
    set arrow from e2x, e2y to e3x, e3y as 1 back
  }
  if(case == 9){
    set arrow from e2x, e2y to e0x, e0y as 1 back
  }
  if(case == 10){
    set arrow from e2x, e2y to e3x, e3y as 1 back
    set arrow from e0x, e0y to e1x, e1y as 1 back
    set arrow from e2x, e2y to e1x, e1y as 2 back
    set arrow from e0x, e0y to e3x, e3y as 2 back
  }
  if(case == 11){
    set arrow from e2x, e2y to e1x, e1y as 1 back
  }
  if(case == 12){
    set arrow from e1x, e1y to e3x, e3y as 1 back
  }
  if(case == 13){
    set arrow from e1x, e1y to e0x, e0y as 1 back
  }
  if(case == 14){
    set arrow from e0x, e0y to e3x, e3y as 1 back
  }

}

plot \
  NaN notitle

