module pimod
double precision dppi
end module

program main

!   This demonstrates numerical error when attempting to find the arc length
!   of an irregular function.
!   David H Bailey    2007-01-09

use pimod
implicit none
integer i, j, k, n
parameter (n = 10000000)
double precision fun, h, s1, t1, t2
external fun

t1 = -1.d0
dppi = acos (t1)
s1 = 0.d0
t1 = 0.d0
h = dppi / dble (n)

do i = 1, n
  t2 = fun (i * h)
  s1 = s1 + sqrt (h**2 + (t2 - t1)**2)
  t1 = t2
enddo

write (6, '(1p,d25.15)') s1
stop
end

function fun (x)
use pimod
implicit none
integer k, nn
parameter (nn = 10)
double precision d1
double precision fun, t1, x

t1 = x

do k = 1, nn
  d1 = 2.d0 ** k
  t1 = t1 + sin (d1 * x) / d1
enddo

fun = t1
return
end
