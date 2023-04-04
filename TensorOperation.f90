! ----------------------------------------------------------------------------------------
!---------------
! Copyright (C)  2023  Federico Pisanò (project leader and scientific supervisor)
!                      Zheng Li (code developer)
!                      Haoyuan Liu (code developer)
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
! USA.
Subroutine GetTrace(v, res)
!
!***********************************************************************
!
!     Function: To return the trace of a 6 vector
!
!***********************************************************************
!
Implicit none

double Precision, intent(in) :: v(6)
double precision, intent(out) :: res

double precision, parameter :: zero  = 0.0D0

res = zero
res = v(1)+v(2)+v(3)

Return
End

Subroutine GetDevPart(v, devP)
!
!***********************************************************************
!
!     Function: To return the deviatoric part of a vector
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v(6)
double precision, intent(out) :: devP(6)
double precision, parameter :: one3 = 1.0d0/3.0d0

double precision res

call GetTrace(v, res)
res = res * one3
!     call COPYRVEC(v, devP, 6)
devP(1) = v(1) - res
devP(2) = v(2) - res
devP(3) = v(3) - res
devP(4) = v(4)
devP(5) = v(5)
devP(6) = v(6)

Return
End

Subroutine SingleDot(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes v1.v2, v1 and v2 should be both in their "contravariant" form
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res(6)

double precision, parameter :: one2  = 0.5D0

res(1) = v1(1)*v2(1) + v1(4)*v2(4) + v1(6)*v2(6)
res(2) = v1(4)*v2(4) + v1(2)*v2(2) + v1(5)*v2(5)
res(3) = v1(6)*v2(6) + v1(5)*v2(5) + v1(3)*v2(3)
res(4) = one2*(v1(1)*v2(4) + v1(4)*v2(1) + v1(4)*v2(2) + v1(2)*v2(4) + v1(6)*v2(5) + v1(5)*v2(6))
res(5) = one2*(v1(4)*v2(6) + v1(6)*v2(4) + v1(2)*v2(5) + v1(5)*v2(2) + v1(5)*v2(3) + v1(3)*v2(5))
res(6) = one2*(v1(1)*v2(6) + v1(6)*v2(1) + v1(4)*v2(5) + v1(5)*v2(4) + v1(6)*v2(3) + v1(3)*v2(6))

Return
End

Subroutine SingleDotCube(v1, res)
!
!***********************************************************************
!
!     Function: computes v1.v2, v1 and v2 should be both in their "contravariant" form
!
!***********************************************************************
!
implicit none

double precision, intent(in) :: v1(6)
double precision, intent(out) :: res(6)

double precision, parameter :: one2  = 0.5D0
double precision, parameter :: two  = 2.0D0

res(1) = v1(1)*v1(1)*v1(1) + two*v1(1)*v1(4)*v1(4) + two*v1(1)*v1(6)*v1(6) + v1(2)*v1(4)*v1(4) + two*v1(5)*v1(4)*v1(6) + v1(3)*v1(6)*v1(6)
res(2) = v1(2)*v1(2)*v1(2) + two*v1(2)*v1(4)*v1(4) + two*v1(2)*v1(5)*v1(5) + v1(1)*v1(4)*v1(4) + two*v1(6)*v1(4)*v1(5) + v1(3)*v1(5)*v1(5)
res(3) = v1(3)*v1(3)*v1(3) + two*v1(3)*v1(5)*v1(5) + two*v1(3)*v1(6)*v1(6) + v1(2)*v1(5)*v1(5) + two*v1(4)*v1(5)*v1(6) + v1(1)*v1(6)*v1(6)
res(4) = v1(2)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(5)*(v1(1)*v1(6) + v1(3)*v1(6) + v1(4)*v1(5)) + v1(4)*(v1(1)*v1(1) + v1(4)*v1(4) + v1(6)*v1(6))
res(5) = v1(3)*(v1(2)*v1(5) + v1(3)*v1(5) + v1(4)*v1(6)) + v1(6)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(5)*(v1(2)*v1(2) + v1(4)*v1(4) + v1(5)*v1(5))
res(6) = v1(5)*(v1(1)*v1(4) + v1(2)*v1(4) + v1(5)*v1(6)) + v1(3)*(v1(1)*v1(6) + v1(3)*v1(6) + v1(4)*v1(5)) + v1(6)*(v1(1)*v1(1) + v1(4)*v1(4) + v1(6)*v1(6))

Return
End

Subroutine DoubleDot2_2_Contr(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, both "contravariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0
double precision, parameter :: two = 2.0D0

integer i

res = zero

do i = 1, 6
    if (i > 3) then
        res = res + two * v1(i) * v2(i)
    else
        res = res + v1(i) * v2(i)
    end if
end do

Return
End

Subroutine DoubleDot2_2_Cov(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, both "covariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one2 = 0.5D0

integer i

res = zero

do i = 1, 6
    if (i > 3) then
        res = res + one2 * v1(i) * v2(i)
    else
        res = res + v1(i) * v2(i)
    end if
end do

Return
End

Subroutine DoubleDot2_2_Mixed(v1, v2, res)
!
!***********************************************************************
!
!     Function: computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
double precision, intent(out) :: res

double precision, parameter :: zero = 0.0D0

integer i

res = zero

do i = 1, 6
    res = res + v1(i)*v2(i)
end do

Return
End

Subroutine GetNorm_Contr(v, res)
!
!***********************************************************************
!
!     Function: computes contravariant (stress-like) norm of input 6x1 tensor
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v(6)
double precision, intent(out) :: res

call DoubleDot2_2_Contr(v, v, res)
res = Dsqrt(res)

Return
End

Subroutine GetNorm_Cov(v, res)
!
!***********************************************************************
!
!     Function: computes covariant (strain-like) norm of input 6x1 tensor
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v(6)
double precision, intent(out) :: res

call DoubleDot2_2_Cov(v, v, res)
res = Dsqrt(res)

Return
End

subroutine Dyadic2_2(v1, v2, n, res)
!
!***********************************************************************
!
!     Function: computes dyadic product for two vector-storage arguments
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6), v2(6)
integer, intent(in) :: n
double precision, intent(out) :: res(6,6)

integer i, j

Do i=1,6
    Do j=1,6
        res(i,j) = v1(i)*v2(j)
    end do
end do

Return
End

subroutine ToContraviant(v1, v2)
!
!***********************************************************************
!
!     Function: Transform to the 'contravariant' form
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6)
double precision, intent(out) :: v2(6)

double precision, parameter :: one2 = 0.5D0

integer i

Do i = 1,6
    if (i > 3) then
        v2(i) = one2 * v1(i)
    else
        v2(i) = v1(i)
    end if
end do

Return
End

subroutine ToCovariant(v1, v2)
!
!***********************************************************************
!
!     Function: Transform to the 'covariant' form
!
!***********************************************************************
!
implicit none
double precision, intent(in) :: v1(6)
double precision, intent(out) :: v2(6)

double precision, parameter :: two = 2.0D0

integer i

Do i = 1,6
    if (i > 3) then
        v2(i) = two * v1(i)
    else
        v2(i) = v1(i)
    end if
end do

Return
End