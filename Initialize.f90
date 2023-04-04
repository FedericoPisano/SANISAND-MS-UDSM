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
!------------------------------------------------------------------------------------------
SUBROUTINE initialize(stress, alpha, alpha_in, alphaM, sMM_plus,    &
           sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p, Params)
!
! Purpose: Assigne the initial values to the local variables 
!           representative of the StVar0
! Arguments:
!            I/O  Type
!  stress     I   R()   : 
!  Params     I   R()  :
!  p          O   R    : 
!  q_Norm_New O   R    :
!  alpha      O   R()  :
!  alpha_in   O   R()  :
!  alphaM     O   R()  :
!  sMM_plus   O   R    : current radius resulting from + first  term in the isotropic hardening  [Memory rotation]
!  sMM_minus  O   R    : second term in the isotropic hardening  [volumetric plastic strains]
!  Void       O   R    :
!  psi        O   R    :
!  Mb         O   R    :
!  Md         O   R    :
!  cos3theta  O   R    :

Use Parameters_position_in_params_iMod1    
implicit none

double precision, intent(in) :: stress(6), Params(25)
double precision, intent(out) :: alpha(6), alpha_in(6), alphaM(6),  &
    sMM_plus, sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p

double precision rr(6), PAR_m, PAR_nb, PAR_nd, rrNorm, small, snrr(6), PAR_Mc, Mbc_max
    

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: zero = 0.0D0

PAR_Mc = Params(i_PAR_Mc)
PAR_m  = Params(i_PAR_m)
PAR_nb = Params(i_PAR_nb)
PAR_nd = Params(i_PAR_nd)

small = 1.0D-7

CALL GetTrace(stress, p)
p = one3 * p
CALL GetDevPart(stress, rr)
CALL GetNorm_Contr(rr, q_Norm_New)
q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
rr = rr / p

alpha = rr
alpha_in = rr
alphaM = rr
sMM_plus = PAR_m 
sMM_minus = zero

Void = Params(i_PAR_e_init)

CALL GetPSI(Void, p, Params, psi)
Mbc_max = Params(i_Mbc_max)
Mb = PAR_Mc * EXP(-PAR_nb * psi)! this is for tx compression not the current 
if (Mb > Mbc_max) then
    Mb = Mbc_max
end if

Md = PAR_Mc * EXP(PAR_nd * psi)

CALL GetNorm_Contr(rr, rrNorm)
IF (rrNorm < small) THEN
    rrNorm = small
END IF
snrr = rr/rrNorm
CALL GetLodeAngle(snrr, cos3theta)

RETURN
END