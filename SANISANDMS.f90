! ----------------------------------------------------------------------------------------
!---------------
! 'SANISAND-MS': a user-defined soil model for PLAXIS
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
Subroutine SANISANDMS ( IDTask, iMod, IsUndr,   &
                        iStep, iTer, iEl, Int,  &
                        X, Y, Z, &
                        Time0, dTime,    &
                        Props, Sig0, Swp0, StVar0,   &
                        dEps, D, BulkW,  &
                        Sig, Swp, StVar, ipl,    &
                        nStat, NonSym, iStrsDep, iTimeDep, iTang,   &
                        iAbort )

! ----------------------------------------------------------------------------------------
! User-defined soil model: SANISAND-MS (iMod=1)
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : RETURN number of state variables
!                       5 : inquire matrix properties
!                           RETURN switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure at start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulk modulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : =1 for Non-Symmetric D-matrix
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iAbort   O   I    : =1 to force arrest of calculation
!   
! ----------------------------------------------------------------------------------------
! Expected content of Props(1..50)
! 
!  1 : G0           Shear modulus coefficient
!  2 : nu           Poisson's ratio
!  3 : Mc           critical state M value (in triaxial compression) 
!  4 : c            Me/Mc ratio between critical stress ratios in triaxial extention and compression
!  5 : lambda_c     slope of critical state line from void ratio-ln_p'-diagram
!  6 : e0           reference void ratio of the critical state line at pc = 0 kPa
!  7 : ksi          critical state line costant
!  8 : m            radious of the yield surface in the stress ratio space
!  9 : h0           hardening parameter
! 10 : ch           hardening parameter
! 11 : nb           bounding surface parameter
! 12 : A0           dilatancy parameter
! 13 : nd           dilatancy surface parameter
! 14 : mu0          ratcheting parameter
! 15 : zeta         memory surface shrinkage parameter
! 16 : beta         dilatancy memory parameter
! 17 : TolF         yield surface tolerance
! 18 : TolR         error tolerance for the RKF23 adaptive substepping scheme
! 19 : e_init       initial void ratio
! 20 : emax         maximum void ratio
! 21 : emin         minimum void ratio
! 22 : Mbc_max      Maximum value of the Bounding stress ratio (in triaxial compression)  
! 23 : El_Scheme    Toggle to choose the integration scheme of the elastic equations: 
!                     [  0  ] = Single step Forward Euler,
!                     [  1  ] = Automatic sub-step RKF23 for elastic trial only, 
!                     [  2  ] = Automatic sub-step RKF23 for elastic trial and path to the intersection stress
! 24 : hM_exc_Eq    Toggle to activate the new hM equation for tackling corner cases
!                     [  0  ] = inactive. 
!                     [  1  ] = active
! 25 : h_max        Toggle to customize the limit for h:
!                     [  0  ] = Default streategy limit separately the exponential term and h
!                     [value] = h is limited to the input value. It could 1.0d8
!
! Note: Toggles 24 and 25 will exclusively be in place during the testing phase
! ----------------------------------------------------------------------------------------

Use Parameters_position_in_params_iMod1
Use State_Var_position_in_StVar_iMod1
implicit none

integer, intent(in) :: IDTask, iMod, IsUndr, iStep, iTer, iEl, Int

double precision, intent(in) :: X, Y, Z, Time0, dTime, Props(*), Swp0, dEps(*) 

double precision, intent(inout) :: Sig0(*), StVar0(*)
double precision, intent(out) :: D(6,6), BulkW,Sig(*),Swp,StVar(*)
integer, intent(out) :: ipl, nStat, NonSym, iStrsDep, iTimeDep, iTang, iAbort
!
!---  Local variables
!
double precision stress(6), strain(6), strainInc(6), Alpha(6), AlphaM(6), alpha_in(6), sNextStress(6),    &
    sNextAlpha(6), sNextAlphaM(6), sNextalpha_in(6), sn(6), TolF, &
    dd(6), bb(6), bM(6), R(6), rr(6),res4(6),res5(6),dDevStrain(6), A, temp4, sNextDGamma,f,PM, b0, B, ZZ, &
    ddsdde(6,6), xNu, PAR_e_init, E_K,E_G, xNu_U, Void, sNextVoid, sNextMM_plus, &
    sNextMM_minus,sMM_plus,  sMM_minus, res, res3, res2, rBtheta, rDtheta, p, psi,  &
    VoidRatio, Fac, dVolStrain, dSwp, dEpsV, DDD, C, Cos3Theta, strainInc1(6), &
    PLNoTension, alfa_trace, emax, emin, PsPsmall, PAR_c, res24, dgdth, Mb, Md, PAR_Mc, PAR_nb, Mc, PAR_nd, &
    q_Norm_New, q_p, Mbc_max

double precision h, hM, fnNext

integer io, i, j, nStatV, itestnan, indicator

double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: two3 = 2.0D0/3.0D0
double precision, parameter :: tol = 1.0D-10
double precision, parameter :: ZERO2D = 1.0D-30


ipl=0

! indicator = 1, open warning output; indicator = 0, close warning output
! 200 presents norm of stress and small stress problem
indicator = 2

io=0
!     IF (iEl==-87 .And. Int==12) io=1 ! possibly write some debug info for a point
!     Number of state variables

nStatV = 31
1060 FORMAT(A, I5, A, I5, A, I5, A, I5)
1070 FORMAT(A, 6F15.10)
1080 FORMAT(A, 1F10.4)
1090 FORMAT(A, 6F35.23)
1110 FORMAT(18F20.15)
1120 FORMAT(42F20.15)

iAbort = 0

xNu         = Props(i_PAR_nu)
PAR_Mc      = Props(i_PAR_Mc)
PAR_c       = Props(i_PAR_c)
PAR_nb      = Props(i_PAR_nb)
PAR_nd      = Props(i_PAR_nd)
TolF        = Props(i_TolF)
PAR_e_init  = Props(i_PAR_e_init)
emax        = props(i_emax)
emin        = props(i_emin)
Mbc_max     = props(i_Mbc_max)

! change the signs (In PLAXIS, compression is negative)	
stress = -Sig0(1:6)
strain = -dEps(7:12)
strainInc = -dEps(1:6)

strainInc1 = strainInc

! ----------------------------------------------------------------------------------------
! get number of state parameters
! ----------------------------------------------------------------------------------------
IF (IDTask == 4) THEN
    nStat = nStatV
    RETURN
End IF

! ----------------------------------------------------------------------------------------
! get matrix attributes
! ----------------------------------------------------------------------------------------
IF (IDTask == 5) THEN ! matrix type
    NonSym   = 0  ! 1 for non-symmetric D-matrix
    iStrsDep = 1  ! 1 for stress dependent D-matrix
    iTang    = 0  ! 1 for tangent D-matrix
    iTimeDep = 0  ! 1 for time dependent D-matrix
End IF

! ----------------------------------------------------------------------------------------
! Initialize state variables
! ----------------------------------------------------------------------------------------
IF (IDTask == 1) THEN
    IF (StVar0(i_resetState) /= 777D0) THEN ! To void reinitialization problem in continuous phases
        CALL initialize(stress, Alpha, alpha_in, AlphaM, sMM_plus,  &
            sMM_minus, Void, Mb, Md, cos3theta, psi, q_Norm_New, p, Props)

        StVar0(i_Alpha_xx)      = Alpha(1)
        StVar0(i_Alpha_yy)      = Alpha(2)
        StVar0(i_Alpha_zz)      = Alpha(3)
        StVar0(i_Alpha_xy)      = Alpha(4)
        StVar0(i_Alpha_yz)      = Alpha(5)
        StVar0(i_Alpha_zx)      = Alpha(6)
        StVar0(i_AlphaM_xx)     = AlphaM(1)
        StVar0(i_AlphaM_yy)     = AlphaM(2)
        StVar0(i_AlphaM_zz)     = AlphaM(3)
        StVar0(i_AlphaM_xy)     = AlphaM(4)
        StVar0(i_AlphaM_yz)     = AlphaM(5)
        StVar0(i_AlphaM_zx)     = AlphaM(6)
        StVar0(i_sMM_plus)      = sMM_plus
        StVar0(i_sMM_minus)     = sMM_minus
        StVar0(i_alpha_in_xx)   = alpha_in(1)
        StVar0(i_alpha_in_yy)   = alpha_in(2)
        StVar0(i_alpha_in_zz)   = alpha_in(3)
        StVar0(i_alpha_in_xy)   = alpha_in(4)
        StVar0(i_alpha_in_yz)   = alpha_in(5)
        StVar0(i_alpha_in_zx)   = alpha_in(6)
        StVar0(i_Void)          = Void
        StVar0(i_resetState)    = 777D0
        StVar0(i_PLNoTension)   = 0.0D0 !deprecated use. Their name is removed from output. They will be repositioned in StVar
        StVar0(i_PsPsmall)      = 0.0D0 !deprecated use. Their name is removed from output. They will be repositioned in StVar
        StVar0(i_Mb)            = Mb
        StVar0(i_Md)            = Md
        StVar0(i_cos3theta)     = cos3theta
        StVar0(i_psi)           = psi
        StVar0(i_q_Norm_New)    = q_Norm_New
        StVar0(i_p)             = p
        StVar0(i_MM)            = sMM_plus
        
        StVar(1:31)             = StVar0(1:31)
    END IF
RETURN
End IF

! ----------------------------------------------------------------------------------------
! calculate stresses
! ----------------------------------------------------------------------------------------
IF (IDTask == 2) THEN
    
    Alpha(1)     = StVar0(i_Alpha_xx)
    Alpha(2)     = StVar0(i_Alpha_yy)
    Alpha(3)     = StVar0(i_Alpha_zz)
    Alpha(4)     = StVar0(i_Alpha_xy)
    Alpha(5)     = StVar0(i_Alpha_yz)
    Alpha(6)     = StVar0(i_Alpha_zx)
    AlphaM(1)    = StVar0(i_AlphaM_xx)
    AlphaM(2)    = StVar0(i_AlphaM_yy)
    AlphaM(3)    = StVar0(i_AlphaM_zz)
    AlphaM(4)    = StVar0(i_AlphaM_xy)
    AlphaM(5)    = StVar0(i_AlphaM_yz)
    AlphaM(6)    = StVar0(i_AlphaM_zx)
    sMM_plus     = StVar0(i_sMM_plus)! current size of the memory as it would be if it was affected only by the expantion due to the hardening
    sMM_minus    = StVar0(i_sMM_minus)! variation of the size of the memory due to the dilatancy
    alpha_in(1)  = StVar0(i_alpha_in_xx)
    alpha_in(2)  = StVar0(i_alpha_in_yy)
    alpha_in(3)  = StVar0(i_alpha_in_zz)
    alpha_in(4)  = StVar0(i_alpha_in_xy)
    alpha_in(5)  = StVar0(i_alpha_in_yz)
    alpha_in(6)  = StVar0(i_alpha_in_zx)
    Void            = StVar0(i_Void)
    PsPsmall        = StVar0(i_PsPsmall)
    Mb              = StVar0(i_Mb)

! 21-12-2021 check input q/p reach Mb or not------------------
    CALL GetTrace(stress, p)
    p = one3 * p
    CALL GetDevPart(stress, rr)
    CALL GetNorm_Contr(rr, q_Norm_New)
    q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
    q_p = q_Norm_New / p
!-------------------------------------------------------------
    !If (iStep == 1 .and. iEl == 1 .and. indicator == 500) then
    !    !WRITE(100, *) 'before integration', iStep, iTer!, iEl, Int
    !    WRITE(100, '(6f25.15)') strainInc
    !    WRITE(100, '(6f25.15)') stress
    !    WRITE(100, '(6f25.15)') Alpha
    !    call flush(100)
    !END IF

    CALL Integration(stress, Void, Alpha, AlphaM, sMM_plus, &
       sMM_minus, alpha_in, sNextStress, strainInc, sNextVoid,   &
       sNextAlpha,sNextAlphaM, sNextMM_plus, sNextMM_minus,    &
       sNextalpha_in, Props, ipl, ddsdde, IDTask, iEl, Int, iStep, &
       indicator, PLNoTension, PsPsmall, iAbort, iTer)

    IF (iAbort == 1) THEN
        RETURN
    END IF
    
    !If (iStep == 1 .and. iEl == 1 .AND. Int = 4.and. indicator == 200) then
        !WRITE(100, *) 'integration succeed', iStep, iTer!, iEl, Int
        !WRITE(100, '(6f25.15)') strainInc
        !WRITE(100, '(6f25.15)') sNextStress
        !WRITE(100, '(6f25.15)') sNextAlpha
        !call flush(100)
    !END IF
    CALL GetNormalToYield2(sNextStress, sNextAlpha, iStep, sn)
    CALL GetTrace(sNextStress, p)
    p = one3 * p
    
    ! Calculate the current stress ratio
    CALL GetDevPart(sNextStress, rr)
    CALL GetNorm_Contr(rr, q_Norm_New)
    q_Norm_New = dsqrt(3.0D0/2.0D0) * q_Norm_New
    
    ! Calculate the current Bounding, Critical state and Dilatancy stress ratios
    CALL GetLodeAngle(sn, cos3Theta)
    CALL g(cos3Theta, PAR_c, res24, dgdth)
    CALL GetPSI(sNextVoid, p, Props, psi)
    Mb = PAR_Mc * EXP(-PAR_nb * psi)
    if (Mb > Mbc_max) then
        Mb = Mbc_max
    endif
    Mb = res24 * Mb
    Mc = res24 * PAR_Mc
    Md = res24 * PAR_Mc * EXP(PAR_nd * psi)

    Sig(1:6)     = -sNextStress
    StVar(i_Alpha_xx)   = sNextAlpha(1)
    StVar(i_Alpha_yy)   = sNextAlpha(2)
    StVar(i_Alpha_zz)   = sNextAlpha(3)
    StVar(i_Alpha_xy)   = sNextAlpha(4)
    StVar(i_Alpha_yz)   = sNextAlpha(5)
    StVar(i_Alpha_zx)   = sNextAlpha(6)
    StVar(i_AlphaM_xx)  = sNextAlphaM(1)
    StVar(i_AlphaM_yy)  = sNextAlphaM(2)
    StVar(i_AlphaM_zz)  = sNextAlphaM(3)
    StVar(i_AlphaM_xy)  = sNextAlphaM(4)
    StVar(i_AlphaM_yz)  = sNextAlphaM(5)
    StVar(i_AlphaM_zx)  = sNextAlphaM(6)
    StVar(i_sMM_plus)    = sNextMM_plus
    StVar(i_sMM_minus)   = sNextMM_minus
    StVar(i_alpha_in_xx) = sNextalpha_in(1)
    StVar(i_alpha_in_yy) = sNextalpha_in(2)
    StVar(i_alpha_in_zz) = sNextalpha_in(3)
    StVar(i_alpha_in_xy) = sNextalpha_in(4)
    StVar(i_alpha_in_yz) = sNextalpha_in(5)
    StVar(i_alpha_in_zx) = sNextalpha_in(6)
    StVar(i_Void)           = sNextVoid
    StVar(i_resetState)     = 777D0
    StVar(i_PLNoTension)    = PLNoTension
    StVar(i_PsPsmall)       = PsPsmall
    StVar(i_Mb)             = Mb
    StVar(i_Md)             = Md
    StVar(i_cos3theta)      = cos3Theta         !cos(3*Lode angle(sn)) 
    StVar(i_psi)            = psi
    StVar(i_q_Norm_New)     = q_Norm_New
    StVar(i_p)              = p
    StVar(i_MM)            = sNextMM_plus + sNextMM_minus
!------------------------------------------------------    
    CALL GetF(sNextStress, sNextAlpha, Props, fnNext)
    alfa_trace = sNextAlpha(1) + sNextAlpha(2) + sNextAlpha(3)
!------------------------------------------------------   
    IF (IsUndr == 1) THEN
        xNu_U = 0.495d0 ! Undrained Poissons' ratio
        Fac = (1.d0+xNu_U)/(1.0d0-2.0d0*xNu_U)-(1.0d0+xNu)/(1.0d0-2.0d0*xNu)
        ! 24-11-2021
        CALL GetElasticModuli1(sNextStress,sNextVoid,Props, E_K, E_G)
        Fac = 2D0*E_G/3D0*Fac
        BulkW = Fac
        dEpsV = dEps(1) + dEps(2) + dEps(3)
        dSwp  = BulkW * dEpsV
        Swp   = Swp0 + dSwp
    ELSE
        Swp = Swp0
    END IF
!
!     Nan check: Sig, sNextAlpha, sNextAlphaM, sNextalpha_in, sNextMM_plus, sNextMM_minus, sNextVoid
    Do i=1,6
        If ( IsNaN( Sig( i))) Then
            Write(100,*), 'Sig(i)', Sig(i)
            Write(100,*), 'iEl', iEl
            call flush(100)
            iabort = 1
            return
        End If
    End Do

    Do i=1,6
        If (IsNaN( sNextAlpha(i))) Then
                Write(100,*), 'sNextAlpha',sNextAlpha
                Write(100,*), 'iEl',iEl
                call flush(100)
                iabort = 1
                return
        End If
    End Do

! 06-29-2021
    Do i = 1, 6
        If ( IsNaN( sNextAlphaM( i))) Then
        Write(100,*), 'sNextAlphaM', sNextAlphaM
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
        endif
    End Do

    Do i = 1,6
        If ( IsNaN( sNextalpha_in(i))) Then
            Write(100,*), 'sNextalpha_in',sNextalpha_in
            Write(100,*), 'iEl',iEl
            call flush(100)
            iabort = 1
            return
        End If
    End Do

    If ( IsNaN( sNextMM_plus)) Then
        Write(100,*), 'sNextMM_plus',sNextMM_plus
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
    End If

    If ( IsNaN( sNextMM_minus)) Then
        Write(100,*), 'sNextMM_minus',sNextMM_minus
        Write(100,*), 'iEl',iEl
        call flush(100)
        iabort = 1
        return
    End If

    If ( IsNaN( sNextVoid)) Then 
        iabort = 1
        return
    End If

! NaN output
! 06-29-2021  
!    IF ((itestnan == 1) .AND. indicator == 2) THEN
!        WRITE(100,*) 'IDTask = ', IDTask
!        WRITE(100, 1070) 'Inistress =', stress(1:6)
!        WRITE(100, 1070) 'IniAlpha =', Alpha
!        WRITE(100, 1070) AlphaM
!        WRITE(100, 1080) sMM_plus
!        WRITE(100, 1080) sMM_minus
!        WRITE(100, 1070) 'Inialpha_in =', alpha_in
!        WRITE(100, 1080) 'IniVoid =', Void
!        WRITE(100, 1090) 'strainInc =', strainInc
!        WRITE(100,1060)'iStep =',iStep,' iTer =',iTer,' iEl = ', iEl, ' Int = ', Int
!        WRITE(100, 1070) 'sNextStress = ', sNextStress
!        WRITE(100, 1070) 'sNextAlpha = ', sNextAlpha
!        WRITE(100, 1070) 'sNextAlphaM = ', sNextAlphaM
!        WRITE(100, 1080) 'sNextMM_plus = ', sNextMM_plus
!        WRITE(100, 1080) 'sNextMM_minus = ', sNextMM_minus
!        WRITE(100, 1070) 'sNextalpha_in = ', sNextalpha_in
!        WRITE(100, 1080) 'sNextVoid = ', sNextVoid
!        call flush(100) 
!        iAbort = 1
!         return
!    END IF

RETURN
End IF

! ----------------------------------------------------------------------------------------
! calculate material stiffness matrix (D-matrix)
! ----------------------------------------------------------------------------------------
IF (IDTask == 6 .or. IDTask == 3) THEN ! matrix type
    !write(100,*), IDTask

    Alpha(1)     = StVar0(i_Alpha_xx)
    Alpha(2)     = StVar0(i_Alpha_yy)
    Alpha(3)     = StVar0(i_Alpha_zz)
    Alpha(4)     = StVar0(i_Alpha_xy)
    Alpha(5)     = StVar0(i_Alpha_yz)
    Alpha(6)     = StVar0(i_Alpha_zx)
    AlphaM(1)    = StVar0(i_AlphaM_xx)
    AlphaM(2)    = StVar0(i_AlphaM_yy)
    AlphaM(3)    = StVar0(i_AlphaM_zz)
    AlphaM(4)    = StVar0(i_AlphaM_xy)
    AlphaM(5)    = StVar0(i_AlphaM_yz)
    AlphaM(6)    = StVar0(i_AlphaM_zx)
    sMM_plus     = StVar0(i_sMM_plus)
    sMM_minus    = StVar0(i_sMM_minus)
    alpha_in(1)  = StVar0(i_alpha_in_xx)
    alpha_in(2)  = StVar0(i_alpha_in_yy)
    alpha_in(3)  = StVar0(i_alpha_in_zz)
    alpha_in(4)  = StVar0(i_alpha_in_xy)
    alpha_in(5)  = StVar0(i_alpha_in_yz)
    alpha_in(6)  = StVar0(i_alpha_in_zx)
!------------------------------------------------------------
    Void            = StVar0(i_Void)    
    stress=-Sig0(1:6)
    
    CALL GetElasticModuli1(stress, Void, Props, E_K, E_G)
    CALL GetStiffness(E_K, E_G, D)

! calculate bulk modulus
!        BulkW = 0  
!    IF (IsUndr == 1) THEN
!        xNu_U = 0.495d0
!        Fac = (1+xNu_U)/(1-2*xNu_U) - (1+xNu)/(1-2*xNu)
!        CALL GetTrace(strainInc, res)
!        VoidRatio = Void - (1.0d0 + Void) * res
!        CALL GetElasticModuli(stress, VoidRatio, Props, E_K, E_G)
!        Fac = 2D0*E_G/3D0*Fac
!        BulkW = Fac
!    END IF
    
    !WRITE(100, *) 'IDTask', IDTask, iStep, iTer
    !WRITE(100,'(6F14.6)') D
    !call flush(100)

End IF

Return

End


