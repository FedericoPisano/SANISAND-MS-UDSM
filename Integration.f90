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
SUBROUTINE Integration(CurStress, CurVoid, CurAlpha, CurAlphaM, &
    CurMM_plus, CurMM_minus, Curalpha_in, sNextStress,  &
    strainInc, sNextVoid, sNextAlpha, sNextAlphaM, sNextMM_plus,   &
    sNextMM_minus, sNextalpha_in, Params, ipl, aCep, IDTask,    &
    iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
!
! Purpose: Update the stress tensor and the state variables
!
! Arguments:
!                  I/O   Type
!  IDTask           I    I()   :
!  iEl              I    I()   :
!  Int              I    I()   :
!  iStep            I    I()   :
!  iTer             I    I()  :
!  indicator        I    I()  :
!  strainInc        I    R() 
!  Params           I    R()  :
!  CurStress        I    R()  :
!  CurVoid          I    R()  :
!  CurAlpha         I    R()  :
!  CurAlphaM        I    R()  : 
!  CurMM_plus       I    R()  : the radius of the memory as affected only by the hardening
!  CurMM_minus      I    R()  : variation of the radius of the memory due to the Dilatancy
!  Curalpha_in      I    R()  :
!  sNextStress      O    R()  :
!  sNextVoid        O    R()  :
!  sNextAlpha       O    R()  :
!  sNextAlphaM      O    R()  :
!  sNextMM_plus     O    R()  : 
!  sNextMM_minus    O    R()  :
!  sNextalpha_in    O    R()  :
!  aCep             O    R()  :
!  PLNoTension      O    R    :
!  iAbort           O    i    :
!  ipl              IO   I    :
!  PsPsmall         IO   R    :
!
use extrainfo, only: onekPa
Use Parameters_position_in_params_iMod1

implicit none
double precision, intent(in) :: CurStress(6), CurAlpha(6), CurAlphaM(6),CurMM_plus, CurMM_minus, Curalpha_in(6),   &
    strainInc(6), Params(25)
integer, intent (in) :: IDTask, iEl, Int, iStep, indicator, iTer
integer, intent (out) :: iAbort
double precision, intent(in) :: CurVoid
double precision, intent(out) :: sNextStress(6), sNextAlpha(6), sNextAlphaM(6),sNextMM_plus, sNextMM_minus, sNextalpha_in(6),   &
    aCep(6,6), PLNoTension, sNextVoid
integer, intent(inout) :: ipl
double precision, intent(inout) :: PsPsmall

double precision dSigma(6), aC(6,6), sn(6), a0, a1, TolF, p, El_Scheme,   &
    sNextDGamma, sNextVoidRatio, res2, res1, pn, f, fn, E_K,E_G,    &
    elasticRatio, strainIncNew(6), InterVoid, res3(6),    &
    xN1(3),xN2(3),xN3(3),sNextStress1,sNextStress2,sNextStress3,PNextStress,QNextStress,sDevNextStress(6),pNext, res4(6)

double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0D0
double precision, parameter :: one3 = 1.0D0/3.0D0
double precision, parameter :: PAR_P_atm = 100.0D0
double precision, parameter :: small = 1.0D-10
double precision, parameter :: tol = 1.0D-6
double precision PAR_Pmin
!
logical, external :: isElIntegrationRKF23, isElIntegrationOneStepEuler,   &
    isElIntegrationRKF23alsoForIntersection, isElIntegrationRKF23NotForIntersection

iAbort = 0

PAR_Pmin = onekPa/100d0
PLNoTension = zero

a0 = zero
a1 = one

El_Scheme = Params(i_El_Scheme)
TolF = Params(i_TolF)
InterVoid = CurVoid

CALL GetElasticModuli1(CurStress, CurVoid, Params, E_K, E_G)
CALL GetStiffness(E_K, E_G, aC)

! Calculate the elastic trial stress with the integration scheme specified in input via the "El_Scheme" parameter 
if (isElIntegrationRKF23(Params)) then
    call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainInc, Params, &
        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
endif

if (iabort == 1 .or. isElIntegrationOneStepEuler(Params)) then    
	dSigma = MATMUL(aC, strainInc)
	sNextStress = CurStress + dSigma
else
    dSigma = sNextStress - CurStress
endif
iAbort = 0

pNext = (sNextStress(1) + sNextStress(2) + sNextStress(3))*one3
call PrnSig(1,sNextStress,xN1,xN2,xN3,sNextStress1,sNextStress2,sNextStress3,PNextStress,QNextStress)
!if (min(sNextStress1,sNextStress2,sNextStress3)<small) then
!if (pNext < tol) then
!    !sDevNextStress = PAR_Pmin * CurAlpha
!    !snextstress(1)  = sdevnextstress(1) + par_pmin
!    !snextstress(2)  = sdevnextstress(2) + par_pmin
!    !snextstress(3)  = sdevnextstress(3) + par_pmin
!    !snextstress(4)  = sdevnextstress(4)
!    !snextstress(5)  = sdevnextstress(5)
!    !snextstress(6)  = sdevnextstress(6)
!    
!    !p = PAR_Pmin
!    !CALL GetF(sNextStress, CurAlpha, Params, f)
!end if

CALL GetF(sNextStress, CurAlpha, Params, f)
CALL GetTrace(sNextStress, p)
p = one3 * p

CALL GetF(CurStress, CurAlpha, Params, fn) ! just for a check
CALL GetNorm_Contr(strainInc, res1)

!Check the trial stress
IF (f <= TolF) THEN
    !Pure elastic loading/unloading (trial inside the yield)'
    sNextAlpha    = CurAlpha
    sNextAlphaM   = CurAlphaM
    sNextMM_plus  = CurMM_plus
    sNextMM_minus = CurMM_minus
    sNextalpha_in = Curalpha_in
    sNextDGamma   = zero
    aCep = aC
    CALL GetTrace(strainInc, res2)
    sNextVoid = InterVoid - (1.0d0 + InterVoid) * res2
    ipl = 0
ELSE
    !Chenck the initial stress
    ipl = 1
    CALL GetF(CurStress, CurAlpha, Params, fn)
    CALL GetTrace(CurStress, pn)
    pn = one3 * pn
!   !Check if the initial stress is illegal
    IF (fn > TolF) THEN
        WRITE(100, *) '    trial stress is outside and initial outside the yield surface!', istep, fn
        WRITE(100, *) 'El, Int, iStep, iTer'
        WRITE(100, *) iEl, Int, iStep, iTer
        call flush(100)
        iAbort = 1
        RETURN
        CALL RungeKutta23(CurStress, CurVoid, CurAlpha, CurAlphaM,  &
            CurMM_plus, CurMM_minus, sNextalpha_in,     &
            strainInc, Params, sNextStress, sNextAlpha, &
            sNextAlphaM, sNextMM_plus, sNextMM_minus, &
            sNextVoidRatio, aCep, IDTask, iEl, Int,  &
            iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
            
            sNextVoid = sNextVoidRatio
        
        IF (iAbort == 1) THEN
            RETURN
        END IF
!   !Check if the initial stress is iniside the yield surface
    ELSEIF (fn < -TolF) THEN
!       !Calculate the intersection with the yield
        CALL IntersectionFactor(CurStress, CurVoid, strainInc,  &
                    CurAlpha, a0, a1, Params, elasticRatio, iAbort)
        IF (iAbort == 1) THEN
                RETURN
        END IF
        res3 = elasticRatio * strainInc
        !dSigma = MATMUL(aC, res3)

		! Calculate the stress at the intersection point with the yield
        if(isElIntegrationRKF23alsoForIntersection(Params))then
	        call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, res3, Params, &
	        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
	        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
		endif
        if (iabort == 1 .or. isElIntegrationRKF23NotForIntersection(Params)) then
            dSigma = MATMUL(aC, res3)
            sNextStress = CurStress + dSigma
            iAbort = 0
  
        else
            dSigma = sNextStress - CurStress
        endif
        
        CALL GetTrace(res3, res2)
        InterVoid = CurVoid - (one + CurVoid) * res2
        strainIncNew = (one - elasticRatio) * strainInc
        res3 = CurStress + dSigma
        
        !Compute the normal to the yield in the intersection point
        CALL GetNormalToYield(res3, CurAlpha, iStep, sn, iAbort)
        IF (iAbort == 1) THEN
            RETURN
        END IF
        !Check for the updating of alfa_initial---------------------------
        res4 = CurAlpha - Curalpha_in
        CALL DoubleDot2_2_Contr(res4, sn, res1)
        IF (res1 < -tol) THEN
            sNextalpha_in = CurAlpha
        ELSE
            sNextalpha_in = Curalpha_in
        END IF
        !-------------------------------------------------------------------
        CALL RungeKutta23(res3, InterVoid,    &
                    CurAlpha, CurAlphaM, CurMM_plus, &
                    CurMM_minus, sNextalpha_in, strainIncNew,  &
                    Params, sNextStress, sNextAlpha,  &
                    sNextAlphaM, sNextMM_plus, sNextMM_minus, &
                    sNextVoidRatio, aCep, IDTask, iEl, Int, &
                    iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
            
        sNextVoid = sNextVoidRatio
        IF (iAbort == 1) THEN
            RETURN
        END IF
    ELSE
        !Check if unloading occurs
        CALL GetNormalToYield(CurStress, CurAlpha, iStep, sn, iAbort)
        IF (iAbort == 1) THEN
            WRITE(100, *) 'failed to compute sn before chcking for unloading'
            call flush(100)
            RETURN
        END IF
        CALL DoubleDot2_2_Contr(sn, dSigma, res1)
        IF (res1 > -DSQRT(TolF)) THEN
            
            !No unloading (pure plastic loading): Check for the updateing of alfa_initial----------------
            res4 = CurAlpha - Curalpha_in
            CALL DoubleDot2_2_Contr(res4, sn, res1)
            IF (res1 < -tol) THEN
                sNextalpha_in = CurAlpha
            ELSE
                sNextalpha_in = Curalpha_in
            END IF
            !--------------------------------------------------------------------------------------------
            CALL RungeKutta23(CurStress, CurVoid, CurAlpha, &
                CurAlphaM,CurMM_plus, CurMM_minus, sNextalpha_in, &
                strainInc, Params, sNextStress, sNextAlpha,    &
                sNextAlphaM, sNextMM_plus, sNextMM_minus,  &
                sNextVoidRatio, aCep, IDTask, iEl, Int,  &
                iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)               
            sNextVoid = sNextVoidRatio
            IF (iAbort == 1) THEN
                RETURN
            END IF
        ELSE

            !Elastic unloading followed by plastic unloading: search the intersection with the yield surface
            CALL IntersectionFactor_Unloading(CurStress,CurVoid,    &
                            strainInc,CurAlpha, Params, elasticRatio, iAbort)
            IF (iAbort == 1) THEN
                RETURN
            END IF
            strainIncNew = elasticRatio * strainInc
            !dSigma = matmul(aC, strainIncNew)
            
            ! Calculate the stress at the intersection point with the yield
			if(isElIntegrationRKF23alsoForIntersection(Params))then
		        call EL_RungeKutta23(CurStress, InterVoid, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, strainIncNew, Params, &
		        sNextStress, sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, sNextVoid, aCep, &
		        IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
			endif            
            
	        if (iabort == 1 .or. isElIntegrationRKF23NotForIntersection(Params)) then
	            dSigma = MATMUL(aC, strainIncNew)
	            sNextStress = CurStress + dSigma
	            iAbort = 0
	        else
	            dSigma = sNextStress - CurStress
	        endif            
                       
            CALL GetTrace(strainIncNew, res2)
                InterVoid = CurVoid - (one + CurVoid) * res2
                strainIncNew = (one - elasticRatio) * strainInc
                res3 = CurStress + dSigma
            !Check for the updateing of alfa_initial
            CALL GetNormalToYield(res3, CurAlpha, iStep, sn, iAbort)
            IF (iAbort == 1) THEN
                RETURN
            END IF
            !Check for the updating of alfa_initial------------------------------------------------------
            res4 = CurAlpha - Curalpha_in
            CALL DoubleDot2_2_Contr(res4, sn, res1)
            IF (res1 < -tol) THEN
                sNextalpha_in = CurAlpha
            ELSE
                sNextalpha_in = Curalpha_in
            END IF
            !--------------------------------------------------------------------------------------------
            CALL RungeKutta23(res3, InterVoid,    &
                    CurAlpha, CurAlphaM, CurMM_plus,    &
                    CurMM_minus, sNextalpha_in, strainIncNew, &
                    Params, sNextStress, &
                    sNextAlpha, sNextAlphaM, sNextMM_plus, sNextMM_minus, &
                    sNextVoidRatio, aCep,    &
                    IDTask, iEl, Int, iStep, indicator, PLNoTension, PsPsmall, iAbort, iTer)
                
                    sNextVoid = sNextVoidRatio
            IF (iAbort == 1) THEN
                 RETURN
            END IF
        END IF
    END IF
END IF

RETURN
END
