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
!----------------------------------------------------------------------------------------
! Subroutines in this file:
!
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iMod , ModelName )
!  Subroutine GetParamCount( iMod , nParam )
!  Subroutine GetParamName ( iMod , iParam, ParamName )
!  Subroutine GetParamUnit ( iMod , iParam, Units )
!  Subroutine GetStateVarCount( iMod , nVar )
!  Subroutine GetStateVarName ( iMod , iVar, Name )
!  Subroutine GetStateVarUnit ( iMod , iVar, Unit )
!
! Local:
!  Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!  Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
!----------------------------------------------------------------------------------------
module extraInfo
      
double precision :: oneAtm = 100d0
double precision :: onekPa =   1d0
      
end module extraInfo
    
Subroutine setDLLExtraInfo( iInfo, rInfo )
use extraInfo, only: oneAtm, onekPa
implicit none      
integer,          intent(in) :: iInfo(*)
double precision, intent(in) :: rInfo(*)
      
!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: setDLLExtraInfo

! Makes it possible to use extra info from kernel
!iInfo ( 1 ) = iVal !Special option from PLAXIS Input Phase Settings
!rInfo ( 1 ) = OneAtm
!iInfo ( 2 ) = TempUnit

oneAtm = rInfo(1)
onekPa = oneAtm / 100.0d0
      
End Subroutine setDLLExtraInfo


Subroutine GetModelCount(nMod)!renamed
!
! Return the maximum model number (iMod) in this DLL
!
implicit none
Integer (Kind=4), intent(out) :: nMod

!DEC$ ATTRIBUTES DLLExport, StdCall, reference ::  GetModelCount

nMod = 1 ! Maximum model number (iMod) in current DLL

Return
End ! GetModelCount

Subroutine GetModelName( iMod , ModelName )
!
! Return the name of the different models
!
implicit none
Integer, intent(in) :: iMod
Character (Len= * ), intent(out) :: ModelName
Character (Len=255) tName

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetModelName
!someone ssigne tName
Select Case (iMod)
Case (1)
    ModelName = 'SANISANDMS'
Case Default
    ModelName = 'not in DLL'
End Select
call Add_Str_Length( ModelName )

Return
End ! GetModelName

Subroutine GetParamCount( iMod , nParam )
!
! Return the number of parameters of the different models
!
implicit none
integer, intent(in) :: iMod
integer, intent(out) :: nParam

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamCount

Select Case (iMod)
Case ( 1 )
    nParam = 25
Case Default
    nParam = 0
End Select
Return
End ! GetParamCount

Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!
! Return the parameters name and units of the different models
!
! Units: use F for force unit
!            L for length unit
!            T for time unit
!
implicit none
integer, intent(in) :: iMod , iParam
Character (Len=255), intent(out) :: ParamName, Units

Select Case (iMod)
Case (1)
    ! ModName = 'SANISANDMS'
    Select Case (iParam)
        Case (1)
            ParamName = 'G_0#'       ; Units     = '-'
        Case (2)
            ParamName = '@n'      ; Units     = '-'
        Case (3)
            ParamName = 'M_c#'   ; Units     = '-'
        Case (4)
            ParamName = 'c'    ; Units     = '-'
        Case (5)
            ParamName = '@l_c#'       ; Units     = '-'
        Case (6)
            ParamName = 'e_0#'       ; Units     = '-'
        Case (7)
            ParamName = '@x'       ; Units     = '-'
        Case (8)
            ParamName = 'm'       ; Units     = '-'
        Case (9)
            ParamName = 'h_0#'       ; Units     = '-'
        Case (10)
            ParamName = 'c_h#'       ; Units     = '-'
        Case (11)
            ParamName = 'n^b#'       ; Units     = '-'
        Case (12)
            ParamName = 'A_0#'       ; Units     = '-'
        Case (13)
            ParamName = 'n^d#'       ; Units     = '-'
        Case (14)
            ParamName = '@m_0#'       ; Units     = '-'
        Case (15)
            ParamName = '@z'       ; Units     = '-'
        Case (16)
            ParamName = '@b'       ; Units     = '-'
        Case (17)
            ParamName = 'TolF'       ; Units     = '-'
        Case (18)
            ParamName = 'TolR'       ; Units     = '-'
        Case (19)
            ParamName = 'e_ini#'       ; Units     = '-'
        Case (20)
            ParamName = 'e_max#'       ; Units     = '-'
        Case (21)
            ParamName = 'e_min#'       ; Units     = '-'  
        Case (22)
            ParamName = 'M_c#^b#max'       ; Units     = '-'
        Case (23)
            ParamName = 'ElScheme#'       ; Units     = '-'
        Case (24)
            ParamName = 'h^M#ExcEqu#'     ; Units     = '-' 
        Case (25)
            ParamName = 'h_Max#'     ; Units     = '-'   
        Case Default
            ParamName = '???'     ; Units     = '???'
    End Select  
Case Default
    ! model not in DLL
    ParamName = ' N/A '     ; Units     = ' N/A '
End Select

Return
    End ! GetParamAndUnit

!
! Interfaces to routines to report parameter names, counts etc.
!
Subroutine Add_Str_Length( aString )
Implicit None
Character*(*) aString
Character *255 tString
Integer Lt              ! length of incoming string
!
! routine should add the length of the string as the first character
!
tString = aString
Lt      = Len_Trim(tString)
aString = Char(Lt) // tString(1:Lt)
!
End Subroutine Add_Str_Length

Subroutine GetParamName( iMod , iParam, ParamName )
!
! Return the parameters name of the different models
!
Implicit None
Integer  iMod, iParam
Character (Len=255) ParamName, Units

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamName

Call GetParamAndUnit(iMod,iParam,ParamName,Units)
Call Add_Str_Length( ParamName )


End ! GetParamName

Subroutine GetParamUnit( iMod , iParam, Units )
!
! Return the units of the different parameters of the different models
!
Implicit None
Integer  iMod, iParam
Character (Len=255) ParamName, Units

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamUnit

Call GetParamAndUnit(iMod,iParam,ParamName,Units)
Call Add_Str_Length( Units )

End ! GetParamUnit

Subroutine GetStateVarCount( iMod , nVar )
!
! Return the number of state variables of the different models
!
!DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarCount
implicit none
integer, intent(in) :: iMod
integer, intent(out) :: nVar

Select Case (iMod)
Case (1)
nVar = 31
Case Default
nVar = 0
End Select

Return
End


Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
!
! Return the name and unit of the different state variables of the different models
!
implicit none
integer, intent(in) :: iMod , iVar
Character (Len=255), intent(out) :: Name, Unit

Select Case (iMod)
Case (1)
    Select Case (iVar)
        Case (1)
            Name = 'alpha_11#'         ; Unit = '-'
        Case (2)
            Name = 'alpha_22#'         ; Unit = '-'
        Case (3)
            Name = 'alpha_33#'         ; Unit = '-'
        Case (4)
            Name = 'alpha_12#'         ; Unit = '-'
        Case (5)
            Name = 'alpha_23#'         ; Unit = '-'
        Case (6)
            Name = 'alpha_31#'         ; Unit = '-'
        Case (7)
            Name = 'alpha^#M_11#'       ; Unit = '-'
        Case (8)
            Name = 'alpha^#M_22#'         ; Unit = '-'
        Case (9)
            Name = 'alpha^#M_33#'         ; Unit = '-'
        Case (10)
            Name = 'alpha^#M_12#'         ; Unit = '-'
        Case (11)
            Name = 'alpha^#M_23#'         ; Unit = '-'
        Case (12)
            Name = 'alpha^#M_31#'         ; Unit = '-'
        Case (13)
            Name = 'MM_+#'         ; Unit = '-'
        Case (14)
            Name = 'MM_-#'         ; Unit = '-'
        Case (15)
            Name = 'alpha^#in_11#'         ; Unit = '-'
        Case (16)
            Name = 'alpha^#in_22#'         ; Unit = '-'
        Case (17)
            Name = 'alpha^#in_33#'         ; Unit = '-'
        Case (18)
            Name = 'alpha^#in_12#'         ; Unit = '-'
        Case (19)
            Name = 'alpha^#in_23#'         ; Unit = '-'
        Case (20)
            Name = 'alpha^#in_31#'         ; Unit = '-'
        Case (21)
            Name = 'void ratio'            ; Unit = '-'
        Case (22)
            Name = 'reinitialize'          ; Unit = '-'
        Case (23)
            Name = 'unused'           ; Unit = '-'
        Case (24)
            Name = 'unused'              ; Unit = '-'
        Case (25)
            Name = 'M^b'                    ; Unit = '-'
        Case (26)
            Name = 'M^d'                    ; Unit = '-'
        Case (27)
            Name = 'cos3@q'             ; Unit = '-'
        Case (28)
            Name = '@y'                   ; Unit = '-'
        Case (29)
            Name = 'q'                     ; Unit = '-'
        Case (30)
            Name = 'p'                     ; Unit = '-'
        Case (31)
            Name = 'm^M'                    ; Unit = '-'
        Case Default
            Name='N/A'              ; Unit = '?'
    End Select
Case Default
Name='N/A'                  ; Unit = '?'
End Select

Return
End
    
Subroutine GetStateVarName( iMod , iVar, Name )
!
! Return the name of the different state variables
! of the different models
!
Implicit None
Integer  iMod, iVar
Character (Len=255) Name, Unit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarName

Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
Call Add_Str_Length( Name )

End ! GetStateVarName

Subroutine GetStateVarUnit( iMod , iVar, Unit )
!
! Return the units of the different state variables of the different models
!
Implicit None
Integer  iMod, iVar
Character (Len=255) Name, Unit

!DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarUnit

Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
Call Add_Str_Length( Unit )

End ! GetStateVarUnit
    
!----------------------------------------------------------------------------------    
module SANISANDMS_constants
!
double precision, parameter :: one3 = 1.0d0/3.0d0
double precision, parameter :: one6 = 1.0d0/6.0d0
double precision, parameter :: two3 = 2.0d0/3.0d0
double precision, parameter :: three2 = 3.0D0 / 2.0D0
double precision, parameter :: root23   = DSQRT(2.0D0/3.0D0)
double precision, parameter :: root32   = DSQRT(3.0D0/2.0D0)
double precision, parameter :: zeroP5 = 0.5D0
double precision, parameter :: zero = 0.0D0
double precision, parameter :: one = 1.0d0
double precision, parameter :: two = 2.0d0
double precision, parameter :: three = 3.0D0
!double precision, parameter :: fou = 4.0d0
!double precision, parameter :: fiv = 5.0d0
!double precision, parameter :: six = 6.0d0
!double precision, parameter :: nin = 9.0d0
!double precision, parameter :: ts  = 27.0d0
!
!double precision, parameter :: tot = 1.5d0
!
!double precision, parameter :: pi = 3.14159265358979d0
!
!double precision, parameter :: tol_for_thrsh = 5.0E-02
!double precision, parameter :: tol_for_small = 1.0E-08
!
! unused in state dependent subroutine
double precision, parameter :: onePM5= 1.0D-5
double precision, parameter :: twoeight = 28.0D0
double precision, parameter :: twoeight625 = 28.0D0/625.0D0
double precision, parameter :: onetwofive = 125.0D0
double precision, parameter :: onetwofive625 = 125.0D0/625.0D0
double precision, parameter :: four6 = 4.0D0/6.0D0
double precision, parameter :: fivefoursix = 546.0D0
double precision, parameter :: fivefoursix625 = 546.0D0/625.0D0
double precision, parameter :: fivefour = 54.0D0
double precision, parameter :: fivefour625 = 54.0D0/625.0D0
double precision, parameter :: threeseveneight = 378.0D0
double precision, parameter :: threeseveneight625 = 378.0D0/625.0D0
double precision, parameter :: sixtwofive = 625.0D0
double precision, parameter :: seven27 = 7.0D0/27.0D0
double precision, parameter :: ten27 = 10.0D0/27.0D0
double precision, parameter :: one27 = 1.0D0/27.0D0
double precision, parameter :: onefour336 = 14.0D0/336.0D0
double precision, parameter :: threefive336 = 35.0D0/336.0D0
double precision, parameter :: onesixtwo336 = 162.0D0/336.0D0
double precision, parameter :: onetwofive336 = 125.0D0/336.0D0
double precision, parameter :: fourtwo336 = 42.0D0/336.0D0
double precision, parameter :: twotwofour336 = 224.0D0/336.0D0
double precision, parameter :: twoone336 = 21.0D0/336.0D0

end module SANISANDMS_constants 
    
Module Parameters_position_in_params_iMod1
!! Purpose: Associate a position of the parameter in the array params
!!
Implicit None
!
        integer, parameter :: i_PAR_G0        = 1
        integer, parameter :: i_PAR_nu        = 2
        integer, parameter :: i_PAR_Mc        = 3
        integer, parameter :: i_PAR_c         = 4
        integer, parameter :: i_PAR_lambda_c  = 5
        integer, parameter :: i_PAR_e0        = 6
        integer, parameter :: i_PAR_ksi       = 7
        integer, parameter :: i_PAR_m         = 8
        integer, parameter :: i_PAR_h0        = 9
        integer, parameter :: i_PAR_ch        = 10
        integer, parameter :: i_PAR_nb        = 11
        integer, parameter :: i_PAR_A0        = 12
        integer, parameter :: i_PAR_nd        = 13
        integer, parameter :: i_PAR_mu0       = 14
        integer, parameter :: i_PAR_zeta      = 15
        integer, parameter :: i_PAR_beta      = 16
        integer, parameter :: i_TolF          = 17   !Yield surface Tolerance
        integer, parameter :: i_TolR          = 18   !RKF23 Error Tolerance
        integer, parameter :: i_PAR_e_init    = 19
        integer, parameter :: i_emax          = 20
        integer, parameter :: i_emin          = 21
        integer, parameter :: i_Mbc_max       = 22   !limit to the opening of the BS in TXCompression
        integer, parameter :: i_El_Scheme     = 23   !integration scheme for the the elastic trial: [0] = Single step Forward Euler, [1] = Automatic sub-step RKF23 for trial only, [2] = Automatic sub-step RKF23 for trial and intersection stress 
        integer, parameter :: i_hM_exc_Eq     = 24
        integer, parameter :: i_Custom_h_max  = 25
    end module Parameters_position_in_params_iMod1

Module State_Var_position_in_StVar_iMod1
!! Purpose: Associate a position of the state variable in the StVar and StVar0 arrays
!!
implicit none
        integer, parameter :: i_Alpha_xx      = 1
        integer, parameter :: i_Alpha_yy      = 2
        integer, parameter :: i_Alpha_zz      = 3
        integer, parameter :: i_Alpha_xy      = 4
        integer, parameter :: i_Alpha_yz      = 5
        integer, parameter :: i_Alpha_zx      = 6
        integer, parameter :: i_AlphaM_xx     = 7
        integer, parameter :: i_AlphaM_yy     = 8
        integer, parameter :: i_AlphaM_zz     = 9
        integer, parameter :: i_AlphaM_xy     = 10
        integer, parameter :: i_AlphaM_yz     = 11
        integer, parameter :: i_AlphaM_zx     = 12
        integer, parameter :: i_sMM_plus      = 13
        integer, parameter :: i_sMM_minus     = 14
        integer, parameter :: i_alpha_in_xx   = 15
        integer, parameter :: i_alpha_in_yy   = 16
        integer, parameter :: i_alpha_in_zz   = 17
        integer, parameter :: i_alpha_in_xy   = 18
        integer, parameter :: i_alpha_in_yz   = 19
        integer, parameter :: i_alpha_in_zx   = 20
        integer, parameter :: i_Void          = 21
        integer, parameter :: i_resetState    = 22 
        integer, parameter :: i_PLNoTension   = 23! deprecated output. plastic strain integrated without accuracy
        integer, parameter :: i_PsPsmall      = 24! deprecated output. indicator of inaccurate solution for low stress levels
        integer, parameter :: i_Mb            = 25
        integer, parameter :: i_Md            = 26
        integer, parameter :: i_cos3theta     = 27
        integer, parameter :: i_psi           = 28
        integer, parameter :: i_q_Norm_New    = 29
        integer, parameter :: i_p             = 30
        integer, parameter :: i_MM            = 31

!    
end  module State_Var_position_in_StVar_iMod1
!---------------------------------------------------------------------------------- 