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
!--------------------------------------------------------------------------------------
Subroutine User_Mod ( IDTask, iMod, IsUndr, &
                      iStep, iTer, iEl, Int,   &
                      X, Y, Z, &
                      Time0, dTime,    &
                      Props, Sig0, Swp0, StVar0,   &
                      dEps, D, BulkW,  &
                      Sig, Swp, StVar, ipl,    &
                      nStat, NonSym, iStrsDep, iTimeDep, iTang,    &
                      iPrjDir, iPrjLen, iAbort )

!Test
! Purpose: User supplied soil model
!          Example: iMod=1 : Drucker-Prager
!                   iMod=2 : Mohr-Coulomb
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetric D-matrix
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
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent matrix
!  iAbort   O   I    : =1 to force stopping of calculation
!
!use MIntegration
implicit none
      
integer, intent(in) :: IDTask, iMod, IsUndr, iStep, iTer, iEl,  &
    Int, iPrjDir(*), iPrjLen
double precision, intent(in) :: X, Y, Z, Time0, dTime, Props(*),    &
    Swp0, dEps(*)
double precision, intent(inout) :: Sig0(*), StVar0(*)
double precision, intent(out) :: D(6,6), BulkW,Sig(*),Swp,StVar(*)
integer, intent(out) :: nStat, NonSym, iStrsDep, iTimeDep, &
    iTang, iAbort
integer, intent(inout) :: ipl

integer iounit

Data iounit / 0 /
Save iounit
!
!---  Local variables
!
Character*100 BaseName, BaseName1, BaseName2, BaseName3

!DEC$ ATTRIBUTES DLLExport, StdCall, reference :: User_Mod

BaseName = 'SANISANDMS'
BaseName1 = 'TensionStrategy'
BaseName2 = 'DifferenceSig1Sig3'
BaseName3 = 'DifferenceSig1Sig3_2'
! Possibly open a file for debugging purposes
If (iounit==0) Then

    Call Open_Dbg_File (iPrjDir, iPrjLen, BaseName)

    Call Open_Dbg_File1 (iPrjDir, iPrjLen, BaseName1)

    Call Open_Dbg_File2 (iPrjDir, iPrjLen, BaseName2)

    Call Open_Dbg_File3 (iPrjDir, iPrjLen, BaseName3)
    
    Write(100,*)'File 1 opened: ', Trim(baseName)
    Write(102,*)'File 2 opened: ', Trim(BaseName1)
    Write(103,*)'File 3 opened: ', Trim(BaseName2)
    Write(104,*)'File 3 opened: ', Trim(BaseName3)

    ! maybe write some more info on version to debug file ?
    !write(100,*)'Compiled : ',__DATE__,' ', __TIME__
    !call flush(100)
    !!DEC$ IF DEFINED(_X86_)
    ! this 32-bit ??
    !Write(100,1050)'IF32',__INTEL_COMPILER
    !!DEC$ ELSE
    ! this 64-bit ??
    !Write(100,1050)'IF64', __INTEL_COMPILER
    !!DEC$ ENDIF

1050   format ( 1X,A,1x,I0 )
    iounit = 1
!    Call WriVec(1,'Props',Props,50)
    Call Flush(100)
End If

Call WriIvl( -1, 'iounit',iounit )
Call WriIvl( -1, 'IDTask',IDTask )
Select Case (iMod)
Case (1)   ! SANISANDMS
    Call SANISANDMS( IDTask, iMod, IsUndr, iStep, iTer, iEl, &
                   Int, X, Y, Z, Time0, dTime, &
                   Props, Sig0, Swp0, StVar0,  &
                   dEps, D, BulkW, Sig, Swp, StVar, ipl,   &
                   nStat, NonSym, iStrsDep, iTimeDep, iTang,   &
                   iAbort)
Case Default
    Write(1,*) 'invalid model number in UsrMod', iMod
    Write(1,*) 'IDTask: ',IDTask
    Stop 'invalid model number in UsrMod'
    iAbort=1
    Return
End Select ! iMod
!If (IDTask == 5.And.iel+int==2) Then
!    Write(1,*)'nStat   : ',nStat
!    Write(1,*)'NonSym  : ',NonSym
!    Write(1,*)'StrsDep : ',iStrsDep
!    Write(1,*)'TimeDep : ',iTimeDep
!    Write(1,*)'Tangent : ',iTang
!End If
!If (IDTask == -333 .And. iel+int == -1234) Then
!!        Write(100,*)'IDTask: ',IDTask,' iStep,iTer',iStep,iTer 
!!        Call Flush(100)
!End If
Call WriIvl( -1, 'IDTask end',IDTask )

Return
End ! User_Mod

! **********************************************************************

Subroutine Open_Dbg_File( iPrjDir, iPrjLen, BaseName )
Implicit None

Integer, intent(in) :: iPrjLen, iPrjDir(*)
Character*(*), intent(in):: BaseName

Character*255 PrjDir, Dbg_Name

Integer i, nErr, ios

PrjDir=' '
Do i=1,iPrjLen
    PrjDir(i:i) = Char( iPrjDir(i) )
End Do

Dbg_Name=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
nErr=0
1 Continue
Open( Unit= 100, File= Dbg_Name,iostat=ios)
If (ios==0) Close(Unit=100,Status='delete',iostat=ios)

If (ios.Ne.0) Then
    !
    ! in case of error try ...udsmex1 or udsmex2 or ..
    !
    nErr=nErr+1
    Dbg_Name=PrjDir(:iPrjLen)//'data.'//    &
      trim(BaseName)//char(48+nErr)//'.rr0'

    If (nErr.Lt.10) Goto 1
End If

Open( Unit= 100, File= Dbg_Name,blocksize=4096)

End Subroutine Open_Dbg_File

! **********************************************************************


!Subroutine GetStateVarCount( iMod , nVar )
!!
!! Return the number of state variables of the different models
!!
!Implicit None
!Integer  iMod, nVar
!
!!DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarCount
!
!Call Get_StateVar_Count( iMod , nVar )
!
!End ! GetStateVarCount 
! **********************************************************************

Subroutine Open_Dbg_File1( iPrjDir, iPrjLen, BaseName )
Implicit None

Integer, intent(in) :: iPrjLen, iPrjDir(*)
Character*(*), intent(in):: BaseName

Character*255 PrjDir, Dbg_Name

Integer i, nErr, ios

PrjDir=' '
Do i=1,iPrjLen
    PrjDir(i:i) = Char( iPrjDir(i) )
End Do

Dbg_Name=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
nErr=0
1 Continue
Open( Unit= 102, File= Dbg_Name,iostat=ios)
If (ios==0) Close(Unit=102,Status='delete',iostat=ios)

If (ios.Ne.0) Then
    !
    ! in case of error try ...udsmex1 or udsmex2 or ..
    !
    nErr=nErr+1
    Dbg_Name=PrjDir(:iPrjLen)//'data.'//    &
      trim(BaseName)//char(48+nErr)//'.rr0'

    If (nErr.Lt.10) Goto 1
End If

Open( Unit= 102, File= Dbg_Name,blocksize=4096)

End Subroutine Open_Dbg_File1
!------------------------------------------------------------------------------------      
Subroutine Open_Dbg_File2( iPrjDir, iPrjLen, BaseName )
Implicit None

Integer, intent(in) :: iPrjLen, iPrjDir(*)
Character*(*), intent(in):: BaseName

Character*255 PrjDir, Dbg_Name

Integer i, nErr, ios

PrjDir=' '
Do i=1,iPrjLen
    PrjDir(i:i) = Char( iPrjDir(i) )
End Do

Dbg_Name=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
nErr=0
1 Continue
Open( Unit= 103, File= Dbg_Name,iostat=ios)
If (ios==0) Close(Unit=103,Status='delete',iostat=ios)

If (ios.Ne.0) Then
    !
    ! in case of error try ...udsmex1 or udsmex2 or ..
    !
    nErr=nErr+1
    Dbg_Name=PrjDir(:iPrjLen)//'data.'//    &
      trim(BaseName)//char(48+nErr)//'.rr0'

    If (nErr.Lt.10) Goto 1
End If

Open( Unit= 103, File= Dbg_Name,blocksize=4096)

End Subroutine Open_Dbg_File2
!------------------------------------------------------------------------------------     
Subroutine Open_Dbg_File3( iPrjDir, iPrjLen, BaseName )
Implicit None

Integer, intent(in) :: iPrjLen, iPrjDir(*)
Character*(*), intent(in):: BaseName

Character*255 PrjDir, Dbg_Name

Integer i, nErr, ios

PrjDir=' '
Do i=1,iPrjLen
    PrjDir(i:i) = Char( iPrjDir(i) )
End Do

Dbg_Name=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
nErr=0
1 Continue
Open( Unit= 104, File= Dbg_Name,iostat=ios)
If (ios==0) Close(Unit=104,Status='delete',iostat=ios)

If (ios.Ne.0) Then
    !
    ! in case of error try ...udsmex1 or udsmex2 or ..
    !
    nErr=nErr+1
    Dbg_Name=PrjDir(:iPrjLen)//'data.'//    &
      trim(BaseName)//char(48+nErr)//'.rr0'

    If (nErr.Lt.10) Goto 1
End If

Open( Unit= 104, File= Dbg_Name,blocksize=4096)

End Subroutine Open_Dbg_File3


! **********************************************************************