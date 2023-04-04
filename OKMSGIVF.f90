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
!
Subroutine OK_MessageBox(t)
use dfwin
use dfwinty
integer (kind=Int_Ptr_Kind()) hWnd

character*(*) t

character(len=256) :: msg, title

!  Display a messagebox with an OK button
!  Note that all strings must be null terminated for C's sake

msg    = Trim(t)  // char(0)
title = 'UDSM_Example' // char(0)

hWnd = 0

iret = MessageBox( hWnd, msg, title, MB_OK )

Return
End