! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine impft0(isor, ift, ibl, fmoy, fetyp, &
                  frms, fmax, fmin)
!
!     IMPRESSION DES FORCES TANGENTIELLES AMV
!
!
!
    implicit none
    integer(kind=8) :: isor, ift
    real(kind=8) :: fmoy, frms, fetyp, fmax, fmin
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ibl
!-----------------------------------------------------------------------
    if (ift .eq. 1) then
        if (ibl .eq. 1) then
            write (isor, *)
            write (isor, *) ' ***** STATISTIQUES FORCE TANGENTE 1 *****'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
            write (isor, *) '!IB! FT1 MOY     ! FT1 E.TYPE  ! FT1 RMS     !',&
     &               ' FT1 MIN     ! FT1 MAX     !'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
        else if (ibl .eq. 0) then
            write (isor, *)
            write (isor, *) ' ***** STATISTIQUES GLOBALES FTANG1 *****'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
            write (isor, *) '!IB! FT1 MOY     ! FT1 E.TYPE  ! FT1 RMS     !',&
     &               ' FT1 MIN     ! FT1 MAX     !'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
        end if
        write (isor, 10) ibl, fmoy, fetyp, frms, fmin, fmax
    else if (ift .eq. 2) then
        if (ibl .eq. 1) then
            write (isor, *)
            write (isor, *) ' ***** STATISTIQUES FORCE TANGENTE 2 *****'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
            write (isor, *) '!IB! FT2 MOY     ! FT2 E.TYPE  ! FT2 RMS     !',&
     &               ' FT2 MIN     ! FT2 MAX     !'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
        else if (ibl .eq. 0) then
            write (isor, *)
            write (isor, *) ' ***** STATISTIQUES GLOBALES FTANG2 *****'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
            write (isor, *) '!IB! FT2 MOY     ! FT2 E.TYPE  ! FT2 RMS     !',&
     &               ' FT2 MIN     ! FT2 MAX     !'
            write (isor, *) '+--+-------------+-------------+-------------+'//&
     &               '-------------+-------------+'
        end if
        write (isor, 10) ibl, fmoy, fetyp, frms, fmin, fmax
    end if
!
10  format(' !', i2, '!', 1pe12.5, ' !', 1pe12.5, ' !', 1pe12.5, ' !',&
      &        1pe12.5, ' !', 1pe12.5, ' !')
!
end subroutine
