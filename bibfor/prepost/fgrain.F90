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
!
subroutine fgrain(pic, npic, itrv, ncyc, sigmin, &
                  sigmax)
!      COMPTAGE DES CYCLES PAR LA METHODE RAINFLOW (POSTDAM)
!       ----------------------------------------------------------------
!      IN  PIC     VECTEUR  DES PICS
!          NPIC    NOMBRE   DE  PICS
!          ITRV    VECTEUR  DE TRAVAIL ENTIER
!      OUT SIGMAX  CONTRAINTES MAXIMALES DES CYCLES
!          SIGMIN  CONTRAINTES MINIMALES DES CYCLES
!      OUT  NCYC    NOMBRE  DE  CYCLE
!       ----------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    real(kind=8) :: pic(*), x, y, e1, e2, e3, sigmax(*), sigmin(*)
    real(kind=8) :: r1, r2, rd, rad
    integer(kind=8) :: npic, ncyc, itrv(*), npicb
    aster_logical :: lresi, cyczer
!       ----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifm, j, k, niv, npicr
!-----------------------------------------------------------------------
    lresi = .false.
    npicb = npic
    cyczer = .true.
!
!     --- RECUPERATION DU NIVEAU D'IMPRESSION ---
!
    call infniv(ifm, niv)
!
    do i = 1, npicb
        itrv(i) = i
    end do
    ncyc = 0
!
    do i = 2, npicb
        if ((pic(i) .gt. pic(1)) .or. (pic(i) .lt. pic(1))) then
            cyczer = .false.
        end if
    end do
!
    if (cyczer) then
        sigmax(1) = pic(1)
        sigmin(1) = pic(1)
        ncyc = 1
!
        call utmess('A', 'FATIGUE1_39')
!
        goto 999
    end if
!
!
1   continue
    i = 1
    j = 1
!
2   continue
    if (i+3 .gt. npicb) then
        goto 100
    end if
    e1 = abs(pic(itrv(i+1))-pic(itrv(i)))
    e2 = abs(pic(itrv(i+2))-pic(itrv(i+1)))
    e3 = abs(pic(itrv(i+3))-pic(itrv(i+2)))
!
    if (e1 .ge. e2 .and. e3 .ge. e2) then
        ncyc = ncyc+1
        if (pic(itrv(i+1)) .ge. pic(itrv(i+2))) then
            sigmax(ncyc) = pic(itrv(i+1))
            sigmin(ncyc) = pic(itrv(i+2))
        else
            sigmax(ncyc) = pic(itrv(i+2))
            sigmin(ncyc) = pic(itrv(i+1))
        end if
        do k = i+2, j+2, -1
            itrv(k) = itrv(k-2)
        end do
        j = j+2
        i = j
        goto 2
    else
        i = i+1
        goto 2
    end if
!
!  --- TRAITEMENT DU RESIDU -------
!
100 continue
    if (.not. lresi) then
        npicr = npicb-2*ncyc
        do i = 1, npicr
            itrv(i) = itrv(2*ncyc+i)
        end do
        r1 = pic(itrv(1))
        r2 = pic(itrv(2))
        rad = pic(itrv(npicr-1))
        rd = pic(itrv(npicr))
        x = (rd-rad)*(r2-r1)
        y = (rd-rad)*(r1-rd)
        if (x .gt. 0.d0 .and. y .lt. 0.d0) then
            do i = 1, npicr
                itrv(i+npicr) = itrv(i)
            end do
            npicb = 2*npicr
        else if (x .gt. 0.d0 .and. y .ge. 0.d0) then
! -- ON ELIMINE  R1 ET RN
            do i = npicr, 2, -1
                itrv(i+npicr-2) = itrv(i)
            end do
            npicb = 2*npicr-2
        else if (x .lt. 0.d0 .and. y .lt. 0.d0) then
! -- ON ELIMINE R1
            do i = npicr, 2, -1
                itrv(i+npicr-1) = itrv(i)
            end do
            npicb = 2*npicr-1
        else if (x .lt. 0.d0 .and. y .ge. 0.d0) then
! -- ON ELIMINE RN
            do i = npicr, 1, -1
                itrv(i+npicr-1) = itrv(i)
            end do
            npicb = 2*npicr-1
        end if
        lresi = .true.
        goto 1
    end if
!
!     --- IMPRESSION DES PICS EXTRAITS DE LA FONCTION ----
    if (niv .eq. 2) then
        write (ifm, *)
        write (ifm, '(1X,A)') 'PICS APRES LE COMPTAGE RAINFLOW'
        write (ifm, *)
        write (6, *) 'NOMBRE DE CYCLES = ', ncyc
        write (ifm, *)
        write (ifm, '(1X,A)') '     CHARGEMENT_MAX     CHARGEMENT_MIN'
        write (ifm, *)
        write (ifm, '(2(1X,E18.6))') (sigmax(i), sigmin(i), i=1, ncyc)
!         DO 106 I=1,NCYC
!             WRITE (IFM,'(2(1X,E18.6))'), SIGMAX(I),SIGMIN(I)
! 106     CONTINUE
!
    end if
!
!
999 continue
!
!
!
end subroutine
