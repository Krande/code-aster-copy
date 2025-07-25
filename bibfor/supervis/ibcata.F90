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
subroutine ibcata(ier)
    implicit none
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterc/jdcget.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/ibcatc.h"
#include "asterfort/lxcadr.h"
#include "asterfort/uldefi.h"
#include "asterfort/utmess.h"
#include "asterfort/utremt.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
    integer(kind=8) :: ier
!     ALLOCATION ET LECTURE DES DIFFERENTS CATALOGUES
!     ------------------------------------------------------------------
!     ROUTINE(S) UTILISEE(S) :
!         GETRES GETFAC GETVXX
!         IBOPER
!     ------------------------------------------------------------------
    real(kind=8) :: temps(7)
    real(kind=8) :: valr
    character(len=16) :: motfac
!     ------------------------------------------------------------------
!     --- DEFAUT POUR LES CATALOGUES D'ELEMENTS
!     --- 3  = CATAELEM
!     --- 4  = BASELEM
!-----------------------------------------------------------------------
    integer(kind=8) :: icata, ieltdf, ier1, iocc, iplace, iun, mxcata
    integer(kind=8) :: mxdfca, nbcata, nbnom, nbocc, nbunit
!-----------------------------------------------------------------------
    parameter(ieltdf=4)
!     ------------------------------------------------------------------
    parameter(mxdfca=4, mxcata=10)
    character(len=32) :: dfnom(mxdfca), nom(mxcata)
    character(len=24) :: valk
    integer(kind=8) :: dfunit(mxdfca), unite(mxcata)
    integer(kind=8) :: i, irestart
!     ------------------------------------------------------------------
!     OPTIONS PAR DEFAUT :
!
    data(dfnom(i), dfunit(i), i=1, mxdfca)/&
     &    'COMMANDE_PRIVEE  ', 03,&
     &    'COMMANDE         ', 02,&
     &    'CATAELEM         ', 04,&
     &    'ELEMBASE         ', 04/
!     ------------------------------------------------------------------
!
    call uldefi(6, ' ', 'MESSAGE', 'A', 'N', &
                'N')
!     --- LA ROUTINE NE S'INTERRESSE QU'AU MOT CLE FACTEUR "CATALOGUE" -
!     --- DANS LA COMMANDE DEBUT
    ier = 0
    motfac = 'CATALOGUE'
!
!     --- NOMBRE DE CATALOGUES SPECIFIES PAR L'UTILISATEUR ---
    call getfac(motfac, nbocc)
!
    if (nbocc .gt. mxcata) then
        ier = ier+1
        call utmess('F', 'SUPERVIS_18')
        nbocc = mxcata
    end if
!
    iun = 1
    do iocc = 1, nbocc
        call getvtx(motfac, 'FICHIER', iocc=iocc, nbval=iun, vect=nom(iocc), &
                    nbret=nbnom)
        call lxcadr(nom(iocc))
        call getvis(motfac, 'UNITE', iocc=iocc, nbval=iun, vect=unite(iocc), &
                    nbret=nbunit)
        if (nbunit .eq. 0) then
            call utremt(nom(iocc), dfnom, mxdfca, iplace)
            if (iplace .gt. 0) unite(iocc) = dfunit(iplace)
        end if
    end do
!
!     --- CATALOGUE DES ELEMENTS ---
    nbcata = 0
    call uttcpu('CPU.IBCATA', 'INIT', ' ')
    call uttcpu('CPU.IBCATA', 'DEBUT', ' ')
    do icata = 1, nbocc
        if (nom(icata) .eq. dfnom(3) .or. nom(icata) .eq. dfnom(4)) then
            if (unite(icata) .gt. 0) then
                call ibcatc(nom(icata), unite(icata), ier1)
                ier = ier+ier1
            end if
            nom(icata) = '  '
            nbcata = nbcata+1
        end if
    end do
    irestart = jdcget('Continue')
    if (nbcata .eq. 0 .and. irestart .eq. 0) then
        call ibcatc(dfnom(ieltdf), dfunit(ieltdf), ier1)
        ier = ier+ier1
    end if
    call uttcpu('CPU.IBCATA', 'FIN', ' ')
    call uttcpr('CPU.IBCATA', 7, temps)
    valr = temps(7)
    valk = ' '
    call utmess('I', 'SUPERVIS_52', sk=valk, sr=valr)
!
!     --- VERIFICATION DE LA COMPLETUDE DE L'EXECUTION ---
    do icata = 1, nbocc
        if (nom(icata) .ne. ' ') then
            call utmess('F', 'SUPERVIS_20', sk=nom(icata))
            ier = ier+1
        end if
    end do
!
    if (ier .gt. 0) then
        call utmess('F', 'SUPERVIS_21')
    end if
!
end subroutine
