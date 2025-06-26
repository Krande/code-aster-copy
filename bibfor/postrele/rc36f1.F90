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
subroutine rc36f1(nbsigr, nocc, saltij, isk, isl, &
                  nk, nl, n0, nbp12, nbp23, &
                  nbp13, sigr, yapass, typass, nsitup)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbsigr, nocc(*), isk, isl, nk, nl, n0, nsitup, nbp12, nbp23
    integer(kind=8) :: nbp13, sigr(*)
    real(kind=8) :: saltij(*)
    aster_logical :: yapass
    character(len=3) :: typass
!
!     CALCUL DU FACTEUR D'USAGE POUR LES SITUATIONS DE PASSAGE
!     DETERMINATION DU CHEMIN DE PASSAGE
! OUT : N0     : NOMBRE D'OCCURRENCE
! OUT : YAPASS : UNE SITUATION DE PASSAGE EXISTE
! OUT : TYPASS : PASSAGE D'UNE SITUATION A UNE AUTRE
!                1_2 : PASSAGE GROUPE 1 A GROUPE 2
!                1_2 : PASSAGE GROUPE 2 A GROUPE 3
!                1_3 : PASSAGE GROUPE 1 A GROUPE 3
! OUT : NSITUP : NUMERO DU CHEMIN DE SITUATION DE PASSAGE
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ig1, ig2, nbsips, jnpass, i, k, i1, nsitu, numg1, numg2
    integer(kind=8) :: sipass, npass, ioc1, ioc2
    real(kind=8) :: salmia, salmib, salt1, salt2, salt3, salt4, saltam, saltbm
    aster_logical :: chemin
    integer(kind=8), pointer :: situ_group(:) => null()
!     ------------------------------------------------------------------
!
    call jeveuo('&&RC32SI.SITU_GROUP', 'L', vi=situ_group)
!
    ioc1 = sigr(isk)
    ioc2 = sigr(isl)
    numg1 = situ_group(1+2*ioc1-2)
    ig1 = situ_group(1+2*ioc1-1)
    numg2 = situ_group(1+2*ioc2-2)
    ig2 = situ_group(1+2*ioc2-1)
!
    yapass = .false.
    typass = '?_?'
    nsitup = 0
!
    if (numg1 .eq. numg2) then
! ------ MEME GROUPE
        n0 = min(nk, nl)
        goto 999
    else if (numg1 .eq. ig2) then
! ------ MEME GROUPE
        n0 = min(nk, nl)
        goto 999
    else if (numg2 .eq. ig1) then
! ------ MEME GROUPE
        n0 = min(nk, nl)
        goto 999
    else if (ig1 .eq. ig2) then
! ------ MEME GROUPE
        n0 = min(nk, nl)
        goto 999
    end if
!
    if ((numg1 .eq. 1 .and. numg2 .eq. 2) .or. (numg1 .eq. 2 .and. numg2 .eq. 1)) then
        if (nbp12 .eq. 0) then
            if ((ig1 .eq. 1 .and. ig2 .eq. 3) .or. (ig1 .eq. 3 .and. ig2 .eq. 1)) then
                typass = '1_3'
                yapass = .true.
            elseif ((ig1 .eq. 2 .and. ig2 .eq. 3) .or. (ig1 .eq. 3 &
                                                        .and. ig2 .eq. 2)) then
                typass = '2_3'
                yapass = .true.
            end if
        else
            typass = '1_2'
            yapass = .true.
        end if
    elseif ((numg1 .eq. 2 .and. numg2 .eq. 3) .or. (numg1 .eq. 3 .and. &
                                                    numg2 .eq. 2)) then
        if (nbp23 .eq. 0) then
            if ((ig1 .eq. 1 .and. ig2 .eq. 2) .or. (ig1 .eq. 2 .and. ig2 .eq. 1)) then
                typass = '1_2'
                yapass = .true.
            elseif ((ig1 .eq. 1 .and. ig2 .eq. 3) .or. (ig1 .eq. 3 &
                                                        .and. ig2 .eq. 1)) then
                typass = '1_3'
                yapass = .true.
            end if
        else
            typass = '2_3'
            yapass = .true.
        end if
    elseif ((numg1 .eq. 1 .and. numg2 .eq. 3) .or. (numg1 .eq. 3 .and. &
                                                    numg2 .eq. 1)) then
        if (nbp13 .eq. 0) then
            if ((ig1 .eq. 1 .and. ig2 .eq. 2) .or. (ig1 .eq. 2 .and. ig2 .eq. 1)) then
                typass = '1_2'
                yapass = .true.
            elseif ((ig1 .eq. 2 .and. ig2 .eq. 3) .or. (ig1 .eq. 3 &
                                                        .and. ig2 .eq. 2)) then
                typass = '2_3'
                yapass = .true.
            end if
        else
            typass = '1_3'
            yapass = .true.
        end if
    end if
!
! --- RECHERCHE DU CHEMIN DE PASSAGE
!
    call jelira('&&RC32SI.PASSAGE_'//typass, 'LONUTI', nbsips)
    call jeveuo('&&RC32SI.PASSAGE_'//typass, 'L', jnpass)
    chemin = .false.
    salmia = 1.d+50
    salmib = 1.d+50
    do i = 1, nbsips
        sipass = zi(jnpass+i-1)
        do k = 1, nbsigr
            if (sigr(k) .eq. sipass) then
                ioc1 = k
                goto 14
            end if
        end do
        call utmess('F', 'POSTRCCM_36')
14      continue
        npass = max(nocc(2*(ioc1-1)+1), nocc(2*(ioc1-1)+2))
        if (npass .eq. 0) goto 10
        chemin = .true.
! --------- ON RECHERCHE LE MIN DES SALT MAX
        saltam = 0.d0
        saltbm = 0.d0
        do k = 1, nbsigr
            i1 = 4*nbsigr*(k-1)
!            COLONNE _A
            salt1 = saltij(i1+4*(ioc1-1)+1)
            salt3 = saltij(i1+4*(ioc1-1)+3)
            if (salt1 .gt. saltam) then
                saltam = salt1
                nsitu = ioc1
            end if
            if (salt3 .gt. saltam) then
                saltam = salt3
                nsitu = ioc1
            end if
!            COLONNE _B
            salt2 = saltij(i1+4*(ioc1-1)+2)
            salt4 = saltij(i1+4*(ioc1-1)+4)
            if (salt2 .gt. saltbm) then
                saltbm = salt2
                nsitu = ioc1
            end if
            if (salt4 .gt. saltbm) then
                saltbm = salt4
                nsitu = ioc1
            end if
        end do
!
        if (saltam .lt. salmia) then
            salmia = saltam
            nsitup = nsitu
        end if
        if (saltbm .lt. salmib) then
            salmib = saltbm
            nsitup = nsitu
        end if
!
10      continue
    end do
    if (chemin) then
        npass = max(nocc(2*(nsitup-1)+1), nocc(2*(nsitup-1)+2))
        n0 = min(nk, nl, npass)
    else
        yapass = .false.
        n0 = min(nk, nl)
    end if
!
999 continue
!
end subroutine
