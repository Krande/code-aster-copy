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
subroutine aceat2(nbtuy, eltuy, notuy, nbpart, noex1, &
                  noex2, nbmap, elpar, nopar, nno)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
    integer(kind=8) :: nno, nbtuy, eltuy(nbtuy), notuy(nno*nbtuy), nbpart, noex1(nbpart)
    integer(kind=8) :: noex2(nbpart), nbmap(nbpart), elpar(nbpart, nbtuy)
    integer(kind=8) :: nopar(nbpart, nno, nbtuy)
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES TUYAUX
! ----------------------------------------------------------------------
! IN  :
! ----------------------------------------------------------------------
!
!     STOCKAGE DES NUMEROS DE NOEUDS EXTREMITES
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iext1, iext2, im1, ima, ipa, jma, kp
    integer(kind=8) :: nbe, nbext1, nbext2, nex1, ni1, ni2, ni3
    integer(kind=8) :: ni4, nj1, nj2, nj3, nj4
!-----------------------------------------------------------------------
    nbext1 = 0
    nbext2 = 0
    do ima = 1, nbtuy
        iext1 = 0
        iext2 = 0
!JMP         NI1 = NOTUY(3*IMA-2)
!JMP         NI2 = NOTUY(3*IMA-1)
        ni1 = notuy(nno*(ima-1)+1)
        ni2 = notuy(nno*(ima-1)+2)
        do jma = 1, nbtuy
            if (jma .ne. ima) then
!JMP               NJ1 = NOTUY(NNO*JMA-2)
!JMP               NJ2 = NOTUY(NNO*JMA-1)
                nj1 = notuy(nno*(jma-1)+1)
                nj2 = notuy(nno*(jma-1)+2)
                if (ni1 .eq. nj2) then
                    iext1 = 1
                end if
                if (ni2 .eq. nj1) then
                    iext2 = 1
                end if
            end if
        end do
        if (iext1 .eq. 0) then
            nbext1 = nbext1+1
            noex1(nbext1) = ni1
        end if
        if (iext2 .eq. 0) then
            nbext2 = nbext2+1
            noex2(nbext2) = ni2
        end if
    end do
    ASSERT(nbext1 .eq. nbext2)
    ASSERT(nbext1 .eq. nbpart)
!
! --- VERIFICATION ET STOCKAGE DES PARTIES CONNEXES
!     HYPOTHESE : LES MAILLES SONT TOUTES ORIENTEES DANS LE MEME SENS
!
    im1 = 0
    do ipa = 1, nbpart
        nex1 = noex1(ipa)
! RECHERCHE DE LA PREMIERE MAILLE
        do ima = 1, nbtuy
!JMP            NI1 = NOTUY(NNO*IMA-2)
!JMP            NI2 = NOTUY(NNO*IMA-1)
!JMP            NI3 = NOTUY(NNO*IMA)
            ni1 = notuy(nno*(ima-1)+1)
            ni2 = notuy(nno*(ima-1)+2)
            ni3 = notuy(nno*(ima-1)+3)
            if (nno .eq. 4) then
                ni4 = notuy(nno*(ima-1)+4)
            end if
            if (nex1 .eq. ni1) then
                nbe = 1
                elpar(ipa, nbe) = eltuy(ima)
                nopar(ipa, 1, nbe) = ni1
                nopar(ipa, 2, nbe) = ni2
                nopar(ipa, 3, nbe) = ni3
                if (nno .eq. 4) then
                    nopar(ipa, 4, nbe) = ni4
                end if
                goto 21
            end if
        end do
21      continue
        im1 = ima
41      continue
! SI NI2 EST UNE EXTREMIE, ON CHANGE DE PARTIE
!JMP         NI2 = NOTUY(3*IM1-1)
        ni2 = notuy(nno*(im1-1)+2)
        do kp = 1, nbpart
            if (ni2 .eq. noex2(kp)) goto 11
        end do
! RECHERCHE DE LA MAILLE ATTENANTE A IM1
        do jma = 1, nbtuy
            if (im1 .eq. jma) goto 42
!JMP            NJ1 = NOTUY(3*JMA-2)
!JMP            NJ2 = NOTUY(3*JMA-1)
!JMP            NJ3 = NOTUY(3*JMA)
            nj1 = notuy(nno*(jma-1)+1)
            nj2 = notuy(nno*(jma-1)+2)
            nj3 = notuy(nno*(jma-1)+3)
            if (nno .eq. 4) then
                nj4 = notuy(nno*(jma-1)+4)
            end if
            if (ni2 .eq. nj1) then
                nbe = nbe+1
                elpar(ipa, nbe) = eltuy(jma)
                nopar(ipa, 1, nbe) = nj1
                nopar(ipa, 2, nbe) = nj2
                nopar(ipa, 3, nbe) = nj3
                if (nno .eq. 4) then
                    nopar(ipa, 4, nbe) = nj4
                end if
                goto 43
            end if
42          continue
        end do
43      continue
        im1 = jma
        goto 41
11      continue
        nbmap(ipa) = nbe
    end do
end subroutine
