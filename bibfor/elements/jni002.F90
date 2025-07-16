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
subroutine jni002(elrefa, nmaxob, liobj, nbobj)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/elraca.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfd2.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/inmat4.h"
#include "asterfort/jeexin.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: elrefa
!
! ======================================================================
! BUT : INITIALISER LES ELREFA
! ======================================================================
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbpg(MT_NBFAMX), iret, ndim, nno, nnos, nbfpg
    integer(kind=8) :: nmaxob, nbobj, lonfam, ifam, lon2, decal, idim
    integer(kind=8) :: ipg, ino, nderiv, jvi, jvr, npg, nno2, jdim
    real(kind=8) :: rvide
    real(kind=8) :: xpg(3*MT_NBPGMX), poipg(MT_NBPGMX)
    real(kind=8) :: ff(MT_NNOMAX), dff(3, MT_NNOMAX), dff2(3, 3, MT_NNOMAX)
    character(len=24) :: liobj(nmaxob)
    character(len=8) :: nofpg(MT_NBFAMX)
!
!     NBPGMX, NBNOMX, NBFAMX SE REFERER A ELRACA
!
! DEB ------------------------------------------------------------------
!
!
    nbobj = 2
    ASSERT(nmaxob .gt. nbobj)
    liobj(1) = '&INEL.'//elrefa//'.ELRA_I'
    liobj(2) = '&INEL.'//elrefa//'.ELRA_R'
!
!
!
    call jeexin('&INEL.'//elrefa//'.ELRA_I', iret)
    if (iret .gt. 0) goto 999

! - Get list of integration schemes of geometric support
    call elraca(elrefa, &
                nbfpg_=nbfpg, fapg_=nofpg, nbpg_=nbpg, &
                ndim_=ndim, nno_=nno, nnos_=nnos)

    ASSERT((ndim .ge. 0) .and. (ndim .le. 3))
    ASSERT((nno .gt. 0) .and. (nno .le. MT_NNOMAX))
    ASSERT((nbfpg .gt. 0) .and. (nbfpg .le. MT_NBFAMX))
!
!
    call wkvect(liobj(1), 'V V I', 4+nbfpg, jvi)
    zi(jvi-1+1) = ndim
    zi(jvi-1+2) = nbfpg
    zi(jvi-1+3) = nno
    zi(jvi-1+4) = nnos
    lon2 = 0
    do ifam = 1, nbfpg
        npg = nbpg(ifam)
        ASSERT((npg .gt. 0) .and. (npg .le. MT_NBPGMX))
        zi(jvi-1+4+ifam) = npg
!
!       ON VEUT STOCKER : W(IPG),GEOM(IDIM,IPG)
!                         FF(INO,IPG) ET
!                         DFF(IDIM,INO,IPG)
!                         DFF2(IDIM,JDIM,INO,IPG)
!                         MAPGNO(INO,IPG)
        lonfam = npg
        lonfam = lonfam+npg*ndim
        lonfam = lonfam+npg*nno
        lonfam = lonfam+npg*nno*ndim
        lonfam = lonfam+npg*nno*ndim*ndim
        lonfam = lonfam+2+npg*nno
        lon2 = lon2+lonfam
    end do
!
    call wkvect(liobj(2), 'V V R', lon2, jvr)
!
    decal = 0
    do ifam = 1, nbfpg
!
!       -- COORDONNEES ET POIDS DES POINTS DE GAUSS :
!       ------------------------------------------------
        call elraga(elrefa, nofpg(ifam), ndim, npg, xpg, poipg)
        do ipg = 1, npg
            decal = decal+1
            zr(jvr-1+decal) = poipg(ipg)
        end do
        do ipg = 1, npg
            do idim = 1, ndim
                decal = decal+1
                zr(jvr-1+decal) = xpg(ndim*(ipg-1)+idim)
            end do
        end do
!
!
!       -- VALEURS DES FONCTIONS DE FORME :
!       ------------------------------------------------
        do ipg = 1, npg
            call elrfvf(elrefa, xpg(ndim*(ipg-1)+1), ff, nno)
            do ino = 1, nno
                decal = decal+1
                zr(jvr-1+decal) = ff(ino)
            end do
        end do
!
!
!       -- DERIVEES 1ERES DES FONCTIONS DE FORME :
!       ------------------------------------------------
        do ipg = 1, npg
            call elrfdf(elrefa, xpg(ndim*(ipg-1)+1), dff, nno, nderiv)
            ASSERT(nderiv .eq. ndim)
            do ino = 1, nno
                do idim = 1, ndim
                    decal = decal+1
                    zr(jvr-1+decal) = dff(idim, ino)
                end do
            end do
        end do
!
!
!       -- DERIVEES 2EMES DES FONCTIONS DE FORME :
!       ------------------------------------------------
        do ipg = 1, npg
            call elrfd2(elrefa, xpg(ndim*(ipg-1)+1), 9*MT_NNOMAX, dff2, nno2, nderiv)
            if (nderiv .eq. 0) then
                ASSERT(nno2 .eq. 0)
                rvide = r8vide()
            else
                ASSERT(nderiv .eq. ndim)
                ASSERT(nno2 .eq. nno)
            end if
            do ino = 1, nno
                do jdim = 1, ndim
                    do idim = 1, ndim
                        decal = decal+1
                        if (nderiv .ne. 0) then
                            zr(jvr-1+decal) = dff2(idim, jdim, ino)
                        else
                            zr(jvr-1+decal) = rvide
                        end if
                    end do
                end do
            end do
        end do

!
!
!       -- MATRICE GAUSS -> NOEUDS : MAPGNO
!       ------------------------------------------------
!       ON STOCKE DANS L'ORDRE : NNO,NPG,MAPGNO
        call inmat4(elrefa, nno, nnos, npg, nofpg(ifam), &
                    zr(jvr-1+decal+1))
        decal = decal+2+npg*nno
!
    end do
!
999 continue
!
end subroutine
