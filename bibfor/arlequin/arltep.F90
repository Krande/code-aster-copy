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
! aslint: disable=W1306
!
subroutine arltep(ndim, coors, npgs, &
                  kpgs, nns, fctfs, &
                  elrefc, nnc, coorc, &
                  fctfc, dfdxc, dfdyc, dfdzc)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/assert.h"
#include "asterfort/reereg.h"
#include "asterfort/arlelr.h"
#include "asterfort/elrfdf.h"
#include "asterfort/arljac.h"
#include "asterfort/jedema.h"
!
    integer(kind=8) :: nns, npgs, kpgs
    character(len=8) :: elrefc
    integer(kind=8) :: nnc, ndim
    real(kind=8) ::  coors(ndim*nns), coorc(ndim*nnc)
    real(kind=8) ::  fctfs(nns*npgs)
    real(kind=8) ::  fctfc(ndim*ndim*nnc)
    real(kind=8) ::  dfdxc(nnc), dfdyc(nnc), dfdzc(nnc)

! ----------------------------------------------------------------------

! CALCUL DES MATRICES DE COUPLAGE ARLEQUIN
! OPTION ARLQ_MATR

! EXTENSION DES FF ET DFF D'UNE MAILLE C DANS UNE MAILLE S PLUS GRANDE

! ----------------------------------------------------------------------

! ON AGIT EN 3 TEMPS:
! 1/ CALCUL DES COORD. DS ESPACE REEL DU PT DE GAUSS CONSIDERE
! 2/ CALCUL DES COORD. ESPACE PARA. DANS LA MAILLE C DE CE POINT
!    PAR METHODE DE NEWTON (REEREG)
! 3/ CALCUL DES FF ET DFF DE LA MAILLE C EN CE POINT

! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NPGS   : NOMBRE DE PT DE GAUSS DE LA MAILLE S
! IN  KPGS   : NUMERO PT DE GAUSS DANS MAILLE S
! IN  NNS    : NOMBRE DE NOEUDS DE LA MAILLE S
! IN  FCTFS  : FCT. FORME DE MAILLE S
! IN  ELREFC : TYPE DE LA MAILLE C
! IN  NNC    : NOMBRE DE NOEUDS DE LA MAILLE C
! IN  COORC  : COORD. NOEUDS DE MAILLE C
! OUT FCTFC  : FCT. FORME DE MAILLE C AU POINT DE GAUSS KPGS
!               DE LA MAILLE S
! OUT DFDXC  : DERIVEE FCT. FORME/X DE MAILLE C AU POINT DE GAUSS KPGS
!              DE LA MAILLE S (SI LCALDF = .TRUE.)
! OUT DFDYC  : DERIVEE FCT. FORME/Y DE MAILLE C AU POINT DE GAUSS KPGS
!              DE LA MAILLE S (SI LCALDF = .TRUE.)
! OUT DFDZC  : DERIVEE FCT. FORME/Z DE MAILLE C AU POINT DE GAUSS KPGS
!              DE LA MAILLE S (SI LCALDF = .TRUE.)

! ----------------------------------------------------------------------

    real(kind=8) ::  coorr(3), coorpc(3)
    real(kind=8) :: zero
    integer(kind=8) ::  idim, jdim, ino, ibid, jdecal
    real(kind=8) ::  dfr(3, nnc)
    real(kind=8) ::  invjac(3, 3)
    integer(kind=8) ::  iret

! ----------------------------------------------------------------------
    call jemarq()

    if (kpgs > npgs) then
        ASSERT(.false.)
    end if

! --- COORDONNEES DS ESPACE REEL DU POINT DE GAUSS IPG DE LA MAILLE S

    jdecal = (kpgs-1)*nns
    zero = 0.d0
    coorr(:) = zero
    coorpc(:) = zero

    do ino = 1, nns
        do idim = 1, ndim
            coorr(idim) = coorr(idim)+ &
                          fctfs(jdecal+ino)* &
                          coors(ndim*(ino-1)+idim)
        end do
    end do

! --- COORDONNEES DU PT DE GAUSS DS ESPACE PARA DE MAILLE C

    call reereg('S', elrefc, nnc, coorc, coorr, ndim, &
                coorpc, iret)

! --- VALEURS DES FONCTIONS DE FORME DE ELT C EN COORPC

    call arlelr(elrefc, coorpc, ndim*ndim*nnc, fctfc, ibid)

! --- DERIVEES DES FONCTIONS DE FORME DE REFERENCE EN COORPC

    call elrfdf(elrefc, coorpc, dfr)

! --- INVERSE DE LA JACOBIENNE

    call arljac(nnc, ndim, dfr, coorc, invjac)

! --- DERIVEES DES FONCTIONS DE FORMES CLASSIQUES EN COORPC

    do ino = 1, nnc
        dfdxc(ino) = invjac(1, 1)*dfr(1, ino)
        dfdyc(ino) = invjac(1, 2)*dfr(1, ino)
        dfdzc(ino) = invjac(1, 3)*dfr(1, ino)
        do jdim = 2, ndim
            dfdxc(ino) = dfdxc(ino)+ &
                         invjac(jdim, 1)*dfr(jdim, ino)
            dfdyc(ino) = dfdyc(ino)+ &
                         invjac(jdim, 2)*dfr(jdim, ino)
            dfdzc(ino) = dfdzc(ino)+ &
                         invjac(jdim, 3)*dfr(jdim, ino)
        end do
    end do

    call jedema()

end subroutine
