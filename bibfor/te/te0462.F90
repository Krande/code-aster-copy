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
subroutine te0462(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fmater.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvlg.h"
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
!     CALCUL DES COORDONNEES DES SOUS POINTS DE GAUSS SUR LES FAMILLE
!     DE LA LISTE MATER
!     POUR LES ELEMENTS : DKT
! ----------------------------------------------------------------------
!     NOMBRE MAX DE FAMILLE DANS MATER
    integer(kind=8) :: nfpgmx
!     NOMBRE DE NIVEAUX PAR COUCHE
    integer(kind=8) :: nbniv
!     DIMENSION
    integer(kind=8) :: ndim
    parameter(nfpgmx=10, nbniv=3, ndim=3)
!
    integer(kind=8) :: ndim1, nno, nnos, npg, jgano, idfde, ipoids, ivf
    integer(kind=8) :: igeom, jtab(7), icopg, inbf, icoq, iret, decpo, iad
    integer(kind=8) :: nbsp, nbcou, nfpg, decfpg
    integer(kind=8) :: ifpg, ig, icou, iniv, ino
    real(kind=8) :: pgl(3, 3), xx, yy, zz
    real(kind=8) :: epais, excen, gm1(3), gm2(3), epc, bas, hh
    aster_logical :: grille
    character(len=8) :: fami(nfpgmx)
    data gm1/0.d0, 0.d0, 1.d0/
! ----------------------------------------------------------------------
!
    if (.not. (lteatt('MODELI', 'DKT') .or. lteatt('MODELI', 'GRC'))) then
        ASSERT(.false.)
    end if
    grille = lteatt('MODELI', 'GRC')
!
!     NOMBRE DE NOEUDS
    if (lteatt('TYPMA', 'QU4')) then
        nno = 4
    else if (lteatt('TYPMA', 'TR3')) then
        nno = 3
    else
        ASSERT(.false.)
    end if
!
!
    call jevech('PGEOMER', 'L', igeom)
!
!     ZR(ICOPG) : COORDONNEES DE SOUS-POINTS DE GAUSS
    call tecach('OOO', 'PCOOPGM', 'E', iret, nval=7, &
                itab=jtab)
    icopg = jtab(1)
    nbsp = jtab(7)
    ASSERT(nbsp .gt. 0)
!
    call jevech('PCACOQU', 'L', icoq)
!
    if (grille) then
        excen = zr(icoq+3)
    else
!       ELEMENTS A SOUS POINTS : DKT
        call jevech('PNBSP_I', 'L', inbf)
        nbcou = zi(inbf)
        epais = zr(icoq)
        excen = zr(icoq+4)
        bas = -epais/2.d0+excen
        epc = epais/nbcou
    end if
!
! ON UTILISE LE VECTEUR NORMAL DE LA PLAQUE
    if (nno .eq. 3) then
        call dxtpgl(zr(igeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(igeom), pgl)
    end if
!
    call utpvlg(1, 3, pgl, gm1, gm2)
!
    call fmater(nfpgmx, nfpg, fami)
    decfpg = 0
    do ifpg = 1, nfpg
!
        call elrefe_info(fami=fami(ifpg), ndim=ndim1, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
        do ig = 1, npg
!
!         CALCUL DES COORDONNEES DES POINTS DE GAUSS
            xx = 0.d0
            yy = 0.d0
            zz = 0.d0
            do ino = 1, nno
                xx = xx+zr(igeom+3*(ino-1)+0)*zr(ivf+(ig-1)*nno+ino-1)
                yy = yy+zr(igeom+3*(ino-1)+1)*zr(ivf+(ig-1)*nno+ino-1)
                zz = zz+zr(igeom+3*(ino-1)+2)*zr(ivf+(ig-1)*nno+ino-1)
            end do
!
            if (grille) then
                decpo = ndim*(decfpg+ig-1)
                iad = icopg+decpo
                zr(iad+0) = xx+excen*gm2(1)
                zr(iad+1) = yy+excen*gm2(2)
                zr(iad+2) = zz+excen*gm2(3)
            else
                decpo = nbcou*nbniv*ndim*(decfpg+ig-1)
                do icou = 1, nbcou
                    do iniv = 1, nbniv
                        hh = bas+dble(icou-1)*epc+dble(iniv-1)*epc/2.d0
                        iad = icopg+decpo+(icou-1)*nbniv*ndim+(iniv-1)*ndim
                        zr(iad+0) = xx+hh*gm2(1)
                        zr(iad+1) = yy+hh*gm2(2)
                        zr(iad+2) = zz+hh*gm2(3)
                    end do
                end do
            end if
        end do
        decfpg = decfpg+npg
    end do
!
!
!
end subroutine
