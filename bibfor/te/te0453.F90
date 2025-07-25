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

subroutine te0453(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmgeom.h"
#include "asterfort/r8inir.h"
#include "asterfort/writeVector.h"
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DES DEFORMATIONS AUX POINTS DE GAUSS
!          DES ELEMENTS INCOMPRESSIBLES 3D
!
!          OPTION : 'EPSI_ELGA'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
!
    aster_logical :: axi, grand
    real(kind=8) :: eps(6), vpg(162), poids, dfdi(60), f(3, 3), rbid, tmp
    integer(kind=8) :: jgano, ndim, ncmp, nno, npg, kpg, ksig, nnos
    integer(kind=8) :: ipoids, ivf, idfde, igeom, idepl
! ......................................................................
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    grand = .false.
    axi = .false.
    ncmp = 6
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
!
    call r8inir(6, 0.d0, eps, 1)
    call r8inir(162, 0.d0, vpg, 1)
!
    do kpg = 1, npg
        call nmgeom(ndim, nno, axi, grand, zr(igeom), &
                    kpg, ipoids, ivf, idfde, zr(idepl), &
                    .true._1, poids, dfdi, f, eps, &
                    rbid)
!       RECUPERATION DE LA DEFORMATION
        do ksig = 1, ncmp
            if (ksig .le. ndim) then
                tmp = 1.d0
            else
                tmp = sqrt(2.d0)
            end if
            vpg(ncmp*(kpg-1)+ksig) = eps(ksig)/tmp
        end do
!
    end do
!
!      AFFECTATION DU VECTEUR EN SORTIE
    call writeVector('PDEFOPG', npg*ncmp, vpg)
!
end subroutine
