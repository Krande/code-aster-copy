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

subroutine dxsith(nomte, mater, sigma)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8vide.h"
#include "asterfort/dmatcp.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
    integer(kind=8) :: mater
    real(kind=8) :: sigma(*)
    character(len=16) :: nomte
!
!     BUT:
!       CALCUL DES CONTRAINTES VRAIES
!        (==SIGMA_MECA - SIGMA_THER).
!
! ----------------------------------------------------------------------
!
!
!
!
    integer(kind=8) :: nbepsg
    parameter(nbepsg=8)
!
    integer(kind=8) :: ndim, nnoel, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: i, j, icou, icpg, igauh, ipg, ipgh, iret, ibid, nbcmp, nbcou
    integer(kind=8) :: npgh
    integer(kind=8) :: jnbspi, itab(8)
!
    real(kind=8) :: d(4, 4), angl_naut(3), inst, zero, epsth(nbepsg)
!
    character(len=4) :: fami
!
    aster_logical :: dkg
!
! ----------------------------------------------------------------------
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nnoel, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    zero = 0.0d0
    angl_naut(:) = zero
!
    dkg = .false.
!
    nbcmp = 6
!
    if ((nomte .eq. 'MEDKTG3') .or. (nomte .eq. 'MEDKQG4')) then
        dkg = .true.
    end if
!
! --- RECUPERATION DE L'INSTANT
!     -------------------------
    call tecach('ONO', 'PINSTR', 'L', iret, nval=8, &
                itab=itab)
    ibid = itab(1)
    if (iret .eq. 0) then
        inst = zr(ibid)
    else
        inst = r8vide()
    end if
!
! --- RECUPERATION DU NOMBRE DE COUCHE ET DE SOUS-POINT
!     -------------------------------------------------
    if (dkg) then
        nbcou = 1
        npgh = 1
    else
        call jevech('PNBSP_I', 'L', jnbspi)
        npgh = 3
        nbcou = zi(jnbspi-1+1)
        if (nbcou .le. 0) then
            call utmess('F', 'ELEMENTS_46')
        end if
    end if
!
! --- BOUCLE SUR LES POINTS DE GAUSS DE LA SURFACE:
!     ---------------------------------------------
    do ipg = 1, npg
        do icou = 1, nbcou
            do igauh = 1, npgh
                icpg = nbcmp*npgh*nbcou*(ipg-1)+nbcmp*npgh*(icou-1)+ &
                       nbcmp*(igauh-1)
!
!         -- INTERPOLATION DE ALPHA EN FONCTION DE LA TEMPERATURE
!         ----------------------------------------------------
                ipgh = npgh*(icou-1)+igauh
                call verift('RIGI', ipg, ipgh, '+', mater, &
                            epsth_=epsth(1))
!
                epsth(2) = epsth(1)
                epsth(3) = zero
                epsth(4) = zero
                epsth(5) = zero
                epsth(6) = zero
!
!           -- CALCUL DE LA MATRICE DE HOOKE
!           --------------------------------
                call dmatcp('RIGI', mater, inst, '+', ipg, &
                            igauh, angl_naut, d)
!
!           -- CALCUL DES CONTRAINTES VRAIES (==SIGMA_MECA - SIGMA_THER)
!           -- AU POINT D'INTEGRATION COURANT
!           ------------------------------------------------------------
                do i = 1, 4
                    do j = 1, 4
                        sigma(icpg+i) = sigma(icpg+i)-epsth(j)*d(i, j)
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
