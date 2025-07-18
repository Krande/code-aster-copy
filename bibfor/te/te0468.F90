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

subroutine te0468(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "asterfort/vff2dn.h"
    character(len=16) :: option, nomte
! ......................................................................
!
!     BUT: CALCUL DU FLUX HYDRAULIQUE NORMAL
!          SUR DES SEG2
!          OPTION : 'FLHN_ELGA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: poids, nx, ny, flx(4), fly(4), flun
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfde, igeom
    integer(kind=8) :: iflux, ivectu, k, i, iad
    character(len=24) :: valkm(3)
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ifl, jgano, nbflux, ndim, nnos
    real(kind=8) :: s, t
!-----------------------------------------------------------------------
    call elrefe_info(elrefe='SE2', fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTR', 'L', iflux)
    call jevech('PFLHN', 'E', ivectu)
!
!  CALCUL DU NBRE DE CMP CALCULEES DU FLUX
    if (nomte .eq. 'HM_DPSE3' .or. nomte .eq. 'THM_DPSE3' .or. nomte .eq. &
        'HM_AXIS_SE3' .or. nomte .eq. 'THM_AXIS_SE3' .or. nomte .eq. 'H_DPSE3') then
        nbflux = 1
    else if (nomte .eq. 'THV_DPSE3' .or. nomte .eq. 'THV_AXIS_SE3') then
        nbflux = 2
    elseif (nomte .eq. 'HHM_DPSE3' .or. nomte .eq. 'THH_DPSE3' &
            .or. nomte .eq. 'THHM_DPSE3' .or. nomte .eq. 'HH_DPSE3' &
            .or. nomte .eq. 'HHM_AXIS_SE3' .or. nomte .eq. 'THH_AXIS_SE3' .or. nomte &
            .eq. 'THHM_AXIS_SE3' .or. nomte .eq. 'HH_AXIS_SE3') then
        nbflux = 3
    elseif (nomte .eq. 'HH2M_DPSE3' .or. nomte .eq. 'THH2_DPSE3' &
            .or. nomte .eq. 'THH2M_DPSE3' .or. nomte .eq. 'HH2_DPSE3' &
            .or. nomte .eq. 'HH2M_AXIS_SE3' .or. nomte .eq. 'THH2_AXIS_SE3' .or. nomte &
            .eq. 'THH2M_AXIS_SE3' .or. nomte .eq. 'HH2_AXIS_SE3') then
        nbflux = 4
    else
        valkm(1) = option
        valkm(2) = nomte
        valkm(3) = 'TE0468'
        call utmess('F', 'CALCULEL7_2', nk=3, valk=valkm)
    end if
!    BOUCLE SUR LES CMP
    do ifl = 1, nbflux
!
!    BOUCLE SUR LES POINTS DE GAUSS
        do kp = 1, npg
            k = (kp-1)*nno
! CALCUL DES FLUX AU POINT DE GAUSS KP A PARTIR DES FLUX AUX NOEUDS
            s = 0.d0
            t = 0.d0
            do i = 1, nno
                iad = iflux+2*(ifl-1)+2*nbflux*(i-1)
                s = s+zr(iad)*zr(ivf+k+i-1)
                t = t+zr(iad+1)*zr(ivf+k+i-1)
            end do
            flx(kp) = s
            fly(kp) = t
! CALCUL DE LA NORMALE
            call vff2dn(ndim, nno, kp, ipoids, idfde, &
                        zr(igeom), nx, ny, poids)
            flun = nx*flx(kp)+ny*fly(kp)
            zr(ivectu+nbflux*(kp-1)+ifl-1) = flun
        end do
    end do
end subroutine
