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
subroutine te0493(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
! ......................................................................
!
!     BUT: CALCUL DU FLUX HYDRAULIQUE NORMAL
!          SUR DES ELEMENTS DE BORD DE 3D (FACE8 ET FACE6)
!          OPTION : 'FLHN_ELGA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: nx, ny, nz, flx, fly, flz, flun, s, t, u, jac
    real(kind=8) :: sx(4, 4), sy(4, 4), sz(4, 4)
    integer(kind=8) :: nno, kp, npg, ipoids, ivf, idfdx, idfdy, igeom
    integer(kind=8) :: iflux, ivectu, k, i, iad
    integer(kind=8) :: idec, jdec, kdec
    character(len=24) :: valkm(3)
    aster_logical :: tria
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ifl, ino, j, jgano, jno, nbflux, ndim
    integer(kind=8) :: nnos
!-----------------------------------------------------------------------
    tria = .false.
!  CALCUL DU NBRE DE CMP CALCULEES DU FLUX
    if (nomte .eq. 'HM_FACE8' .or. nomte .eq. 'THM_FACE8' .or. nomte .eq. 'H_FACE8') then
        nbflux = 1
    elseif (nomte .eq. 'HM_FACE6' .or. nomte .eq. 'THM_FACE6' .or.&
 & nomte .eq. 'H_FACE6') then
        nbflux = 1
        tria = .true.
    else if (nomte .eq. 'THV_FACE8') then
        nbflux = 2
    else if (nomte .eq. 'THV_FACE6') then
        nbflux = 2
        tria = .true.
    elseif (nomte .eq. 'HHM_FACE8' .or. nomte .eq. 'THH_FACE8'&
 &.or. nomte .eq. 'THHM_FACE8' .or. nomte .eq. 'HH_FACE8') then
        nbflux = 3
    elseif (nomte .eq. 'HHM_FACE6' .or. nomte .eq. 'THH_FACE6'&
 &.or. nomte .eq. 'THHM_FACE6' .or. nomte .eq. 'HH_FACE6') then
        nbflux = 3
        tria = .true.
    elseif (nomte .eq. 'HH2M_FACE8' .or. nomte .eq. 'THH2_FACE8' .or. nomte &
            .eq. 'THH2M_FACE8' .or. nomte .eq. 'HH2_FACE8') then
        nbflux = 4
    elseif (nomte .eq. 'HH2M_FACE6' .or. nomte .eq. 'THH2_FACE6' .or. nomte &
            .eq. 'THH2M_FACE6' .or. nomte .eq. 'HH2_FACE6') then
        nbflux = 4
        tria = .true.
    else
        valkm(1) = option
        valkm(2) = nomte
        valkm(3) = 'TE0493'
        call utmess('F', 'CALCULEL7_2', nk=3, valk=valkm)
    end if
!
    if (tria) then
        call elrefe_info(elrefe='TR3', fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    else
        call elrefe_info(elrefe='QU4', fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                         npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    end if
    idfdy = idfdx+1
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTR', 'L', iflux)
    call jevech('PFLHN', 'E', ivectu)
!
! --- CALCUL DES PRODUITS VECTORIELS OMI X OMJ ---
!
    do ino = 1, nno
        i = igeom+3*(ino-1)-1
        do jno = 1, nno
            j = igeom+3*(jno-1)-1
            sx(ino, jno) = zr(i+2)*zr(j+3)-zr(i+3)*zr(j+2)
            sy(ino, jno) = zr(i+3)*zr(j+1)-zr(i+1)*zr(j+3)
            sz(ino, jno) = zr(i+1)*zr(j+2)-zr(i+2)*zr(j+1)
        end do
    end do
!
!    BOUCLE SUR LES CMP
    do ifl = 1, nbflux
!
!    BOUCLE SUR LES POINTS DE GAUSS
        do kp = 1, npg
            k = (kp-1)*nno
! CALCUL DES FLUX AU POINT DE GAUSS KP A PARTIR DES FLUX AUX NOEUDS
            s = 0.d0
            t = 0.d0
            u = 0.d0
            do i = 1, nno
                iad = iflux+3*(ifl-1)+3*nbflux*(i-1)
                s = s+zr(iad)*zr(ivf+k+i-1)
                t = t+zr(iad+1)*zr(ivf+k+i-1)
                u = u+zr(iad+2)*zr(ivf+k+i-1)
            end do
            flx = s
            fly = t
            flz = u
! --- CALCUL DE LA NORMALE AU POINT DE GAUSS KP ---
            nx = 0.0d0
            ny = 0.0d0
            nz = 0.0d0
            kdec = (kp-1)*nno*ndim
            do i = 1, nno
                idec = (i-1)*ndim
                do j = 1, nno
                    jdec = (j-1)*ndim
                    nx = nx+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sx(i, j)
                    ny = ny+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sy(i, j)
                    nz = nz+zr(idfdx+kdec+idec)*zr(idfdy+kdec+jdec)*sz(i, j)
                end do
            end do
            jac = sqrt(nx*nx+ny*ny+nz*nz)
            flun = (nx*flx+ny*fly+nz*flz)/jac
            zr(ivectu+nbflux*(kp-1)+ifl-1) = flun
        end do
    end do
!
end subroutine
