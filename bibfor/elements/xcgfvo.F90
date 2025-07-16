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
subroutine xcgfvo(option, ndim, nnop, fno)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
    character(len=16) :: option
    integer(kind=8) :: ndim, nnop
    real(kind=8) :: fno(ndim*nnop)
!
! person_in_charge: samuel.geniaut at edf.fr
!
!    BUT : CALCUL DES CHARGES VOLUMIQUES AUX NOEUD DE L'ELEM PARENT
!         POUR LES OPTIONS CALC_G, CALC_G_F, CALC_K_G, CALC_K_G_F
!
!
! IN  OPTION : OPTION DE CALCUL
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNOP   : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! OUT FNO    : FORCES NODALES CORRESPONDANT AUX CHARGES VOLUMIQUES
!
!
!
    integer(kind=8) :: igeom, imate, iforc, iforf, itemps, ipesa, irota
    integer(kind=8) :: iret, ino, j, kk
    integer(kind=8), parameter :: mxstac = 1000
    aster_logical :: fonc
    real(kind=8) :: valpar(4), rbid, om, omo, val(1), rhocst
    integer(kind=8) :: icodre(1)
    character(len=8) :: nompar(4)
    character(len=16) :: phenom
!
    rbid = 0.d0
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(ndim*nnop .le. mxstac)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
!     PARAMETRES DES FORCES VOLUMIQUES
    if (option .eq. 'CALC_G_XFEM' .or. option .eq. 'CALC_K_G_XFEM') then
        fonc = .false.
        call jevech('PFRVOLU', 'L', iforc)
    else if (option .eq. 'CALC_G_XFEM_F' .or. option .eq. 'CALC_K_G_XFEM_F') then
        fonc = .true.
        call jevech('PFFVOLU', 'L', iforf)
        call jevech('PINSTR', 'L', itemps)
    else
        ASSERT(.false.)
    end if
!
    call tecach('ONO', 'PPESANR', 'L', iret, iad=ipesa)
    call tecach('ONO', 'PROTATR', 'L', iret, iad=irota)
!
!     INITIALISATION DE FNO
    fno(:) = 0.d0
!
!     ------------------------------------------------------------------
!                     TRAITEMENT DES FORCES VOLUMIQUES
!     ------------------------------------------------------------------
!
!     FORCES VOLUMIQUES FONCTION
    if (fonc) then
!
        nompar(1) = 'X'
        nompar(2) = 'Y'
        valpar(ndim+1) = zr(itemps)
        if (ndim .eq. 2) then
            nompar(3) = 'INST'
        else if (ndim .eq. 3) then
            nompar(3) = 'Z'
            nompar(4) = 'INST'
        end if
!
!       INTERPOLATION DE LA FORCE (FONCTION PAR ELEMENT) AUX NOEUDS
        do ino = 1, nnop
            do j = 1, ndim
                valpar(j) = zr(igeom+ndim*(ino-1)+j-1)
            end do
            do j = 1, ndim
                kk = ndim*(ino-1)+j
                call fointe('FM', zk8(iforf+j-1), ndim+1, nompar, valpar, &
                            fno(kk), iret)
            end do
        end do
!
!     FORCES VOLUMIQUES CONSTANTES (AUX NOEUDS)
    else
!
        do ino = 1, nnop
            do j = 1, ndim
                fno(ndim*(ino-1)+j) = zr(iforc+ndim*(ino-1)+j-1)
            end do
        end do
!
    end if
!
!     ------------------------------------------------------------------
!            TRAITEMENT DES FORCES DE PESANTEUR OU DE ROTATION
!     ------------------------------------------------------------------
!
    if ((ipesa .ne. 0) .or. (irota .ne. 0)) then
!
!       on est sur de la presence de RHO suite a l'appel a cgverho
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
        call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                    ' ', phenom, 1, ' ', [rbid], &
                    1, 'RHO', val, icodre(1), 1)
        rhocst = val(1)
!
        if (ipesa .ne. 0) then
            do ino = 1, nnop
                do j = 1, ndim
                    kk = ndim*(ino-1)+j
                    fno(kk) = fno(kk)+rhocst*zr(ipesa)*zr(ipesa+j)
                end do
            end do
        end if
!
        if (irota .ne. 0) then
            om = zr(irota)
            do ino = 1, nnop
                omo = 0.d0
                do j = 1, ndim
                    omo = omo+zr(irota+j)*zr(igeom+ndim*(ino-1)+j-1)
                end do
                do j = 1, ndim
                    kk = ndim*(ino-1)+j
                    fno(kk) = fno(kk)+rhocst*om*om*(zr(igeom+kk-1)-omo*zr( &
                                                    irota+j))
                end do
            end do
        end if
!
    end if
!
!     ------------------------------------------------------------------
!
end subroutine
