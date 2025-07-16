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

subroutine te0126(option, nomte)
!
! ------------------------------------------------------------------------------
!
!     CALCUL DE L'OPTION TEMP_ELGA
!
! ------------------------------------------------------------------------------
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvarc.h"
!
    character(len=16) :: nomte, option
!
! ------------------------------------------------------------------------------
!
    integer(kind=8)         :: ndim, nno, npg
    integer(kind=8)         :: nbcou, ipg, icou, sp_couche, iptc, ksp, iret4, indx
    integer(kind=8)         :: jnbspi, jtherm
!
    real(kind=8)    :: tp
!
    character(len=8)    :: fami
!
    aster_logical   :: elem_dkt, elem_gri
!
! ------------------------------------------------------------------------------
!
    if (option .ne. 'TEMP_ELGA') then
        ASSERT(.false.)
    end if
    !
    ! Les éléments traités
    elem_dkt = (nomte .eq. 'MEDKTR3') .or. (nomte .eq. 'MEDKQU4') .or. &
               (nomte .eq. 'MEDSQU4') .or. (nomte .eq. 'MEDSTR3') .or. &
               (nomte .eq. 'MEQ4QU4') .or. (nomte .eq. 'MET3TR3') .or. &
               (nomte .eq. 'MEC3QU9H') .or. (nomte .eq. 'MEC3TR7H')
    !
    elem_gri = (nomte .eq. 'MEGCQU4') .or. (nomte .eq. 'MEGCTR3')
    ASSERT(elem_dkt .or. elem_gri)
    !
    fami = 'RIGI'
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg)
    !
    ! Les couches
    !   nbcou       : Nombre de couche
    !   sp_couche   : sous-points par couche
    sp_couche = 0; nbcou = 0
    if (elem_dkt) then
        call jevech('PNBSP_I', 'L', jnbspi)
        nbcou = zi(jnbspi-1+1)
        ASSERT(nbcou .gt. 0)
        sp_couche = 3
    else if (elem_gri) then
        nbcou = 1
        sp_couche = 1
    end if
    !
    call jevech('PTEMP_R', 'E', jtherm)
    !
    ! Boucle sur les points de Gauss
    do ipg = 1, npg
        ! Boucle sur les couches
        do icou = 1, nbcou
            ! Boucle sur les points dans une couche
            do iptc = 1, sp_couche
                ! sous-point dans la couche
                ksp = (icou-1)*sp_couche+iptc
                call rcvarc(' ', 'TEMP', '+', fami, ipg, ksp, tp, iret4)
                if (iret4 .ne. 0) tp = r8nnem()
                ! Adresse du sous-point sur l'élément support
                indx = jtherm+sp_couche*(nbcou*(ipg-1)+icou-1)+iptc-1
                zr(indx) = tp
            end do
        end do
    end do

end subroutine
