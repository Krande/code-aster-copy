! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgTableG(cgField, cgTheta, cgStudy)
!
use calcG_type
!
    implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utimsd.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajvi.h"
#include "asterfort/tbajvk.h"
#include "asterfort/tbajvr.h"
#include "asterfort/titre.h"
#include "asterfort/tableGinit.h"
#include "asterfort/getvis.h"

!
    type(CalcG_field), intent(inout) :: cgField
    type(CalcG_theta), intent(in) :: cgTheta
    type(CalcG_study), intent(in) :: cgStudy
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!     Create table
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nxpara = 15, nbmxpa = 20
    integer            :: nbpara
    integer            :: livi(nbmxpa)
    real(kind=8)       :: livr(nbmxpa), g_irwin
    complex(kind=8)    :: livc(nbmxpa)
    character(len=8)   :: litypa(nxpara), k8bid
    character(len=24)  :: livk(nbmxpa), basloc
    character(len=16)  :: linopa(nxpara)
!
!   A EFFACER
    real(kind=8)       :: temp, coor_x, coor_y, coor_z, abs_curv
!
    call jemarq()
!
    call tableGinit(cgField%table_g, cgStudy%option, cgField%ndim, nxpara, &
                    cgStudy%l_modal, nbpara, linopa, litypa)
!
!   Récupération des coord des pts du fond
    basloc=cgTheta%crack//'.BASLOC'

!   Application de la formule d'Irwin
    g_irwin= cgStudy%gth(5)*cgStudy%gth(5)+cgStudy%gth(6)*cgStudy%gth(6)&
            + cgStudy%gth(7)*cgStudy%gth(7)
!
!====== VALEURS BIDON A EFFACER ==========
    k8bid = 'K8_BIDON'
    temp = 0.d0
    coor_x = 0.d0
    coor_y = 0.d0
    coor_z = 0.d0
    abs_curv = 0.d0
!   Pour coor_x, coor_y, coor_z , il faut
!   récupérer les valeurs de basloc en 2D
!   et identifier les pts du fond de basloc
!   pour le cas 3D : A FAIRE
!=========================================
!
    call tbajvr(cgField%table_g, nbpara, 'TEMP', temp, livr)
    call tbajvk(cgField%table_g, nbpara, 'COMPORTEMENT', k8bid, livk)
!
    call tbajvi(cgField%table_g, nbpara, 'NUME_ORDRE', cgStudy%nume_ordre, livi)
    call tbajvr(cgField%table_g, nbpara, 'INST', cgStudy%time, livr)
!
    call tbajvr(cgField%table_g, nbpara, 'COOR_X', coor_x, livr)
    call tbajvr(cgField%table_g, nbpara, 'COOR_Y', coor_y, livr)

    if (cgField%ndim.eq.3) then
        call tbajvr(cgField%table_g, nbpara, 'COOR_Z', coor_z, livr)
        call tbajvr(cgField%table_g, nbpara, 'ABSC_CURV_NORM', abs_curv, livr)
    endif

    call tbajvr(cgField%table_g, nbpara, 'G', cgStudy%gth(1), livr)
!
    if (cgStudy%option .eq. "K") then
        call tbajvr(cgField%table_g, nbpara, 'K1', cgStudy%gth(2), livr)
        call tbajvr(cgField%table_g, nbpara, 'K2', cgStudy%gth(3), livr)
        if (cgField%ndim .eq. 3) then
            call tbajvr(cgField%table_g, nbpara, 'K3', cgStudy%gth(4), livr)
        endif
        call tbajvr(cgField%table_g, nbpara, 'G_IRWIN', g_irwin, livr)        
    endif
!
    call tbajli(cgField%table_g, nbpara, linopa, livi, livr, livc, livk, 0)
!
    call jedema()
!
end subroutine
