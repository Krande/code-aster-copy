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
subroutine mfront_get_mater_value(extern_addr, BEHinteg, rela_comp, fami, kpg, &
                                  ksp, imate, props, nprops)
!
    use Behaviour_type
!
    implicit none
!
#include "asterc/mgis_get_props.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/mat_proto.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: extern_addr
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=16), intent(in) :: rela_comp
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp, imate, nprops
    real(kind=8), intent(out) :: props(nprops)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (MFront)
!
! Get material properties
!
! --------------------------------------------------------------------------------------------------
!
! In  extern_addr       : address of the MGISBehaviour object
! In  BEHinteg         : parameters for integration of behaviour
! In  rela_comp        : RELATION comportment
! In  fami             : Gauss family for integration point rule
! In  imate            : coded material address
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npropmax = 197
    integer(kind=8) :: codrel(npropmax), i
    character(len=16) :: nomres(npropmax)
    real(kind=8) :: zalpha
    integer(kind=8)      :: nb_phasis, meta_type, nbcoef
    integer(kind=8), parameter :: nb_para = 3
    real(kind=8) :: para_vale(nb_para)
    character(len=16), parameter :: para_name(nb_para) = (/'X', 'Y', 'Z'/)
!
! --------------------------------------------------------------------------------------------------
!
! - Coordinates of current Gauss point
!
    para_vale = BEHinteg%behavESVA%behavESVAGeom%coorElga(kpg, :)
!
    ASSERT(nprops <= npropmax)
    call mgis_get_props(extern_addr, nomres)

! - Get parameters
!
    if (rela_comp .eq. 'MFRONT') then
        nbcoef = nprops
        call mat_proto(BEHinteg, &
                       fami, kpg, ksp, '+', imate, rela_comp, &
                       nbcoef, props)
        ASSERT(nbcoef == nprops)
    else
! ----- Get the properties values (enter under 'rela_comp' in DEFI_MATERIAU)
        props(1:nprops) = r8nnem()
        if (BEHinteg%behavESVA%tabcod(ZFERRITE) .eq. 1) then
            meta_type = 1
            nb_phasis = 5
            call metaGetPhase(fami, '+', kpg, ksp, meta_type, &
                              nb_phasis, zcold_=zalpha)
            do i = 1, nprops
                if (nomres(i) (1:4) .eq. 'meta') then
                    call rcvalb(fami, 1, 1, '+', imate, &
                                ' ', rela_comp, 1, 'META', [zalpha], &
                                1, nomres(i), props(i), codrel(i), 2)
                else
                    call rcvalb(fami, kpg, ksp, '+', imate, &
                                ' ', rela_comp, 0, ' ', [0.d0], &
                                1, nomres(i), props(i), codrel(i), 1)
                end if
            end do
        elseif (BEHinteg%behavESVA%tabcod(ZALPHPUR) .eq. 1) then
            meta_type = 2
            nb_phasis = 3
            call utmess('F', 'COMPOR4_24')
        else
            if (BEHinteg%behavESVA%lGeomInESVA) then
                call rcvalb(fami, kpg, ksp, '+', imate, &
                            ' ', rela_comp, 0, ' ', [0.d0], &
                            nprops, nomres, props, codrel, 1)
            else
                call rcvalb(fami, kpg, ksp, '+', imate, &
                            ' ', rela_comp, nb_para, para_name, para_vale, &
                            nprops, nomres, props, codrel, 1)
            end if
        end if
    end if

end subroutine
