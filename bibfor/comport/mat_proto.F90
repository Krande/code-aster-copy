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
subroutine mat_proto(BEHinteg, &
                     fami, kpg, ksp, poum, imate, itface, &
                     nprops, props)
!
    use Behaviour_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterfort/rcadlv.h"
#include "asterfort/assert.h"
#include "asterfort/rccoma.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in)          :: kpg
    integer(kind=8), intent(in)          :: ksp
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in)          :: imate
    character(len=*), intent(in) :: itface
    integer(kind=8), intent(inout)       :: nprops
    real(kind=8), intent(out)    :: props(*)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (MFront/UMAT)
!
! Calculer les coef. materiau pour l'interface umat et mfront
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHinteg         : parameters for integration of behaviour
!       in   fami      : famille de point de gauss (rigi,mass,...)
!       in   kpg,ksp   : numero du (sous)point de gauss
!       in   poum      : '+' /'-'
!       in   nprops    : en entree : dimension du tableau props
!       in   itface    : nom de l'interface de prototypage
!       out  nprops    : en sortie : nombre de coefficients recuperes dans props
!            props(*)  : coeficients du materiau
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)      :: i, jadr, icodre, ncoef
    real(kind=8) :: phase(5), zalpha
    integer(kind=8)      :: meta_type, nb_phasis
    character(len=16) :: elas_keyword
    integer(kind=8), parameter :: nb_para = 3
    real(kind=8) :: para_vale(nb_para)
    character(len=16), parameter :: para_name(nb_para) = (/'X', 'Y', 'Z'/)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(itface .eq. 'UMAT' .or. itface .eq. 'MFRONT')
    ASSERT(nprops .gt. 0 .and. nprops .le. 197)
    props(1:nprops) = r8nnem()
!
! - Coordinates of current Gauss point
!
    para_vale = BEHinteg%behavESVA%behavESVAGeom%coorElga(kpg, :)
!
! - Get material properties
!
    call rccoma(imate, 'ELAS_META', 0, elas_keyword, icodre)
    if (icodre .eq. 0) then
        call metaGetType(meta_type, nb_phasis)
        call metaGetPhase(fami, poum, kpg, ksp, meta_type, &
                          nb_phasis, phase, zcold_=zalpha)
        call rcadlv(fami, kpg, ksp, poum, imate, ' ', itface, &
                    'LISTE_COEF', 1, ['META'], [zalpha], jadr, ncoef, icodre, 1)
    else
        if (BEHinteg%behavESVA%lGeomInESVA) then
            call rcadlv(fami, kpg, ksp, poum, imate, ' ', itface, &
                        'LISTE_COEF', 0, [' '], [0.d0], jadr, ncoef, icodre, 1)
        else
            call rcadlv(fami, kpg, ksp, poum, imate, ' ', itface, &
                        'LISTE_COEF', nb_para, para_name, para_vale, jadr, ncoef, icodre, 1)
        end if
    end if
    ASSERT(icodre .eq. 0)
!
! - Copy properties
!
    ASSERT(ncoef .le. nprops)
    do i = 1, ncoef
        props(i) = zr(jadr-1+i)
    end do
    nprops = ncoef

end subroutine
