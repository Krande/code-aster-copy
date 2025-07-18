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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine te0353(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/meta_vpta_coef.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/Metallurgy_type.h"
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN, C_PLAN, AXIS
! Option: CHAR_MECA_META_Z
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  i, k, lgpg, iret, ispg
    real(kind=8) :: sigmo
    character(len=4) :: fami
    character(len=16) :: metaPhasName, rela_comp, valk(2)
    integer(kind=8) :: j_mate, j_mater
    integer(kind=8) :: meta_type, nb_phasis
    real(kind=8) :: young, nu, deuxmu
    integer(kind=8) :: j_sigm
    integer(kind=8) :: nb_sigm, elas_id
    integer(kind=8) :: j_poids, j_vf, j_dfde, j_geom
    integer(kind=8) :: nno, ipg, npg, jtab(7)
    integer(kind=8) :: j_vectu
    integer(kind=8) :: j_vari
    real(kind=8) :: sig(4), sigdv(4)
    real(kind=8) :: dfdx(9), dfdy(9), poids, r, co_axis
    real(kind=8) :: coef, trans
    real(kind=8) :: zcold_curr
    real(kind=8) :: phas_prev(5), phas_curr(5), temp
    aster_logical :: l_axi, l_temp
    character(len=16) :: elas_keyword
    character(len=16) :: metaRela, metaGlob
    character(len=16), pointer :: compor(:) => null()
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!
    ispg = 1
    fami = 'RIGI'
!
! - Element reference
!
    call elrefe_info(fami=fami, nno=nno, &
                     npg=npg, jpoids=j_poids, jvf=j_vf, jdfde=j_dfde)
    call tecach('OOO', 'PVARIPR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
    l_axi = .false.
    if (lteatt('AXIS', 'OUI')) then
        l_axi = .true.
    end if
!
! - Geometry
!
    call jevech('PGEOMER', 'L', j_geom)
!
! - Comportement
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_comp = compor(RELA_NAME)
    metaRela = compor(META_RELA)
    metaGlob = compor(META_GLOB)
!
! - Cannot evaluate command variables effect for Mfront behaviors
!
    if ((metaRela .eq. 'MFRONT') .or. (metaRela .eq. 'AnisoLemaitre') &
        .or. (metaRela(1:4) .eq. 'Meta')) then
        goto 99
    end if
!
! - Get type of phasis
!
    metaPhasName = compor(META_PHAS)
    call metaGetType(meta_type, nb_phasis)
    ASSERT(nb_phasis .le. 5)
    if ((meta_type .eq. META_NONE) .or. (rela_comp .eq. 'META_LEMA_ANI')) then
        goto 99
    end if
!
! - Check type of phasis
!
    valk(1) = metaPhasName
    if (metaPhasName .eq. 'ACIER') then
        if (meta_type .ne. META_STEEL) then
            valk(2) = 'ZIRC'
            call utmess('F', 'COMPOR3_8', nk=2, valk=valk)
        end if
    elseif (metaPhasName .eq. 'ZIRC') then
        if (meta_type .ne. META_ZIRC) then
            valk(2) = 'ACIER'
            call utmess('F', 'COMPOR3_8', nk=2, valk=valk)
        end if
    else
        goto 99
    end if
!
! - Stresses
!
    nb_sigm = 4
    call jevech('PCONTMR', 'L', j_sigm)
!
! - Internal variables
!
    call jevech('PVARIPR', 'L', j_vari)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', j_mate)
    j_mater = zi(j_mate)
!
! - Output vector
!
    call jevech('PVECTUR', 'E', j_vectu)
    do i = 1, nno
        zr(j_vectu+2*i-2) = 0.d0
        zr(j_vectu+2*i-1) = 0.d0
    end do
!
    do ipg = 1, npg
        k = (ipg-1)*nno
!
! ----- Derived of shape functions
!
        call dfdm2d(nno, ipg, j_poids, j_dfde, zr(j_geom), &
                    poids, dfdx, dfdy)
!
! ----- Radius for axi-symmetric
!
        r = 0.d0
        do i = 1, nno
            r = r+zr(j_geom+2*(i-1))*zr(j_vf+k+i-1)
        end do
!
! ----- Get current temperature
!
        call rcvarc(' ', 'TEMP', '+', fami, ipg, &
                    1, temp, iret)
        l_temp = iret .eq. 0
!
! ----- Get phasis
!
        phas_prev(:) = 0.d0
        phas_curr(:) = 0.d0
        call metaGetPhase(fami, '-', ipg, ispg, meta_type, &
                          nb_phasis, phas_prev)
        call metaGetPhase(fami, '+', ipg, ispg, meta_type, &
                          nb_phasis, phas_curr, zcold_=zcold_curr)
!
! ----- Get elastic parameters
!
        call get_elas_id(j_mater, elas_id, elas_keyword)
        call get_elas_para(fami, j_mater, '+', ipg, ispg, &
                           elas_id, elas_keyword, &
                           e_=young, nu_=nu)
        ASSERT(elas_id .eq. 1)
        deuxmu = young/(1.d0+nu)
!
! ----- Compute coefficients for second member
!
        call meta_vpta_coef(metaRela, metaGlob, &
                            lgpg, fami, ipg, j_mater, &
                            l_temp, temp, meta_type, nb_phasis, phas_prev, &
                            phas_curr, zcold_curr, young, deuxmu, coef, &
                            trans)
!
! ----- Compute geometric coefficient for axisymmetric
!
        if (l_axi) then
            poids = poids*r
            co_axis = 1.d0/r
        else
            co_axis = 0.d0
        end if
!
! ----- Compute stresses
!
        sigmo = 0.d0
        do i = 1, 3
            sigmo = sigmo+zr(j_sigm+(ipg-1)*nb_sigm+i-1)
        end do
        sigmo = sigmo/3.d0
!
        do i = 1, nb_sigm
            sigdv(i) = zr(j_sigm+(ipg-1)*nb_sigm+i-1)-sigmo*kron(i)
            sig(i) = coef*(1.5d0*trans*sigdv(i))
            sig(i) = deuxmu*sig(i)
        end do
!
! ----- Second member
!
        do i = 1, nno
            zr(j_vectu+2*i-2) = zr(j_vectu+2*i-2)+ &
                                poids*(sig(1)*dfdx(i)+sig(3)*zr(j_vf+k+i-1)*co_axis+sig(4)*dfdy(i))
            zr(j_vectu+2*i-1) = zr(j_vectu+2*i-1)+ &
                                poids*(sig(2)*dfdy(i)+sig(4)*dfdx(i))
        end do
    end do
!
99  continue
end subroutine
