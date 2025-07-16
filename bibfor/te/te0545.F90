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
subroutine te0545(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefv.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ngfint.h"
#include "asterfort/nglgic.h"
#include "asterfort/ngvlog.h"
#include "asterfort/nmgvmb.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_GRAD_INCO, 3D_GRAD_VARI
!           D_PLAN_GRAD_INCO, D_PLAN_GRAD_VARI
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: typmod(2)
    character(len=16) :: defo_comp
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: axi, matsym
    integer(kind=8) :: nnoQ, nnoL, npg, ndim, nddl, neps, lgpg
    integer(kind=8) :: jv_poids, jv_vfQ, jv_dfdeQ, jv_vfL, jv_dfdeL
    integer(kind=8) :: imate, icontm, ivarim, iinstm, iinstp, ideplm, ideplp
    integer(kind=8) :: ivectu, icontp, ivarip, imatns, imatuu, icarcr, ivarix, igeom, icoret
    integer(kind=8) :: iret, nnos, jv_ganoQ, jv_ganoL, itab(7)
    integer(kind=8) :: codret
    real(kind=8) :: angmas(3)
    real(kind=8), allocatable :: b(:, :, :), w(:, :), ni2ldc(:, :)
    aster_logical :: lMatr, lVect, lSigm, lVari
!
! --------------------------------------------------------------------------------------------------
!
    matsym = ASTER_FALSE
    ivectu = 1
    icontp = 1
    ivarip = 1
    icoret = 1
    imatns = 1
!
! - Type of modelling
!
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))
    axi = typmod(1) .eq. 'AXIS'
!
! - Get parameters of element
!
    call elrefv('RIGI', ndim, nnoL, nnoQ, nnos, &
                npg, jv_poids, jv_vfL, jv_vfQ, jv_dfdeL, &
                jv_dfdeQ, jv_ganoL, jv_ganoQ)
!
! - PARAMETRES EN ENTREE ET DIMENSION
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call tecach('OOO', 'PDEPLPR', 'L', iret, nval=2, &
                itab=itab)
    nddl = itab(2)
!
! - Behaviour
!
    defo_comp = compor(DEFO)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, lMatr, lVect, lVari, &
                         lSigm, codret)
!
!    NOMBRE DE VARIABLES INTERNES
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=itab)
    lgpg = max(itab(6), 1)*itab(7)
!
! - PARAMETRES EN SORTIE
!
    if (lMatr) then
        matsym = .false.
        call jevech('PMATUNS', 'E', imatns)
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', icoret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        zr(ivarip:ivarip-1+npg*lgpg) = zr(ivarix:ivarix-1+npg*lgpg)
    end if
!
!    ORIENTATION DU MASSIF
    call getElemOrientation(ndim, nnoQ, igeom, angmas)
!
! - CALCUL DES FORCES INTERIEURES ET MATRICES TANGENTES
!
    if (defo_comp .eq. 'GDEF_LOG') then
        if (lteatt('INCO', 'C5GV')) then
            call nglgic('RIGI', option, typmod, ndim, nnoQ, &
                        nnoL, npg, nddl, jv_poids, zr(jv_vfQ), &
                        zr(jv_vfL), jv_dfdeQ, jv_dfdeL, zr(igeom), compor, &
                        zi(imate), lgpg, zr(icarcr), angmas, zr(iinstm), &
                        zr(iinstp), matsym, zr(ideplm), zr(ideplp), zr(icontm), &
                        zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatns), &
                        lMatr, lVect, lSigm, lVari, codret)
        else
            call ngvlog('RIGI', option, typmod, ndim, nnoQ, &
                        nnoL, npg, nddl, jv_poids, zr(jv_vfQ), &
                        zr(jv_vfL), jv_dfdeQ, jv_dfdeL, zr(igeom), compor, &
                        zi(imate), lgpg, zr(icarcr), angmas, zr(iinstm), &
                        zr(iinstp), matsym, zr(ideplm), zr(ideplp), zr(icontm), &
                        zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatns), &
                        lMatr, lVect, lSigm, lVari, codret)
        end if
    else if (defo_comp(1:5) .eq. 'PETIT') then
        call nmgvmb(ndim, nnoQ, nnoL, npg, axi, &
                    zr(igeom), zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                    jv_poids, nddl, neps, b, w, &
                    ni2ldc)
        matsym = .not. (nint(zr(icarcr-1+CARCRI_MATRSYME)) .gt. 0)
        call ngfint(option, typmod, ndim, nddl, neps, &
                    npg, w, b, compor, 'RIGI', &
                    zi(imate), angmas, lgpg, zr(icarcr), zr(iinstm), &
                    zr(iinstp), zr(ideplm), zr(ideplp), ni2ldc, zr(icontm), &
                    zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), matsym, &
                    zr(imatuu), zr(imatns), lMatr, lVect, lSigm, &
                    codret)
        deallocate (b)
        deallocate (w)
        deallocate (ni2ldc)
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (lSigm) then
        zi(icoret) = codret
    end if
!
end subroutine
