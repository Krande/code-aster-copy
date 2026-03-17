! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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

subroutine te0054(option, nomte)

    use Behaviour_module, only: behaviourOption
!
    implicit none

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/nmsfin.h"
#include "asterfort/nmsfon.h"

    character(len=16), intent(in) :: option, nomte

! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN_MIX_STA et 3D_MIX_STA
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA, FORC_NODA et REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: typmod(2), elrefe
    aster_logical :: axi
    integer(kind=8) :: nno, nnos, npg, ndim, lgpg, nddl, neps, i
    integer(kind=8) :: jv_poids, jv_vf, jv_dfde
    integer(kind=8) :: imate, icontm, ivarim, iinstm, iinstp, ideplm, ideplp, icompo
    integer(kind=8) :: ivectu, icontp, ivarip, imatuu, icarcr, ivarix, igeom, icoret
    integer(kind=8) :: icont
    integer(kind=8) :: iret, itab(7)
    integer(kind=8) :: codret
    real(kind=8) :: angmas(3), sigref, lagref
    real(kind=8), allocatable:: sref(:)
    aster_logical :: lMatr, lVect, lSigm, lVari, refe

!
! --------------------------------------------------------------------------------------------------
!
    ivectu = 1
    icontp = 1
    ivarip = 1
    icoret = 1
    imatuu = 1
!
! - Type of modelling
!
    call teattr('S', 'TYPMOD', typmod(1))
    typmod(2) = ' '
    axi = typmod(1) .eq. 'AXIS'
    refe = ASTER_FALSE
!
! - Get parameters of element
!
    call elref1(elrefe)
    call elrefe_info(elrefe=elrefe, fami='RIGI', &
                     ndim=ndim, nno=nno, &
                     npg=npg, jpoids=jv_poids, &
                     jvf=jv_vf, jdfde=jv_dfde)
    nddl = 3*nno*ndim

! - PARAMETRES EN ENTREE ET DIMENSION
!
    call jevech('PGEOMER', 'L', igeom)
    if (option .eq. "FORC_NODA") then
        call jevech('PCOMPOR', 'L', icompo)
        call jevech('PSIEFR', 'L', icont)
    else if (option .eq. "REFE_FORC_NODA") then
        allocate (sref(4*ndim))
    else
        call jevech('PMATERC', 'L', imate)
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PCOMPOR', 'L', icompo)
        call jevech('PCARCRI', 'L', icarcr)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
        call tecach('OOO', 'PDEPLPR', 'L', iret, nval=2, &
                    itab=itab)
!
!    NOMBRE DE VARIABLES INTERNES
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                    itab=itab)
        lgpg = max(itab(6), 1)*itab(7)

    end if

! - Select objects to construct from option name
!
    call behaviourOption(option, zk16(icompo), lMatr, lVect, lVari, &
                         lSigm, codret)

    if (lMatr) then
        call jevech('PMATUNS', 'E', imatuu)
    end if
    if ((lVect) .or. (option .eq. "REFE_FORC_NODA") .or. (option .eq. "FORC_NODA")) then
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
    call getElemOrientation(ndim, nno, igeom, angmas)

!    OPTION FORC_NODA
    if (option .eq. "FORC_NODA") then
        call nmsfon(refe, ndim, nno, npg, nddl, &
                    zr(igeom), zr(jv_vf), jv_dfde, &
                    jv_poids, zr(icont), zr(ivectu))
!    OPTION REFE_FORC_NODA
    else if (option .eq. "REFE_FORC_NODA") then
        refe = ASTER_TRUE
        call terefe('SIGM_REFE', 'MECA_MIXSTA', sigref)
        call terefe('LAGR_REFE', 'MECA_MIXSTA', lagref)
        sref(1:2*ndim) = sigref
        sref(2*ndim+1:4*ndim) = lagref

        call nmsfon(refe, ndim, nno, npg, nddl, &
                    zr(igeom), zr(jv_vf), jv_dfde, &
                    jv_poids, transpose(spread(sref, 1, npg)), zr(ivectu))
!    OPTIONS RAPH_MECA, FULL_MECA_*, RIGI_MECA_*
    else
        call nmsfin('RIGI', option, typmod, ndim, nno, &
                    npg, nddl, jv_poids, zr(jv_vf), jv_dfde, &
                    zr(igeom), zk16(icompo), &
                    zi(imate), lgpg, zr(icarcr), angmas, zr(iinstm), &
                    zr(iinstp), zr(ideplm), zr(ideplp), zr(icontm), &
                    zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatuu), &
                    lMatr, lVect, lSigm, lVari, codret)
    end if

    if (refe) then
        deallocate (sref)
    end if

    if (lSigm) then
        zi(icoret) = codret
    end if

end subroutine
