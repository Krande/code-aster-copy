! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine te0100(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmdlog.h"
#include "asterfort/nmgpfi.h"
#include "asterfort/nmgrla.h"
#include "asterfort/nmplxd.h"
#include "asterfort/nmtstm.h"
#include "asterfort/rcangm.h"
#include "asterfort/tecach.h"
#include "asterfort/tgveri.h"
#include "asterfort/Behaviour_type.h"
#include "blas/dcopy.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: D_PLAN, C_PLAN, AXIS
!           AXIS_SI, C_PLAN_SI (QUAD8), D_PLAN_SI (QUAD8)
!           AXIS_ELDI, PLAN_ELDI
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
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
!
    character(len=8) :: typmod(2)
    character(len=4) :: fami
    integer, parameter :: sz_tens = 4, ndim = 2
    integer :: nno, npg, imatuu, lgpg, iret
    integer :: igeom, imate
    integer :: icontm, ivarim
    integer :: iinstm, iinstp, ideplm, ideplp, icompo, icarcr
    integer :: ivectu, icontp, ivarip
    integer :: ivarix, jv_mult_comp
    integer :: jtab(7)
    real(kind=8) :: angl_naut(7)
    aster_logical :: matsym
    character(len=16) :: mult_comp, defo_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer :: codret
    integer :: jv_codret
!     POUR TGVERI
    real(kind=8) :: sdepl(3*9), svect(3*9), scont(6*9), smatr(3*9*3*9)
    real(kind=8) :: epsilo
    real(kind=8) :: varia(2*3*9*3*9)
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    ivarip = 1
    imatuu = 1
    ivectu = 1
    ivarix = 1
    jv_codret = 1
    fami = 'RIGI'
    codret = 0
!
! - Type of finite element
!
    typmod = " "
    if (lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (lteatt('C_PLAN', 'OUI')) then
        typmod(1) = 'C_PLAN'
    else if (lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - Get input fields
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', icompo)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PMULCOM', 'L', jv_mult_comp)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
!
! - Properties of behaviour
!
    mult_comp = zk16(jv_mult_comp-1+1)
    defo_comp = zk16(icompo-1+DEFO)
!
    call FECell%init()
    ! if (defo_comp == "PETIT_REAC") then
    !     call FECell%updateCoordinates(dispPrev+dispIncr)
    ! end if
    call FEQuad%initCell(FECell, fami)
    call FEBasis%initCell(FECell)
!
    nno = FECell%nbnodes
    ASSERT(nno .le. 9)
    npg = FEQuad%nbQuadPoints
!
! - Get orientation
!
    call rcangm(ndim, FECell%barycenter(), angl_naut)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, zk16(icompo), &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)
!
! - Get output fields
!
    if (lMatr) then
        call nmtstm(zr(icarcr), imatuu, matsym)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        call dcopy(npg*lgpg, zr(ivarix), 1, zr(ivarip), 1)
    end if
    if (option .eq. 'RIGI_MECA_IMPLEX') then
        call jevech('PCONTXR', 'E', icontp)
        call dcopy(npg*sz_tens, zr(icontm), 1, zr(icontp), 1)
    end if

500 continue

    if (defo_comp .eq. 'PETIT') then

        call nmplxd(FECell, FEBasis, FEQuad, nno, npg, ndim, &
                    typmod, option, zi(imate), &
                    zk16(icompo), mult_comp, lgpg, zr(icarcr), &
                    zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), &
                    angl_naut, zr(icontm), zr(ivarim), &
                    matsym, zr(icontp), zr(ivarip), &
                    zr(imatuu), zr(ivectu), codret)
        if (codret .ne. 0) goto 999

    else if (defo_comp .eq. 'PETIT_REAC') then
        call nmplxd(FECell, FEBasis, FEQuad, nno, npg, ndim, &
                    typmod, option, zi(imate), &
                    zk16(icompo), mult_comp, lgpg, zr(icarcr), &
                    zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), &
                    angl_naut, zr(icontm), zr(ivarim), &
                    matsym, zr(icontp), zr(ivarip), &
                    zr(imatuu), zr(ivectu), codret)
        if (codret .ne. 0) goto 999

    else if (defo_comp .eq. 'SIMO_MIEHE') then
        call nmgpfi(fami, option, typmod, ndim, nno, &
                    npg, zr(igeom), &
                    zk16(icompo), zi(imate), mult_comp, lgpg, zr(icarcr), &
                    angl_naut, zr(iinstm), zr(iinstp), zr(ideplm), zr(ideplp), &
                    zr(icontm), zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), &
                    zr(imatuu), codret)
        if (codret .ne. 0) goto 999

    else if (defo_comp .eq. 'GREEN_LAGRANGE') then
        call nmgrla(FECell, FEBasis, FEQuad, &
                    option, typmod, zi(imate), &
                    ndim, nno, npg, lgpg, &
                    zk16(icompo), zr(icarcr), mult_comp, &
                    zr(iinstm), zr(iinstp), &
                    zr(ideplm), &
                    zr(ideplp), angl_naut, &
                    zr(icontm), zr(icontp), &
                    zr(ivarim), zr(ivarip), &
                    matsym, zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999

    else if (defo_comp .eq. 'GDEF_LOG') then
        call nmdlog(FECell, FEBasis, FEQuad, option, typmod, ndim, nno, npg, &
                    zk16(icompo), mult_comp, zi(imate), lgpg, &
                    zr(icarcr), angl_naut, zr(iinstm), zr(iinstp), matsym, &
                    zr(ideplm), zr(ideplp), zr(icontm), zr(ivarim), zr(icontp), &
                    zr(ivarip), zr(ivectu), zr(imatuu), codret)
        if (codret .ne. 0) goto 999

    else
        ASSERT(ASTER_FALSE)
    end if

! ----- Calcul eventuel de la matrice TGTE par PERTURBATION

    call tgveri(option, zr(icarcr), zk16(icompo), nno, zr(igeom), &
                ndim, ndim*nno, zr(ideplp), sdepl, zr(ivectu), &
                svect, sz_tens*npg, zr(icontp), scont, npg*lgpg, &
                zr(ivarip), zr(ivarix), zr(imatuu), smatr, matsym, &
                epsilo, varia, iret)
    if (iret .ne. 0) goto 500

999 continue
!
! - Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', jv_codret)
        zi(jv_codret) = codret
    end if
!
end subroutine
