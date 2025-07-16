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
subroutine te0139(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
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
#include "asterfort/tgveri_use.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "jeveux.h"
!
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
!           3D_SI (HEXA20)
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
    integer(kind=8) :: sz_tens, ndim
    integer(kind=8) :: nno, npg, imatuu, lgpg, iret
    integer(kind=8) :: igeom, imate, iuse
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: iinstm, iinstp, ideplm, ideplp, icarcr
    integer(kind=8) :: ivectu, icontp, ivarip
    integer(kind=8) :: ivarix
    integer(kind=8) :: jtab(7)
    real(kind=8) :: angl_naut(7)
    aster_logical :: matsym
    character(len=16), pointer :: compor(:) => null(), v_mult_comp(:) => null()
    character(len=16) :: mult_comp, defo_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: codret
    integer(kind=8) :: jv_codret
!     POUR TGVERI
    real(kind=8) :: sdepl(3*27), svect(3*27), scont(6*27)
    real(kind=8) :: epsilo, disp_curr(MAX_BV)
    real(kind=8), pointer :: varia(:) => null(), smatr(:) => null()
    blas_int :: b_incx, b_incy, b_n
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
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PMULCOM', 'L', vk16=v_mult_comp)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
!
! - Properties of behaviour
!
    mult_comp = v_mult_comp(1)
    defo_comp = compor(DEFO)
!
    call FECell%init()
    nno = FECell%nbnodes
    ASSERT(nno .le. 27)
    ndim = FECell%ndim
    sz_tens = 2*ndim
!
    call tgveri_use(option, zr(icarcr), compor, iuse)
    if (iuse == 1) then
        allocate (varia(2*3*27*3*27))
        allocate (smatr(3*27*3*27))
    end if
!
    if (defo_comp == "PETIT_REAC") then
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ideplm), b_incx, disp_curr, b_incy)
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, zr(ideplp), b_incx, disp_curr, &
                   b_incy)
        call FECell%updateCoordinates(disp_curr)
    end if
!
    call FEQuad%initCell(FECell, fami)
    npg = FEQuad%nbQuadPoints
!
    call FEBasis%initCell(FECell)
!
!
! - Type of finite element
!
    if (ndim == 3) then
        typmod(1) = '3D'
    else
        if (lteatt('AXIS', 'OUI')) then
            typmod(1) = 'AXIS'
        else if (lteatt('C_PLAN', 'OUI')) then
            typmod(1) = 'C_PLAN'
        else if (lteatt('D_PLAN', 'OUI')) then
            typmod(1) = 'D_PLAN'
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
    typmod(2) = ' '
!
!
! - Get orientation
!
    call rcangm(ndim, FECell%barycenter(), angl_naut)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, lMatr, lVect, lVari, &
                         lSigm, codret)
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
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if
    if (option .eq. 'RIGI_MECA_IMPLEX') then
        call jevech('PCONTXR', 'E', icontp)
        b_n = to_blas_int(npg*sz_tens)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(icontm), b_incx, zr(icontp), b_incy)
    end if
!
500 continue
!
    if (defo_comp .eq. 'PETIT') then
        call nmplxd(FECell, FEBasis, FEQuad, nno, npg, &
                    ndim, typmod, option, zi(imate), compor, &
                    mult_comp, lgpg, zr(icarcr), zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), angl_naut, zr(icontm), zr(ivarim), &
                    matsym, zr(icontp), zr(ivarip), zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999
!
    else if (defo_comp .eq. 'PETIT_REAC') then
        call nmplxd(FECell, FEBasis, FEQuad, nno, npg, &
                    ndim, typmod, option, zi(imate), compor, &
                    mult_comp, lgpg, zr(icarcr), zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), angl_naut, zr(icontm), zr(ivarim), &
                    matsym, zr(icontp), zr(ivarip), zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999
!
    else if (defo_comp .eq. 'SIMO_MIEHE') then
        call nmgpfi(fami, option, typmod, ndim, nno, &
                    npg, zr(igeom), compor, zi(imate), mult_comp, &
                    lgpg, zr(icarcr), angl_naut, zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), zr(icontm), zr(ivarim), zr(icontp), &
                    zr(ivarip), zr(ivectu), zr(imatuu), codret)
        if (codret .ne. 0) goto 999
!
    else if (defo_comp .eq. 'GREEN_LAGRANGE') then
        call nmgrla(FECell, FEBasis, FEQuad, option, typmod, &
                    zi(imate), ndim, nno, npg, lgpg, &
                    compor, zr(icarcr), mult_comp, zr(iinstm), zr(iinstp), &
                    zr(ideplm), zr(ideplp), angl_naut, zr(icontm), zr(icontp), &
                    zr(ivarim), zr(ivarip), matsym, zr(imatuu), zr(ivectu), &
                    codret)
        if (codret .ne. 0) goto 999
!
    else if (defo_comp .eq. 'GDEF_LOG') then
        call nmdlog(FECell, FEBasis, FEQuad, option, typmod, &
                    ndim, nno, npg, compor, mult_comp, &
                    zi(imate), lgpg, zr(icarcr), angl_naut, zr(iinstm), &
                    zr(iinstp), matsym, zr(ideplm), zr(ideplp), zr(icontm), &
                    zr(ivarim), zr(icontp), zr(ivarip), zr(ivectu), zr(imatuu), &
                    codret)
        if (codret .ne. 0) goto 999
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ----- Calcul eventuel de la matrice TGTE par PERTURBATION
    call tgveri(option, zr(icarcr), compor, nno, zr(igeom), &
                ndim, ndim*nno, zr(ideplp), sdepl, zr(ivectu), &
                svect, sz_tens*npg, zr(icontp), scont, npg*lgpg, &
                zr(ivarip), zr(ivarix), zr(imatuu), smatr, matsym, &
                epsilo, varia, iret)
    if (iret .ne. 0) then
        goto 500
    end if
!
999 continue
!
! - Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', jv_codret)
        zi(jv_codret) = codret
    end if
!
! - Free large arrays
    if (iuse == 1) then
        deallocate (smatr)
        deallocate (varia)
    end if
!
end subroutine
