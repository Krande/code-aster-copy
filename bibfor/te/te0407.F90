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
subroutine te0407(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmas3d.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_SI (for HEXA8)
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
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: mult_comp, defo_comp, rela_comp, type_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    character(len=4), parameter :: fami = 'RIGI'
    integer(kind=8) :: nno, npg, i, imatuu, lgpg, ndim
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: icontm, ivarim, jv_mult_comp
    integer(kind=8) :: iinstm, iinstp, ideplm, ideplp, icarcr
    integer(kind=8) :: ivectu, icontp, ivarip
    integer(kind=8) :: ivarix, iret
    integer(kind=8) :: jtab(7), jcret, codret
    real(kind=8) :: def(6, 3, 8), dfdi(8, 3)
    real(kind=8) :: angl_naut(3)
    character(len=8), parameter :: typmod(2) = (/'3D', '  '/)
    type(Behaviour_Integ) :: BEHinteg
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    imatuu = 1
    ivectu = 1
    icontp = 1
    ivarip = 1
    codret = 0

! - Get element parameters
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .eq. 8)

! - Get input fields
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PMULCOM', 'L', jv_mult_comp)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
! - Get orientation
    call getElemOrientation(ndim, nno, igeom, angl_naut)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    mult_comp = zk16(jv_mult_comp-1+1)
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    type_comp = compor(INCRELAS)

! - Get output fields
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        call jevech('PVARIMP', 'L', ivarix)
        b_n = to_blas_int(npg*lgpg)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
    end if

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              fami, zi(imate), &
                              BEHinteg)

! - HYPER-ELASTICITE
    if (type_comp .eq. 'COMP_ELAS') then
        call utmess('F', 'ELEMENTSSI_3')
    end if

! - HYPO-ELASTICITE
    if (defo_comp(6:10) .eq. '_REAC') then
        do i = 1, ndim*nno
            zr(igeom+i-1) = zr(igeom+i-1)+zr(ideplm+i-1)+zr(ideplp+i-1)
        end do
    end if
!
    if (defo_comp(1:5) .eq. 'PETIT') then
        call nmas3d(BEHinteg, &
                    fami, nno, npg, ipoids, ivf, &
                    idfde, zr(igeom), typmod, option, zi(imate), &
                    compor, mult_comp, lgpg, zr(icarcr), zr(iinstm), &
                    zr(iinstp), zr(ideplm), zr(ideplp), angl_naut, zr(icontm), &
                    zr(ivarim), dfdi, def, zr(icontp), zr(ivarip), &
                    zr(imatuu), zr(ivectu), codret)
    else
        call utmess('F', 'ELEMENTSSI_1', sk=defo_comp)
    end if

! - Save return code
    if (lSigm) then
        zi(jcret) = codret
    end if
!
end subroutine
