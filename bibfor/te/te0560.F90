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
!
subroutine te0560(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefv.h"
#include "asterfort/jevech.h"
#include "asterfort/massup.h"
#include "asterfort/nmgvno.h"
#include "asterfort/nmtstm.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D_GVNO
!           D_PLAN_GVNO
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!          MASS_MECA*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: fami
    aster_logical :: matsym
    integer(kind=8) :: nb_DOF
    integer(kind=8) :: nnoQ, npg, imatuu, lgpg, ndim
    integer(kind=8) :: jv_poids, jv_vfQ, jv_dfdeQ, jvGeom, jvMaterc
    integer(kind=8) :: nnoL, jv_vfL, jv_dfdeL, jv_ganoL
    integer(kind=8) :: icontm, ivarim
    integer(kind=8) :: jvInstmr, jvInstpr, ideplm, ideplp, jvCarcri
    integer(kind=8) :: ivectu, icontp, ivarip, nnos, jv_ganoQ
    integer(kind=8) :: ivarix, iret
    integer(kind=8) :: jtab(7), jcret, codret
    integer(kind=8) :: icodr1(1)
    character(len=8) :: typmod(2)
    character(len=16) :: defoComp
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: lVect, lMatr, lVari, lSigm
    blas_int :: b_incx, b_incy, b_n
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0

! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('S', 'TYPMOD2', typmod(2))

! - Get parameters of element
    fami = 'RIGI'
    if (option .eq. 'MASS_MECA') then
        fami = 'MASS'
    end if
    call elrefv(fami, ndim, nnoL, nnoQ, nnos, &
                npg, jv_poids, jv_vfL, jv_vfQ, jv_dfdeL, &
                jv_dfdeQ, jv_ganoL, jv_ganoQ)
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)

! - Input fields
    call jevech('PGEOMER', 'L', jvGeom)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nnoQ, jvGeom, materPara%lcsPara)
!
    if (option(1:9) .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', imatuu)
! ----- nb_DOF: displacements (2 or 3) + VARI
        nb_DOF = ndim+2

        call massup(zi(jvMaterc), &
                    option, ndim, nb_DOF, nnoQ, nnoL, &
                    npg, jv_poids, jv_dfdeQ, &
                    zr(jvGeom), zr(jv_vfQ), imatuu, icodr1, jvGeom, &
                    jv_vfQ)

    else
        if (option .eq. 'RIGI_MECA_ELAS' .or. &
            option .eq. 'FULL_MECA_ELAS' .or. &
            option .eq. 'RAPH_MECA') then
            fami = 'ELAS'
        else
            fami = 'RIGI'
        end if

! ----- Get input fields
        call jevech('PINSTMR', 'L', jvInstmr)
        call jevech('PINSTPR', 'L', jvInstpr)
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
        lgpg = max(jtab(6), 1)*jtab(7)

! ----- Get fields for behaviour
        call jevech('PCOMPOR', 'L', vk16=compor)
        call jevech('PCARCRI', 'L', jvCarcri)

! ----- Get parameters for behaviour
        defoComp = compor(DEFO)
        if (defoComp .ne. 'PETIT') then
            call utmess('F', 'ELEMENTS3_16', sk=defoComp)
        end if

! ----- Initialisation of behaviour datastructure
        call behaviourInit(BEHInteg)

! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  compor, zr(jvCarcri), &
                                  zr(jvInstmr), zr(jvInstpr), &
                                  materPara, BEHInteg)

! ----- PARAMETRES EN SORTIE
        ivectu = 1
        icontp = 1
        ivarip = 1
        if (lMatr) then
            call nmtstm(zr(jvCarcri), imatuu, matsym)
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
        if (lSigm) then
            call jevech('PCODRET', 'E', jcret)
!
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
!
        call nmgvno(BEHInteg, &
                    ndim, nnoQ, nnoL, npg, &
                    jv_poids, zr(jv_vfQ), zr(jv_vfL), jv_dfdeQ, jv_dfdeL, &
                    zr(jvGeom), typmod, option, compor, &
                    lgpg, zr(jvCarcri), zr(jvInstmr), zr(jvInstpr), zr(ideplm), &
                    zr(ideplp), zr(icontm), zr(ivarim), zr(icontp), &
                    zr(ivarip), zr(imatuu), zr(ivectu), codret)
!
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
end subroutine
