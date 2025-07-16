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
subroutine thmCompNonLin(option, ds_thm)
!
    use THM_type
    use Behaviour_module
    use Behaviour_type
!
    implicit none
!
#include "asterc/ismaem.h"
#include "asterf_types.h"
#include "asterfort/assthm.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/jevech.h"
#include "asterfort/thmGetElemPara.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! Non-linear options
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codret
    real(kind=8) :: angl_naut(3)
    integer(kind=8) :: jv_geom, jv_matr, jv_vect, jv_sigmm, jv_varim, jv_cret
    integer(kind=8) :: jv_mater, jv_instm, jv_instp, jv_dispm
    integer(kind=8) :: jv_dispp, jv_carcri, jv_varip, jv_sigmp
    aster_logical :: l_axi
    character(len=3) :: inte_type
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    integer(kind=8) :: dimdep, dimdef, dimcon, dimuel
    integer(kind=8) :: nddls, nddlm, nddl_p1, nddl_p2, nddl_meca, nddl_2nd
    integer(kind=8) :: ndim, nno, nnos
    integer(kind=8) :: npi, npg, nbvari
    integer(kind=8) :: jv_poids, jv_func, jv_dfunc, jv_poids2, jv_func2, jv_dfunc2, jv_gano
    character(len=8) :: typmod(2)
    integer(kind=8):: lg_vi, lg_sig
    real(kind=8), allocatable:: varip(:), sigp(:), deplp(:)
    aster_logical :: lVect, lMatr, lVari, lSigm, lMatrPred
    character(len=16) :: compor_copy(COMPOR_SIZE), rela_meca
    integer(kind=8) :: iCompor
    type(Behaviour_Integ) :: BEHinteg
    character(len=4), parameter :: fami = 'FPG1'
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0
    lMatrPred = option(1:9) .eq. 'RIGI_MECA'

! - Get all parameters for current element
    call thmGetElemPara(ds_thm, l_axi, &
                        typmod, inte_type, ndim, &
                        mecani, press1, press2, tempe, second, &
                        dimdep, dimdef, dimcon, dimuel, &
                        nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                        nno, nnos, &
                        npi, npg, &
                        jv_poids, jv_func, jv_dfunc, &
                        jv_poids2, jv_func2, jv_dfunc2, &
                        jv_gano)

! - Input fields
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mater)
    call jevech('PINSTMR', 'L', jv_instm)
    call jevech('PINSTPR', 'L', jv_instp)
    call jevech('PDEPLMR', 'L', jv_dispm)
    call jevech('PDEPLPR', 'L', jv_dispp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jv_carcri)
    call jevech('PVARIMR', 'L', jv_varim)
    call jevech('PCONTMR', 'L', jv_sigmm)

! - Make copy of COMPOR map
    do iCompor = 1, COMPOR_SIZE
        compor_copy(iCompor) = compor(iCompor)
    end do

! - Force DEFO_LDC="MECANIQUE" for THM
    rela_meca = compor_copy(MECA_NAME)

! - Something (very) strange with Hujeux => glute
    if (rela_meca .ne. "HUJEUX") then
        compor_copy(DEFO_LDC) = "MECANIQUE"
    end if

! - Select objects to construct from option name
    call behaviourOption(option, compor_copy, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Output fields
    if (lMatr) then
        call jevech('PMATUNS', 'E', jv_matr)
    else
        jv_matr = ismaem()
    end if

    if (lVect) then
        call jevech('PVECTUR', 'E', jv_vect)
    else
        jv_vect = ismaem()
    end if

! - Get frame orientation for anisotropy
    call getElemOrientation(ndim, nno, jv_geom, angl_naut)

! - Number of (total) internal variables
    read (compor_copy(NVAR), '(I16)') nbvari

! - Intermediate arrays to be safe when the addresses do not exist
    lg_sig = dimcon*npi
    allocate (sigp(lg_sig))
    if (lMatrPred) then
        sigp(1:lg_sig) = zr(jv_sigmm:jv_sigmm+lg_sig-1)
    else
        sigp(1:lg_sig) = 0
    end if

    lg_vi = nbvari*npi
    allocate (varip(lg_vi))
    if (lMatrPred) then
        varip(1:lg_vi) = zr(jv_varim:jv_varim+lg_vi-1)
    else
        varip(1:lg_vi) = 0
    end if

! - Prepare reference configuration
    allocate (deplp(dimuel))
    deplp(1:dimuel) = zr(jv_dispm:jv_dispm+dimuel-1)+zr(jv_dispp:jv_dispp+dimuel-1)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, zr(jv_carcri), &
                              zr(jv_instm), zr(jv_instp), &
                              fami, zi(jv_mater), &
                              BEHinteg)

! - Save
    ds_thm%ds_behaviour%BEHinteg = BEHinteg

! - Compute
    call assthm(ds_thm, option, zi(jv_mater), &
                lMatr, lSigm, lVect, &
                lVari, lMatrPred, l_axi, &
                typmod, inte_type, angl_naut, &
                ndim, nbvari, nno, nnos, &
                npg, npi, &
                nddls, nddlm, nddl_meca, &
                nddl_p1, nddl_p2, nddl_2nd, &
                dimdef, dimcon, dimuel, &
                mecani, press1, press2, tempe, second, &
                compor_copy, zr(jv_carcri), &
                jv_poids, jv_poids2, &
                jv_func, jv_func2, &
                jv_dfunc, jv_dfunc2, &
                zr(jv_geom), zr(jv_dispm), deplp, &
                zr(jv_sigmm), sigp, &
                zr(jv_varim), varip, &
                zr(jv_instm), zr(jv_instp), &
                zr(jv_matr), zr(jv_vect), codret)

! - Copy fields if required
    if (lSigm) then
        call jevech('PCONTPR', 'E', jv_sigmp)
        zr(jv_sigmp:jv_sigmp+lg_sig-1) = sigp(1:lg_sig)

        call jevech('PCODRET', 'E', jv_cret)
        zi(jv_cret) = codret
    end if

    if (lVari) then
        call jevech('PVARIPR', 'E', jv_varip)
        zr(jv_varip:jv_varip+lg_vi-1) = varip(1:lg_vi)
    end if

! Memory management
    deallocate (deplp)
    deallocate (sigp)
    deallocate (varip)

!
end subroutine
