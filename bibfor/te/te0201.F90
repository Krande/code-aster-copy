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
subroutine te0201(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmfi2d.h"
#include "asterfort/tecach.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: PLAN_JOINT
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
    integer(kind=8), parameter :: ndim = 2
    character(len=4), parameter :: fami = "RIGI"
    integer(kind=8) :: kk, i, j, npg
    integer(kind=8) :: igeom, imater, icarcr, idepm, iddep, icoret
    integer(kind=8) :: icontm, icontp, ivect, imatr
    integer(kind=8) :: ivarim, ivarip, jtab(7), iret, iinstm, iinstp
    integer(kind=8) :: lgpg, codret
    real(kind=8) :: mat(8, 8), fint(8), sigmo(6, 2), sigma(6, 2)
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: typmod(2)
    character(len=16) :: rela_comp
    aster_logical :: matsym
    aster_logical :: lVect, lMatr, lVari, lSigm
    blas_int :: b_incx, b_incy, b_n
    type(Behaviour_Integ) :: BEHinteg
!
! --------------------------------------------------------------------------------------------------
!
    npg = 2
    ivarip = 1
    icoret = 1
    icontp = 1
    ivect = 1
    icoret = 1

! - Type of finite element
    if (lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else
        typmod(1) = 'PLAN'
    end if
    typmod(2) = 'ELEMJOIN'

! - Get input fields
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imater)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PDEPLMR', 'L', idepm)
    call jevech('PDEPLPR', 'L', iddep)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              fami, zi(imater), &
                              BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    rela_comp = compor(RELA_NAME)

! - Get output fields
    if (lMatr) then
        matsym = .true.
        if (rela_comp .eq. 'JOINT_MECA_RUPT') matsym = .false.
        if (rela_comp .eq. 'JOINT_MECA_FROT') matsym = .false.
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivect)
    end if
!
!     CONTRAINTE -, RANGEE DANS UN TABLEAU (6,NPG)
    sigmo = 0.d0
    sigmo(1, 1) = zr(icontm)
    sigmo(2, 1) = zr(icontm+1)
    sigmo(1, 2) = zr(icontm+2)
    sigmo(2, 2) = zr(icontm+3)

! CALCUL DES CONTRAINTES, VIP, FORCES INTERNES ET MATR TANG ELEMENTAIRES
    call nmfi2d(BEHInteg, &
                npg, lgpg, zi(imater), option, zr(igeom), &
                zr(idepm), zr(iddep), sigmo, sigma, fint, &
                mat, zr(ivarim), zr(ivarip), zr(iinstm), zr(iinstp), &
                zr(icarcr), compor, typmod, lMatr, lVect, lSigm, &
                codret)

! - Save matrix
    if (lMatr) then
        if (matsym) then
            call jevech('PMATUUR', 'E', imatr)
            kk = 0
            do i = 1, 8
                do j = 1, i
                    zr(imatr+kk) = mat(i, j)
                    kk = kk+1
                end do
            end do
        else
            call jevech('PMATUNS', 'E', imatr)
            kk = 0
            do i = 1, 8
                do j = 1, 8
                    zr(imatr+kk) = mat(i, j)
                    kk = kk+1
                end do
            end do
        end if
    end if

! - Save stresses
    if (lSigm) then
        zr(icontp) = sigma(1, 1)
        zr(icontp+1) = sigma(2, 1)
        zr(icontp+2) = sigma(1, 2)
        zr(icontp+3) = sigma(2, 2)
    end if

! - Save internal forces
    if (lVect) then
        b_n = to_blas_int(8)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, fint, b_incx, zr(ivect), b_incy)
    end if

! - Save return code
    if (lSigm) then
        call jevech('PCODRET', 'E', icoret)
        zi(icoret) = codret
    end if
!
end subroutine
