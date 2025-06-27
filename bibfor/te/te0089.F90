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
subroutine te0089(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/vff2dn.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/Behaviour_type.h"
#include "asterc/r8prem.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: AXIS_FLUI_STRU, 2D_FLUI_STRU
!
! Option: RIGI_MECA/FORC_NODA/FULL_MECA/RAPH_MECA/RIGI_MECA_HYST/RIGI_MECA_TANG
! RIGI_MECA_HYST actuellement pas possible car matrice complexe non symétrique
! n'est pas prevue
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: mmat(9, 9)
    real(kind=8) :: nx, ny, norm(2), ul(9)
    real(kind=8) :: poids, celer
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: jv_geom, jv_mate, jv_matr
    integer(kind=8) :: jvVect, jvDisp, jv_codret, jvDispm, jvDispp
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: nno, npg, ndim
    integer(kind=8) :: ij, i
    integer(kind=8) :: n1, n2
    integer(kind=8) :: ino1, ino2, ipg, ind1, ind2, jdim
    integer(kind=8) :: ldec
    integer(kind=8) :: j_mater, iret, codret
    character(len=16) :: FEForm, rela_comp
    aster_logical :: l_axis
    real(kind=8) :: r
    aster_logical :: lVect, lMatr, lVari, lSigm
!
! --------------------------------------------------------------------------------------------------
!
    lVect = ASTER_FALSE
    lMatr = ASTER_FALSE
    lVari = ASTER_FALSE
    lSigm = ASTER_FALSE
!
! - Check behaviour
!
    if (option(1:9) .eq. 'FULL_MECA' .or. &
        option .eq. 'RAPH_MECA' .or. &
        option .eq. 'RIGI_MECA_TANG') then
        call jevech('PCOMPOR', 'L', vk16=compor)
! ----- Select objects to construct from option name
        call behaviourOption(option, compor, &
                             lMatr, lVect, &
                             lVari, lSigm, &
                             codret)
        rela_comp = compor(RELA_NAME)
        if (rela_comp .ne. 'ELAS') then
            call utmess('F', 'FLUID1_1')
        end if
    end if

    mmat = 0.d0
!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
!
! - Get element parameters
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    l_axis = (lteatt('AXIS', 'OUI'))
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. 3)
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, cele_r_=celer)
!
! - Loop on Gauss points
!
    do ipg = 1, npg
        ldec = (ipg-1)*nno
        call vff2dn(ndim, nno, ipg, ipoids, idfde, &
                    zr(jv_geom), nx, ny, poids)
        norm(1) = nx
        norm(2) = ny
        if (l_axis) then
            r = 0.d0
            do ino1 = 1, nno
                r = r+zr(jv_geom+2*(ino1-1))*zr(ivf+ldec+ino1-1)
            end do
            poids = poids*r
        end if
        if (FEForm .eq. 'U_P') then
            do ino1 = 1, nno
                do ino2 = 1, nno
                    do jdim = 1, 2
                        ind1 = 3*(ino1-1)+jdim
                        ind2 = 3*(ino2-1)+3
                        if (celer .le. r8prem()) then
                            mmat(ind1, ind2) = 0.d0
                        else
                            mmat(ind1, ind2) = mmat(ind1, ind2)- &
                                               poids*norm(jdim)* &
                                               zr(ivf+ldec+ino1-1)*zr(ivf+ldec+ino2-1)
                        end if
                    end do
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end do
!
! - Output field
!
    if (option(1:9) .eq. 'RIGI_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        if (FEForm .eq. 'U_P') then
            call jevech('PMATUNS', 'E', jv_matr)
            do ino2 = 1, 3*nno
                do ino1 = 1, 3*nno
                    ij = ino2+3*nno*(ino1-1)
                    zr(jv_matr+ij-1) = mmat(ino1, ino2)
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if
!
! - Save vector
!
    if (lVect .or. option .eq. 'FORC_NODA') then
        call jevech('PVECTUR', 'E', jvVect)
        if (FEForm .eq. 'U_P') then
            if (lVect) then
                call jevech('PDEPLMR', 'L', jvDispm)
                call jevech('PDEPLPR', 'L', jvDispp)
                do i = 1, 3*nno
                    zr(jvVect+i-1) = 0.d0
                    ul(i) = zr(jvDispm+i-1)+zr(jvDispp+i-1)
                end do
            elseif (option .eq. "FORC_NODA") then
                call jevech('PDEPLAR', 'L', jvDisp)
                do i = 1, 3*nno
                    zr(jvVect+i-1) = 0.d0
                    ul(i) = zr(jvDisp+i-1)
                end do
            else
                ASSERT(ASTER_FALSE)
            end if

            do n1 = 1, 3*nno
                do n2 = 1, 3*nno
                    zr(jvVect+n1-1) = zr(jvVect+n1-1)+mmat(n1, n2)*ul(n2)
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if

! - Save return code
    if (lSigm) then
        if (FEForm .eq. 'U_P') then
            call jevech('PCODRET', 'E', jv_codret)
            zi(jv_codret) = 0
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end if
!
end subroutine
