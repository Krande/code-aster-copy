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

subroutine te0333(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcgr.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/epsvmc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/granvi.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/Behaviour_type.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: EPSP_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: mxcmel = 162
    integer(kind=8), parameter :: nbsgm = 6
    real(kind=8) :: epsi_meca(mxcmel), epsi_plas(mxcmel)
    real(kind=8) :: sigma(nbsgm)
    real(kind=8) :: epsi_creep(nbsgm)
    integer(kind=8) :: i, ndim, nno, nbsig, idsig
    integer(kind=8) :: npg, ipoids, ivf, idfde, igau, isig, igeom, idepl, itemps, imate
    integer(kind=8) :: idefp, nbvari, ivari, nvi, nvif, ibid, jtab(7), iret, ibid2
    real(kind=8) :: c1, c2, trsig
    real(kind=8) :: angl_naut(3), nharm, e, nu, zero, un, tempg, time
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: mod3d
    integer(kind=8) :: elas_id
    character(len=16) :: optio2, kit_comp_1, kit_comp_2, rela_comp, elas_keyword
    aster_logical :: l_creep, l_temp
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    un = 1.d0
    nharm = zero
    mod3d = '3D'
!
! - Finite element informations
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Number of stress components
!
    nbsig = nbsigm()
    ASSERT(nbsig .eq. nbsgm)
!
! - Geometry
!
    call jevech('PGEOMER', 'L', igeom)
!
! - Material parameters
!
    call jevech('PMATERC', 'L', imate)
!
! - Orthotropic parameters
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! - Current time
!
    call jevech('PINSTR', 'L', itemps)
    time = zr(itemps)
!
! - Current displacements (nodes)
!
    call jevech('PDEPLAR', 'L', idepl)
!
! - Current stresses (gauss points)
!
    call jevech('PCONTRR', 'L', idsig)
!
! - Comportment
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_comp = compor(RELA_NAME)
    kit_comp_1 = compor(CREEP_NAME)
    kit_comp_2 = compor(PLAS_NAME)
!
! - Internal variables
!
    call jevech('PVARIGR', 'L', ivari)
    call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, &
                itab=jtab)
    nbvari = max(jtab(6), 1)*jtab(7)
!
! - Elasticity: only isotropic !
!
    call get_elas_id(zi(imate), elas_id, elas_keyword)
    if (elas_id .ne. 1) then
        call utmess('F', 'ELEMENTS6_2')
    end if
!
! - Compute mechanical strains epsi_meca = epsi_tota - epsi_varc
! -- Command variables strains: epsi_varc (contains thermics, drying, ...)
! -- Total strains: epsi_tota
!
    optio2 = 'EPME_ELGA'
    call epsvmc('RIGI', nno, ndim, nbsig, npg, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                time, angl_naut, nharm, optio2, epsi_meca)
!
! - Creep strains: epsi_creep
!
    if (rela_comp(1:13) .ne. 'BETON_GRANGER' .and. &
        (rela_comp(1:7) .ne. 'KIT_DDI' .or. kit_comp_1(1:13) .ne. 'BETON_GRANGER')) then
        l_creep = .false.
        do i = 1, mxcmel
            epsi_plas(i) = zero
        end do
        do i = 1, nbsig
            epsi_creep(i) = zero
        end do
    else
        call granvi(mod3d, ibid, ibid2, nvif)
        l_creep = .true.
    end if
!
! - Materials parameters depend on temperature ?
!
    l_temp = .false.
    if (rela_comp(1:15) .eq. 'BETON_DOUBLE_DP') then
        nvi = 3
        l_temp = .true.
    else if (rela_comp(1:7) .eq. 'KIT_DDI') then
        if (kit_comp_2(1:15) .eq. 'BETON_DOUBLE_DP') then
            if (kit_comp_1(1:13) .eq. 'BETON_GRANGER') then
                nvi = nvif+3
                l_temp = .true.
            else
                call utmess('F', 'ELEMENTS3_76')
            end if
        end if
    end if
!
! - Loop on Gauss points
!
    do igau = 1, npg
!
! ----- Get current temperature
!
        call rcvarc(' ', 'TEMP', '+', 'RIGI', igau, &
                    1, tempg, iret)
!
! ----- Change temperature from internal variable (maximum) for BETON_DOUBLE_DP/BETON_GRANGER
!
        if (l_temp) then
            if (tempg .lt. zr(ivari+(igau-1)*nbvari+nvi-1)) then
                tempg = zr(ivari+(igau-1)*nbvari+nvi-1)
            end if
        end if
!
! ----- Get elastic parameters (only isotropic elasticity)
!
        call get_elas_id(zi(imate), elas_id, elas_keyword)
        call get_elas_para('RIGI', zi(imate), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           time=time, temp=tempg, e_=e, nu_=nu)
        ASSERT(elas_id .eq. 1)
!
! ----- Compute creep strains (current Gauss point)
!
        if (l_creep) then
            call calcgr(igau, nbsig, nbvari, zr(ivari), nu, &
                        epsi_creep)
        end if
!
! ----- Compute stresses (current Gauss point)
!
        do i = 1, nbsig
            sigma(i) = zr(idsig+(igau-1)*nbsig+i-1)
        end do
        trsig = sigma(1)+sigma(2)+sigma(3)
!
! ----- Compute plastic strains (current Gauss point) epsi_plas = epsi_tota - epsi_elas - epsi_creep
! -- Creep strains: epsi_creep
! -- Elastic strains: epsi_elas
!
        c1 = (un+nu)/e
        c2 = nu/e
        epsi_plas(nbsig*(igau-1)+1) = epsi_meca(nbsig*(igau-1)+1)-(c1*sigma(1)-c2*trsig)- &
                                      epsi_creep(1)
        epsi_plas(nbsig*(igau-1)+2) = epsi_meca(nbsig*(igau-1)+2)-(c1*sigma(2)-c2*trsig)- &
                                      epsi_creep(2)
        epsi_plas(nbsig*(igau-1)+3) = epsi_meca(nbsig*(igau-1)+3)-(c1*sigma(3)-c2*trsig)- &
                                      epsi_creep(3)
        epsi_plas(nbsig*(igau-1)+4) = epsi_meca(nbsig*(igau-1)+4)-c1*sigma(4)-epsi_creep(4)
        epsi_plas(nbsig*(igau-1)+5) = epsi_meca(nbsig*(igau-1)+5)-c1*sigma(5)-epsi_creep(5)
        epsi_plas(nbsig*(igau-1)+6) = epsi_meca(nbsig*(igau-1)+6)-c1*sigma(6)-epsi_creep(6)
    end do
!
! - Plastic strain output
!
    call jevech('PDEFOPG', 'E', idefp)
    do igau = 1, npg
        do isig = 1, nbsig
            zr(idefp+nbsig*(igau-1)+isig-1) = epsi_plas(nbsig*(igau-1)+isig)
        end do
    end do
!
end subroutine
