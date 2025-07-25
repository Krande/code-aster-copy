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

subroutine te0334(option, nomte)
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
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
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
! Elements: 2D
! Option: EPSP_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: mxcmel = 54
    integer(kind=8), parameter :: nbsgm = 4
    real(kind=8) :: epsi_meca(mxcmel), epsi_plas(mxcmel)
    real(kind=8) :: sigma(nbsgm)
    real(kind=8) :: epsi_creep(nbsgm)
    integer(kind=8) :: i, ndim, nno, nbsig, idsig
    integer(kind=8) :: npg, ipoids, ivf, idfde, igau, isig, igeom, idepl, idefp, itemps
    integer(kind=8) :: imate, nbvari, ivari, jtab(7), iret
    real(kind=8) :: c1, c2, trsig
    real(kind=8) :: angl_naut(3), nharm, e, nu, zero, un, time
    integer(kind=8) :: elas_id
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: optio2, kit_comp_2, rela_comp, elas_keyword
    aster_logical :: l_creep
!
! --------------------------------------------------------------------------------------------------
!
    zero = 0.d0
    un = 1.d0
    nharm = zero
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
! - Orthotropic parameters: cannot use => zero
!
    angl_naut(:) = 0.d0
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
    kit_comp_2 = compor(CREEP_NAME)
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
! - Stress plane warning
!
    if (lteatt('C_PLAN', 'OUI')) then
        if (rela_comp .ne. 'VMIS_ISOT_LINE' .and. &
            rela_comp(1:4) .ne. 'ELAS' .and. &
            rela_comp .ne. 'VMIS_ISOT_TRAC') then
            call utmess('A', 'ELEMENTS6_3', sk=rela_comp)
        end if
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
        (rela_comp(1:7) .ne. 'KIT_DDI' .or. kit_comp_2(1:13) .ne. 'BETON_GRANGER')) then
        l_creep = .false.
        do i = 1, mxcmel
            epsi_plas(i) = zero
        end do
        do i = 1, nbsig
            epsi_creep(i) = zero
        end do
    else
        l_creep = .true.
    end if
!
! - Loop on Gauss points
!
    do igau = 1, npg
!
! ----- Get elastic parameters (only isotropic elasticity)
!
        call get_elas_id(zi(imate), elas_id, elas_keyword)
        call get_elas_para('RIGI', zi(imate), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           time=time, e_=e, nu_=nu)
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
!
        if (lteatt('C_PLAN', 'OUI')) then
            trsig = sigma(1)+sigma(2)
        else
            trsig = sigma(1)+sigma(2)+sigma(3)
        end if
!
! ----- Compute plastic strains (current Gauss point) epsi_plas = epsi_tota - epsi_elas - epsi_creep
! -- Creep strains: epsi_creep
! -- Elastic strains: epsi_elas
!
        c1 = (un+nu)/e
        c2 = nu/e
        epsi_plas(nbsig*(igau-1)+1) = epsi_meca(nbsig*(igau-1)+1)- &
                                      (c1*sigma(1)-c2*trsig)- &
                                      epsi_creep(1)
        epsi_plas(nbsig*(igau-1)+2) = epsi_meca(nbsig*(igau-1)+2)- &
                                      (c1*sigma(2)-c2*trsig)- &
                                      epsi_creep(2)
        if (lteatt('C_PLAN', 'OUI')) then
            epsi_plas(nbsig*(igau-1)+3) = -(epsi_plas(nbsig*(igau-1)+1)+ &
                                            epsi_plas(nbsig*(igau-1)+2))
        else
            epsi_plas(nbsig*(igau-1)+3) = epsi_meca(nbsig*(igau-1)+3)- &
                                          (c1*sigma(3)-c2*trsig)- &
                                          epsi_creep(3)
        end if
        epsi_plas(nbsig*(igau-1)+4) = epsi_meca(nbsig*(igau-1)+4)- &
                                      c1*sigma(4)- &
                                      epsi_creep(4)
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
end subroutine
