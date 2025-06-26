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

subroutine te0312(option, nomte)
!
! --------------------------------------------------------------------------------------------------
!   Les éléments qui appellent ce "Te" ne peuvent pas calculer les chargements listés ci-dessous
!
!   On vérifie que :
!       les caractéristiques matériaux donnés par l'utilisateur sont nulles de cette façon
!       même si le champ est donné, il ne sera pas pris en compte.
!
!   CHARGEMENTS :
!       SÉCHAGE     : CHAR_MECA_SECH_R
!       HYDRATATION : CHAR_MECA_HYDR_R
!       TEMPÉRATURE : CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
!
    implicit none

#include "jeveux.h"
#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lteatt.h"
#include "asterfort/jeexin.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/utmess.h"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lmater, iret, ig, icp
    integer(kind=8) :: icodre(2), kpg, spt
    integer(kind=8) :: ndim, nno, nnos, npg1, ipoids, ivf, idfde, jgano, igau
    real(kind=8)      :: bendog(1), kdessi(1), alpha(1), epsa(6)
    character(len=8)  :: fami, poum, materi
    character(len=16) :: mult_comp
!
    aster_logical     :: IsPmf
!
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
!
    integer(kind=8), pointer :: cpri(:) => null()
    character(len=16), pointer :: compor(:) => null()
    character(len=24), pointer :: cprk24(:) => null()
!.......................................................................
!
    call jevech('PMATERC', 'L', lmater)
    fami = 'FPG1'; kpg = 1; spt = 1; poum = '+'
    epsa(:) = 0.0d0
!
    IsPmf = ASTER_FALSE
!   Si c'est une PMF : Il peut y avoir plusieurs matériaux sur la maille
    if (lteatt('TYPMOD2', 'PMF')) then
        call jevech('PCOMPOR', 'L', vk16=compor)
        mult_comp = compor(MULTCOMP)
        ! Récupération de la SD_COMPOR où le comportement des groupes de fibres est stocké
        call jeexin(mult_comp, iret)
        ASSERT(iret .ne. 0)
        call jeveuo(mult_comp, 'L', vk24=cprk24)
        call jeveuo(mult_comp(1:8)//'.CPRI', 'L', vi=cpri)
        ! Ceinture et bretelle : seulement pour les multifibres
        ASSERT(cpri(MULTI_FIBER_TYPE) .eq. 3)
        ! Récupération des caractéristiques des fibres
        call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug)
        IsPmf = ASTER_TRUE
    else
        materi = ' '
        nbgrfi = 1
    end if
    !
    ! On boucle sur les groupes de fibre pour aller chercher le matériau du groupe
    do ig = 1, nbgrfi
        if (IsPmf) then
            icp = (nug(ig)-1)*MULTI_FIBER_SIZEK+MULTI_FIBER_MATER
            materi = cprk24(icp) (1:8)
        end if
        !
        if (option .eq. 'CHAR_MECA_HYDR_R') then
            call rcvalb(fami, kpg, spt, poum, zi(lmater), materi, 'ELAS', &
                        0, ' ', [0.d0], 1, 'B_ENDOGE', bendog, icodre, 0)
            ! BENDOGE ABSENT => CHARGEMENT NUL
            if ((icodre(1) .eq. 0) .and. (bendog(1) .ne. 0.d0)) then
                call utmess('F', 'ELEMENTS_22', sk=nomte)
            end if
        else if (option .eq. 'CHAR_MECA_SECH_R') then
            call rcvalb(fami, kpg, spt, poum, zi(lmater), materi, 'ELAS', &
                        0, ' ', [0.d0], 1, 'K_DESSIC', kdessi, icodre, 0)
            ! KDESSI ABSENT => CHARGEMENT NUL
            if ((icodre(1) .eq. 0) .and. (kdessi(1) .ne. 0.d0)) then
                call utmess('F', 'ELEMENTS_23', sk=nomte)
            end if
        else if (option .eq. 'CHAR_MECA_TEMP_R') then
            call rcvalb(fami, kpg, spt, poum, zi(lmater), materi, 'ELAS', &
                        0, ' ', [0.d0], 1, 'ALPHA', alpha, icodre, 0)
            ! ALPHA ABSENT => CHARGEMENT NUL
            if ((icodre(1) .eq. 0) .and. (alpha(1) .ne. 0.d0)) then
                call utmess('F', 'ELEMENTS_19', sk=nomte)
            end if
        else if (option .eq. 'CHAR_MECA_EPSA_R') then
            ! EPSA DOIT ETRE NUL SUR LES ELEMENTS DE STRUCTURE
            call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                             jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)

            do igau = 1, npg1
                epsa(:) = 0.0d0
                call rcvarc(' ', 'EPSAXX', '+', 'RIGI', igau, &
                            1, epsa(1), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(1) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
                call rcvarc(' ', 'EPSAYY', '+', 'RIGI', igau, &
                            1, epsa(2), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(2) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
                call rcvarc(' ', 'EPSAZZ', '+', 'RIGI', igau, &
                            1, epsa(3), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(3) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
                call rcvarc(' ', 'EPSAXY', '+', 'RIGI', igau, &
                            1, epsa(4), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(4) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
                call rcvarc(' ', 'EPSAXZ', '+', 'RIGI', igau, &
                            1, epsa(5), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(5) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
                call rcvarc(' ', 'EPSAYZ', '+', 'RIGI', igau, &
                            1, epsa(6), icodre(1))
                if ((icodre(1) .eq. 0) .and. (epsa(6) .ne. 0.d0)) then
                    call utmess('F', 'ELEMENTS_93', sk=nomte)
                end if
            end do
        else
            ASSERT(.false.)
        end if
    end do
!
end subroutine
