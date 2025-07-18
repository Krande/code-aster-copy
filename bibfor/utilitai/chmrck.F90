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
subroutine chmrck(chmat, nomrc, nommat, nbmtrc)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=8) :: chmat, nommat(*)
    character(len=16) :: nomrc
    integer(kind=8) :: nbmtrc
!
!     ==================================================================
!     ! UTILITAIRE - RECHERCHE DES MATERIAUX D'UN CHAM_MATER QUI       !
!     ! UTILISENT UNE RELATION DE COMPORTEMENT DONNEE               RM !
!     ==================================================================
!     !                                                                !
!     !    ETANT DONNES UN CHAM_MATER ET UNE RELATION DE COMPORTEMENT  !
!     !    CHERCHER LES MATERIAUX DU CHAM_MATER QUI UTILISE CETTE RC   !
!     !                                                                !
!     !    DANS UN MATERIAU, IL N'Y A QU'UNE RC D'UN TYPE DONNE        !
!     !                                                                !
!     ==================================================================
! IN  ! CHMAT  ! K8  ! NOM DU CONCEPT CHAM_MATER                       !
! IN  ! NOMRC  ! K8  ! NOM DE LA RC CHERCHEE                           !
! IN  ! NBMTCH ! IS  ! LONGUEUR DU .VALE DU CHAM_MATER                 !
! OUT ! NBMTRC ! IS  ! NOMBRE DE MAT QUI UTILISE LA RC                 !
! OUT ! NOMMAT ! K8  !LISTE DES MAT DU CHAM_MATER QUI UTILISENT LA RC  !
!     ==================================================================
!
!
! --- VARIABLES LOCALES ---
    character(len=8) :: kmat, kbid
    character(len=24) :: krc
    integer(kind=8) :: arc, imat, nbrc, ipos, ncmpmx, izone, i, nbzone
    integer(kind=8) :: l1, nbzmax, k
    integer(kind=8), pointer :: desc(:) => null()
    character(len=8), pointer :: vale(:) => null()
    parameter(ncmpmx=30)
!
! ====================== DEBUT DU PROGRAMME ============================
!
    call jemarq()
    kbid = ' '
    call jeveuo(chmat//'.CHAMP_MAT .VALE', 'L', vk8=vale)
    call jelira(chmat//'.CHAMP_MAT .VALE', 'LONMAX', l1)
    call jeveuo(chmat//'.CHAMP_MAT .DESC', 'L', vi=desc)
    nbzmax = desc(2)
    nbzone = desc(3)
!     ON VERIFIE QUE LA TAILLE DE LA CARTE EST BIEN TOUJOURS DE 30
!     PAR ZONE
    ASSERT(l1 .eq. (nbzmax*ncmpmx))
!
!
    nbmtrc = 0
    do izone = 1, nbzone
        do i = 1, ncmpmx
            imat = (izone-1)*ncmpmx+i
            kmat = vale(imat)
            if (kmat .eq. ' ') goto 50
            if (kmat .eq. 'TREF=>') goto 50
            krc = kmat//'.MATERIAU.NOMRC'
            call jeveuo(krc, 'L', arc)
            call jelira(krc, 'LONMAX', nbrc)
            ipos = 0
            do k = 1, nbrc
                if (zk32(arc+k-1) .eq. nomrc) then
                    ipos = k
                end if
            end do
            if (ipos .gt. 0) then
                nbmtrc = nbmtrc+1
                nommat(nbmtrc) = kmat
            end if
        end do
50      continue
    end do
!
!
    call jedema()
end subroutine
