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

subroutine xfisco(noma, modelx)
! person_in_charge: patrick.massin at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=8) :: noma, modelx
!
! CREATION D'UN CHAMP ÉLÉMENTAIRE DE CONNECTIVITÉ DES FISSURES BRANCHÉES
!
!
!
!
!
    integer(kind=8) ::  jcesd, jcesl, ibid, iret, nncp
    integer(kind=8) :: jcesd2, jcesl2, iad, iad3
    integer(kind=8) :: ima, nbma, ifiss, ifis2, ifis3, nfiss, nfis2
    character(len=19) :: ces, ces2, ligrel, chglo
    character(len=8) :: nomfis, nomfi3, licmp(2), valk(3)
    character(len=8), pointer :: jonfiss(:) => null()
    integer(kind=8), pointer :: joncoef(:) => null()
    integer(kind=8), pointer :: nbsp(:) => null()
    character(len=8), pointer :: cesv2(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- RECUPERATION DU NOMBRE DE SOUS POINT (NBRE DE FISSURES VUES)
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
!
! --- CONSTRUCTION DU CHAMP SIMPLE TEMPORAIRE
!
    ces = '&&XFISCO.FISCO'
    licmp(1) = 'X1'
    licmp(2) = 'X2'
!
    call cescre('V', ces, 'ELEM', noma, 'NEUT_I', &
                2, licmp, [ibid], nbsp, [-2])
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESV', 'E', vi=cesv)
    call jeveuo(ces//'.CESL', 'E', jcesl)
!
! --- RECUPERATION DU CHAMP ELEM S CONTENANT LE NOM DES FISSURES VUES
!
    ces2 = '&&XCONNO.CES2'
    call jeveuo(ces2//'.CESD', 'L', jcesd2)
    call jeveuo(ces2//'.CESV', 'E', vk8=cesv2)
    call jeveuo(ces2//'.CESL', 'E', jcesl2)
!
! --- RECUPERATION NOMBRE DE MAILLES
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
!
! --- BOUCLE SUR LES MAILLES
!
    do ima = 1, nbma
        nfiss = nbsp(ima)
!
! --- BOUCLE SUR LES FISSURE DE LA MAILLE
!
        do ifiss = 1, nfiss
            call cesexi('S', jcesd2, jcesl2, ima, 1, &
                        ifiss, 1, iad)
            ASSERT(iad .gt. 0)
            nomfis = cesv2(iad)
            call jeexin(nomfis//'.JONFISS', iret)
            if (iret .ne. 0) then
                call jeveuo(nomfis//'.JONFISS', 'L', vk8=jonfiss)
                call jeveuo(nomfis//'.JONCOEF', 'L', vi=joncoef)
                call jelira(nomfis//'.JONFISS', 'LONMAX', nfis2)
!
! --- BOUCLE SUR LES FISSURES DE LA MAILLE IFIS3, AUTRE QUE IFISS
!
                do ifis3 = 1, nfiss
                    if (ifis3 .eq. ifiss) goto 130
!
! --- RECUPERATION DU NOM GLOBALE  NOMFI3 DE IFIS3
!
                    call cesexi('S', jcesd2, jcesl2, ima, 1, &
                                ifis3, 1, iad)
                    nomfi3 = cesv2(iad)
!
! --- ON REGARDE SI LA FISSURE NOMFI3 EST CONNECTÉ À NOMFIS
!
                    do ifis2 = 1, nfis2
                        if (jonfiss(ifis2) .eq. nomfi3) then
                            if (ifis3 .gt. ifiss) then
                                valk(1) = nomfis
                                valk(2) = nomfi3
                                call utmess('F', 'XFEM_46', nk=2, valk=valk)
                            end if
                            call cesexi('S', jcesd, jcesl, ima, 1, &
                                        ifiss, 1, iad)
                            if (iad .gt. 0) then
                                call cesexi('S', jcesd, jcesl, ima, 1, &
                                            ifis3, 1, iad3)
                                if (cesv(iad) .eq. cesv(iad3)) then
                                    iad = -iad
                                else
                                    valk(1) = nomfis
                                    valk(3) = nomfi3
                                    call utmess('F', 'XFEM_47', nk=3, valk=valk)
                                end if
                            end if
                            valk(2) = nomfi3
                            zl(jcesl-1-iad) = .true.
                            cesv(1-1-iad) = ifis3
                            call cesexi('S', jcesd, jcesl, ima, 1, &
                                        ifiss, 2, iad)
                            if (iad .gt. 0) iad = -iad
                            zl(jcesl-1-iad) = .true.
                            cesv(1-1-iad) = joncoef(ifis2)
                        end if
                    end do
130                 continue
                end do
            end if
!
! --- SI ON A RIEN TROUVER
!
            call cesexi('S', jcesd, jcesl, ima, 1, &
                        ifiss, 1, iad)
            if (iad .lt. 0) then
                zl(jcesl-1-iad) = .true.
                cesv(1-1-iad) = 0
                call cesexi('S', jcesd, jcesl, ima, 1, &
                            ifiss, 2, iad)
                ASSERT(iad .lt. 0)
                zl(jcesl-1-iad) = .true.
                cesv(1-1-iad) = 0
            end if
        end do
!
    end do
!
! --- CONVERSION CHAM_ELEM_S -> CHAM_ELEM
!
    chglo = modelx(1:8)//'.FISSCO'
    ligrel = modelx(1:8)//'.MODELE'
    call cescel(ces, ligrel, 'TOPOSE', 'PFISCO', 'OUI', &
                nncp, 'G', chglo, 'F', ibid)
    call detrsd('CHAM_ELEM_S', ces)
!
    call jedema()
end subroutine
