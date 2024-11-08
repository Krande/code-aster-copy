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
subroutine ace_get_node_reparti(nbocc, mclef, infdonn, grplmax)
!
! --------------------------------------------------------------------------------------------------
!     AFFE_CARA_ELEM
!
!       Évaluation du nombre de noeud à répartir
!
! --------------------------------------------------------------------------------------------------
!   Entrée :
!       nbocc   : nombre d’occurrence du mot clef
!       mclef   : le mot clef
!       infdonn : information sur le concept affe_cara_elem
!       grplmax : tampon déjà dimensionné au nombre maxi de groupe de maille
!
!   Sortie
!       infdonn : mise à jour des informations
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_info_type
!
    implicit none
    integer(kind=8)         :: nbocc
    character(len=16)       :: mclef
    type(cara_elem_info)    :: infdonn
    character(len=24)       :: grplmax(*)
!

#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/getvtx.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!

    character(len=8)  :: noma
    integer(kind=8)   :: jdme, GroupeMaxOccur
!
    character(len=24) :: magrma, manoma
!-----------------------------------------------------------------------
    integer(kind=8) :: ii, in
    integer(kind=8) :: ioc, ldgm, nbnoe, nbgr
    integer(kind=8) :: nbmail, NoeudMaxMaille, MailleMaxOccur, MailleOccur
!-----------------------------------------------------------------------
!
    if (nbocc .eq. 0) goto 999
!
    call jemarq()
!   Récupère les informations
    noma = infdonn%maillage
    jdme = infdonn%jmodmail
    GroupeMaxOccur = infdonn%GroupeMaxOccur
!
    magrma = noma//'.GROUPEMA'
    manoma = noma//'.CONNEX'
!
    NoeudMaxMaille = 0
    MailleMaxOccur = 0
!
!   On parcours les occurences
!       On éclate les GROUP_MA
!           MailleOccur     : nombre de maille affectée dans une occurence
!           NoeudMaxMaille  : nombre maximum de noeud dans une maille
!
!   MailleMaxOccur = max( MailleOccur ) : nombre maximum de maille dans une occurence
!   NoeudMaxMaille                      : nombre maximum de noeud dans une maille
!
!   Le nombre maximum des noeuds concernés par la répartition :
!       MailleMaxOccur*NoeudMaxMaille
!
!   ACE_RIGI_MISS_3D
    if (mclef .eq. 'RIGI_MISS_3D') then
        do ioc = 1, nbocc
!           Dans le catalogue pour GROUP_MA_POI1 : "o", max=1
            MailleOccur = 0
            call getvtx(mclef, 'GROUP_MA_POI1', iocc=ioc, nbval=GroupeMaxOccur, &
                        vect=grplmax, nbret=nbgr)
            if (nbgr .ne. 0) then
                do ii = 1, nbgr
                    call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nbmail)
                    MailleOccur = MailleOccur+nbmail
!                   Il y a 1 noeud par maille
                    NoeudMaxMaille = max(NoeudMaxMaille, 1)
                end do
            end if
            MailleMaxOccur = max(MailleMaxOccur, MailleOccur)
!           Dans le catalogue pour GROUP_MA_SEG2 : "f", max=1
            MailleOccur = 0
            call getvtx(mclef, 'GROUP_MA_SEG2', iocc=ioc, nbval=GroupeMaxOccur, &
                        vect=grplmax, nbret=nbgr)
            if (nbgr .ne. 0) then
                do ii = 1, nbgr
                    call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nbmail)
                    MailleOccur = MailleOccur+nbmail
!                   Il y a 2 noeuds par maille
                    NoeudMaxMaille = max(NoeudMaxMaille, 2)
                end do
            end if
            MailleMaxOccur = max(MailleMaxOccur, MailleOccur)
        end do
!
!   ACE_RIGI_PARASOL, ACE_MASS_AJOU
    else
        do ioc = 1, nbocc
            MailleOccur = 0
            call getvtx(mclef, 'GROUP_MA', iocc=ioc, nbval=GroupeMaxOccur, &
                        vect=grplmax, nbret=nbgr)
            do ii = 1, nbgr
                call jelira(jexnom(magrma, grplmax(ii)), 'LONUTI', nbmail)
                call jeveuo(jexnom(magrma, grplmax(ii)), 'L', ldgm)
                MailleOccur = MailleOccur+nbmail
                do in = 0, nbmail-1
                    call jelira(jexnum(manoma, zi(ldgm+in)), 'LONMAX', nbnoe)
                    NoeudMaxMaille = max(NoeudMaxMaille, nbnoe)
                end do
            end do
            MailleMaxOccur = max(MailleMaxOccur, MailleOccur)
        end do
    end if
!
    infdonn%NoeudMaxMaille = max(NoeudMaxMaille, infdonn%NoeudMaxMaille)
    infdonn%MailleMaxOccur = max(MailleMaxOccur, infdonn%MailleMaxOccur)
!
!   Les mailles surfaciques avec le plus grand nombre de noeud : QUA9
!   On pourrait donc fixer NoeudMaxMaille à 9 sans vérifier
!   * Pour ACE_RIGI_PARASOL, ACE_MASS_AJOU ==> supprimer la boucle "in"
!   * Pour ACE_RIGI_MISS_3D. Plus de max(NoeudMaxMaille, 1 | 2)
    ASSERT(infdonn%NoeudMaxMaille .le. 9)
!
    call jedema()
999 continue
end subroutine
