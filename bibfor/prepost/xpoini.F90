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

subroutine xpoini(maxfem, modele, malini, modvis, licham, &
                  resuco, resux, prefno, nogrfi)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=2) :: prefno(4)
    character(len=8) :: maxfem, modele, malini, resuco, resux, modvis
    character(len=24) :: licham, nogrfi
!
!
!               RECUPERATION DES ENTREES SORTIES
!               POUR LES OPERATEURS DE POST-TRAITEMENT X-FEM
!
!
!   OUT
!       MAXFEM : MAILLAGE X-FEM
!       MODELE : MODELE FISSURE
!       MALINI : MAILLAGE SAIN
!       MODVIS : MODELE DE VISU (X-FEM)
!       LICHAM : LISTE DES CHAMPS A POST-TRAITER
!       RESUCO : NOM DU CONCEPT RESULTAT DONT ON EXTRAIT LES CHAMPS
!       RESUX  : NOM DU CONCEPT RESULTAT A CREER
!       PREFNO : PREFERENCES POUR LE NOMAGE DES NOUVELLES ENTITES
!       NOGRFI : NOM DU GROUPE DES NOEUDS SITUES SUR LA FISSURE
!
    integer(kind=8) :: iret, ibid, jlicha, jxc, i
    integer(kind=8) :: nbcham, nchmax
!     NOMBRE MAX DE CHAMPS A POST-TRAITER
    parameter(nchmax=3)
    character(len=8) :: k8b
    character(len=16) :: k16b, nomcmd, tysd, linom(nchmax)
    character(len=19) :: k19bid
!
!
    call jemarq()
!
!     NOM DE LA COMMANDE (POST_MAIL_XFEM OU POST_CHAM_XFEM)
    call getres(k8b, k16b, nomcmd)
!
!     ------------------------------------------------------------------
    if (nomcmd .eq. 'POST_MAIL_XFEM') then
!     ------------------------------------------------------------------
!
!       NOM DU MAILLAGE DE SORTIE : MAXFEM
        call getres(maxfem, k16b, k16b)
!
!       MODELE ENRICHI : MODELE
        call getvid(' ', 'MODELE', scal=modele, nbret=iret)
        call exixfe(modele, iret)
        if (iret .eq. 0) then
            call utmess('F', 'XFEM_3', sk=modele)
        end if
!
!
!       PREFERENCES POUR LE NOMAGE DES NOUVELLES ENTITES
        call getvtx(' ', 'PREF_NOEUD_X', scal=prefno(1), nbret=ibid)
        call getvtx(' ', 'PREF_NOEUD_M', scal=prefno(2), nbret=ibid)
        call getvtx(' ', 'PREF_NOEUD_P', scal=prefno(3), nbret=ibid)
        call getvtx(' ', 'PREF_MAILLE_X', scal=prefno(4), nbret=ibid)
        call getvtx(' ', 'PREF_GROUP_CO', scal=nogrfi, nbret=ibid)
!
!     ------------------------------------------------------------------
    else if (nomcmd .eq. 'POST_CHAM_XFEM') then
!     ------------------------------------------------------------------
!
!       NOM DE LA SD RESULTAT A CREER : RESUX
        call getres(resux, k16b, k16b)
!
!       MODELE DE VISU ET MAILLAGE DE VISU (X-FEM)
        call getvid(' ', 'MODELE_VISU', scal=modvis, nbret=iret)
        call dismoi('NOM_MAILLA', modvis, 'MODELE', repk=maxfem)
!
!       VERIFICATION QUE LE MODELE de VISU N'EST PAS X-FEM
        call exixfe(modvis, ibid)
        if (ibid .eq. 1) then
            call utmess('F', 'XFEM2_30')
        end if
!
!       NOM ET TYPE DE LA SD RESULTAT EN ENTREE : RESUCO
        call getvid(' ', 'RESULTAT', scal=resuco, nbret=ibid)
!
!       MODELE ENRICHI ASSOCIE AU RESULTAT EN ENTREE
        call dismoi('NOM_MODELE', resuco, 'RESULTAT', repk=modele)
!
!       NOM DES CHAMPS A POST-TRAITER
        call gettco(resuco, tysd)
!
        if (tysd(1:9) .eq. 'MODE_MECA') then
            nbcham = 1
            linom(1) = 'DEPL'
        else if (tysd(1:9) .eq. 'EVOL_NOLI') then
!         A CORRIGER SUITE FICHE 15408
!         PB POST-TRAITEMENT VARIABLES INTERNES SI CONTACT P2P1 (GLUTE)
            call jeveuo(modele//'.XFEM_CONT', 'L', jxc)
            if (zi(jxc-1+1) .eq. 3) then
                write (6, *) 'ON NE PEUT PAS POST-TRAITER LE CHAMP VARI_ELGA'
                nbcham = 2
                linom(1) = 'DEPL'
                linom(2) = 'SIEF_ELGA'
            else
                nbcham = 3
                linom(1) = 'DEPL'
                linom(2) = 'SIEF_ELGA'
                linom(3) = 'VARI_ELGA'
            end if
        else if (tysd(1:9) .eq. 'EVOL_ELAS') then
            call rsexch(' ', resuco, 'SIEF_ELGA', 1, k19bid, &
                        iret)
            if (iret .eq. 0) then
                nbcham = 2
                linom(1) = 'DEPL'
                linom(2) = 'SIEF_ELGA'
            else
                nbcham = 1
                linom(1) = 'DEPL'
            end if
        else if (tysd(1:9) .eq. 'EVOL_THER') then
            nbcham = 1
            linom(1) = 'TEMP'
        end if
!
        call wkvect(licham, 'V V K16', nbcham, jlicha)
        do i = 1, nbcham
            zk16(jlicha-1+i) = linom(i)
        end do
!
!     ----------------------------------------------------------------
    end if
!     ------------------------------------------------------------------
!
!     MAILLAGE INITIAL : MALINI
!     CE MAILLAGE EST CELUI ASSOCIE AU MODELE ENRICHI
!     SAUF DANS LE CAS DU CONTACT AU ARETE 'P1P1A'
    call jeveuo(modele//'.XFEM_CONT', 'L', jxc)
!
!     MAILLAGE_SAIN NE SERT A RIEN :
!     RECUPERATION DU MAILLAGE ASSOCIE AU MODELE
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=malini)
!
    call jedema()
end subroutine
