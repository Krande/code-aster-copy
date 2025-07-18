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
subroutine mpmod2(basemo, nommes, nbmesu, nbmtot, basepr, &
                  vnoeud, vrange, vcham)
!
!
!     PROJ_MESU_MODAL : PROJECTION DE LA MATRICE MODALE SUR LES CMP
!                       DES NOEUDS MESURE
!
!     IN  : BASEMO : NOM DE LA BASE DE PROJECTION
!     IN  : NOMMES : NOM DE LA MESURE
!     IN  : NBMESU : NOMBRE DE MESURE (DATASET 58)
!     IN  : NBMTOT : NOMBRE DE VECTEURS DE BASE
!     OUT  : BASEPR : NOM BASE PROJETEE SUIVANT DIRECTION MESURE
!     IN  : VNOEUD : NOM RANGEMENT NOEUD MESURE
!     IN  : VRANGE : NOM CORRESPONDANCE CMP SUIVANT VNOEUD
!     IN  : VCHAM : NOM CORRESPONDANCE NOMCHAMP SUIVANT VNOEUD
!
    implicit none
!     ------------------------------------------------------------------
!
!
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8prem.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsprj.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mpjeft.h"
#include "asterfort/mpmod3.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: basemo, nommes
    character(len=24) :: vnoeud, vrange, basepr, vcham
    integer(kind=8) :: nbmesu, nbmtot
!
    character(len=8) :: nomres
    character(len=16) :: nomcha, corres, nomch, typres, k16bid, nomchm
    character(len=19) :: chamno, ch1s, ch2s
    character(len=24) :: vorien
!
    integer(kind=8) :: lred, lori, lrange
    integer(kind=8) :: imesu, ii, imode, iret
    integer(kind=8) :: iposd, icmp, ino
    integer(kind=8) :: lnoeud, nnoema, ncmpma
    integer(kind=8) :: jcnsv, jcnsl, jcnsk
    integer(kind=8) :: ibid, nbcmpi, nbcham, lch, ich, lcham
!
    real(kind=8) :: vori(3)
    real(kind=8) :: val, vect(3)
    integer(kind=8), pointer :: ordr(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! RECUPERATION DU NOM DU CONCEPT RESULTAT
    call getres(nomres, typres, k16bid)
!
! =============================================
! RECUPERATION DES OBJETS LIES A LA MESURE
! =============================================
!
    call mpmod3(basemo, nommes, nbmesu, nbmtot, vcham, &
                vnoeud, vrange, vorien, nnoema, ncmpma)
!
    call jeveuo(vnoeud, 'L', lnoeud)
    call jeveuo(vrange, 'L', lrange)
    call jeveuo(vcham, 'L', lcham)
    call jeveuo(vorien, 'L', lori)
!
! RECUPERATION DES NOMS DU CHAMP MESURE
!
    call getvtx('MODELE_MESURE', 'NOM_CHAM', iocc=1, nbval=0, nbret=nbcham)
    if (nbcham .ne. 0) then
        nbcham = -nbcham
    else
        call utmess('A', 'ALGORITH10_93')
    end if
    call wkvect('&&LISTE_CHAMP', 'V V K16', nbcham, lch)
    call getvtx('MODELE_MESURE', 'NOM_CHAM', iocc=1, nbval=nbcham, vect=zk16(lch), &
                nbret=ibid)
!
!     -> OBJET MATRICE MODALE REDUITE SUIVANT DIRECTION DE MESURE
!
    basepr = nomres//'.PROJM    .PJMBP'
!
    call wkvect(basepr, 'G V R', nnoema*ncmpma*nbmtot, lred)
!
!     INITIALISATION DE BASEPR
    do ii = 1, nnoema*ncmpma*nbmtot
        zr(lred-1+ii) = 0.d0
    end do
!
! RECUPERATION SD CORRESPONDANCE ENTRE MAILLAGE MODELE/MESURE
!
    corres = '&&PJEFTE.CORRESP'
    call mpjeft(corres)
!
! BOUCLE SUR LES CHAMPS MESURES
! POUR L'INSTANT ON NE RAJOUTE PAS DE PONDERATION SUR LES MESURE
! A FAIRE EVENTUELLEMENT : EN FONCTION DU TYPE DE CHAMP
!
    do ich = 1, nbcham
        nomcha = zk16(lch-1+ich)
        nomch = nomcha
! MEMES VECTEURS DE BASE POUR : DEPL, VITE ET ACCE
        if (nomch .eq. 'VITE' .or. nomch .eq. 'ACCE') nomch = 'DEPL'
! MEMES VECTEURS DE BASE POUR LES CONTRAINTES
        if (nomch(1:4) .eq. 'SIEF') nomch = 'SIGM_NOEU'
! MEMES VECTEURS DE BASE POUR LES DEFORMATIONS
        if (nomch(1:4) .eq. 'EPSI') nomch = 'EPSI_NOEU'
!
        call jeveuo(basemo//'           .ORDR', 'L', vi=ordr)
!
        ch1s = '&&PJEFPR.CH1S'
        ch2s = '&&PJEFPR.CH2S'
!
! BOUCLE SUR TOUS LES MODES
!
        do imode = 1, nbmtot
!
            call rsexch('F', basemo, nomch, ordr(imode), chamno, &
                        iret)
!
!       2-1 : TRANSFORMATION DE CHAMNO EN CHAM_NO_S : CH1S
            call detrsd('CHAM_NO_S', ch1s)
            call cnocns(chamno, 'V', ch1s)
!
!       2-2 : PROJECTION DU CHAM_NO_S : CH1S -> CH2S
            call detrsd('CHAM_NO_S', ch2s)
            call cnsprj(ch1s, corres, 'V', ch2s, iret)
            if (iret .gt. 0) then
                call utmess('F', 'ALGORITH6_25')
            end if
!
            call jeveuo(ch2s//'.CNSK', 'L', jcnsk)
            call jeveuo(ch2s//'.CNSD', 'L', vi=cnsd)
            call jeveuo(ch2s//'.CNSC', 'L', vk8=cnsc)
            call jeveuo(ch2s//'.CNSV', 'L', jcnsv)
            call jeveuo(ch2s//'.CNSL', 'L', jcnsl)
!
            nbcmpi = cnsd(2)
!
! BOUCLE SUR LES POINTS DE MESURE
!
            do imesu = 1, nbmesu
!
                nomchm = zk16(lcham-1+imesu)
! MEMES VECTEURS DE BASE POUR : DEPL, VITE ET ACCE
                if (nomchm .eq. 'VITE' .or. nomchm .eq. 'ACCE') nomchm = 'DEPL'
! MEMES VECTEURS DE BASE POUR LES CONTRAINTES
                if (nomchm(1:4) .eq. 'SIEF') nomchm = 'SIGM_NOEU'
! MEMES VECTEURS DE BASE POUR LES DEFORMATIONS
                if (nomchm(1:4) .eq. 'EPSI') nomchm = 'EPSI_NOEU'
!
! NUMERO DU NOEUD ASSOCIE A IMESU : INO
                ino = zi(lnoeud-1+imesu)
!
                if ((nomch(1:4) .eq. 'DEPL') .and. (nomchm(1:4) .eq. 'DEPL')) then
!
! RECUPERATION DIRECTION DE MESURE (VECTEUR DIRECTEUR)
                    do ii = 1, 3
                        vori(ii) = zr(lori-1+(imesu-1)*3+ii)
                    end do
!
! NORMALISATION DU VECTEUR DIRECTEUR
                    val = 0.d0
                    do ii = 1, 3
                        val = val+vori(ii)*vori(ii)
                    end do
                    val = sqrt(val)
                    if (val .lt. r8prem()) then
                        call utmess('F', 'ALGORITH6_26')
                    end if
                    do ii = 1, 3
                        vori(ii) = vori(ii)/val
                    end do
!
! RECUPERATION DU CHAMP AU NOEUD (BASE)
!
                    do icmp = 1, nbcmpi
                        if (cnsc(icmp) .eq. 'DX') vect(1) = zr(jcnsv-1+(ino-1)*nbcmpi+icmp)
                        if (cnsc(icmp) .eq. 'DY') vect(2) = zr(jcnsv-1+(ino-1)*nbcmpi+icmp)
                        if (cnsc(icmp) .eq. 'DZ') vect(3) = zr(jcnsv-1+(ino-1)*nbcmpi+icmp)
                    end do
!
! CALCUL DE LA BASE RESTREINTE
!
                    iposd = (imode-1)*nbmesu+imesu
                    zr(lred-1+iposd) = 0.d0
!
                    do ii = 1, 3
                        zr(lred-1+iposd) = zr(lred-1+iposd)+vect(ii)*vori(ii)
                    end do
!
                else if ((nomch(1:14) .eq. 'EPSI_NOEU') .and. &
                         (nomchm(1:14) .eq. 'EPSI_NOEU')) then
!
                    iposd = (imode-1)*nbmesu+imesu
                    do icmp = 1, nbcmpi
                        if (cnsc(icmp) .eq. zk8(lrange-1+imesu)) zr(lred-1+iposd) = &
                            zr( &
                            jcnsv-1+(ino-1)*nbcmpi+icmp)
                    end do
!
                else if ((nomch(1:14) .eq. 'SIGM_NOEU') .and. &
                         (nomchm(1:14) .eq. 'SIGM_NOEU')) then
!
                    iposd = (imode-1)*nbmesu+imesu
                    do icmp = 1, nbcmpi
                        if (cnsc(icmp) .eq. zk8(lrange-1+imesu)) zr(lred-1+iposd) = &
                            zr( &
                            jcnsv-1+(ino-1)*nbcmpi+icmp)
                    end do
!
                end if
!
! FIN DE LA BOUCLE SUR LES POINTS DE MESURE
!
            end do
!
! FIN DE LA BOUCLE SUR LES MODES
!
        end do
!
! FIN BOUCLE SUR LES NOMCHA
!
    end do
!
    call jeecra(basepr, 'LONUTI', nbmesu*nbmtot)
!
! DESTRUCTION DES VECTEURS DE TRAVAIL
!
    call detrsd('CHAM_NO_S', ch1s)
    call detrsd('CHAM_NO_S', ch2s)
    call detrsd('CORRESP_2_MAILLA', corres)
!
    call jedema()
!
end subroutine
