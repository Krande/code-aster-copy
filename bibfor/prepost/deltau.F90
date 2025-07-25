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
subroutine deltau(jrwork, jnbpg, nbpgt, nbordr, ordini, &
                  nmaini, nbmap, numpaq, tspaq, nommet, &
                  nomcri, nomfor, grdvie, forvie, forcri, &
                  cesr)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/acgrdo.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcpare.h"
#include "asterfort/rnomat.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: jrwork, jnbpg, nbpgt, nbordr, nmaini, numpaq, nbmap
    integer(kind=8) :: tspaq, ordini
    character(len=16) :: nomcri, nommet, nomfor, forvie, forcri, grdvie
    character(len=19) :: cesr
! ---------------------------------------------------------------------
! BUT: DETERMINER LE PLAN INCLINE POUR LEQUEL DELTA_TAU EST MAXIMUM
!      POUR CHAQUE POINT DE GAUSS D'UN <<PAQUET>> DE MAILLES.
! ---------------------------------------------------------------------
! ARGUMENTS:
! JRWORK     IN    I  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                       L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                       ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                       DU <<PAQUET>> DE MAILLES.
! JNBPG      IN    I  : ADRESSE DU VECTEUR CONTENANT LE NOMBRE DE
!                       POINT DE GAUSS DE CHAQUE MAILLE DU MAILLAGE.
! NBPGT      IN    I  : NOMBRE TOTAL DE POINTS DE GAUSS A TRAITER.
! NBORDR     IN    I  : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                       STRUCTURE DE DONNEES RESULTAT.
! ORDINI     IN    I  : ORDRE INITIAL POUR LE CHARGEMENT CYCLIQUE
! NMAINI     IN    I  : NUMERO DE LA 1ERE MAILLE DU <<PAQUET>> DE
!                       MAILLES COURANT.
! NBMAP      IN    I  : NOMBRE DE MAILLES DANS LE <<PAQUET>> DE
!                       MAILLES COURANT.
! NUMPAQ     IN    I  : NUMERO DU PAQUET DE MAILLES COURANT.
! TSPAQ      IN    I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                       COURANT.
! NOMMET     IN    K16: NOM DE LA METHODE DE CALCUL DU CERCLE
!                       CIRCONSCRIT.
! NOMCRI     IN    K16: NOM DU CRITERE AVEC PLANS CRITIQUES.
! CESR       IN    K19: NOM DU CHAMP SIMPLE DESTINE A RECEVOIR LES
!                       RESULTATS.
!
! REMARQUE :
!  - LA TAILLE DU SOUS-PAQUET EST EGALE A LA TAILLE DU <<PAQUET>> DE
!    MAILLES DIVISEE PAR LE NOMBRE DE NUMERO D'ORDRE (NBORDR).
!-----------------------------------------------------------------------
!
    integer(kind=8) :: kwork, jcerd, jcerl, jad
    integer(kind=8) :: iret, imap, icesd, icesl, icesv, ibid
    integer(kind=8) :: ipg
    integer(kind=8) :: nbpg, sompgw, nbpgp, l
    integer(kind=8) :: icmp
    real(kind=8) :: vala, valb, coefpa, vresu2(24), valpar(35)
    integer(kind=8) :: icodwo
    character(len=8) :: chmat1, nommat
    character(len=10) :: optio
    character(len=19) :: chmat, cesmat
    real(kind=8), pointer :: cerv(:) => null()
!
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    call jemarq()
!
! !
! OBTENTION DES ADRESSES '.CESD', '.CESL' ET '.CESV' DU CHAMP SIMPLE
! DESTINE A RECEVOIR LES RESULTATS : DTAUM, ....
!
    call jeveuo(cesr//'.CESD', 'L', jcerd)
    call jeveuo(cesr//'.CESL', 'E', jcerl)
    call jeveuo(cesr//'.CESV', 'E', vr=cerv)
! !
! !
! RECUPERATION MAILLE PAR MAILLE DU MATERIAU DONNE PAR L'UTILISATEUR
!
    call getvid(' ', 'CHAM_MATER', scal=chmat1, nbret=iret)
    chmat = chmat1//'.CHAMP_MAT'
    cesmat = '&&DELTAU.CESMAT'
    call carces(chmat, 'ELEM', ' ', 'V', cesmat, &
                'A', iret)
    call jeveuo(cesmat//'.CESD', 'L', icesd)
    call jeveuo(cesmat//'.CESL', 'L', icesl)
    call jeveuo(cesmat//'.CESV', 'L', icesv)
!
!
! CONSTRUCTION DU VECTEUR : CONTRAINTE = F(NUMERO D'ORDRE) EN CHAQUE
! POINT DE GAUSS DU PAQUET DE MAILLES.
    l = 1
    nbpg = 0
    nbpgp = 0
    kwork = 0
    sompgw = 0
!
    do imap = nmaini, nmaini+(nbmap-1)
        if (imap .gt. nmaini) then
            kwork = 1
            sompgw = sompgw+zi(jnbpg+imap-2)
        end if
        nbpg = zi(jnbpg+imap-1)
! SI LA MAILLE COURANTE N'A PAS DE POINTS DE GAUSS, LE PROGRAMME
! PASSE DIRECTEMENT A LA MAILLE SUIVANTE.
        if (nbpg .eq. 0) then
            goto 400
        end if
!
        nbpgp = nbpgp+nbpg
        if ((l*int(nbpgt/10.0d0)) .lt. nbpgp) then
            write (6, *) numpaq, '   ', (nbpgp-nbpg)
            l = l+1
        end if
!
! RECUPERATION DU NOM DU MATERIAU AFFECTE A LA MAILLE COURANTE
! ET DES PARAMETRES ASSOCIES AU CRITERE CHOISI POUR LA MAILLE COURANTE.
!
        optio = 'DOMA_ELGA'
        call rnomat(icesd, icesl, icesv, imap, nomcri, &
                    ibid, ibid, ibid, optio, vala, &
                    valb, coefpa, nommat)
!
        call rcpare(nommat, 'FATIGUE', 'WOHLER', icodwo)
        if (icodwo .eq. 1) then
            call utmess('F', 'FATIGUE1_90', sk=nomcri(1:16))
        end if
!
!
        do ipg = 1, nbpg
!
!            call jerazo('&&DELTAU.VECTPG', tneces, 1)
!
! REMPACER PAR ACMATA
            call acgrdo(nbordr, ordini, kwork, sompgw, jrwork, &
                        tspaq, ipg, nommet, nommat, nomcri, &
                        vala, coefpa, nomfor, grdvie, forvie, &
                        forcri, valpar, vresu2)
!
!
! C AFFECTATION DES RESULTATS DANS UN CHAM_ELEM SIMPLE
!
            do icmp = 1, 24
                call cesexi('C', jcerd, jcerl, imap, ipg, &
                            1, icmp, jad)
!
!              -- TOUTES LES MAILLES NE SAVENT PAS CALCULER LA FATIGUE :
                if (jad .eq. 0) then
                    ASSERT(icmp .eq. 1)
                    ASSERT(ipg .eq. 1)
                    goto 400
                end if
                jad = abs(jad)
                zl(jcerl-1+jad) = .true.
                cerv(jad) = vresu2(icmp)
!
            end do
!
        end do
400     continue
    end do
!
! MENAGE
!
    call detrsd('CHAM_ELEM_S', cesmat)
! !
!
    call jedema()
end subroutine
