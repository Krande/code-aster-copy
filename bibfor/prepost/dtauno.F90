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

subroutine dtauno(jrwork, lisnoe, nbnot, nbordr, ordini, &
                  nnoini, nbnop, tspaq, nommet, nomcri, &
                  nomfor, grdvie, forvie, forcri, nommai, &
                  cnsr, nommap, post, valpar, vresu)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/acgrdo.h"
#include "asterfort/carces.h"
#include "asterfort/cncinv.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rcpare.h"
#include "asterfort/recofa.h"
#include "asterfort/rnomat.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: jrwork, nbnot, lisnoe(nbnot), nbordr, nnoini, nbnop
    integer(kind=8) :: tspaq, ordini
    aster_logical :: post
    real(kind=8) :: vresu(24), valpar(35)
    character(len=8) :: nommai, nommap
    character(len=16) :: nomcri, nommet, nomfor, forvie, forcri, grdvie
    character(len=19) :: cnsr
! ---------------------------------------------------------------------
! BUT: DETERMINER LE PLAN INCLINE POUR LEQUEL DELTA_TAU EST MAXIMUM
!      POUR CHAQUE NOEUD D'UN <<PAQUET>> DE NOEUDS.
! ---------------------------------------------------------------------
! ARGUMENTS:
! JRWORK     IN    I  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                       L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                       ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                       DU <<PAQUET>> DE MAILLES.
! LISNOE     IN    I  : LISTE COMPLETE DES NOEUDS A TRAITER.
! NBNOT      IN    I  : NOMBRE TOTAL DE NOEUDS A TRAITER.
! NBORDR     IN    I  : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                       STRUCTURE DE DONNEES RESULTAT.
! ORDINI     IN    I  : ORDRE INITIAL POUR LE CHARGEMENT CYCLIQUE
! NNOINI     IN    I  : NUMERO DU 1ER NOEUD DU <<PAQUET>> DE
!                       NOEUDS COURANT.
! NBNOP      IN    I  : NOMBRE DE NOEUDS DANS LE <<PAQUET>> DE
!                       NOEUDS COURANT.
! TSPAQ      IN    I  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE NOEUDS
!                       COURANT.
! NOMMET     IN    K16: NOM DE LA METHODE DE CALCUL DU CERCLE
!                       CIRCONSCRIT.
! NOMCRI     IN    K16: NOM DU CRITERE AVEC PLANS CRITIQUES.
! NOMMAI     IN    K8 : NOM UTILISATEUR DU MAILLAGE.
! CNSR       IN    K19: NOM DU CHAMP SIMPLE DESTINE A RECEVOIR LES
!                       RESULTATS.
!
! REMARQUE :
!  - LA TAILLE DU SOUS-PAQUET EST EGALE A LA TAILLE DU <<PAQUET>> DE
!    NOEUDS DIVISEE PAR LE NOMBRE DE NUMERO D'ORDRE (NBORDR).
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ki, l, jcnrd, jcnrl, ibidno
    integer(kind=8) :: iret, nbma, adrma, icesd, icesl, icesv
    integer(kind=8) :: inop, nunoe
    integer(kind=8) :: jtypma
    integer(kind=8) :: icmp, kwork, somnow, cnbno
    integer(kind=8) :: vali(2), jad, ima
    integer(kind=8) :: icodwo
!
    real(kind=8) :: coepre, vala, valb
    real(kind=8) :: coefpa
!
    character(len=8) :: chmat1, nommat
    character(len=10) :: optio
    character(len=19) :: chmat, cesmat, ncncin
    character(len=24) :: typma
    real(kind=8), pointer :: cnsv(:) => null()
!
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
!
    call jemarq()
!
!
! OBTENTION DES ADRESSES '.CESD', '.CESL' ET '.CESV' DU CHAMP SIMPLE
! DESTINE A RECEVOIR LES RESULTATS : DTAUM, ....
    if (.not. post) then
        call jeveuo(cnsr//'.CNSD', 'L', jcnrd)
        call jeveuo(cnsr//'.CNSL', 'E', jcnrl)
        call jeveuo(cnsr//'.CNSV', 'E', vr=cnsv)
!
! RECUPERATION DU COEFFICIENT DE PRE-ECROUISSAGE DONNE PAR L'UTILISATEUR
!
        call getvr8(' ', 'COEF_PREECROU', scal=coepre, nbret=iret)
!
! RECUPERATION MAILLE PAR MAILLE DU MATERIAU DONNE PAR L'UTILISATEUR
!
        call getvid(' ', 'CHAM_MATER', scal=chmat1, nbret=iret)
        chmat = chmat1//'.CHAMP_MAT'
        cesmat = '&&DTAUNO.CESMAT'
        call carces(chmat, 'ELEM', ' ', 'V', cesmat, &
                    'A', iret)
        call jeveuo(cesmat//'.CESD', 'L', icesd)
        call jeveuo(cesmat//'.CESL', 'L', icesl)
        call jeveuo(cesmat//'.CESV', 'L', icesv)
!
    end if
!
!  CONSTRUCTION DU VECTEUR : CONTRAINTE = F(NUMERO D'ORDRE) EN CHAQUE
! NOEUDS DU PAQUET DE MAILLES.
    l = 1
    cnbno = 0
    kwork = 0
    somnow = 0
!
    ncncin = '&&DTAUNO.CNCINV'
!
    if (.not. post) then
!
        call cncinv(nommai, [0], 0, 'V', ncncin)
        typma = nommai//'.TYPMAIL'
        call jeveuo(typma, 'L', jtypma)
!
    end if
!
    do inop = nnoini, nnoini+(nbnop-1)
!
        if (inop .gt. nnoini) then
            kwork = 1
            somnow = somnow+1
        end if
!
        cnbno = cnbno+1
        if ((l*int(nbnot/10.0d0)) .lt. cnbno) then
            l = l+1
        end if
!
! RECUPERATION DU NOM DU MATERIAU AFFECTE A LA MAILLE OU AUX MAILLES
! QUI PORTENT LE NOEUD COURANT.
        if (.not. post) then
            nunoe = lisnoe(inop)
            call jelira(jexnum(ncncin, nunoe), 'LONMAX', nbma)
            call jeveuo(jexnum(ncncin, nunoe), 'L', adrma)
!
            ki = 0
            optio = 'DOMA_NOEUD'
            do ima = 1, nbma
                call rnomat(icesd, icesl, icesv, ima, nomcri, &
                            adrma, jtypma, ki, optio, vala, &
                            valb, coefpa, nommat)
            end do
!
            call rcpare(nommat, 'FATIGUE', 'WOHLER', icodwo)
            if (icodwo .eq. 1) then
                call utmess('F', 'FATIGUE1_90', sk=nomcri(1:16))
            end if
!
            if (ki .eq. 0) then
                vali(1) = nunoe
                vali(2) = nbma
                call utmess('A', 'PREPOST5_10', ni=2, vali=vali)
            end if
        end if
!
!        call jerazo('&&DTAUNO.VECTNO', tneces, 1)
!
! C  IBIDNO JOUE LE ROLE DE IPG DNAS TAURLO
!
        if (post) then
            nommat = nommap
!
! RECUPERER LES COEEF DE CRITERES
!
            call recofa(nomcri, nommat, vala, valb, coefpa)
!
        end if
!
        ibidno = 1
!
!
! REMPACER PAR ACMATA
        call acgrdo(nbordr, ordini, kwork, somnow, jrwork, &
                    tspaq, ibidno, nommet, nommat, nomcri, &
                    vala, coefpa, nomfor, grdvie, forvie, &
                    forcri, valpar, vresu)
!
! AFFECTATION DES RESULTATS DANS UN CHAM_ELEM SIMPLE
        if (.not. post) then
            do icmp = 1, 24
                jad = 24*(nunoe-1)+icmp
                zl(jcnrl-1+jad) = .true.
                cnsv(jad) = vresu(icmp)
!
            end do
        end if
!
    end do
!
! MENAGE
!
    if (.not. post) then
        call detrsd('CHAM_ELEM_S', cesmat)
    end if
!
!
    call jedetr('&&DTAUNO.CNCINV')
!
    call jedema()
end subroutine
