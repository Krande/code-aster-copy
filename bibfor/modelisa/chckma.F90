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

subroutine chckma(nomu, dtol)
    implicit none
!-----------------------------------------------------------------------
!
!       ROUTINE DE VERIFICATION DU MAILLAGE :
!       1- RECHERCHE (ET ELIMINATION SI DEMANDEE) DES NOEUDS ORPHELINS
!       2- RECHERCHE (ET ELIMINATION SI DEMANDEE) DES MAILLES DOUBLES
!       3- RECHERCHE DES MAILLES APLATIES
!
!       IN,OUT : NOMU   NOM DU CONCEPT MAILLAGE PRODUIT PAR LA COMMANDE
!       IN     : DTOL   TOLERANCE POUR TESTER L APPLATISST DES MAILLES
!
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterfort/cncinv.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: nomu
    real(kind=8) :: dtol
!
!
! ----- DECLARATIONS
!
    integer(kind=8) :: iaconx, ilconx, ima, nbnm, nbnm2, it, jcoor, ifm, niv
    integer(kind=8) :: jdrvlc, jcncin, numail, imail
    integer(kind=8) :: iadr, iadr0, nbm, nbm0, iadtyp, nb200
    integer(kind=8) :: ja, jb, tabma(200), i, j, k1, k2, knso, kmdb, l
    character(len=8) :: noxa, noxb, tyma, mess_k8(3)
    character(len=24) :: ncncin
    real(kind=8) :: dm, dp, aplat, drap, mess_r(1)
    real(kind=8) :: xa, xb, ya, yb, za, zb
    character(len=24) :: cooval, connex, nsolo, mdoubl
    integer(kind=8) :: nbmail, nbnoeu
    integer(kind=8) :: insolo, imdoub, iatyma, nmdoub
    aster_logical :: indic, alarme, erreur
!
    call jemarq()
    call infniv(ifm, niv)
!
    cooval = nomu//'.COORDO    .VALE'
    connex = nomu//'.CONNEX         '
    call dismoi('NB_MA_MAILLA', nomu, 'MAILLAGE', repi=nbmail)
    call dismoi('NB_NO_MAILLA', nomu, 'MAILLAGE', repi=nbnoeu)
    if (nbnoeu == 0) goto 100
    call jeveuo(nomu//'.TYPMAIL', 'L', iatyma)
!
    call jeveuo(connex, 'L', iaconx)
    call jeveuo(jexatr(connex, 'LONCUM'), 'L', ilconx)
    call jeveuo(cooval, 'L', jcoor)
!
    ncncin = '&&CHCKMA.CONNECINVERSE  '
    call cncinv(nomu, [0], 0, 'V', ncncin)
    call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', jdrvlc)
    call jeveuo(jexnum(ncncin, 1), 'L', jcncin)
!
!     -----------------------------------------------------------
!     RECHERCHE DES NOEUDS ORPHELINS (ATTACHES A AUCUNE MAILLE)
!     A PARTIR DE LA CONNECTIVITE INVERSE RENVOYEE PAR CNCINV
!     -----------------------------------------------------------
    nsolo = '&&CHCKMA.NSOLO          '
    call wkvect(nsolo, 'V V I', max(nbnoeu, 1), insolo)
!
    it = 0
    knso = 0
    alarme = .false.
    nb200 = 0
    call utmess('I', 'MODELISA8_9')
    do ja = 1, nbnoeu
        iadr = zi(jdrvlc+ja-1)
        nbm = zi(jdrvlc+ja+1-1)-zi(jdrvlc+ja-1)
        if (nbm .gt. 200) then
            nb200 = 1
            noxa = int_to_char8(ja)
            mess_k8(1) = noxa
            call utmess('A', 'MODELISA8_10', nk=1, valk=mess_k8)
        end if
        do imail = 1, nbm
            numail = zi(jcncin+iadr-1+imail-1)
            if (numail .eq. 0) then
                knso = knso+1
                zi(insolo-1+knso) = ja
                noxa = int_to_char8(ja)
                mess_k8(1) = noxa
                call utmess('A', 'MODELISA8_11', nk=1, valk=mess_k8)
                alarme = .true.
            end if
        end do
    end do
    if (alarme) then
        call utmess('A', 'MODELISA4_6')
    end if
    if (nb200 .eq. 1) then
        call utmess('A', 'MODELISA4_7')
    end if
!
!
!     -----------------------------------------------------------
!     RECHERCHE DES MAILLES DOUBLES : C EST A DIRE LES MAILLES
!     DE NUMEROS DIFFERENTS QUI ONT LES MEMES NOEUDS EN SUPPORT :
!     POUR CHAQUE PREMIER NOEUD DE CHAQUE MAILLE, ON
!     REGARDE LES AUTRES MAILLES POSSEDANT CE NOEUD DANS LA
!     CONNECTIVITE INVERSE : LES CANDIDATS AU DOUBLON Y SONT
!     FORCEMMENT. CA EVITE UN ALGO EN N2.
!     -----------------------------------------------------------
!
    mdoubl = '&&CHCKMA.MDOUBLE'
    nmdoub = nbmail
    call wkvect(mdoubl, 'V V I', max(nmdoub, 1), imdoub)
!
!     BOUCLE SUR TOUTES LES MAILLES DU MAILLAGE
!
    it = 0
    kmdb = 0
    alarme = .false.
    erreur = .false.
    do ima = 1, nbmail
        nbnm = zi(ilconx-1+ima+1)-zi(ilconx+ima-1)
        iadr0 = zi(jdrvlc+zi(iaconx+1+it-1)-1)
        nbm0 = zi(jdrvlc+zi(iaconx+1+it-1)+1-1)-zi(jdrvlc+zi(iaconx+1+it-1)-1)
        i = 1
        do ja = 1, nbm0
            if (zi(jcncin+iadr0-1+ja-1) .ne. ima) then
                tabma(i) = zi(jcncin+iadr0-1+ja-1)
                i = i+1
            end if
!         -- POUR NE PAS DEBORDER DE  TABMA :
            if (i .gt. 199) goto 99
        end do
!
!     SI NBM0 DIFFERENT DE I : UN NOEUD DE LA MAILLE EST PRESENT
!     PLUSIEURS FOIS DANS LA CONNECTIVITE DE CELLE CI
        if (nbm0 .ne. i) then
            erreur = .true.
            noxa = int_to_char8(ima)
            iadtyp = iatyma-1+ima
            call jenuno(jexnum('&CATA.TM.NOMTM', zi(iadtyp)), tyma)
            mess_k8(1) = noxa
            mess_k8(2) = tyma
            call utmess('A', 'MODELISA8_12', nk=2, valk=mess_k8)
        end if
        nbm0 = i
!
!     TABMA CONTIENT LA LISTE DES MAILLES (HORS IMA) QUI
!     CONTIENNENT LE PREMIER NOEUD DE IMA
!
        if (nbnm .gt. 1) then
            do i = 1, nbm0-1
                nbnm2 = zi(ilconx-1+tabma(i)+1)-zi(ilconx+tabma(i)-1)
!
!     COMPARAISON DES NOEUDS DE IMA AVEC CEUX DES MAILLES DE TABMA
!     SI LES CARDINAUX SONT DEJA DIFFERENTS (NBNM) : ON SAUTE
!
                if ((nbnm2 .eq. nbnm) .and. (tabma(i) .lt. ima)) then
                    do j = 1, nbnm2
                        k1 = zi(iaconx-1+zi(ilconx+tabma(i)-1)+j-1)
                        indic = .false.
                        do l = 1, nbnm
                            k2 = zi(iaconx-1+zi(ilconx+ima-1)+l-1)
                            if (k1 .eq. k2) indic = .true.
                        end do
                        if (.not. indic) goto 102
                    end do
                    kmdb = kmdb+1
                    if (kmdb .gt. nmdoub) then
                        nmdoub = 2*nmdoub
                        call juveca(mdoubl, nmdoub)
                    end if
                    call jeveuo(mdoubl, 'E', imdoub)
                    zi(imdoub-1+kmdb) = tabma(i)
                    noxa = int_to_char8(ima)
                    noxb = int_to_char8(tabma(i))
                    iadtyp = iatyma-1+ima
                    call jenuno(jexnum('&CATA.TM.NOMTM', zi(iadtyp)), tyma)
                    mess_k8(1) = noxa
                    mess_k8(2) = noxb
                    mess_k8(3) = tyma
                    call utmess('A', 'MODELISA8_13', nk=3, valk=mess_k8)
                    alarme = .true.
                end if
!
102             continue
            end do
!
        else if (nbm0 .gt. 1) then
            noxa = int_to_char8(ima)
            mess_k8(1) = noxa
            call utmess('A', 'MODELISA8_14', nk=1, valk=mess_k8)
        end if
!
99      continue
        it = it+nbnm
    end do
    if (alarme) then
        call utmess('A', 'MODELISA4_8')
    end if
!
!     -----------------------------------------------------------
!     CALCUL POUR CHAQUE MAILLE DU RAPPORT MINIMUM ENTRE LA PLUS
!     PETITE ARRETE ET LA PLUS GRANDE POUR REPERER LES ELEMENTS
!     TRES APLATIS VOIRE DEGENERES. LE RAPPORT MIN TOLERE EST :
!     DTOL = 1 POURCENT
!     -----------------------------------------------------------
!
    it = 0
    alarme = .false.
    do ima = 1, nbmail
        nbnm = zi(ilconx-1+ima+1)-zi(ilconx+ima-1)
        dm = r8maem()
        dp = r8miem()
        if (nbnm .gt. 1) then
!
            do ja = 1, nbnm-1
                do jb = ja+1, nbnm
                    xa = zr(jcoor-1+3*(zi(iaconx+ja+it-1)-1)+1)
                    ya = zr(jcoor-1+3*(zi(iaconx+ja+it-1)-1)+2)
                    za = zr(jcoor-1+3*(zi(iaconx+ja+it-1)-1)+3)
                    xb = zr(jcoor-1+3*(zi(iaconx+jb+it-1)-1)+1)
                    yb = zr(jcoor-1+3*(zi(iaconx+jb+it-1)-1)+2)
                    zb = zr(jcoor-1+3*(zi(iaconx+jb+it-1)-1)+3)
                    aplat = (xa-xb)**2+(ya-yb)**2+(za-zb)**2
                    if (aplat .lt. dm) dm = aplat
                    if (aplat .gt. dp) dp = aplat
                end do
            end do
            if (dp .gt. 0.d0) then
                drap = sqrt(dm/dp)
                if (drap .lt. dtol) then
                    alarme = .true.
                    noxa = int_to_char8(ima)
                    iadtyp = iatyma-1+ima
                    call jenuno(jexnum('&CATA.TM.NOMTM', zi(iadtyp)), tyma)
                    mess_k8(1) = noxa
                    mess_k8(2) = tyma
                    mess_r(1) = drap
                    call utmess('A', 'MODELISA8_15', nk=2, valk=mess_k8, &
                                nr=1, valr=mess_r)
                end if
            end if
!
        end if
        it = it+nbnm
    end do
    if (alarme) then
        call utmess('A', 'MODELISA4_9', sr=dtol)
    end if
!     ON ARRETE EN ERREUR SUR MAILLE DEGENEREE
    if (erreur) then
        call utmess('F', 'MODELISA4_10')
    end if
!
!     -----------------------------------------------------------
!     MENAGE DANS LE MAILLAGE : ON DETRUIT NOEUDS ORPHELINS ET
!                               MAILLES DOUBLES
!     -----------------------------------------------------------
!
!     CA RESTE A FAIRE ...
!     LES NOEUDS ORPHELINS SONT RANGES DANS &&CHCKMA.NSOLO(1:KNSO)
!     LES MAILLES DOUBLES SONT RANGEES DANS &&CHCKMA.MDOUBLE(1:KMDB)
!
100 continue
    call jedetr('&&CHCKMA.NSOLO          ')
    call jedetr('&&CHCKMA.MDOUBLE')
    call jedetr('&&CHCKMA.CONNECINVERSE  ')
    call jedema()
end subroutine
