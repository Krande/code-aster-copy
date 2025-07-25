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

subroutine xajuls(noma, nbma, cnslt, cnsln, jconx1, &
                  jconx2, clsm, typdis, critlst)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/conare.h"
#include "asterfort/dismoi.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/xajuls_stop.h"
!
    integer(kind=8) :: nbma, jconx1, jconx2, clsm
    character(len=8) :: noma
    character(len=16) :: typdis
    character(len=19) :: cnslt, cnsln
    real(kind=8), optional :: critlst
!
! person_in_charge: patrick.massin at edf.fr
!
!     ------------------------------------------------------------------
!     XFEM : REAJUSTEMENT DES LEVEL SETS (BOOK III 06/02/04)
!     -        ---
!     BUT : ON MODIFIE LES VALEURS DE LS AUX NOEUDS SI TROP
!           PROCHES DE 0 POUR EVITER LES ERREURS D'INTEGRATION
!
!    ENTREE :
!              IFM    :   FICHIER D'IMPRESSION
!              NOMA   :   OBJET MAILLAGE
!              NBMA   :   NOMBRE DE MAILLES DU MAILLAGE
!              CNSLN  :   LEVEL-SET NORMALE  (PLAN DE LA FISSURE)
!              CNSLT  :   LEVEL-SET TANGENTE (TRACE DE LA FISSURE)
!       JCONX1,JCONX2 :   INDICES DE LA CONNECTIVITE
!
!    SORTIE :
!              CNSLN  :   LEVEL-SET NORMALE
!              CNSLT  :   LEVEL-SET TANGENTE
!              CLSM   :   NOMBRE DE LEVEL SETS MODIFIEES
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: jma, ima, itypma, ar(12, 3), nbar, ia
    integer(kind=8) :: na, nb, nm, nunoa, nunob, nunom
    integer(kind=8) :: nmaabs, ndime, ndim
    real(kind=8) :: d1, lsna, lsnb, crilsn, lsta, lstb, crilst, d2, r8pre, crlst
    real(kind=8) :: lsnm, lstm, lsnmax, lstmax, penal, d3, fit_to_vertex(2)
    character(len=19) :: mai
    character(len=8) :: typma
    real(kind=8), pointer :: lnsv(:) => null()
    real(kind=8), pointer :: ltsv(:) => null()
    aster_logical :: ajust
!
    parameter(fit_to_vertex=(/1.d-6, 1d-6/), crlst=1.d-6, penal=0.05)
!
!-----------------------------------------------------------------------
!     DEBUT
!-----------------------------------------------------------------------
    call jemarq()
!
    call jeveuo(cnsln//'.CNSV', 'E', vr=lnsv)
    call jeveuo(cnslt//'.CNSV', 'E', vr=ltsv)
!
    r8pre = r8prem()
    d2 = 999.d0
    if (present(critlst)) then
        crilst = critlst
    else
        crilst = crlst
    end if
!
!     COMPTEUR DES LSN ET LST MODIFIÉES
    clsm = 0
    mai = noma//'.TYPMAIL'
    call jeveuo(mai, 'L', jma)
!
!     RECUPERATION DE LA DIMENSION DE L'ESPACE
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndime)
!
!     BOUCLE SUR TOUTES LES MAILLES DU MAILLAGE
!
!     BOOLEEN POUR SAVOIR SI UN AJSTEMENT A ETE EFFECTUE DANS LE CAS QUADRATIQUE
100 continue
    ajust = .false.
!
    do ima = 1, nbma
        nmaabs = ima
        itypma = zi(jma-1+ima)
        call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
!
!       UTILISATION DE FIT-TO-VERTX POUR LE MOMENT:
        if (ismali(typma)) then
            crilsn = fit_to_vertex(1)
        else
            crilsn = fit_to_vertex(2)
        end if
!       RECUPERATION DE LA DIMENSION TOPOLOGIQUE DE L'ELEMENT
        call dismoi('DIM_TOPO', typma, 'TYPE_MAILLE', repi=ndim)
!
!       LES ELEMENTS DE BORD NE SONT PAS TRAITES
        if (ndim .lt. ndime) goto 200
!
!       BOUCLE SUR LES ARETES DE LA MAILLE VOLUMIQUE
        call conare(typma, ar, nbar)
        lsnmax = 0.d0
        lstmax = 0.d0
!
!       BOUCLE POUR RECUPERER LE MAX DE ABS(LSN) SUR UNE MAILLE
!       (UNIQUEMENT LES NOEUDS SOMMETS) ET FIT TO VERTEX LES NOEUDS
!       SOMMETS
!
        do ia = 1, nbar
            na = ar(ia, 1)
            nb = ar(ia, 2)
            nunoa = zi(jconx1-1+zi(jconx2+nmaabs-1)+na-1)
            nunob = zi(jconx1-1+zi(jconx2+nmaabs-1)+nb-1)
            lsna = lnsv((nunoa-1)+1)
            lsnb = lnsv((nunob-1)+1)
            lsta = ltsv((nunoa-1)+1)
            lstb = ltsv((nunob-1)+1)
            if (abs(lsna) .gt. lsnmax) lsnmax = abs(lsna)
            if (abs(lsnb) .gt. lsnmax) lsnmax = abs(lsnb)
            if (abs(lsta) .gt. lstmax) lstmax = abs(lsta)
            if (abs(lstb) .gt. lstmax) lstmax = abs(lstb)
!
!         REAJUSTEMENT DE LA LEVEL SET NORMALE AUX NOEUDS SOMMETS, QUAND
!         LA VALEUR D'UN LSN DIVISE PAR LA DIFFERENCE DES VALEURS AUX
!         DEUX EXTREMITES SONT INFERIEURES D'UN CERTAINE NOMBRE, ON MET
!         LES LSN A ZERO.
          if (abs(lsna-lsnb) .gt. r8pre .and. (typdis .ne. 'COHESIF' .or. lsna*lsnb .lt. 0.d0)) then
                d1 = lsna/(lsna-lsnb)
                if (abs(d1) .le. crilsn) then
!              REAJUSTEMENT DE LSNA
                    lnsv((nunoa-1)+1) = 0.d0
                    clsm = clsm+1
                end if
                if (abs(d1-1.d0) .le. crilsn) then
!              REAJUSTEMENT DE LSNB
                    lnsv((nunob-1)+1) = 0.d0
                    clsm = clsm+1
                end if
                if (.not. ismali(typma)) then
                    nm = ar(ia, 3)
                    nunom = zi(jconx1-1+zi(jconx2+nmaabs-1)+nm-1)
                    lsnm = lnsv((nunom-1)+1)
                    if (abs(lsnm/(lsna-lsnb)) .le. crilsn) then
                        lnsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                    end if
                end if
            end if
!
!         REAJUSTEMENT DE LA LEVEL SET TANGENTE AUX NOEUDS SOMMETS,
!         QUAND LA VALEUR D'UN LST DIVISE PAR LA DIFFERENCE DES VALEURS
!         AUX DEUX EXTREMITES SONT INFERIEURES D'UN CERTAINE NOMBRE, ON
!         MET LES LST A ZERO.
!
            if (abs(lsta-lstb) .gt. r8pre) then
                d1 = lsta/(lsta-lstb)
                if (abs(d1) .le. crilst) then
!              REAJUSTEMENT DE LSTA
                    ltsv((nunoa-1)+1) = 0.d0
                    clsm = clsm+1
                end if
                if (abs(d1-1.d0) .le. (crilst)) then
!              REAJUSTEMENT DE LSTB
                    ltsv((nunob-1)+1) = 0.d0
                    clsm = clsm+1
                end if
                if (.not. ismali(typma)) then
                    nm = ar(ia, 3)
                    nunom = zi(jconx1-1+zi(jconx2+nmaabs-1)+nm-1)
                    lstm = ltsv((nunom-1)+1)
                    if (abs(lstm/(lsta-lstb)) .le. crilst) then
                        ltsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                    end if
                end if
            end if
!
        end do
!
        if (.not. ismali(typma)) then
!
            do ia = 1, nbar
                na = ar(ia, 1)
                nb = ar(ia, 2)
                nunoa = zi(jconx1-1+zi(jconx2+nmaabs-1)+na-1)
                nunob = zi(jconx1-1+zi(jconx2+nmaabs-1)+nb-1)
!
                nm = ar(ia, 3)
                nunom = zi(jconx1-1+zi(jconx2+nmaabs-1)+nm-1)
!
                lsna = lnsv((nunoa-1)+1)
                lsnb = lnsv((nunob-1)+1)
                lsta = ltsv((nunoa-1)+1)
                lstb = ltsv((nunob-1)+1)
!
                lsnm = lnsv((nunom-1)+1)
                lstm = ltsv((nunom-1)+1)
!
!!!!!!!!!!!! TRAITEMENT DE LSN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            REAJUSTEMENT DES CONFIGURATIONS RENTRANTES
                if (lsna .eq. 0.d0 .and. lsnb .eq. 0.d0 .and. lsnm .ne. 0.d0) then
                    d1 = lsnm/lsnmax
                    if (abs(d1) .le. penal) then
!             REAJUSTEMENT A ZERO DE LSNM AUX NOEUDS MILIEUX,QUAND LA
!             VALEUR DE LSNM EST INFERIEURE A 1% LSNMAX ET QUE LES
!             EXTREMITES DE L'ARETE ON LSN=0
                        lnsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                        ajust = .true.
                    end if
                else if ((lsna*lsnm) .lt. 0.d0 .and. (lsnb*lsnm) .lt. 0.d0) then
                    d1 = lsna/lsnmax
                    d2 = lsnb/lsnmax
                    d3 = lsnm/lsnmax
                    if ((abs(d1) .le. penal) .and. (abs(d2) .le. penal) .and. &
                        (abs(d3) .le. penal)) then
!             REAJUSTEMENT A ZERO DE LSNM ET D'UNE LSNM D'EXTREMITE LORSQUE
!             LA LSN CHANGE DEUX FOIS DE SIGNE SUR L'ARETE ET QUE LES LSN
!             SONT INFERIEURES A 1% LSNMAX
                        lnsv((nunom-1)+1) = 0.d0
                        if (d1 .gt. d2) then
                            lnsv((nunob-1)+1) = 0.d0
                        else
                            lnsv((nunoa-1)+1) = 0.d0
                        end if
                        clsm = clsm+2
                        ajust = .true.
                    end if
                else if (lsna .eq. 0.d0 .and. (lsnb*lsnm) .lt. 0.d0) then
                    d3 = lsnm/lsnmax
                    if (abs(d3) .le. penal) then
!             REAJUSTEMENT A ZERO DE LSNM LORSQUE LA LSN EST NULLE EN A ET
!             CHANGE DE SIGNE SUR L'ARETE ET QUE LES LSN SONT INFERIEURES
!             A 1% LSNMAX
                        lnsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                        ajust = .true.
                    end if
                else if ((lsna*lsnm) .lt. 0.d0 .and. lsnb .eq. 0.d0) then
                    d3 = lsnm/lsnmax
                    if (abs(d3) .le. penal) then
!             REAJUSTEMENT A ZERO DE LSNM LORSQUE LA LSN EST NULLE EN B ET
!             CHANGE DE SIGNE SUR L'ARETE ET QUE LES LSN SONT INFERIEURES
!             A 1% LSNMAX
                        lnsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                        ajust = .true.
                    end if
                else if ((lsna*lsnb) .gt. 0.d0 .and. lsnm .eq. 0.d0) then
                    d1 = lsna/lsnmax
                    d2 = lsnb/lsnmax
                    if ((abs(d1) .le. penal) .and. (abs(d2) .le. penal)) then
!             REAJUSTEMENT A ZERO D'UNE LSN D'EXTREMITE LORSQUE LA LSN
!             EST NULLE EN M ET DE MEME SIGNE AUX EXTREMITES DE L'ARETE ET
!             QUE LES LSN SONT INFERIEURES A 1% LSNMAX
                        if (d1 .gt. d2) then
                            lnsv((nunob-1)+1) = 0.d0
                        else
                            lnsv((nunoa-1)+1) = 0.d0
                        end if
                        clsm = clsm+1
                        ajust = .true.
                    end if
                end if
!
!!!!!!!!!!!!! TRAITEMENT DE LST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!            REAJUSTEMENT DES CONFIGURATIONS RENTRANTES
                if (lsta .eq. 0.d0 .and. lstb .eq. 0.d0 .and. lstm .ne. 0.d0) then
!                   dans le cas ou lstmax est nul, on court-circuite
!                   (sous certaines conditions tres particulieres)
!                   l'iteration de la boucle sur les aretes pour eviter
!                   une division par zero
                    if (abs(lstmax) .lt. r8pre) then
                        call xajuls_stop(noma, cnslt, jconx1, jconx2, ima)
                        cycle
                    end if
                    d1 = lstm/lstmax
                    if (abs(d1) .le. penal) then
!             REAJUSTEMENT A ZERO DE LSNM AUX NOEUDS MILIEUX,QUAND LA
!             VALEUR DE LSNM EST INFERIEURE A 1% LSNMAX ET QUE LES
!             EXTREMITES DE L'ARETE ON LSN=0
                        ltsv((nunom-1)+1) = 0.d0
                        clsm = clsm+1
                        ajust = .true.
                    end if
                else if ((lsta*lstm) .lt. 0.d0 .and. (lstb*lstm) .lt. 0.d0) then
                    d1 = lsta/lstmax
                    d2 = lstb/lstmax
                    d3 = lstm/lstmax
                    if ((abs(d1) .le. penal) .and. (abs(d2) .le. penal) .and. &
                        (abs(d3) .le. penal)) then
!             REAJUSTEMENT A ZERO DE LSNM ET D'UNE LSNM D'EXTREMITE LORSQUE
!             LA LSN CHANGE DEUX FOIS DE SIGNE SUR L'ARETE ET QUE LES LSN
!             SONT INFERIEURES A 1% LSNMAX
                        ltsv((nunom-1)+1) = 0.d0
                        ltsv((nunob-1)+1) = 0.d0
                        ltsv((nunoa-1)+1) = 0.d0
                        clsm = clsm+3
                        ajust = .true.
                    end if
                else if (lsta .eq. 0.d0 .and. (lstb*lstm) .lt. 0.d0) then
                    d2 = lstb/lstmax
                    d3 = lstm/lstmax
                    if ((abs(d3) .le. penal) .and. (abs(d2) .le. penal)) then
!             REAJUSTEMENT A ZERO DE LSNM LORSQUE LA LSN EST NULLE EN A ET
!             CHANGE DE SIGNE SUR L'ARETE ET QUE LES LSN SONT INFERIEURES
!             A 1% LSNMAX
                        ltsv((nunom-1)+1) = 0.d0
                        ltsv((nunob-1)+1) = 0.d0
                        clsm = clsm+2
                        ajust = .true.
                    end if
                else if ((lsta*lstm) .lt. 0.d0 .and. lstb .eq. 0.d0) then
                    d1 = lstm/lstmax
                    d3 = lstm/lstmax
                    if ((abs(d1) .le. penal) .and. (abs(d3) .le. penal)) then
!             REAJUSTEMENT A ZERO DE LSNM LORSQUE LA LSN EST NULLE EN B ET
!             CHANGE DE SIGNE SUR L'ARETE ET QUE LES LSN SONT INFERIEURES
!             A 1% LSNMAX
                        ltsv((nunom-1)+1) = 0.d0
                        ltsv((nunoa-1)+1) = 0.d0
                        clsm = clsm+2
                        ajust = .true.
                    end if
                else if ((lsta*lstb) .gt. 0.d0 .and. lstm .eq. 0.d0) then
                    d1 = lsta/lstmax
                    d2 = lstb/lstmax
                    if ((abs(d1) .le. penal) .and. (abs(d2) .le. penal)) then
!             REAJUSTEMENT A ZERO D'UNE LSN D'EXTREMITE LORSQUE LA LSN
!             EST NULLE EN M ET DE MEME SIGNE AUX EXTREMITES DE L'ARETE ET
!             QUE LES LSN SONT INFERIEURES A 1% LSNMAX
                        ltsv((nunob-1)+1) = 0.d0
                        ltsv((nunoa-1)+1) = 0.d0
                        clsm = clsm+2
                        ajust = .true.
                    end if
                end if
!
            end do
!
        end if
!
200     continue
    end do
!
!     TANT QUE DES AJUSTEMENTS SONT EFFECTUES POUR LES LEVEL SET QUADRATIQUES,
!     ON REITERE LA PROCEDURE AFIN DE N'AVOIR AUCUN PROBLEME DE COHERENCE SUITE
!     AUX AJUSTEMENTS
    if (ajust) go to 100
!
    if (clsm .gt. 0) then
        call utmess('A', 'XFEM_63', si=clsm)
    end if
!
!-----------------------------------------------------------------------
!     FIN
!-----------------------------------------------------------------------
    call jedema()
end subroutine
