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

subroutine gveri3(chfond, taillr, config, lnoff, liss, &
                  ndeg, trav1, trav2, trav3, &
                  typdis)
!
!
    implicit none
!
!     ------------------------------------------------------------------
!
! FONCTION REALISEE:     DANS LE CADRE DE X-FEM et FEM
!
!     - METHODES THETA_LAGRANGE
!
!         POUR CHAQUE NOEUD DU FOND DE FISSURE GAMM0 ON RECUPERE
!         LE DOUBLET (RINF, RSUP )
!
!     - METHODE THETA_LEGENDRE
!
!         POUR CHAQUE NOEUD DU FOND DE FISSURE GAMM0 ON RECUPERE
!         LE TRIPLET ( DEGRE DES POLYNOMES DE LEGENDRE, RINF, RSUP )
!
!     ------------------------------------------------------------------
! ENTREE:
!        CHFOND : NOMS DES NOEUDS
!        TAILLR : TAILLES DE MAILLES CONNECTEES AUX NOEUDS
!        CONFIG : CONFIGURATION DE LA FISSURE EN FEM
!        LNOFF  : NOMBRE DE NOEUD DE GAMM0
!        LISS   : TYPE DE LISSAGE
!        NDEG   : DEGRE DES POLYNOMES DE LEGENDRE
!
! SORTIE:
!        RINF          ( OBJET TRAV1 )
!        RSUP          ( OBJET TRAV2 )
!        MODULE(THETA) ( OBJET TRAV3 )
!     ------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/glegen.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=16) :: typdis
    character(len=24) :: trav0, trav1, trav2, trav3, chfond, absgam, taillr, liss
    character(len=8) :: config, nompar(1), rinff, rsupf
!
    integer(kind=8) :: lnoff, ndeg, nbre, nr, nrf, nbpar, i, j
    integer(kind=8) :: iadrt0, iadrt1, iadrt2, iadrt3, ifon, iadabs, ier
    integer(kind=8) :: iatmno
!
    real(kind=8) :: maxtai, mintai, rinf, rsup, xl, valpar(1), valres
    real(kind=8) :: valr(2)
!
!
    call jemarq()
!
    if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') .or. (liss .eq. 'MIXTE')) then
        nbre = lnoff-1
    else
        nbre = ndeg
    end if
!
! ALLOCATION DE 3 OBJETS DE TRAVAIL
!
    if (typdis .ne. 'COHESIF') then
        trav0 = '&&VERIFG.GAM0'//'           '
        trav1 = '&&VERIFG.RINF'//'           '
        trav2 = '&&VERIFG.RSUP'//'           '
        call wkvect(trav0, 'V V K8', lnoff, iadrt0)
        call wkvect(trav1, 'V V R', lnoff, iadrt1)
        call wkvect(trav2, 'V V R', lnoff, iadrt2)
    end if
    trav3 = '&&VERIFG.THET'//'           '
    call wkvect(trav3, 'V V R', (nbre+1)*lnoff, iadrt3)
    if (typdis .eq. 'COHESIF') goto 98
!
    call getvr8('THETA', 'R_INF', iocc=1, scal=rinf, nbret=nr)
    call getvr8('THETA', 'R_SUP', iocc=1, scal=rsup, nbret=nr)
    if (nr .ne. 0 .and. rsup .le. rinf) then
        call utmess('F', 'RUPTURE1_6')
    end if
    call getvid('THETA', 'R_INF_FO', iocc=1, scal=rinff, nbret=nrf)
    call getvid('THETA', 'R_SUP_FO', iocc=1, scal=rsupf, nbret=nrf)
!     RECUPERATION DE RINF ET DE RSUP DANS LA SD FOND_FISS
    if (nr .eq. 0 .and. nrf .eq. 0) then
        if (config .eq. 'DECOLLEE') then
            call utmess('F', 'RUPTURE1_7')
        end if
        call jeveuo(taillr, 'L', iatmno)
        maxtai = 0.d0
        mintai = zr(iatmno)
        do j = 1, lnoff
            maxtai = max(maxtai, zr(iatmno-1+j))
            mintai = min(mintai, zr(iatmno-1+j))
        end do
        rinf = 2*maxtai
        rsup = 4*maxtai
        valr(1) = rinf
        valr(2) = rsup
        call utmess('I', 'RUPTURE1_5', nr=2, valr=valr)
        valr(1) = mintai
        valr(2) = maxtai
        if (maxtai .gt. 2*mintai) then
            call utmess('A', 'RUPTURE1_16', nr=2, valr=valr)
        end if
    end if
!
    call jeveuo(chfond, 'L', ifon)
    absgam = '&&GVERI3.TEMP     .ABSCU'
    call wkvect(absgam, 'V V R', lnoff, iadabs)
    do i = 1, lnoff
        zr(iadabs-1+(i-1)+1) = zr(ifon-1+4*(i-1)+4)
    end do
    xl = zr(iadabs-1+(lnoff-1)+1)
!
    if ((liss .ne. 'LAGRANGE') .and. (liss .ne. 'LAGRANGE_NO_NO') .and. (liss .ne. 'MIXTE')) then
!
! METHODE THETA_LEGENDRE
!
        do j = 1, lnoff
            zk8(iadrt0+j-1) = 'PTFONFIS'
            if (nrf .ne. 0) then
                nbpar = 1
                nompar(1) = 'ABSC'
                valpar(1) = zr(iadabs+j-1)
                call fointe('FM', rinff, nbpar, nompar, valpar, &
                            valres, ier)
                zr(iadrt1+j-1) = valres
                call fointe('FM', rsupf, nbpar, nompar, valpar, &
                            valres, ier)
                zr(iadrt2+j-1) = valres
                if (zr(iadrt2+j-1) .le. zr(iadrt1+j-1)) then
                    call utmess('F', 'RUPTURE1_6')
                end if
            else
                zr(iadrt1+j-1) = rinf
                zr(iadrt2+j-1) = rsup
            end if
        end do
!
        call glegen(nbre, lnoff, xl, absgam, zr(iadrt3))
!
    else if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') &
             .or. (liss .eq. 'MIXTE')) then
!
! METHODES THETA_LAGRANGE
!
        do j = 1, lnoff
            zk8(iadrt0+j-1) = 'PTFONFIS'
            if (nrf .ne. 0) then
                nbpar = 1
                nompar(1) = 'ABSC'
                valpar(1) = zr(iadabs+j-1)
                call fointe('FM', rinff, nbpar, nompar, valpar, &
                            valres, ier)
                zr(iadrt1+j-1) = valres
                call fointe('FM', rsupf, nbpar, nompar, valpar, &
                            valres, ier)
                zr(iadrt2+j-1) = valres
                if (zr(iadrt2+j-1) .le. zr(iadrt1+j-1)) then
                    call utmess('F', 'RUPTURE1_6')
                end if
            else
                zr(iadrt1+j-1) = rinf
                zr(iadrt2+j-1) = rsup
            end if
        end do
!
    end if
!
    call jedetr(absgam)
    call jedetr(trav0)
98  continue
!
    call jedema()
end subroutine
