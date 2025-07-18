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
subroutine renrfa(nomfor, valgrd, nrupt, icodre)
!
    implicit none
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/fonbpa.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: nomfor
    integer(kind=8) :: icodre
    real(kind=8) :: valgrd, nrupt
!     OBTENTION DU NOMBRE DE CYLE A LA RUPTURE A PARTIR D'UN GRANDEUR
!     EQUIVALENT ET UNE COURBE DE GRANDEUR_DUREE DE VIE
!     ARGUMENTS D'ENTREE:
!        NOMMAT : NOM UTILISATEUR DU MATERIAU
!        PHENOM : NOM DU PHENOMENE
!        NBPAR  : NOMBRE DE PARAMETRES DANS NOMPAR ET VALPAR
!        NOMPAR : NOMS DES PARAMETRES(EX: TEMPERATURE )
!        VALPAR : VALEURS DES PARAMETRES
!        NBRES  : NOMBRE DE RESULTATS
!        NOMRES : NOM DES RESULTATS (EX: E,NU,... )
!                 TELS QU'IL FIGURENT DANS LA COMMANDE MATERIAU
!     ARGUMENTS DE SORTIE:
!     VALRES : VALEURS DES RESULTATS APRES RECUPERATION ET INTERPOLATION
!     ICODRE : POUR CHAQUE RESULTAT, 0 SI ON A TROUVE, 1 SINON
!
!
!
!
    integer(kind=8) :: jprol, np, nbmx, ndat, i
    character(len=24) :: chnom, cbid
    character(len=8) :: nompf, typfon, nompfi
    real(kind=8) :: lnf(20), grd(20), lnrup, grdmax
    data lnf/1.d7, 5.d6, 1.d6, 5.d5, 2.d5, 1.d5, 5.d4, 2.d4,&
     &           1.2d4, 1.d4, 5.d3, 2.d3, 1.d3, 5.d2,&
     &           2.d2, 1.d2, 5.d1, 2.d1, 1.d1, 1.d0/
!
    call jemarq()
!
! RECUPERER LE NOM DE PARAMETRE GRANDEUR_EQUIVALEN
!
    nbmx = 35
    chnom(20:24) = '.PROL'
    chnom(1:19) = nomfor
    nrupt = 0.d0
!
    call jeveuo(chnom, 'L', jprol)
    typfon = zk24(jprol-1+1) (1:8)
!
    call fonbpa(nomfor, zk24(jprol), cbid, nbmx, np, &
                nompf)
!
    if (typfon .eq. 'FONCTION') then
!
! A PARTIR DE GRANDEUR EQUI, ON RESOUT LE NOMBRE DE CYCLE AVEC
! LA FONCTION DU TYPE: NBRUPT = F(VALGRD)
        call fointe('F', nomfor, np, [nompf], [valgrd], &
                    nrupt, icodre)
!
    else
! A PARTIR DE GRANDEUR EQUI, ON RESOUT LE NOMBRE DE CYCLE AVEC
! LA FONCTION DU TYPE:  VALGRD = F(NBRUPT)
! ON DEFINIT  FONCT = FONCT(NBRUPT) =  F(NBRUPT)- VALGRD
! ON VA RESOUDRE L'EQUATION POUR TROUVER NBRUPT TEL QUE 1<NBRUPT<10^7
!
        ndat = 20
        call fonbpa(nomfor, zk24(jprol), cbid, nbmx, np, &
                    nompfi)
!
! VERIFIER QUE LA FORMULE A LA VARIABLE NRUPT = N_F
        if ((nompfi .ne. 'NBRUP') .or. (np .ne. 1)) then
            call utmess('F', 'FATIGUE1_93')
        end if
!
! VERIFIER QUE  VALGRD < GRDMAX
        call fointe('F', nomfor, np, [nompfi], [lnf(20)], &
                    grdmax, icodre)
!
!
        if (valgrd .gt. grdmax) then
            call utmess('F', 'FATIGUE1_94')
        end if
!
        call fointe('F', nomfor, np, [nompfi], [lnf(1)], &
                    grd(1), icodre)
        do i = 2, ndat
            call fointe('F', nomfor, np, [nompfi], [lnf(i)], &
                        grd(i), icodre)
            if (grd(i) .gt. valgrd) then
!                   NRUPT2 = LNF(I)+
!      &           (LNF(I-1)-LNF(I))/(GRD(I-1)-GRD(I))*(VALGRD-GRD(I))
!
                lnrup = log( &
                        lnf(i))+(log(lnf(i-1))-log(lnf(i)))/(log(grd(i-1))-log(grd(i)))*(log(&
                        &valgrd)-log(grd(i)) &
                        )
                nrupt = exp(lnrup)
!
                goto 20
!
            end if
        end do
    end if
!
20  continue
!
    call jedema()
end subroutine
