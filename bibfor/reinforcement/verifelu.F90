! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine verifelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                    gammas, gammac, clacier, eys, typdiag, uc, &
                    dnsinf, dnssup, effm, effn, verif)
!______________________________________________________________________
!
!      VERIFELU

!      VERIFICATION D'UN TORSEUR D'EFFORTS (N,M)
!      SOLLICITANT UNE SECTION DE FERRAILLAGE CONNUE
!      PAR LA MÉTHODE DU DIAGRAMME D'INTERACTION
!      CRITERE = LIMITATION DES DEFORMATIONS

!      I TYPCO     CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I ALPHACC   COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON EN SUPRESSION
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I FACIER    LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I FBETON    RESISTANCE EN SUPRESSION DU BETON (CONTRAINTE)
!      I GAMMAS    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DES ACIERS
!      I GAMMAC    COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                  DE CALCUL DU BETON
!      I CLACIER   CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                     CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                     CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                     CLACIER = 3 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I EYS       MODULE D'YOUNG DE L'ACIER
!      I TYPDIAG   TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                     TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                     TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I DNSINF    DENSITE DE L'ACIER INFERIEUR
!      I DNSSUP    DENSITE DE L'ACIER SUPERIEUR
!      I EFFM      MOMENT FLECHISSANT A VERIFIER
!      I EFFN      EFFORT NORMAL A VERIFIER
!
!      O VERIF     VERIFICATION DU FERRAILLAGE
!                  = 0 --OK (TORSEUR A L'INTERIEUR DU DIAGRAMME D'INTERACTION)
!                  = 1 --PAS OK (TORSEUR A L'EXTERIEUR DU DIAGRAMME D'INTERACTION)
!
!______________________________________________________________________
!
!
    implicit none
!
#include "asterfort/dintelu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
!
    integer :: typco
    real(kind=8) :: alphacc
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: facier
    real(kind=8) :: fbeton
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    integer :: clacier
    real(kind=8) :: eys
    integer :: typdiag
    integer :: uc
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    real(kind=8) :: effm
    real(kind=8) :: effn
    integer :: verif

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: d, d0, dneg, d0neg
    real(kind=8) :: Esu, Xsup, Calc
    real(kind=8) :: piv_a, piv_b, piv_c
    real(kind=8) :: unite_pa
    integer :: N_ET, N_PC, N_EC, N_PCN, s
    integer :: ntot, ndemi
    logical :: COND_OK
    real(kind=8) :: nrd0, nrd1, mrd0, mrd1
    character(24) :: pnrd, pmrd
    real(kind=8), pointer :: mrd(:) => null()
    real(kind=8), pointer :: nrd(:) => null()

    !Dimensionnement des vecteurs

    pnrd = 'POINT_NRD'
    pmrd = 'POINT_MRD'

    if (typco .eq. 1) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'BAEL91'
        piv_a = 10.0E-3
        piv_b = 3.5E-3
        piv_c = 2.0E-3

    else if (typco .eq. 2) then
!       CALCUL DES PARAMETRES POUR CODIFICATION = 'EC2'

        if (uc .eq. 0) then
            unite_pa = 1.e-6
        elseif (uc .eq. 1) then
            unite_pa = 1.
        end if
        if (clacier .eq. 0) then
            piv_a = 0.9*2.5e-2
        else if (clacier .eq. 1) then
            piv_a = 0.9*5.e-2
        else
            piv_a = 0.9*7.5e-2
        end if
        piv_b = min(3.5E-3, 0.26*0.01+3.5*0.01*(((90.d0-fbeton*unite_pa)/100.d0)**4))
        piv_c = 2.0E-3
        if ((fbeton*unite_pa) .ge. (50.d0)) then
            piv_c = 0.2*0.01+0.0085*0.01*((fbeton*unite_pa-50.d0)**(0.53))
        end if

    end if

    Xsup = piv_b/piv_c
    Esu = piv_a
    N_ET = floor(Esu*1000)+1
    N_EC = ceiling(Xsup*100)+1

    d = ht-enrobi
    d0 = enrobs
    dneg = ht-enrobs
    d0neg = enrobi

    N_PC = ceiling((ht/d)*100)+1
    N_PCN = ceiling((ht/dneg)*100)+1
    ntot = N_ET+N_PC+N_EC+N_EC+N_PCN+N_ET
    ndemi = N_ET+N_PC+N_EC

    call wkvect(pnrd, ' V V R ', ntot, vr=nrd)
    call wkvect(pmrd, ' V V R ', ntot, vr=mrd)

    call dintelu(typco, alphacc, ht, bw, enrobi, enrobs, facier, fbeton, &
                 gammas, gammac, clacier, eys, typdiag, uc, &
                 dnsinf, dnssup, ntot, nrd, mrd)

    nrd0 = nrd(1)
    nrd1 = nrd(ndemi)
    if ((effn .ge. nrd0) .and. (effn .le. nrd1)) then
        COND_OK = .true.
    else
        COND_OK = .false.
        goto 998
    end if

    s = 1
    nrd0 = nrd(s)
    do while (nrd0 .lt. effn)
        s = s+1
        nrd0 = nrd(s)
    end do

    if (s .eq. 1) then
        mrd0 = mrd(1)
    else
        Calc = nrd(s)-nrd(s-1)
        if (abs(Calc) .gt. epsilon(Calc)) then
            mrd0 = ((mrd(s)-mrd(s-1))/(nrd(s)-nrd(s-1)))*(effn-nrd(s-1))+mrd(s-1)
        else
            mrd0 = 0.5*(mrd(s-1)+mrd(s))
        end if
    end if

    s = ndemi+1
    nrd1 = nrd(s)
    do while (nrd1 .gt. effn)
        s = s+1
        nrd1 = nrd(s)
    end do

    if (s .eq. ndemi) then
        mrd1 = mrd(ndemi)
    else
        Calc = nrd(s)-nrd(s-1)
        if (abs(Calc) .gt. epsilon(Calc)) then
            mrd1 = ((mrd(s)-mrd(s-1))/(nrd(s)-nrd(s-1)))*(effn-nrd(s-1))+mrd(s-1)
        else
            mrd1 = 0.5*(mrd(s-1)+mrd(s))
        end if
    end if

    if ((effm .le. mrd0) .and. (effm .ge. mrd1)) then
        COND_OK = .true.
    else
        COND_OK = .false.
    end if

998 continue

    if (COND_OK .eqv. (.true.)) then
        verif = 0
    else
        verif = 1
    end if

    call jedetr(pnrd)
    call jedetr(pmrd)

end subroutine
