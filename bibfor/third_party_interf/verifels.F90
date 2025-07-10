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

subroutine verifels(cequi, ht, bw, enrobi, enrobs, &
                    scmaxi, scmaxs, ssmax, uc, &
                    dnsinf, dnssup, effm, effn, verif)
!______________________________________________________________________
!
!      VERIFELS

!      VERIFICATION D'UN TORSEUR D'EFFORTS (N,M)
!      SOLLICITANT UNE SECTION DE FERRAILLAGE CONNUE
!      PAR LA MÉTHODE DU DIAGRAMME D'INTERACTION
!      CRITERE = LIMITATION DES CONTRAINTES

!      I CEQUI     COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I SCMAXI    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE INF
!      I SCMAXS    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE SUP
!      I SSMAX     CONTRAINTE MAXI DE L'ACIER DE FLEXION
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
#include "extern/dintels.h"
#include "asterfort/wkvect.h"
#include "asterfort/jedetr.h"
!
    real(kind=8) :: cequi
    real(kind=8) :: ht
    real(kind=8) :: bw
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: scmaxi
    real(kind=8) :: scmaxs
    real(kind=8) :: ssmax
    integer(kind=8) :: uc
    real(kind=8) :: dnsinf
    real(kind=8) :: dnssup
    real(kind=8) :: effm
    real(kind=8) :: effn
    integer(kind=8) :: verif

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: Calc
    integer(kind=8) :: s, ntot, ndemi
    logical :: COND_OK
    real(kind=8) :: nrd0, nrd1, mrd0, mrd1
    character(24) :: pnrd, pmrd
    real(kind=8), pointer :: mrd(:) => null()
    real(kind=8), pointer :: nrd(:) => null()

    !Dimensionnement des vecteurs

    pnrd = 'POINT_NRD'
    pmrd = 'POINT_MRD'

    ntot = -1
    call dintels(cequi, ht, bw, enrobi, enrobs, &
                 scmaxi, scmaxs, ssmax, uc, &
                 ntot, ndemi=ndemi)
    call wkvect(pnrd, ' V V R ', ntot, vr=nrd)
    call wkvect(pmrd, ' V V R ', ntot, vr=mrd)

    call dintels(cequi, ht, bw, enrobi, enrobs, &
                 scmaxi, scmaxs, ssmax, uc, &
                 ntot, dnsinf, dnssup, nrd, mrd)

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
            mrd0 = 0.5d0*(mrd(s-1)+mrd(s))
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
            mrd1 = 0.5d0*(mrd(s-1)+mrd(s))
        end if
    end if

    if ((effm .le. mrd0) .and. (effm .ge. mrd1)) then
        COND_OK = .true.
    else
        COND_OK = .false.
    end if

998 continue

    if (COND_OK) then
        verif = 0
    else
        verif = 1
    end if

    call jedetr(pnrd)
    call jedetr(pmrd)

end subroutine
