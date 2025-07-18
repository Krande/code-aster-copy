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
subroutine asmcyc(cmass, ndim, soumat, beta, ni, &
                  nj, na, axok, liax, nbliax, &
                  libid)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 14/03/91
!-----------------------------------------------------------------------
!  BUT:  ASSEMBLER LES MATRICES MASSE COMPLEXES RELATIVES
!   A UN PROBLEME CYCLIQUE SOUS MATRICES ET DEPHASAGE BETA
!
!
! METHODE:
!         TOUTES LES SOUS-MATRICES POSSIBLES SONT PASSEE  EN REVUE
!         PAR APPEL A DES ROUTINES D'ASSEMBLAGE SPECIALISEES
!         QUI TESTENT L'EXISTENCE DE LA SOUS-MATRICE ET NE FONT
!         RIEN SI ELLE N'EXISTE PAS
!         CHAQUE METHODE: CRAIG-BAMPTON (CB), MAC-NEAL (MN), ET
!                         CRAIG-BAMPTON HARMONIQUE (CBH) Y
!                         TROUVE SON COMPTE
!
!-----------------------------------------------------------------------
!
! CMASS    /O/: MATRICE DE MASSE COMPLEXE
! NDIM     /O/: DIMENSION DES MATRICES COMPLEXES
! SOUMAT   /I/: NOM K24 DE LA FAMILLE DES SOUS-MATRICES
! BETA     /I/: ANGLE DE PEPHASAGE EN RADIANS
! NI       /I/: NOMBRE DE MODE PROPRES DE LA BASE
! NJ       /I/: NOMBRE DE DDL GENERALISES DE DROITE
! NA       /I/: NOMBRE DE DDL GENERALISES AXE
! AXOK     /I/: FLAG ASSEMBLAGE TERMES RELATIFS DDL AXE
! LIAX     /I/: LISTE DES NUMERO DDL AXE A ASSEMBLER
! NBLIAX   /I/: NOMBRE DES DDL AXE A ASSEMBLER (<=NA)
! LIBID    /I/: LISTE BIDON LIBID(I)=I DE DIM >=MAX(NI,NJ)
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/acyel1.h"
#include "asterfort/acyel2.h"
#include "asterfort/acyel4.h"
#include "asterfort/acyelt.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ia, id, na, nbliax, ndim, ni
    integer(kind=8) :: nj
    real(kind=8) :: beta
    character(len=24) :: soumat
    complex(kind=8) :: cmass(*)
    integer(kind=8) :: libid(*), liax(nbliax)
    aster_logical :: axok, vrai, faux
!-----------------------------------------------------------------------
    data vrai, faux/.true._1, .false./
!-----------------------------------------------------------------------
!
!
!-----------------MISE A ZERO DES MATRICES COMPLEXES--------------------
!
    do i = 1, ndim*(ndim+1)/2
        cmass(i) = dcmplx(0.d0, 0.d0)
    end do
!
!
!--------------DETERMINATION DES DIMENSIONS DES SOUS-MATRICES-----------
!
!
! DEBUT DES DDL GENERALISES DE DROITE
!
    id = ni+1
!
! DEBUT DES DDL GENERALISES DE AXE
!
    ia = ni+nj+1
!
!--------------ASSEMBLAGE DES TERMES STANDARDS--------------------------
!
!  POUR TOUTES LES METHODES
!
    call acyelt(soumat, 'M0II', ni, cmass, ndim, &
                1, 1, 1.d0)
!
!  POUR CRAIG-BAMPTON ET CRAIG-BAMPTON HARMONIQUE
!
    call acyelt(soumat, 'M0JJ', nj, cmass, ndim, &
                id, id, 1.d0)
!
    call acyel2(soumat, 'M0IJ', ni, nj, faux, &
                libid, ni, libid, nj, cmass, &
                ndim, 1, id, 1.d0)
!
    call acyel4(soumat, 'MPLUSJJ', nj, nj, faux, &
                libid, nj, libid, nj, cmass, &
                ndim, id, id, beta)
!
    call acyel4(soumat, 'MPLUSIJ', ni, nj, faux, &
                libid, ni, libid, nj, cmass, &
                ndim, 1, id, beta)
!
!  POUR MAC-NEAL
!    RIEN
!
!
!-------------------------CAS DE PRESENCE DE DDL AXE--------------------
!
    if (axok) then
!
!
!  POUR TOUTES LES METHODES
!       RIEN
!
!
!  POUR CRAIG-BAMPTON ET CRAIG-BAMPTON HARMONIQUE
!
        call acyel2(soumat, 'M0IA', ni, na, vrai, &
                    libid, ni, liax, nbliax, cmass, &
                    ndim, 1, ia, 1.d0)
        call acyel2(soumat, 'M0AJ', na, nj, vrai, &
                    liax, nbliax, libid, nj, cmass, &
                    ndim, ia, id, 1.d0)
        call acyel1(soumat, 'M0AA', na, na, vrai, &
                    liax, nbliax, liax, nbliax, cmass, &
                    ndim, ia, ia, 1.d0)
        call acyel4(soumat, 'MPLUSIA', ni, na, vrai, &
                    libid, ni, liax, nbliax, cmass, &
                    ndim, 1, ia, beta)
        call acyel4(soumat, 'MPLUSAJ', na, nj, vrai, &
                    liax, nbliax, libid, nj, cmass, &
                    ndim, ia, id, beta)
        call acyel4(soumat, 'MPLUSJA', nj, na, vrai, &
                    libid, nj, liax, nbliax, cmass, &
                    ndim, id, ia, beta)
        call acyel4(soumat, 'MPLUSAA', na, na, vrai, &
                    liax, nbliax, liax, nbliax, cmass, &
                    ndim, ia, ia, beta)
!
!  POUR MAC-NEAL
!
!    RIEN
!
    end if
!
!
end subroutine
