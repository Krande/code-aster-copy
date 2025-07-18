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
subroutine lcummd(vari, nvari, cmat, nmat, sigm, &
                  nstrs, isph, tdt, hini, hfin, &
                  an, bn, cn, cfps, cfpd)
!
!
!
! ROUTINE APPELE DANS NCUMLVFP
! LCUMMD     SOURCE    BENBOU   01/03/26
!-----------------------------------------------------------------------
!_______________________________________________________________________
!
! ROUTINE QUI CALCUL LA MATRICE DE DEFORMATION
!   DE FLUAGE DE DESSICCATION ET DE FLUAGE PROPRE
!
! LA VARIATION DE DEFORMATION DE FLUAGE EST CALCULEE
! AVEC LA RELATION SUIVANTE
!      EPS_N+1 - EPS_N = AN + BN:SIG_N+1 + CN:SIG_N
!
!     EQUATION (3.9-1)
!
! IN  VARI     : VARIABLES INTERNES INITIALES
! IN  NVARI    : DIMENSION DES VECTEURS VARIABLES INTERNES
! IN  CMAT     : VECTEUR DE PARAMETRES (MATERIAU ET AUTRE)
! IN  NMAT     : DIMENSION DE CMAT
! IN  NSTRS    : DIMENSION DES VECTEURS CONTRAINTE ET DEFORMATION
! IN  ISPH     : MODE DE CALCUL DE LA PARTIE SPHERIQUE
! IN  TDT      : PAS DE TEMPS
! IN  HINI     : HUMIDITE INITIALE
! IN  HFIN     : HUMIDITE FINALE
! OUT AN       : CF EQUATION CI-DESSUS
! OUT BN       : CF EQUATION CI-DESSUS
! OUT CN       : CF EQUATION CI-DESSUS
! OUT CFPS     : COEFFICICIENT DE CN SPHERIQUE    (MATRICE TANGENTE)
! OUT CFPD     : COEFFICICIENT DE CN DEVIATORIQUE (MATRICE TANGENTE)
!_______________________________________________________________________
!
    implicit none
#include "asterfort/lcumfb.h"
#include "asterfort/lcumsd.h"
    integer(kind=8) :: nvari, nmat, nstrs, ifpo, isph, ides, i, j
! MODIFI DU 6 JANVIER 2003 - YLP SUPPRESSION DES DECLARATIONS
! IMPLICITES DES TABLEAUX
!      REAL*8 VARI(NVARI),CMAT(NMAT)
    real(kind=8) :: vari(nvari), cmat(nmat)
! MODIFI DU 6 JANVIER 2003 - YLP SUPPRESSION DES DECLARATIONS
! IMPLICITES DES TABLEAUX
!      REAL*8 AN(NSTRS),BN(NSTRS,NSTRS),CN(NSTRS,NSTRS)
!      REAL*8 AFD(NSTRS),BFD(NSTRS,NSTRS),CFD(NSTRS,NSTRS)
!      REAL*8 AFP(NSTRS),BFP(NSTRS,NSTRS),CFP(NSTRS,NSTRS)
    real(kind=8) :: sigm(6)
    real(kind=8) :: an(6), bn(6, 6), cn(6, 6)
    real(kind=8) :: afd(6), bfd(6, 6), cfd(6, 6)
    real(kind=8) :: afp(6), bfp(6, 6), cfp(6, 6)
    real(kind=8) :: cfpd, cfps, hfin, hini, tdt
!
! RECUPERATION DES VALEURS DES PARAMETRES MATERIAU
!
    ifpo = nint(cmat(13))
    ides = nint(cmat(14))
!
! INITIALISATION DES VARIABLES
!
! MODIFI DU 6 JANVIER 2003 6 YLP NSTRS --> 6
!      DO 11 I=1,NSTRS
    do i = 1, 6
        afd(i) = 0.d0
        afp(i) = 0.d0
!        DO 12 J=1,NSTRS
        do j = 1, 6
            bfp(i, j) = 0d0
            bfd(i, j) = 0d0
            cfp(i, j) = 0d0
            cfd(i, j) = 0d0
        end do
    end do
!
!_______________________________________________________________________
!
!  FLUAGE DE DESSICCATION
!
! --> DEBRANCHE le 11 octobre 2002 - ylp.
!
!      IF ((IDES.EQ.1).OR.(IDES.EQ.2)) THEN
!        CALL FLUDES(SIGI,NSTRS,VARI,NVARI,CMAT,NMAT,TDT,HINI,HFIN,
!     $              AFD,BFD,CFD)
!      ENDIF
! --> REBRANCHE le 25 aout 2004 - ylp.
    if ((ides .eq. 1) .or. (ides .eq. 2)) then
        call lcumfb(sigm, nstrs, vari, nvari, cmat, &
                    nmat, tdt, hini, hfin, afd, &
                    bfd, cfd)
    end if
!
!  FLUAGE PROPRE
!
    if (ifpo .ne. 0) then
        call lcumsd(vari, nvari, cmat, nmat, nstrs, &
                    isph, tdt, hini, hfin, afp, &
                    bfp, cfp, cfps, cfpd)
    end if
!
!  CONSTRUCTION DE LA MATRICE DE FLUAGE TOTAL : EQUATION (3.9-2)
!
    do i = 1, nstrs
        an(i) = afd(i)+afp(i)
        do j = 1, nstrs
            bn(i, j) = bfd(i, j)+bfp(i, j)
            cn(i, j) = cfd(i, j)+cfp(i, j)
        end do
    end do
!
end subroutine
