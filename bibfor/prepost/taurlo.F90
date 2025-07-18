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
subroutine taurlo(nbvec, jvectn, jvectu, jvectv, nbordr, &
                  kwork, sompgw, jrwork, tspaq, ipg, &
                  dectau, jvecpg, jnorma)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: nbvec, jvectn, jvectu, jvectv, nbordr, kwork
    integer(kind=8) :: sompgw, jrwork, tspaq, ipg, jvecpg, jnorma, dectau
! person_in_charge: van-xuan.tran at edf.fr
! ---------------------------------------------------------------------
! BUT: CONSTRUIRE LES COMPOSANTES u ET v DU VECTEUR DE CISAILLEMENT TAU
!      DANS LE REPERE LOCAL PERPENDICULAIRE AU VECTEUR NORMAL, POUR
!      TOUS LES VECTEURS NORMAUX A TOUS LES NUMEROS D'ORDRE.
! ----------------------------------------------------------------------
! ARGUMENTS :
!     NBVEC   : IN  : NOMBRE DE VECTEURS NORMAUX.
!     JVECTN  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS NORMAUX.
!     JVECTU  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS u DU PLAN DE CISAILLEMENT.
!     JVECTV  : IN  : ADRESSE DU VECTEUR CONTENANT LES COMPOSANTES DES
!                     VECTEURS v DU PLAN DE CISAILLEMENT.
!     NBORDR  : IN  : NOMBRE DE NUMEROS D'ORDRE.
!     KWORK   : IN  : KWORK = 0 ON TRAITE LA 1ERE MAILLE DU PAQUET DE
!                               MAILLES ;
!                     KWORK = 1 ON TRAITE LA IEME (I>1) MAILLE DU PAQUET
!                               MAILLES.
!     SOMPGW  : IN  : SOMME DES POINTS DE GAUSS DES N MAILLES PRECEDANT
!                     LA MAILLE COURANTE.
!     JRWORK  : IN  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                     L'HISTORIQUE DES TENSEURS DES CONTRAINTES
!                     ATTACHES A CHAQUE POINT DE GAUSS DES MAILLES
!                     DU <<PAQUET>> DE MAILLES.
!     TSPAQ   : IN  : TAILLE DU SOUS-PAQUET DU <<PAQUET>> DE MAILLES
!                     COURANT.
!     IPG     : IN  : IEME POINT DE GAUSS.
!     DECTAU  : IN  : DECALER ZRWORK DECTAU = 0, STRESS
!                                    DECTAU = 6, STRAIN
!                                    DECTAU = 12, PLASTIC STRAIN
!     JVECPG  : IN  : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                     LES COMPOSANTES u ET v DU VECTEUR TAU
!                     (CISAILLEMENT), POUR TOUS LES NUMEROS
!                     D'ORDRE DE CHAQUE VECTEUR NORMAL.
! ----------------------------------------------------------------------
    integer(kind=8) :: ivect, iordr, n, adrs, decal
    real(kind=8) :: nx, ny, nz, ux, uy, uz, vx, vy, vz
    real(kind=8) :: sixx, siyy, sizz, sixy, sixz, siyz, fx, fy, fz
    real(kind=8) :: norm, taux, tauy, tauz, cutau, cvtau
!     ------------------------------------------------------------------
!
!234567                                                              012
!
    call jemarq()
!
    n = 0
!
!
    do ivect = 1, nbvec
        nx = zr(jvectn+(ivect-1)*3)
        ny = zr(jvectn+(ivect-1)*3+1)
        nz = zr(jvectn+(ivect-1)*3+2)
!
        ux = zr(jvectu+(ivect-1)*3)
        uy = zr(jvectu+(ivect-1)*3+1)
        uz = zr(jvectu+(ivect-1)*3+2)
!
        vx = zr(jvectv+(ivect-1)*3)
        vy = zr(jvectv+(ivect-1)*3+1)
        vz = zr(jvectv+(ivect-1)*3+2)
!
        do iordr = 1, nbordr
            decal = 18
            adrs = (iordr-1)*tspaq+kwork*sompgw*decal+(ipg-1)*decal
!
            sixx = zr(jrwork+adrs+0+dectau)
            siyy = zr(jrwork+adrs+1+dectau)
            sizz = zr(jrwork+adrs+2+dectau)
            sixy = zr(jrwork+adrs+3+dectau)
            sixz = zr(jrwork+adrs+4+dectau)
            siyz = zr(jrwork+adrs+5+dectau)
!
! CALCUL DE vect_F = [SIG].vect_n
            fx = sixx*nx+sixy*ny+sixz*nz
            fy = sixy*nx+siyy*ny+siyz*nz
            fz = sixz*nx+siyz*ny+sizz*nz
!
! CALCUL DE NORM = vect_F.vect_n
            norm = fx*nx+fy*ny+fz*nz
!
! CALCUL DE vect_TAU = vect_F - NORM vect_n
            taux = fx-norm*nx
            tauy = fy-norm*ny
            tauz = fz-norm*nz
!
! PROJECTION DU vect_TAU SUR LES VECTEURS u ET v DU REPERE LOCAL
            cutau = ux*taux+uy*tauy+uz*tauz
            cvtau = vx*taux+vy*tauy+vz*tauz
            n = n+1
            zr(jvecpg+(n-1)*2) = cutau
            zr(jvecpg+(n-1)*2+1) = cvtau
            zr(jnorma+(n-1)) = norm
        end do
!
    end do
!
    call jedema()
end subroutine
