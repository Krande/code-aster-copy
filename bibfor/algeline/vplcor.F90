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
subroutine vplcor(ldynam, neq, nbvect, nborto, prorto, &
                  signes, vect, ivecp, pkx, plx)
    implicit none
#include "jeveux.h"
#include "asterfort/mrmult.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ldynam, neq, nborto, nbvect, ivecp
    real(kind=8) :: prorto
    real(kind=8) :: signes(nbvect), vect(neq, nbvect), pkx(neq, nbvect)
    real(kind=8) :: plx(neq)
!     ORTHOGONALISATION D'UN VECTEUR PAR RAPPORT A UNE FAMILLE DE
!     VECTEURS DEJA ORTHOGONAUX
!     LES VECTEURS OBTENUS SONT K-ORTHOGONAUX .
!     ------------------------------------------------------------------
! IN  LDYNAM : IS : DESCRIPTEUR MATRICE DE "RAIDEUR"
! IN  NEQ    : IS : DIMENSION DES VECTEURS
! IN  NBVECT : IS : NOMBRE TOTAL DE VECTEURS (SERT A DIMENSIONER)
! IN  PRORTO : R8 : PRECISON DEMANDEE POUR L'ORTHOGONALISATION
! IN  NBORTO : IS : NOMBRE MAXIMUM D'ORTHOGONALISATION PERMISE.
! IN  PLX    : R8 : VECTEUR DE TRAVAIL
!     SIGNES : R8 : (+/- 1)  SIGNE DU PSEUDO PRODUIT SCALAIRE ENTRE LE
!                   I EME VECTEUR ET LE I-1 EME VECTEUR
!     VECT   : R8 : VECT(1..NEQ,1..NBVECT) PRODUIT DE LA MATRICE DE
!                   RAIDEUR PAR LES VECTEURS DE LA BASE
! IN  IVECP  : IS : NUMERO DU VECTEUR A ORTHOGONALISER AVEC LES IVECP-1
!                   AUTRES PREMIERS VECTEURS
!     ------------------------------------------------------------------
!
!
!     -----------------------------------------------------------------
    integer(kind=8) :: ieq, ito
    real(kind=8) :: coef, xikxi, xjkxi, xjkxis
!     -----------------------------------------------------------------
!
!         --- K-REORTHOGONALISATION COMPLETE DU VECTEUR IVECP
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ior, iortho, jvec
!-----------------------------------------------------------------------
    ior = 0
!
    do jvec = 1, ivecp-1
!
        xjkxi = 0.d0
        do ieq = 1, neq
            xjkxi = xjkxi+pkx(ieq, jvec)*vect(ieq, ivecp)
        end do
!
        iortho = 0
        if (abs(xjkxi) .gt. prorto) then
            ior = 1
            do ito = 1, nborto
                iortho = ito
!
                do ieq = 1, neq
                    plx(ieq) = vect(ieq, ivecp)-xjkxi*signes(jvec)*vect(ieq, jvec)
                end do
!
                xjkxis = 0.d0
                do ieq = 1, neq
                    xjkxis = xjkxis+pkx(ieq, jvec)*plx(ieq)
                end do
!
                if (abs(xjkxis) .lt. prorto) then
                    do ieq = 1, neq
                        vect(ieq, ivecp) = plx(ieq)
                    end do
                    xjkxi = xjkxis
                    goto 100
                else if (abs(xjkxis) .lt. abs(xjkxi)) then
                    do ieq = 1, neq
                        vect(ieq, ivecp) = plx(ieq)
                    end do
                    xjkxi = xjkxis
                else
                    call utmess('A', 'ALGELINE4_76', si=iortho)
                    goto 100
                end if
!
            end do
            iortho = -nborto
!
100         continue
            call mrmult('ZERO', ldynam, vect(1, ivecp), pkx(1, ivecp), 1, &
                        .false._1)
        end if
!
    end do
!
!         --- SI LE VECTEUR IVECP A ETE MODIFIE (IOR=1) ALORS ---
!         ---               ON LE RENORMALISE                 ---
!
    if (ior .eq. 1) then
        xikxi = 0.d0
        do ieq = 1, neq
            xikxi = xikxi+vect(ieq, ivecp)*pkx(ieq, ivecp)
        end do
        signes(ivecp) = sign(1.d0, xikxi)
        coef = 1.d0/sqrt(abs(xikxi))
        do ieq = 1, neq
            vect(ieq, ivecp) = coef*vect(ieq, ivecp)
            pkx(ieq, ivecp) = coef*pkx(ieq, ivecp)
        end do
    end if
!
end subroutine
