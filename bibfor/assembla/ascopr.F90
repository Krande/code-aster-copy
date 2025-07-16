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
subroutine ascopr(lmasym, lmesym, tt, jtmp2, nrmax, &
                  jresl, rcoef, jvalm)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
    character(len=2) :: tt
    integer(kind=8) :: jtmp2, nrmax, jresl
    integer(kind=8) :: jvalm(2)
    aster_logical :: lmasym, lmesym
    real(kind=8) :: rcoef
!     ROUTINE QUI ACCUMULE LES TERMES ELEMENTAIRES DANS LES BLOCS DE LA
!     MATRICE ASSEMBLEE POUR UNE MATRICE REELLE
!-----------------------------------------------------------------------
! IN  K2  TT    : /'RR' : COPIE R->R
!                 /'CC' : COPIE C->C
!                 /'RC' : COPIE R->C
! IN  I   JTMP2 : ADRESSE JEVEUX DE L'OBJET ".TMP2"
! IN  I   NRMAX  : NOMBRE DE REELS A CHARGER.
! IN  I   JRESL : ADRESSE JEVEUX DE L'OBJET ".RESL(IEL)".
! IN  L   LMASYM   : .TRUE. => LA MATR_ASSE EST SYMETRIQUE (2 BLOCS)
! IN  L   LMESYM   : .TRUE. => LA MATRICE ELEMANTAIRE EST SYMETRIQUE
! IN  R   RCOEF   : COEFFICIENT REEL MULTIPLICATEUR.
! IN  I   JVALM  : LISTE DES ADRESSES DES BLOCS DE .VALM
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibloc, iadloc, j, jvalb, kfois, nbfois, permbl(2)
!-----------------------------------------------------------------------
!
!     -- SI ON ASSEMBLE UNE MATRICE ELEMENTAIRE SYMETRIQUE DANS
!        UNE MATRICE ASSEMBLEE NON SYMETRIQUE, IL FAUT COPIER 2 FOIS
    nbfois = 1
    if (lmasym) then
        ASSERT(lmesym)
    else
        if (lmesym) nbfois = 2
    end if
!
!
    do kfois = 1, nbfois
!
        if (kfois .eq. 1) then
            permbl(1) = 1
            permbl(2) = 2
        else
            permbl(1) = 2
            permbl(2) = 1
        end if
!
!
        if (tt .eq. 'RR') then
!       ---------------------------
            do j = 1, nrmax
                ibloc = zi(jtmp2-1+2*(j-1)+1)
                jvalb = jvalm(permbl(ibloc))
                iadloc = zi(jtmp2-1+2*(j-1)+2)
                zr(jvalb+iadloc-1) = zr(jvalb+iadloc-1)+rcoef*zr(jresl- &
                                                                 1+j)
            end do
!
!
        else if (tt .eq. 'CC') then
!       ---------------------------
            do j = 1, nrmax
                ibloc = zi(jtmp2-1+2*(j-1)+1)
                jvalb = jvalm(permbl(ibloc))
                iadloc = zi(jtmp2-1+2*(j-1)+2)
                zc(jvalb+iadloc-1) = zc(jvalb+iadloc-1)+rcoef*zc(jresl- &
                                                                 1+j)
            end do
!
!
        else if (tt .eq. 'RC') then
!       ---------------------------
            do j = 1, nrmax
                ibloc = zi(jtmp2-1+2*(j-1)+1)
                jvalb = jvalm(permbl(ibloc))
                iadloc = zi(jtmp2-1+2*(j-1)+2)
                zc(jvalb+iadloc-1) = zc(jvalb+iadloc-1)+dcmplx(rcoef* &
                                                               zr(jresl-1+j), 0.d0)
            end do
!
!
        else
            ASSERT(.false.)
        end if
!
    end do
!
end subroutine
