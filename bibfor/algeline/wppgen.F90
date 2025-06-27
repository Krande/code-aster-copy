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
subroutine wppgen(lmasse, lamor, lraide, masseg, amorg, &
                  raideg, vect, neq, nbvect, iddl)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mcmult.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: lmasse, lamor, lraide, neq, nbvect, iddl(*)
    real(kind=8) :: masseg(*), amorg(*), raideg(*)
    complex(kind=8) :: vect(neq, *)
!     CALCUL DES PARAMETRES MODAUX :
!            MASSE, AMORTISSEMENT ET RAIDEUR GENERALISES
!     ------------------------------------------------------------------
! IN  LMASSE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE MASSE
!                   = 0  ON NE CALCULE PAS LA MASSE GENERALISEE
! IN  LAMOR  : IS : DESCRIPTEUR NORMALISE DE LA MATRICE D'AMORTISSEMENT
!                   = 0  ON NE CALCULE PAS L'AMORTISSEMENT GENERALISE
! IN  LRAIDE : IS : DESCRIPTEUR NORMALISE DE LA MATRICE DE RIGIDITE
!                   = 0  ON NE CALCULE PAS LA RIGIDITE GENERALISEE
!     ------------------------------------------------------------------
!     REMARQUE : ON FAIT LES CALCULS VECTEURS APRES VECTEURS
!              : C'EST PLUS LONG MAIS PAS DE PB DE TAILLE MEMOIRE
!     ------------------------------------------------------------------
!
!
    complex(kind=8) :: rval, zero
    character(len=24) :: vecaux, vecau1
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ieq, ivect, laux, laux1
!-----------------------------------------------------------------------
    data vecaux/'&&VPPGEN.VECTEUR.AUX0'/
    data vecau1/'&&VPPGEN.VECTEUR.AUX1'/
!     ------------------------------------------------------------------
    call jemarq()
    zero = 0.d0
    call wkvect(vecaux, 'V V C', neq, laux)
    call wkvect(vecau1, 'V V C', neq, laux1)
    laux = laux-1
    laux1 = laux1-1
!
!     --- CALCUL DE LA MASSE GENERALISEE ---
    if (lmasse .ne. 0) then
        do ivect = 1, nbvect
            call mcmult('ZERO', lmasse, vect(1, ivect), zc(laux+1), 1, &
                        .false._1)
            rval = zero
            do ieq = 1, neq
                rval = rval+dconjg(vect(ieq, ivect))*zc(laux+ieq)
            end do
            masseg(ivect) = dble(rval)
        end do
    end if
!
!     --- CALCUL DE L'AMORTISSEMENT GENERALISE ---
    if (lamor .ne. 0) then
        do ivect = 1, nbvect
            call mcmult('ZERO', lamor, vect(1, ivect), zc(laux+1), 1, &
                        .false._1)
            rval = zero
            do ieq = 1, neq
                rval = rval+dconjg(vect(ieq, ivect))*zc(laux+ieq)
            end do
            amorg(ivect) = dble(rval)
        end do
    else
        do ivect = 1, nbvect
            amorg(ivect) = 0.d0
        end do
    end if
!
!     --- CALCUL DE LA RAIDEUR GENERALISEE ---
    if (lraide .ne. 0) then
        do ivect = 1, nbvect
            do ieq = 1, neq
                zc(laux1+ieq) = vect(ieq, ivect)*iddl(ieq)
            end do
            call mcmult('ZERO', lraide, zc(laux1+1), zc(laux+1), 1, &
                        .false._1)
            rval = zero
            do ieq = 1, neq
                rval = rval+dconjg(vect(ieq, ivect))*zc(laux+ieq)*iddl(ieq)
            end do
            raideg(ivect) = dble(rval)
        end do
    end if
!
    call jedetr(vecaux)
    call jedetr(vecau1)
!
    call jedema()
end subroutine
