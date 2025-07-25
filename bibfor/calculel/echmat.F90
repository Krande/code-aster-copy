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

subroutine echmat(matz, ldist, lmhpc, rmin, rmax)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
!
    character(len=*) :: matz
    real(kind=8) :: rmin, rmax
    aster_logical :: ldist, lmhpc
! ---------------------------------------------------------------------
! BUT: DONNER LES VALEURS EXTREMES DES VALEURS ABSOLUES
!      DES TERMES NON NULS DE LA DIAGONALE D'UNE MATR_ASSE
! ---------------------------------------------------------------------
!
!     ARGUMENTS:
!
! IN   MATZ  (K19)     : MATR_ASSE A ANALYSER
! IN   LDIST (LOGICAL) : INDIQUE SI LE CALCUL EST DISTRIBUE AU SENS
!                        DONNEE INCOMPLETE PAR PROC
!
! OUT  RMIN  (R8)      : PLUS PETIT TERME NON NUL (EN VALEUR ABSOLUE)
!                        SUR LA DIAGONALE DE MATZ
! OUT  RMAX  (R8)      : PLUS GRAND TERME (EN VALEUR ABSOLUE)
!                        SUR LA DIAGONALE DE MATZ
! ATTENTION : SI LA MATRICE EST IDENTIQUEMENT NULLE, LA ROUTINE
!             RETOURNE :
!               RMAX=0.D0
!               RMIN=R8MAEM ~1.8E308   (RMIN > RMAX !)
! ---------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: nsmhc, jdelgg, jdelgl, jsmhc, ng, nz, n, imatd
    integer(kind=8) :: jcol, nlong, jvalm1, jcolg
    character(len=1) :: ktyp, base1
    character(len=14) :: nonu
    character(len=19) :: mat19
    character(len=24), pointer :: refa(:) => null()
    integer(kind=8), pointer :: smdi(:) => null()
    real(kind=8), pointer :: rdiag(:) => null()
    complex(kind=8), pointer :: zdiag(:) => null()
    integer(kind=8), pointer :: nlgp(:) => null()
!=================================================================
    call jemarq()
!
    mat19 = matz
    call jeveuo(mat19//'.REFA', 'L', vk24=refa)
    nonu = refa(2) (1:14)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', n)
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    nz = smdi(n)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    ASSERT(nz .le. nsmhc)
!
!
    call jeveuo(nonu//'.NUME.DELG', 'L', jdelgg)
    call jelira(nonu//'.NUME.DELG', 'LONMAX', ng)
    call jeexin(nonu//'.NUML.DELG', imatd)
    if (imatd .ne. 0) then
        call jeveuo(nonu//'.NUML.DELG', 'L', jdelgl)
    else
        jdelgl = jdelgg
        ASSERT(ng .eq. n)
    end if
    !
    call jeexin(nonu//'.NUML.NLGP', imatd)
    if (imatd .ne. 0) then
        call jeveuo(nonu//'.NUML.NLGP', 'L', vi=nlgp)
    end if
!
    call jelira(mat19//'.VALM', 'TYPE', cval=ktyp)
    call jelira(mat19//'.VALM', 'CLAS', cval=base1)
    call jeveuo(jexnum(mat19//'.VALM', 1), 'L', jvalm1)
    call jelira(jexnum(mat19//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
!
!   Le vecteur rdiag/zdiag contient la diagonale de la matrice globale
!
    if (ktyp .eq. 'R') then
        AS_ALLOCATE(vr=rdiag, size=ng)
        rdiag(:) = 0.d0
    else
        AS_ALLOCATE(vc=zdiag, size=ng)
        zdiag(:) = cmplx(0.d0, 0.d0)
    end if
!
!
!     --CALCUL DE RMIN ET RMAX :
!     -----------------------------
    rmin = r8maem()
    rmax = -1.d0
!     CALCUL DE RMIN : PLUS PETIT TERME NON NUL DE LA DIAGONALE
!     CALCUL DE RMAX : PLUS GRAND TERME DE LA DIAGONALE
    do jcol = 1, n
! si le dl est un multiplicateur de Lagrange on passe
        if (zi(jdelgl-1+jcol) .lt. 0) then
            cycle
        end if
!   Indice dans la matrice globale Aster de la colonne locale jcol
        if (imatd .ne. 0) then
            jcolg = nlgp(jcol)
        else
            jcolg = jcol
        end if
!   Lecture des valeurs de rdiag/zdiag
        if (ktyp .eq. 'R') then
            rdiag(jcolg) = rdiag(jcolg)+zr(jvalm1-1+smdi(jcol))
        else
            zdiag(jcolg) = zdiag(jcolg)+zc(jvalm1-1+smdi(jcol))
        end if
!
    end do
!
!     -- SI EXECUTION PARALLELE, IL FAUT COMMUNIQUER
!
    if (ldist) then
        if (ktyp .eq. 'R') then
            call asmpi_comm_vect('MPI_SUM', 'R', nbval=ng, vr=rdiag)
        else
            call asmpi_comm_vect('MPI_SUM', 'C', nbval=ng, vc=zdiag)
        end if
    end if
!
!   Tous les procs possèdent les termes qui correspondent à des ddls physiques.
!   On calcule les valeurs min et max des modules des termes de la diagonale.
!
    if (ktyp .eq. 'R') then
        rmax = maxval(abs(rdiag))
        rmin = minval(abs(rdiag), mask=abs(rdiag) > r8prem())
    else
        rmax = maxval(abs(zdiag))
        rmin = minval(abs(zdiag), mask=abs(zdiag) > r8prem())
    end if
!
    if (lmhpc) then
        if (ktyp .eq. 'R') then
            call asmpi_comm_vect('MPI_MIN', 'R', scr=rmin)
            call asmpi_comm_vect('MPI_MAX', 'R', scr=rmax)
        else
            ASSERT(.false.)
        end if
    end if
!
!   Libération de la mémoire
!
    AS_DEALLOCATE(vr=rdiag)
    AS_DEALLOCATE(vc=zdiag)
    !
    call jedema()
end subroutine
