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
subroutine mrmult(cumul, lmat, vect, xsol, nbvect, &
                  prepos, lrom)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmmvr.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/mtmchc.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    character(len=*) :: cumul
    integer(kind=8) :: lmat, nbvect
    real(kind=8) :: vect(*), xsol(*)
    aster_logical :: prepos
    aster_logical, optional :: lrom
!    EFFECTUE LE PRODUIT D'UNE MATRICE PAR N VECTEURS REELS. LE RESULTAT
!    EST STOCKE DANS N VECTEURS REELS
!     ATTENTION:
!       - MATRICE SYMETRIQUE OU NON, REELLE.
!       - LES VECTEURS INPUT ET OUTPUT REELS DOIVENT ETRE DISTINCTS
!       - POUR LES DDLS ELIMINES PAR AFFE_CHAR_CINE, ON NE PEUT PAS
!         CALCULER XSOL. CES DDLS SONT MIS A ZERO.
!     ------------------------------------------------------------------
! IN  CUMUL  : K4 :
!              / 'ZERO' : XSOL =        MAT*VECT
!              / 'CUMU' : XSOL = XSOL + MAT*VECT
!
! IN  LMAT  :I:  DESCRIPTEUR DE LA MATRICE
! IN  VECT  :R: VECTEUR(S) A MULTIPLIER PAR LA MATRICE
! VAR XSOL  :R: VECTEUR(S) SOLUTION(S)
!               SI CUMUL = 'ZERO' ALORS XSOL EST EN MODE OUT
! IN  NBVECT: I : NOMBRE DE VECTEURS A MULTIPLIER (ET DONC DE SOLUTIONS)
!     ------------------------------------------------------------------
    character(len=3) :: kmpic, kmatd
    character(len=19) :: matas
    integer(kind=8) :: neq, neql, jsmhc, jsmdi
    aster_logical :: lmatd, prepo2
    real(kind=8), pointer :: vectmp(:) => null()
    real(kind=8), pointer :: xtemp(:) => null()
    character(len=24), pointer :: refa(:) => null()
    blas_int :: b_incx, b_incy, b_n
!     ---------------------------------------------------------------
!
    prepo2 = prepos
    call jemarq()
    ASSERT(cumul .eq. 'ZERO' .or. cumul .eq. 'CUMU')
    matas = zk24(zi(lmat+1)) (1:19)
    ASSERT(zi(lmat+3) .eq. 1)
    call jeveuo(matas//'.REFA', 'L', vk24=refa)
    if (.not. present(lrom)) then
        if (refa(3) .eq. 'ELIMF') call mtmchc(matas, 'ELIML')
    end if
    neq = zi(lmat+2)
    AS_ALLOCATE(vr=vectmp, size=neq)
!
    call jeveuo(refa(2) (1:14)//'.SMOS.SMHC', 'L', jsmhc)
    call mtdsc2(zk24(zi(lmat+1)), 'SMDI', 'L', jsmdi)
    call dismoi('MPI_COMPLET', matas, 'MATR_ASSE', repk=kmpic)
!
!
!     1.  MATRICE MPI_INCOMPLET :
!     ----------------------------
    if (kmpic .eq. 'NON') then
        if (cumul .eq. 'CUMU') then
            AS_ALLOCATE(vr=xtemp, size=nbvect*neq)
            b_n = to_blas_int(nbvect*neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, xsol, b_incx, xtemp, b_incy)
        end if
!
        call dismoi('MATR_DISTR', matas, 'MATR_ASSE', repk=kmatd)
        if (kmatd .eq. 'OUI') then
            lmatd = .true.
            neql = zi(lmat+5)
        else
            lmatd = .false.
            neql = 0
        end if
        call mrmmvr('ZERO', lmat, zi(jsmdi), zi4(jsmhc), lmatd, &
                    neq, neql, vect, xsol, nbvect, &
                    vectmp, prepo2)
!       ON DOIT COMMUNIQUER POUR OBTENIR LE PRODUIT MAT-VEC 'COMPLET'
        call asmpi_comm_vect('MPI_SUM', 'R', nbval=nbvect*neq, vr=xsol)
!
        if (cumul .eq. 'CUMU') then
            b_n = to_blas_int(nbvect*neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, xtemp, b_incx, xsol, &
                       b_incy)
            AS_DEALLOCATE(vr=xtemp)
        end if
!
!
!     2.  MATRICE MPI_COMPLET :
!     ----------------------------
    else
        lmatd = .false.
        neql = 0
        call mrmmvr(cumul, lmat, zi(jsmdi), zi4(jsmhc), lmatd, &
                    neq, neql, vect, xsol, nbvect, &
                    vectmp, prepo2)
    end if
!
!
    AS_DEALLOCATE(vr=vectmp)
    call jedema()
end subroutine
