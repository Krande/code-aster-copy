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
! aslint: disable=W0413
!
subroutine assma1(matas, ldist, lmhpc)
!
    implicit none
!
!--------------------------------------------------------------
! BUT : METTRE A L'ECHELLE LES LIGNES ET COLONNES D'UNE MATR_ASSE
!       CORRESPONDANT AUX DDLS DE LAGRANGE
!
! IN/JXVAR : MATAS (K19) : SD_MATR_ASSE  :
!    -- CREATION DE L'OBJET .CONL
!    -- MODIFICATION DE L'OBJET .VALM
! IN LDIST (LOGICAL): INDIQUE SI LE CALCUL EST DISTRIBUE AU SENS
!                     DONNEE INCOMPLETE PAR PROC
!---------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/echmat.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: matas
!---------------------------------------------------------------
    aster_logical :: lmnsy, exilag, ldist, lmhpc
    integer(kind=8) :: nsmhc, jdelgg, jdelgl, jsmhc, ng, nz, n, imatd, comm(1)
    integer(kind=8) :: ilig, jcol, kterm, nlong, nvale, jvalm1, jvalm2, jconl
    character(len=1) :: ktyp, base1
    character(len=14) :: nonu
    character(len=19) :: mat19
    real(kind=8) :: rmin, rmax, rcoef
    integer(kind=8), pointer :: smdi(:) => null()
    character(len=24), pointer :: refa(:) => null()
!=================================================================
    call jemarq()
!
!
!
! 1. *  MISE EN MEMOIRE DES OBJETS JEVEUX
!    *  CALCUL DE  :
!        N  : NOMBRE D'EQUATIONS
!        NZ : NOMBRE DE TERMES NON NULS DANS LA MOITIE SUPERIEURE
!        LMNSY : .TRUE.  : LA MATRICE EST NON SYMETRIQUE
!                .FALSE. : LA MATRICE EST SYMETRIQUE
!        KTYP  : 'R'/'C'
!        BASE1 : 'G'/'V'
!    *  QUELQUES VERIFICATIONS DE COHERENCE
! ---------------------------------------------------------------
    mat19 = matas
    call jeveuo(mat19//'.REFA', 'L', vk24=refa)
    nonu = refa(2) (1:14)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', n)
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    nz = smdi(n)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    ASSERT(nz .le. nsmhc)
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
    call jelira(mat19//'.VALM', 'TYPE', cval=ktyp)
    call jelira(mat19//'.VALM', 'CLAS', cval=base1)
    call jeveuo(jexnum(mat19//'.VALM', 1), 'E', jvalm1)
    call jelira(jexnum(mat19//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
!
    lmnsy = .false.
    call jelira(mat19//'.VALM', 'NMAXOC', nvale)
    if (nvale .eq. 2) lmnsy = .true.
!
    if (lmnsy) then
        call jeveuo(jexnum(mat19//'.VALM', 2), 'E', jvalm2)
        call jelira(jexnum(mat19//'.VALM', 2), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
!
!
!     CALCUL DE EXILAG : .TRUE. : IL EXISTE DES DDLS DE LAGRANGE
    exilag = .false.
    do jcol = 1, n
        if (zi(jdelgl-1+jcol) .lt. 0) then
            exilag = .true.
            goto 10
        end if
10      continue
    end do
    if (imatd .ne. 0) then
        exilag = .true.
    end if
!
    if (lmhpc) then
        comm(1) = 0
        if (exilag) comm(1) = 1
        call asmpi_comm_vect('MPI_MAX', 'I', nbval=1, vi=comm)
        if (comm(1) .eq. 1) exilag = .true.
    end if
!
!     -- S'IL N'Y A PAS DE LAGRANGE, IL N'Y A RIEN A FAIRE :
    if (.not. exilag) goto 40
!
! 2.  CALCUL DU COEFFICIENT DE CONDITIONNEMENT DES LAGRANGES (RCOEF)
! -------------------------------------------------------------------
    call echmat(mat19, ldist, lmhpc, rmin, rmax)
!     -- PARFOIS, LA MATRICE EST == 0.
    if (rmax .eq. 0.d0 .or. rmax .lt. rmin) then
        rcoef = 1.d0
    else
        rcoef = 0.5d0*(rmin+rmax)
    end if
!   RCOEF=1.D0
!
! ---------------------------------------------------------------
    call wkvect(mat19//'.CONL', base1//' V R', ng, jconl)
    do jcol = 1, ng
        if (zi(jdelgg-1+jcol) .eq. 0) then
            zr(jconl-1+jcol) = 1.d0
        else
            zr(jconl-1+jcol) = rcoef
        end if
    end do
!
!
! 4.  MISE A L'ECHELLE DE LA MATRICE
! ---------------------------------------------------------------
    jcol = 1
    do kterm = 1, nz
        if (smdi(jcol) .lt. kterm) jcol = jcol+1
        ilig = zi4(jsmhc-1+kterm)
        if (zi(jdelgl-1+jcol)+zi(jdelgl-1+ilig) .lt. 0) then
            if (ktyp .eq. 'R') then
                zr(jvalm1-1+kterm) = rcoef*zr(jvalm1-1+kterm)
            else
                zc(jvalm1-1+kterm) = rcoef*zc(jvalm1-1+kterm)
            end if
            if (lmnsy) then
                if (ktyp .eq. 'R') then
                    zr(jvalm2-1+kterm) = rcoef*zr(jvalm2-1+kterm)
                else
                    zc(jvalm2-1+kterm) = rcoef*zc(jvalm2-1+kterm)
                end if
            end if
        end if
    end do
    ASSERT(jcol .eq. n)
!
!
40  continue
    call jedema()
end subroutine
