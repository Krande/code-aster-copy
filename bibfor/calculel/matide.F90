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
subroutine matide(matz, nbcmp, licmp, modlag, tdiag, &
                  vdiag)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    character(len=*) :: matz
    integer(kind=8) :: nbcmp
    real(kind=8) :: vdiag
    character(len=8) :: licmp(nbcmp), tdiag
    character(len=16) :: modlag
! ---------------------------------------------------------------------
! BUT: METTRE "A L'IDENTITE"  UNE MATR_ASSE SUR CERTAINS DDLS
! ---------------------------------------------------------------------
!     ARGUMENTS:
! MATZ   IN/JXVAR K19  : MATR_ASSE A MODIFIER
! NBCMP  IN       I    : NOMBRE DE COMPOSANTES DE LICMP
! LICMP  IN       K8   : LISTE DES NOMS DES COMPOSANTES
! MODLAG IN       K16  : MODIFICATION OU PAS DES TERMES DE LAGRANGE
!                        VALEURS : MODI_LAGR_OUI OU MODI_LAGR_NON
! TDIAG  IN       K8   : NORMALISATION DES TERMES MODIFIEES
!                        DE LA DIAGONALE
!                        VALEURS : MAX_ABS, MIN_ABS ou IMPOSE
! VDIAG  IN       R8   : SI TDIAG VAUT IMPOSE ALORS VDIAG PERMET
!                        DE DONNER LA vALEUR A IMPOSER
! ---------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
    integer(kind=8) :: ilig, jcol, kterm, n, nz, nsmdi, jsmhc, nsmhc
    integer(kind=8) :: jdelg, n1, nvale, jvale, nlong, jval2, nucmp, k, jcmp
    integer(kind=8) :: kcmp
    character(len=8) :: nomgd, nocmp
    character(len=14) :: nonu
    character(len=1) :: ktyp
    character(len=19) :: mat19
    aster_logical :: ltypr, lsym, eliml, elimc
    real(kind=8) :: kmax
    complex(kind=8) :: ckmax
    integer(kind=8), pointer :: smdi(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=24), pointer :: refn(:) => null()
    integer(kind=8), pointer :: deeq(:) => null()
    integer(kind=8), pointer :: lddlelim(:) => null()
    integer(kind=8), pointer :: llag(:) => null()
!
!     ------------------------------------------------------------------
    call jemarq()
!
    mat19 = matz
!
    call jeveuo(mat19//'.REFA', 'L', vk24=refa)
    nonu = refa(2)
!
    call jeveuo(nonu//'.SMOS.SMDI', 'L', vi=smdi)
    call jelira(nonu//'.SMOS.SMDI', 'LONMAX', nsmdi)
    call jeveuo(nonu//'.SMOS.SMHC', 'L', jsmhc)
    call jelira(nonu//'.SMOS.SMHC', 'LONMAX', nsmhc)
    call jeveuo(nonu//'.NUME.DELG', 'L', jdelg)
    call jelira(nonu//'.NUME.DELG', 'LONMAX', n1)
    ASSERT(n1 .eq. nsmdi)
!     --- CALCUL DE N
    n = nsmdi
!     --- CALCUL DE NZ
    nz = smdi(n)
!
    ASSERT(nz .le. nsmhc)
    call jelira(mat19//'.VALM', 'NMAXOC', nvale)
    if (nvale .eq. 1) then
        lsym = .true.
    else if (nvale .eq. 2) then
        lsym = .false.
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(jexnum(mat19//'.VALM', 1), 'E', jvale)
    call jelira(jexnum(mat19//'.VALM', 1), 'LONMAX', nlong)
    ASSERT(nlong .eq. nz)
    if (.not. lsym) then
        call jeveuo(jexnum(mat19//'.VALM', 2), 'E', jval2)
        call jelira(jexnum(mat19//'.VALM', 2), 'LONMAX', nlong)
        ASSERT(nlong .eq. nz)
    end if
!
    call jelira(jexnum(mat19//'.VALM', 1), 'TYPE', cval=ktyp)
    ltypr = (ktyp .eq. 'R')
!
!     -- CALCUL DE LA LISTE DES DDLS A ELIMINER :
!     -------------------------------------------
    AS_ALLOCATE(vi=lddlelim, size=n)
    AS_ALLOCATE(vi=llag, size=n)
!
    call jeveuo(nonu//'.NUME.DEEQ', 'L', vi=deeq)
    call jeveuo(nonu//'.NUME.REFN', 'L', vk24=refn)
    call jelira(nonu//'.NUME.DEEQ', 'LONMAX', n1)
    nomgd = refn(2)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jcmp)
    ASSERT(n1 .eq. 2*n)
    do k = 1, n
        nucmp = deeq(2*(k-1)+2)
        if (nucmp .gt. 0) then
            nocmp = zk8(jcmp-1+nucmp)
            do kcmp = 1, nbcmp
                if (nocmp .eq. licmp(kcmp)) then
                    lddlelim(k) = 1
                end if
            end do
        else if (modlag(1:13) .eq. 'MODI_LAGR_OUI') then
            llag(k) = 1
        end if
    end do
!
!
!
!     ------------------------------------------------
!     PARCOURS DES TERMES DE LA MATRICE
!     ------------------------------------------------
    jcol = 1
    kmax = 0.d0
    ckmax = dcmplx(0.d0, 0.d0)
    if (tdiag(1:7) .eq. 'MAX_ABS') then
        if (ltypr) then
            do kterm = 1, nz
                kmax = max(abs(zr(jvale-1+kterm)), abs(kmax))
            end do
            kmax = kmax*vdiag
        else
            do kterm = 1, nz
                ckmax = max(abs(zc(jvale-1+kterm)), abs(ckmax))
            end do
            ckmax = ckmax*vdiag
        end if
    else if (tdiag(1:7) .eq. 'MIN_ABS') then
        if (ltypr) then
            kmax = abs(zr(jvale))
            do kterm = 1, nz
                kmax = max(abs(zr(jvale-1+kterm)), abs(kmax))
            end do
            kmax = kmax*vdiag
        else
            ckmax = abs(zc(jvale))
            do kterm = 1, nz
                ckmax = max(abs(zc(jvale-1+kterm)), abs(ckmax))
            end do
            ckmax = ckmax*vdiag
        end if
    else if (tdiag(1:6) .eq. 'IMPOSE') then
        kmax = vdiag
    else
        ASSERT(.false.)
    end if
!
    do kterm = 1, nz
        if (smdi(jcol) .lt. kterm) jcol = jcol+1
        ilig = zi4(jsmhc-1+kterm)
        elimc = .false.
        eliml = .false.
     if ((lddlelim(jcol) .eq. 1) .and. (llag(jcol) .eq. 0) .and. (llag(ilig) .eq. 0)) elimc = .true.
     if ((lddlelim(ilig) .eq. 1) .and. (llag(ilig) .eq. 0) .and. (llag(jcol) .eq. 0)) eliml = .true.
!
        if (elimc .or. eliml) then
!
!         -- PARTIE TRIANGULAIRE SUPERIEURE :
            if (jcol .eq. ilig) then
                if (ltypr) then
                    zr(jvale-1+kterm) = kmax
                else
                    zc(jvale-1+kterm) = ckmax
                end if
            else
                if (ltypr) then
                    zr(jvale-1+kterm) = 0.d0
                else
                    zc(jvale-1+kterm) = dcmplx(0.d0, 0.d0)
                end if
            end if
!
!         -- PARTIE TRIANGULAIRE INFERIEURE (SI NON-SYMETRIQUE):
            if (.not. lsym) then
                if (jcol .eq. ilig) then
                    if (ltypr) then
                        zr(jval2-1+kterm) = kmax
                    else
                        zc(jvale-1+kterm) = ckmax
                    end if
                else
                    if (ltypr) then
                        zr(jval2-1+kterm) = 0.d0
                    else
                        zc(jval2-1+kterm) = dcmplx(0.d0, 0.d0)
                    end if
                end if
            end if
        end if
!
    end do
!
    AS_DEALLOCATE(vi=lddlelim)
    AS_DEALLOCATE(vi=llag)
    call jedema()
end subroutine
