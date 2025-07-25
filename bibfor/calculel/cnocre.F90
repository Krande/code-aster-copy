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

subroutine cnocre(maz, nomgdz, nbnoz, linoe, ncmpz, &
                  licmp, cnocmp, basez, prof, cnoz)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: maz, nomgdz, cnoz, basez, prof
    integer(kind=8) :: ncmpz, nbnoz, linoe(nbnoz), cnocmp(nbnoz*ncmpz)
    character(len=*) :: licmp(ncmpz)
! ------------------------------------------------------------------
! BUT : CREER UN CHAM_NO A VALEURS NULLES SUR UN PROFIL DEJA EXISTANT
! OU NON (SI LE PROFIL N'EXISTE PAS -- BLANCS ' ' -- IL EST CREE)
! ------------------------------------------------------------------
!     ARGUMENTS:
! MAZ     IN/JXIN  K8  : MAILLAGE DE CNOZ
! NOMGDZ  IN       K8  : NOM DE LA GRANDEUR DE CNOZ
! NBNOZ   IN       I   : NOMBRE DE NOEUDS VOULUES DANS CNOZ
! LINOE   IN       L_I : NOMS DES NOEUDS VOULUES DANS CNOZ
! NCMPZ   IN       I   : NOMBRE DE CMPS VOULUES DANS CNOZ
! LICMP   IN       L_K8: NOMS DES CMPS VOULUES DANS CNOZ
! BASEZ   IN       K1  : BASE DE CREATION POUR CNOZ : G/V/L
! PROF    IN/JXVAR K19 : SD NUME_EQUA SUR LAQUELLE LE CHAM_NO EST CREE
! CNOZ    IN/JXOUT K19 : SD CHAM_NO A CREER
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: ibid, nbno, ino
    integer(kind=8) :: i, k, jcnsl, jcnsv, ncmp
    character(len=3) :: tsca
    character(len=8) :: nomgd
    character(len=19) :: cns
    character(len=8), pointer :: cnsk(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    cns = '&&CNOCRE.CNS'
    call cnscre(maz, nomgdz, ncmpz, licmp, 'V', &
                cns)
!
    call jeveuo(cns//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cns//'.CNSV', 'E', jcnsv)
    call jeveuo(cns//'.CNSL', 'E', jcnsl)
!
    nomgd = cnsk(2)
    nbno = cnsd(1)
    ncmp = cnsd(2)
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
    ASSERT((tsca .eq. 'R') .or. (tsca .eq. 'C'))
    if (tsca .eq. 'R') then
!         -----------
        if (nbnoz .eq. 0) then
            do k = 1, ncmp
                do ino = 1, nbno
                    zl(jcnsl-1+(ino-1)*ncmp+k) = .true.
                    zr(jcnsv-1+(ino-1)*ncmp+k) = 0.0d0
                end do
            end do
!
        else
            do i = 1, nbnoz
                ino = linoe(i)
                do k = 1, ncmp
                    if (cnocmp((i-1)*ncmp+k) .eq. 1) then
                        zl(jcnsl-1+(ino-1)*ncmp+k) = .true.
                        zr(jcnsv-1+(ino-1)*ncmp+k) = 0.0d0
                    end if
                end do
            end do
        end if
!
    else if (tsca .eq. 'C') then
!             -----------
        if (nbnoz .eq. 0) then
            do k = 1, ncmp
                do ino = 1, nbno
                    zl(jcnsl-1+(ino-1)*ncmp+k) = .true.
                    zc(jcnsv-1+(ino-1)*ncmp+k) = (0.0d0, 0.0d0)
                end do
            end do
!
        else
            do i = 1, nbnoz
                ino = linoe(i)
                do k = 1, ncmp
                    if (cnocmp((i-1)*ncmp+k) .eq. 1) then
                        zl(jcnsl-1+(ino-1)*ncmp+k) = .true.
                        zc(jcnsv-1+(ino-1)*ncmp+k) = (0.0d0, 0.0d0)
                    end if
                end do
            end do
        end if
    end if
!
!
    call cnscno(cns, prof, 'NON', basez, cnoz, &
                'F', ibid)
    call detrsd('CHAM_NO_S', cns)
!
    call jedema()
end subroutine
