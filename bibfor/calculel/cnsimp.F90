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

subroutine cnsimp(cnsz, unite)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: cnsz
    integer(kind=8) :: unite
! ---------------------------------------------------------------------
! BUT: IMPRIMER UN CHAM_NO_S
! ---------------------------------------------------------------------
!     ARGUMENTS:
! CNSZ   IN/JXIN  K19 : SD CHAM_NO_S A IMPRIMER
! UNITE  IN       I   : NUMERO DE L'UNITE LOGIQUE D'IMPRESSION
!
! REMARQUE :
!  - POUR L'INSTANT ON IMPRIME AU FORMAT "RESULTAT" LES CHAMPS DE R8
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: jcnsv, jcnsl
    integer(kind=8) :: nbno, k, ino, ncmp, ncmpu, jlval, ik, licmpu(997)
    character(len=8) :: ma, nomgd, nomno
    character(len=3) :: tsca
    character(len=19) :: cns
    character(len=40) :: fmt1, fmt2
    aster_logical :: exicmp
    character(len=8), pointer :: cnsk(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
!
    cns = cnsz
!
    call jeveuo(cns//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns//'.CNSD', 'L', vi=cnsd)
    call jeveuo(cns//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(cns//'.CNSV', 'L', jcnsv)
    call jeveuo(cns//'.CNSL', 'L', jcnsl)
!
    ma = cnsk(1)
    nomgd = cnsk(2)
    nbno = cnsd(1)
    ncmp = cnsd(2)
!
!
!     1- CALCUL DE NCMPU  : NB CMPS UTILISEES DANS LE CHAMP
!            ET DE LICMPU : NUMEROS DES CMPS UTILISEES
!     ------------------------------------------------------------
    ncmpu = 0
    do k = 1, ncmp
        do ino = 1, nbno
            if (zl(jcnsl-1+(ino-1)*ncmp+k)) goto 20
        end do
        goto 30
20      continue
        ncmpu = ncmpu+1
        ASSERT(ncmpu .le. 997)
        licmpu(ncmpu) = k
30      continue
    end do
!
!     -- LE CHAMP EST VIDE : ON SORT
    if (ncmpu .eq. 0) goto 999
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT((tsca .eq. 'R') .or. (tsca .eq. 'K8') .or. (tsca .eq. 'I') .or. (tsca .eq. 'C'))
!
!
!     1- ALLOCATION D'UN TABLEAU DE K16 QUI CONTIENDRA LES VALEURS
!         D'UNE LIGNE A ECRIRE
!     ------------------------------------------------------------
    if (tsca .ne. 'C') then
        call wkvect('&&CNSIMP.LVALEURS', 'V V K16', ncmpu, jlval)
    else
        call wkvect('&&CNSIMP.LVALEURS', 'V V K16', 2*ncmpu, jlval)
    end if
!
!
!     2- FORMAT DES LIGNES :
!     ----------------------
    if (tsca .ne. 'C') then
        fmt1 = '(A12,XXX(''|'',A12))'
        fmt2 = '(A12,XXX(''|'',A12))'
    else
        fmt1 = '(A12,XXX(''|'',A25))'
        fmt2 = '(A12,XXX(''|'',A12,'' '',A12))'
    end if
    call codent(ncmpu, 'D', fmt1(6:8))
    call codent(ncmpu, 'D', fmt2(6:8))
!
!
!     3- ECRITURE DE L'ENTETE DU CHAMP :
!     ---------------------------------------
    write (unite, *) ' '
    write (unite, *) ' GRANDEUR: ', nomgd
    write (unite, *) ' '
    write (unite, fmt1) 'NOEUD', (cnsc(licmpu(ik)), ik=1, ncmpu)
!
!
!     4- ECRITURE DES VALEURS :
!     ---------------------------------------
    do ino = 1, nbno
        nomno = int_to_char8(ino)
!
!       -- ON N'ECRIT UN NOEUD QUE S'IL EXISTE AU MOINS 1 CMP :
        exicmp = .false.
        do ik = 1, ncmpu
            k = licmpu(ik)
            if (zl(jcnsl-1+(ino-1)*ncmp+k)) then
                exicmp = .true.
                goto 50
            end if
        end do
50      continue
        if (.not. exicmp) goto 70
!
!
!
!       -- ON MET LES VALEURS NON AFFECTEES A " " :
        do ik = 1, ncmpu
            k = licmpu(ik)
            if (zl(jcnsl-1+(ino-1)*ncmp+k)) then
                if (tsca .eq. 'R') then
                    write (zk16(jlval-1+ik), '(E12.5,A4)') zr(jcnsv-1+ &
                                                              (ino-1)*ncmp+k), ' '
                else if (tsca .eq. 'K8') then
                    write (zk16(jlval-1+ik), '(A8,A8)') zk8(jcnsv-1+ &
                                                            (ino-1)*ncmp+k), ' '
                else if (tsca .eq. 'C') then
                    write (zk16(jlval-1+2*(ik-1)+1), '(E12.5,A4)') &
                        dble(zc(jcnsv-1+(ino-1)*ncmp+k)), ' '
                    write (zk16(jlval-1+2*(ik-1)+2), '(E12.5,A4)') &
                        dimag(zc(jcnsv-1+(ino-1)*ncmp+k)), ' '
                else if (tsca .eq. 'I') then
                    write (zk16(jlval-1+ik), '(I12,A4)') zi(jcnsv-1+ &
                                                            (ino-1)*ncmp+k), ' '
                end if
            else
                if (tsca .ne. 'C') then
                    write (zk16(jlval-1+ik), '(A16)') ' '
                else
                    write (zk16(jlval-1+2*(ik-1)+1), '(A16)') ' '
                    write (zk16(jlval-1+2*(ik-1)+2), '(A16)') ' '
                end if
            end if
        end do
        if (tsca .ne. 'C') then
            write (unite, fmt2) nomno, (zk16(jlval-1+ik), ik=1, ncmpu)
        else
            write (unite, fmt2) nomno, (zk16(jlval-1+2*(ik-1)+1), &
                                        zk16(jlval-1+2*(ik-1)+2), ik=1, ncmpu)
        end if
!
70      continue
    end do
!
999 continue
!
    call jedetr('&&CNSIMP.LVALEURS')
    call jedema()
end subroutine
