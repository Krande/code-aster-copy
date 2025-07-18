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

subroutine cescar(cesz, cartz, basz)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=*) :: cartz, cesz, basz
! ------------------------------------------------------------------
! BUT: TRANSFORMER UN CHAM_ELEM_S (DE TYPE ELEM)  EN CARTE
! ATTENTION : CETTE ROUTINE EST COUTEUSE POUR LES GROS MAILLAGES
!             JACQUES DEVRA L'AMELIORER PLUTARD
! ------------------------------------------------------------------
!     ARGUMENTS:
! CESZ   IN/JXOUT K19 : SD CHAM_ELEM_S A TRANSFORMER
! CARTZ  IN/JXIN  K19 : SD CARTE A CREER
! BASZ   IN       K1  : BASE DE CREATION POUR CARTZ : G/V
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: jce1d, jce1l, jce1v, nbmam, ncmp, ncmpmx
    integer(kind=8) :: jvalv, iad1, kcmp, ncmpma, nbpt, nbsp, ima
    integer(kind=8) :: k, jvals, nbpaqu, nbcmps, vali(3)
    aster_logical :: idprec, premie
    character(len=1) :: base
    character(len=8) :: ma, nomgd
    character(len=3) :: tsca
    character(len=19) :: cart, ces1
    character(len=24) :: valk(3)
    character(len=8), pointer :: vncmp(:) => null()
    character(len=8), pointer :: cesk(:) => null()
    character(len=8), pointer :: cesc(:) => null()
    integer(kind=8), pointer :: lima(:) => null()
    character(len=8), pointer :: noms(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!      CALL IMPRSD('CHAMP',CESZ,6,'AJOCOT CESCAR IN')
!
    ces1 = cesz
    cart = cartz
    base = basz
!
!
!
!     1- RECUPERATION D'INFORMATIONS DANS CES1 :
!     ------------------------------------------
    call jeveuo(ces1//'.CESK', 'L', vk8=cesk)
    call jeveuo(ces1//'.CESD', 'L', jce1d)
    call jeveuo(ces1//'.CESC', 'L', vk8=cesc)
    call jeveuo(ces1//'.CESV', 'L', jce1v)
    call jeveuo(ces1//'.CESL', 'L', jce1l)
!
    ma = cesk(1)
    nomgd = cesk(2)
    nbmam = zi(jce1d-1+1)
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=ncmpmx)
!
!
    call alcart(base, cart, ma, nomgd)
    call jeveuo(cart//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(cart//'.VALV', 'E', jvalv)
!
    AS_ALLOCATE(vi=lima, size=nbmam)
    AS_ALLOCATE(vk8=noms, size=ncmpmx)
    call wkvect('&&CESCAR.VALS', 'V V '//tsca, ncmpmx, jvals)
!
!     -- POUR ECONOMISER L'ESPACE ET LE TEMPS, ON REGROUPE
!        LES MAILLES SUCCESSIVES QUI PORTENT LES MEMES VALEURS :
!
!     -- IDPREC : .TRUE. -> LA MAILLE EST IDENTIQUE A LA PRECEDENTE
    idprec = .false.
!     -- NBPAQU : NOMBRE DE MAILLES DU "PAQUET" DE MAILLES IDENTIQUES
    nbpaqu = 0
!     -- NBCMPS : NOMBRE DE CMPS DU PAQUET
    nbcmps = 0
!
    do ima = 1, nbmam
        nbpt = zi(jce1d-1+5+4*(ima-1)+1)
        nbsp = zi(jce1d-1+5+4*(ima-1)+2)
        ncmp = zi(jce1d-1+5+4*(ima-1)+3)
        if ((nbpt .gt. 1) .or. (nbsp .gt. 1)) then
            valk(1) = cesz
            valk(2) = cartz
            vali(1) = nbpt
            vali(2) = nbsp
            vali(3) = ima
            call utmess('F', 'MODELISA9_8', nk=2, valk=valk, ni=3, &
                        vali=vali)
        end if
        if (nbpt*nbsp .eq. 0) goto 80
!
!       -- NCMPMA : NBRE DE CMPS SUR LA MAILLE :
        ncmpma = 0
        do kcmp = 1, ncmp
            call cesexi('C', jce1d, jce1l, ima, 1, &
                        1, kcmp, iad1)
            ASSERT(iad1 .ne. 0)
            if (iad1 .gt. 0) then
                ncmpma = ncmpma+1
                vncmp(ncmpma) = cesc(kcmp)
!
                if (tsca .eq. 'R') then
                    zr(jvalv-1+ncmpma) = zr(jce1v-1+iad1)
                else if (tsca .eq. 'C') then
                    zc(jvalv-1+ncmpma) = zc(jce1v-1+iad1)
                else if (tsca .eq. 'I') then
                    zi(jvalv-1+ncmpma) = zi(jce1v-1+iad1)
                else if (tsca .eq. 'K8') then
                    zk8(jvalv-1+ncmpma) = zk8(jce1v-1+iad1)
                else if (tsca .eq. 'K16') then
                    zk16(jvalv-1+ncmpma) = zk16(jce1v-1+iad1)
                else if (tsca .eq. 'K24') then
                    zk24(jvalv-1+ncmpma) = zk24(jce1v-1+iad1)
                else if (tsca .eq. 'K32') then
                    zk32(jvalv-1+ncmpma) = zk32(jce1v-1+iad1)
                else if (tsca .eq. 'K80') then
                    zk80(jvalv-1+ncmpma) = zk80(jce1v-1+iad1)
                else
                    ASSERT(.false.)
                end if
            end if
        end do
        if (ncmpma .eq. 0) goto 80
!
!
        if (nbcmps .eq. 0) then
!         -- C'EST LE 1ER PAQUET QUI COMMENCE
            ASSERT(.not. idprec)
            nbcmps = ncmpma
            premie = .true.
!
        else
!         -- LA MAILLE EST-ELLE COMME LA MAILLE SAUVEGARDEE ?
            if (ncmpma .ne. nbcmps) goto 30
            do k = 1, nbcmps
                if (noms(k) .ne. vncmp(k)) goto 30
                if (tsca .eq. 'R') then
                    if (zr(jvals-1+k) .ne. zr(jvalv-1+k)) goto 30
                else if (tsca .eq. 'C') then
                    if (zc(jvals-1+k) .ne. zc(jvalv-1+k)) goto 30
                else if (tsca .eq. 'I') then
                    if (zi(jvals-1+k) .ne. zi(jvalv-1+k)) goto 30
                else if (tsca .eq. 'K8') then
                    if (zk8(jvals-1+k) .ne. zk8(jvalv-1+k)) goto 30
                else if (tsca .eq. 'K16') then
                    if (zk16(jvals-1+k) .ne. zk16(jvalv-1+k)) goto 30
                else if (tsca .eq. 'K24') then
                    if (zk24(jvals-1+k) .ne. zk24(jvalv-1+k)) goto 30
                else if (tsca .eq. 'K32') then
                    if (zk32(jvals-1+k) .ne. zk32(jvalv-1+k)) goto 30
                else if (tsca .eq. 'K80') then
                    if (zk80(jvals-1+k) .ne. zk80(jvalv-1+k)) goto 30
                end if
            end do
            idprec = .true.
            goto 40
!
30          continue
            idprec = .false.
40          continue
        end if
!
!
        if (.not. idprec) then
!          -- SI LA MAILLE EST DIFFERENTE :
!            - IL FAUT STOCKER LE PAQUET PRECEDENT
!            - PUIS IL FAUT SAUVEGARDER LA NOUVELLE MAILLE
!          -----------------------------------------------------
            if (.not. premie) then
                do k = 1, nbcmps
                    vncmp(k) = noms(k)
                    if (tsca .eq. 'R') then
                        zr(jvalv-1+k) = zr(jvals-1+k)
                    else if (tsca .eq. 'C') then
                        zc(jvalv-1+k) = zc(jvals-1+k)
                    else if (tsca .eq. 'I') then
                        zi(jvalv-1+k) = zi(jvals-1+k)
                    else if (tsca .eq. 'K8') then
                        zk8(jvalv-1+k) = zk8(jvals-1+k)
                    else if (tsca .eq. 'K16') then
                        zk16(jvalv-1+k) = zk16(jvals-1+k)
                    else if (tsca .eq. 'K24') then
                        zk24(jvalv-1+k) = zk24(jvals-1+k)
                    else if (tsca .eq. 'K32') then
                        zk32(jvalv-1+k) = zk32(jvals-1+k)
                    else if (tsca .eq. 'K80') then
                        zk80(jvalv-1+k) = zk80(jvals-1+k)
                    end if
                end do
                call nocart(cart, 3, nbcmps, mode='NUM', nma=nbpaqu, &
                            limanu=lima)
!
!           -- POUR FAIRE LE NOCART, ON A DU ECRASER JVALV.
!           -- IL FAUT LE RETABLIR :
                ncmpma = 0
                do kcmp = 1, ncmp
                    call cesexi('C', jce1d, jce1l, ima, 1, &
                                1, kcmp, iad1)
                    ASSERT(iad1 .ne. 0)
                    if (iad1 .gt. 0) then
                        ncmpma = ncmpma+1
                        vncmp(ncmpma) = cesc(kcmp)
!
                        if (tsca .eq. 'R') then
                            zr(jvalv-1+ncmpma) = zr(jce1v-1+iad1)
                        else if (tsca .eq. 'C') then
                            zc(jvalv-1+ncmpma) = zc(jce1v-1+iad1)
                        else if (tsca .eq. 'I') then
                            zi(jvalv-1+ncmpma) = zi(jce1v-1+iad1)
                        else if (tsca .eq. 'K8') then
                            zk8(jvalv-1+ncmpma) = zk8(jce1v-1+iad1)
                        else if (tsca .eq. 'K16') then
                            zk16(jvalv-1+ncmpma) = zk16(jce1v-1+iad1)
                        else if (tsca .eq. 'K24') then
                            zk24(jvalv-1+ncmpma) = zk24(jce1v-1+iad1)
                        else if (tsca .eq. 'K32') then
                            zk32(jvalv-1+ncmpma) = zk32(jce1v-1+iad1)
                        else if (tsca .eq. 'K80') then
                            zk80(jvalv-1+ncmpma) = zk80(jce1v-1+iad1)
                        end if
                    end if
                end do
            end if
!
            premie = .false.
            nbcmps = ncmpma
            do k = 1, nbcmps
                noms(k) = vncmp(k)
                if (tsca .eq. 'R') then
                    zr(jvals-1+k) = zr(jvalv-1+k)
                else if (tsca .eq. 'C') then
                    zc(jvals-1+k) = zc(jvalv-1+k)
                else if (tsca .eq. 'I') then
                    zi(jvals-1+k) = zi(jvalv-1+k)
                else if (tsca .eq. 'K8') then
                    zk8(jvals-1+k) = zk8(jvalv-1+k)
                else if (tsca .eq. 'K16') then
                    zk16(jvals-1+k) = zk16(jvalv-1+k)
                else if (tsca .eq. 'K24') then
                    zk24(jvals-1+k) = zk24(jvalv-1+k)
                else if (tsca .eq. 'K32') then
                    zk32(jvals-1+k) = zk32(jvalv-1+k)
                else if (tsca .eq. 'K80') then
                    zk80(jvals-1+k) = zk80(jvalv-1+k)
                end if
            end do
            nbpaqu = 1
            lima(nbpaqu) = ima
!
!
        else
!         -- SI LA MAILLE EST IDENTIQUE :
!         --------------------------------
            nbpaqu = nbpaqu+1
            lima(nbpaqu) = ima
        end if
!
80      continue
    end do
!
!     -- IL NE FAUT PAS OUBLIER LE DERNIER PAQUET :
    do k = 1, nbcmps
        vncmp(k) = noms(k)
        if (tsca .eq. 'R') then
            zr(jvalv-1+k) = zr(jvals-1+k)
        else if (tsca .eq. 'C') then
            zc(jvalv-1+k) = zc(jvals-1+k)
        else if (tsca .eq. 'I') then
            zi(jvalv-1+k) = zi(jvals-1+k)
        else if (tsca .eq. 'K8') then
            zk8(jvalv-1+k) = zk8(jvals-1+k)
        else if (tsca .eq. 'K16') then
            zk16(jvalv-1+k) = zk16(jvals-1+k)
        else if (tsca .eq. 'K24') then
            zk24(jvalv-1+k) = zk24(jvals-1+k)
        else if (tsca .eq. 'K32') then
            zk32(jvalv-1+k) = zk32(jvals-1+k)
        else if (tsca .eq. 'K80') then
            zk80(jvalv-1+k) = zk80(jvals-1+k)
        end if
    end do
    call nocart(cart, 3, nbcmps, mode='NUM', nma=nbpaqu, &
                limanu=lima)
!
!
    AS_DEALLOCATE(vi=lima)
    AS_DEALLOCATE(vk8=noms)
    call jedetr('&&CESCAR.VALS')
!
    call jedema()
end subroutine
