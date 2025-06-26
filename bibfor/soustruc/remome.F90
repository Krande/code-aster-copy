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
subroutine remome(promes, modmes, nommac)
!
!     RECUPERATION DES MODES MESURES ET CREATION DE  .PJMMM
!
!     IN  : PROMES : NOM DU CONCEPT PROJ_MESU_MODAL ASSOCIE A LA MESURE
!     IN  : MODMES : NOM DU CONCEPT MODES PROPRES IDENTIFIES
!     IN  : NOMMAC : NOM DU CONCEPT MACR_ELEM_STAT CONCERNE
!
    implicit none
!     ------------------------------------------------------------------
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cnocns.h"
#include "asterfort/dismoi.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/scalai.h"
#include "asterfort/wkvect.h"
    character(len=8) :: promes, modmes, nommac
!
    aster_logical :: zcmplx
!
    character(len=1) :: typval
    character(len=8) :: k8bid, scal
    character(len=16) :: nomcha
    character(len=19) :: chamno, chs
    character(len=24) :: vnoeud, vrange, vmes, vorien, vref
!
    integer(kind=8) :: nbmesu, nbmtot, numord, lmesu, imes, lrange, lori, ii
    integer(kind=8) :: iret, icmp, ino, lnoeud, gd, nbcmp, ibid, lref
    integer(kind=8) :: jcnsv, jcnsl, jcnsk, tord(1)
!
    real(kind=8) :: rbid, vori(3), val, vect(3)
!
    complex(kind=8) :: cbid, valc, vectc(3)
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
    vmes = nommac//'.PROJM    .PJMMM'
!
! RECUPERATION ORDRE DE RANGEMENT MESURE SELON VRANGE ET VNOEUD
    vnoeud = promes//'.PROJM    .PJMNO'
    vrange = promes//'.PROJM    .PJMRG'
    vorien = promes//'.PROJM    .PJMOR'
!
    call jelira(vnoeud, 'LONUTI', nbmesu)
!
    vref = nommac//'.PROJM    .PJMRF'
    call jeveuo(vref, 'L', lref)
    nomcha = zk16(lref-1+2)
!
    call jeveuo(vrange, 'L', lrange)
    call jeveuo(vnoeud, 'L', lnoeud)
    call jeveuo(vorien, 'L', lori)
!
! RECUPERATION ADRESSE DES NUMEROS D'ORDRE ET DU NOM SYMBOLIQUE
!
    call jeveuo(modmes//'           .ORDR', 'L', vi=ordr)
!
    chs = '&&MESURE.CHS'
!
! RECUPERATION DU NB DE VECTEURS PROPRES IDENTIFIES : NBMTOT
    call rsorac(modmes, 'LONUTI', 0, rbid, k8bid, &
                cbid, rbid, 'ABSOLU', tord, 1, &
                ibid)
    nbmtot = tord(1)
!
!
! BOUCLE SUR LES NUMEROS ORDRE
!
    do numord = 1, nbmtot
!        -> EXISTENCE DES CHAMPS DANS LA STRUCTURE DE DONNEES MESURE
        call rsexch('F', modmes, nomcha, ordr(numord), chamno, &
                    iret)
        if (numord .le. 1) then
            call dismoi("NUM_GD", chamno, "CHAM_NO", repi=gd)
            scal = scalai(gd)
            typval = scal(1:1)
            if (typval .eq. 'C') then
                zcmplx = .true.
                call wkvect(vmes, 'G V C', nbmesu*nbmtot, lmesu)
            else
                zcmplx = .false.
                call wkvect(vmes, 'G V R', nbmesu*nbmtot, lmesu)
            end if
            call jeecra(vmes, 'LONUTI', nbmesu*nbmtot)
        end if
!
! TRANSFORMATION DE CHAMNO EN CHAM_NO_S : CHS
        call detrsd('CHAM_NO_S', chs)
        call cnocns(chamno, 'V', chs)
        call jeveuo(chs//'.CNSK', 'L', jcnsk)
        call jeveuo(chs//'.CNSD', 'L', vi=cnsd)
        call jeveuo(chs//'.CNSC', 'L', vk8=cnsc)
        call jeveuo(chs//'.CNSV', 'L', jcnsv)
        call jeveuo(chs//'.CNSL', 'L', jcnsl)
!
        nbcmp = cnsd(2)
!
        do imes = 1, nbmesu
            ino = zi(lnoeud-1+imes)
!
! DIRECTION DE MESURE (VECTEUR DIRECTEUR)
            do ii = 1, 3
                vori(ii) = zr(lori-1+(imes-1)*3+ii)
            end do
!
! NORMALISATION DU VECTEUR DIRECTEUR
            val = 0.d0
            do ii = 1, 3
                val = val+vori(ii)*vori(ii)
            end do
            val = sqrt(val)
            do ii = 1, 3
                vori(ii) = vori(ii)/val
            end do
!
            if (zcmplx) then
                do icmp = 1, nbcmp
                    if (cnsc(icmp) .eq. 'DX') vectc(1) = zc(jcnsv-1+(ino-1)*nbcmp+icmp)
                    if (cnsc(icmp) .eq. 'DY') vectc(2) = zc(jcnsv-1+(ino-1)*nbcmp+icmp)
                    if (cnsc(icmp) .eq. 'DZ') vectc(3) = zc(jcnsv-1+(ino-1)*nbcmp+icmp)
                end do
!
                valc = dcmplx(0.d0, 0.d0)
!
                do ii = 1, 3
                    valc = valc+vectc(ii)*vori(ii)
                end do
                zc(lmesu-1+(numord-1)*nbmesu+imes) = valc
            else
                do icmp = 1, nbcmp
                    if (cnsc(icmp) .eq. 'DX') vect(1) = zr(jcnsv-1+(ino-1)*nbcmp+icmp)
                    if (cnsc(icmp) .eq. 'DY') vect(2) = zr(jcnsv-1+(ino-1)*nbcmp+icmp)
                    if (cnsc(icmp) .eq. 'DZ') vect(3) = zr(jcnsv-1+(ino-1)*nbcmp+icmp)
                end do
                val = 0.d0
                do ii = 1, 3
                    val = val+vect(ii)*vori(ii)
                end do
                zr(lmesu-1+(numord-1)*nbmesu+imes) = val
            end if
        end do
!
! FIN BOUCLE SUR NUMERO ORDRE
    end do
!
    call detrsd('CHAM_NO_S', chs)
!
    call jedema()
!
end subroutine
