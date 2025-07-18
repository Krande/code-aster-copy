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

subroutine irvgen(genein, ifi, nbcmpg, cmpg, lhist)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mgutdm.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: cmpg(*)
    character(len=*) :: genein
    aster_logical :: lhist
!     IMPRESSION D'UN "VECT_ASSE_GENE"
!
!     ------------------------------------------------------------------
    character(len=3) :: typval
    character(len=8) :: k8b, mode, noeu, cmp, mogene, blan, dynsta
    character(len=9) :: typmod
    character(len=14) :: nugene
    character(len=16) :: typrem
    character(len=19) :: gene, basmod
    character(len=24) :: typeba
    aster_logical :: lbase
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, ib, ifi, im, imode
    integer(kind=8) :: ir, istru, j, jdesc, jfreq, jnume
    integer(kind=8) :: jordr, jpara, jvale, nbcmpg, nbmode
!
    real(kind=8) :: ximag, xreal
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=24), pointer :: refe(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    blan = ' '
    gene = genein
!
    call jeveuo(gene//'.DESC', 'L', jdesc)
    call jeveuo(gene//'.REFE', 'L', vk24=refe)
    call jeveuo(gene//'.VALE', 'L', jvale)
    call jelira(gene//'.VALE', 'TYPE', cval=typval)
!
    mode = refe(1) (1:8)
    typrem = ' '
    if (mode .ne. blan) call gettco(mode, typrem)
!
!        --- CALCUL PAR SOUS-STRUCTURATION ---
!
    if ((mode .eq. blan) .or. (typrem .eq. 'MODELE_GENE')) then
        nugene = refe(2) (1:14)
        call jeveuo(nugene//'.NUME.REFN', 'L', jnume)
        mogene = zk24(jnume) (1:8)
        call jeveuo(nugene//'.NUME.NEQU', 'L', jnume)
        nbmode = zi(jnume)
        call jeveuo(nugene//'.NUME.DEEQ', 'L', vi=deeq)
        if (lhist) then
            if (typval(1:1) .eq. 'R') then
                write (ifi, 1010)
            else if (typval(1:1) .eq. 'C') then
                write (ifi, 1040)
            end if
        else
            if (typval(1:1) .eq. 'R') then
                write (ifi, 1020)
            else if (typval(1:1) .eq. 'C') then
                write (ifi, 1030)
            end if
        end if
        ir = 0
        im = 0
        do i = 1, nbmode
            imode = deeq(1+2*(i-1)+1-1)
            istru = deeq(1+2*(i-1)+2-1)
            if (istru .lt. 0) goto 10
            im = im+1
            if (nbcmpg .gt. 0) then
                do j = 1, nbcmpg
                    if (im .eq. cmpg(j)) goto 14
                end do
                goto 10
14              continue
            end if
            call mgutdm(mogene, blan, istru, 'NOM_BASE_MODALE', ib, &
                        mode)
            call rsadpa(mode, 'L', 1, 'FREQ', imode, &
                        0, sjv=jfreq, styp=k8b)
            if (lhist) then
                basmod = mode
                call jeveuo(basmod//'.ORDR', 'L', jordr)
                call rsadpa(basmod, 'L', 1, 'TYPE_DEFO', zi(jordr-1+imode), &
                            0, sjv=jpara, styp=k8b)
                typmod = zk16(jpara) (1:9)
                call rsadpa(basmod, 'L', 1, 'NOEUD_CMP', zi(jordr-1+imode), &
                            0, sjv=jpara, styp=k8b)
                noeu = zk16(jpara) (1:8)
                cmp = zk16(jpara) (9:16)
                if (typval(1:1) .eq. 'R') then
                    write (ifi, 1012) im, zr(jvale+i-1), mode, typmod, &
                        zr(jfreq), noeu, cmp
                else if (typval(1:1) .eq. 'C') then
                    xreal = dble(zc(jvale+i-1))
                    ximag = dimag(zc(jvale+i-1))
                    write (ifi, 1042) im, xreal, ximag, mode, typmod, &
                        zr(jfreq), noeu, cmp
                end if
            else
                if (typval(1:1) .eq. 'R') then
                    write (ifi, 1022) im, zr(jvale+i-1)
                else if (typval(1:1) .eq. 'C') then
                    xreal = dble(zc(jvale+i-1))
                    ximag = dimag(zc(jvale+i-1))
                    write (ifi, 1032) im, xreal, ximag
                end if
            end if
            ir = ir+1
            if (ir .eq. nbcmpg) goto 999
10          continue
        end do
    else
!
!      --- CALCUL TRADITIONNEL ---
!
!---------ON RECUPERE LE TYPE DE BASE MODALE S'IL S'AGIT D'UNE BASE
        call dismoi('TYPE_BASE', mode, 'RESU_DYNA', repk=typeba, arret='C', &
                    ier=ir)
!---------ON RECUPERE LE TYPE DE MODES STAT/DYN
        call rsadpa(mode, 'L', 1, 'TYPE_MODE', 1, &
                    0, sjv=iad, styp=k8b)
        dynsta = zk16(iad) (1:8)
        if (typrem .eq. 'MODE_MECA' .or. typrem .eq. 'MODE_GENE') then
            typmod = '  PROPRE'
            noeu = ' '
            cmp = ' '
            lbase = .false.
            call jelira(gene//'.VALE', 'LONMAX', nbmode)
!             NBMODE = ZI(JDESC+1)
        else if (typeba(1:1) .ne. ' ') then
!          ELSEIF ( TYPREM .EQ. 'BASE_MODALE' ) THEN
            lbase = .true.
            basmod = mode
            call jeveuo(basmod//'.ORDR', 'L', jordr)
            call jelira(gene//'.VALE', 'LONMAX', nbmode)
        else if (dynsta .eq. 'MODE_STA') then
!          ELSEIF ( TYPREM(1:9) .EQ. 'MODE_STAT' ) THEN
            typmod = '  PROPRE'
            noeu = ' '
            cmp = ' '
            lbase = .false.
            call jelira(gene//'.VALE', 'LONMAX', nbmode)
        else
            call utmess('A', 'PREPOST3_9', sk=typrem)
            typmod = '  PROPRE'
            noeu = ' '
            cmp = ' '
            lbase = .false.
            call jelira(gene//'.VALE', 'LONMAX', nbmode)
        end if
        if (lhist) then
            if (typval(1:1) .eq. 'R') then
                write (ifi, 1010)
            else if (typval(1:1) .eq. 'C') then
                write (ifi, 1040)
            end if
        else
            if (typval(1:1) .eq. 'R') then
                write (ifi, 1020)
            else if (typval(1:1) .eq. 'C') then
                write (ifi, 1030)
            end if
        end if
        ir = 0
        do i = 1, nbmode
            call rsadpa(mode, 'L', 1, 'TYPE_MODE', i, &
                        0, sjv=iad, styp=k8b)
            dynsta = zk16(iad) (1:8)
            if (nbcmpg .gt. 0) then
                do j = 1, nbcmpg
                    if (i .eq. cmpg(j)) goto 24
                end do
                goto 20
24              continue
            end if
!             IF (TYPREM(1:9) .EQ. 'MODE_STAT') THEN
            if (dynsta .eq. 'MODE_STA') then
                call rsadpa(mode, 'L', 1, 'NOEUD_CMP', i, &
                            0, sjv=jfreq, styp=k8b)
            else
                call rsadpa(mode, 'L', 1, 'FREQ', i, &
                            0, sjv=jfreq, styp=k8b)
            end if
            if (lhist) then
                if (lbase) then
                    call rsadpa(basmod, 'L', 1, 'TYPE_DEFO', zi(jordr-1+i), &
                                0, sjv=jpara, styp=k8b)
                    typmod = zk16(jpara) (1:9)
                    call rsadpa(basmod, 'L', 1, 'NOEUD_CMP', zi(jordr-1+i), &
                                0, sjv=jpara, styp=k8b)
                    noeu = zk16(jpara) (1:8)
                    cmp = zk16(jpara) (9:16)
                end if
                if (dynsta .eq. 'MODE_STA') then
!               IF (TYPREM(1:9) .EQ. 'MODE_STAT') THEN
                    if (typval(1:1) .eq. 'R') then
                        write (ifi, 1013) i, zr(jvale+i-1), mode, &
                            typmod, zk16(jfreq), noeu, cmp
                    else if (typval(1:1) .eq. 'C') then
                        xreal = dble(zc(jvale+i-1))
                        ximag = dimag(zc(jvale+i-1))
                        write (ifi, 1043) i, xreal, ximag, mode, &
                            typmod, zk16(jfreq), noeu, cmp
                    end if
                else
                    if (typval(1:1) .eq. 'R') then
                        write (ifi, 1012) i, zr(jvale+i-1), mode, &
                            typmod, zr(jfreq), noeu, cmp
                    else if (typval(1:1) .eq. 'C') then
                        xreal = dble(zc(jvale+i-1))
                        ximag = dimag(zc(jvale+i-1))
                        write (ifi, 1042) i, xreal, ximag, mode, &
                            typmod, zr(jfreq), noeu, cmp
                    end if
                end if
            else
                if (typval(1:1) .eq. 'R') then
                    write (ifi, 1022) i, zr(jvale+i-1)
                else if (typval(1:1) .eq. 'C') then
                    xreal = dble(zc(jvale+i-1))
                    ximag = dimag(zc(jvale+i-1))
                    write (ifi, 1032) i, xreal, ximag
                end if
            end if
            ir = ir+1
            if (ir .eq. nbcmpg) goto 999
20          continue
        end do
    end if
!
999 continue
!
1010 format(/, ' NUME_CMP   VALEUR        BASE_MODALE  ',&
   &         'TYPE_MODE     FREQUENCE    APPLICATION')
1040 format(/, ' NUME_CMP          VALEUR              BASE_MODALE  ',&
   &         'TYPE_MODE     FREQUENCE    APPLICATION')
1012 format(1p, 3x, i5, 3x, d12.5, 4x, a8, 4x, a9, 3x, d12.5, 3x, a8, a8)
1013 format(1p, 3x, i5, 3x, d12.5, 4x, a8, 4x, a9, 3x, 3x, a8, 3x, a8, a8)
1042 format(1p, 3x, i5, 3x, d12.5, 1x, d12.5, 4x, a8, 4x, a9, 3x, d12.5, 3x, a8, a8)
1043 format(1p, 3x, i5, 3x, d12.5, 1x, d12.5, 4x, a8, 4x, a9, 3x, a8, 3x, a8, a8)
1020 format(/, ' NUME_CMP   VALEUR')
1030 format(/, ' NUME_CMP          VALEUR')
1022 format(1p, 3x, i5, 3x, d12.5)
1032 format(1p, 3x, i5, 3x, d12.5, 1x, d12.5)
!
    call jedema()
end subroutine
