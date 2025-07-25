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

subroutine pemain(resu, modele, mate, mateco, cara, nh, &
                  nbocc, deform)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlim3.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecham.h"
#include "asterfort/pemica.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/umalma.h"
#include "asterfort/vtgpld.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nh, nbocc
    character(len=*) :: resu, modele, mate, mateco, cara, deform
!     OPERATEUR   POST_ELEM
!     TRAITEMENT DU MOT CLE-FACTEUR "MASS_INER"
!     ------------------------------------------------------------------
!
    integer(kind=8) :: mxvale, nbparr, ibid, iret, iocc, nt, ng, nr, nm, nbgrma, jgr, ig, nbma, jad
    integer(kind=8) :: nbmail, jma, im, nume, nb, ifm, niv, mxval1, nbpar1, mxval2, nbpar2, iorig
    integer(kind=8) :: icage, nbtot, nbMaiT, nre
    parameter(mxval1=16, nbpar1=18)
    parameter(mxval2=25, nbpar2=27)
    real(kind=8) :: zero, orig(3), r8b
    character(len=8) :: k8b, noma, lpain(16), lpaout(5), typarr(nbpar2), valk(2)
    character(len=16) :: noparr(nbpar2)
    character(len=19) :: chelem, chdef
    character(len=24) :: lchin(16), lchout(1), mlggma, valk2(2)
    character(len=24) :: chgeom, chgeo2, chcara(18), chharm, ligrel
    complex(kind=8) :: c16b
    real(kind=8), pointer :: trav1(:) => null()
    integer(kind=8), pointer :: v_allma(:) => null()
!
    data noparr/'LIEU', 'ENTITE', 'MASSE', 'CDG_X', 'CDG_Y', 'CDG_Z', &
        'IX_G', 'IY_G', 'IZ_G', 'IXY_G', 'IXZ_G', 'IYZ_G', 'IX_PRIN_G', &
        'IY_PRIN_G', 'IZ_PRIN_G', 'ALPHA', 'BETA', 'GAMMA', 'X_P', 'Y_P', &
        'Z_P', 'IX_P', 'IY_P', 'IZ_P', 'IXY_P', 'IXZ_P', 'IYZ_P'/
    data typarr/'K24', 'K8', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', &
        'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R'/
!     ------------------------------------------------------------------
!
    call jemarq()
    ibid = 0
    c16b = (0.d0, 0.d0)
!
! --- RECUPERATION DU NIVEAU D'IMPRESSION
    call infniv(ifm, niv)
!
    icage = 0
    zero = 0.0d0
    r8b = 0.0d0
    chdef = deform
    call mecham('MASS_INER', modele, cara, nh, chgeom, &
                chcara, chharm, iret)
    if (iret .ne. 0) goto 60
    noma = chgeom(1:8)
    mlggma = noma//'.GROUPEMA'
!
    call exlim3('MASS_INER', 'V', modele, ligrel)
!
!     --- CALCUL DE L'OPTION ---
    chelem = '&&PEMAIN.MASS_INER'
    lpain(1) = 'PGEOMER'
    if (chdef .ne. ' ') then
        chgeo2 = '&&PEMAIN.CH_GEOMER'
        call vtgpld('CUMU', 1.d0, chgeom, chdef, 'V', &
                    chgeo2)
        lchin(1) = chgeo2
    else
        lchin(1) = chgeom
    end if
    lpain(2) = 'PMATERC'
    lchin(2) = mateco
    lpain(3) = 'PCAORIE'
    lchin(3) = chcara(1)
    lpain(4) = 'PCADISM'
    lchin(4) = chcara(3)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6)
    lpain(6) = 'PCACOQU'
    lchin(6) = chcara(7)
    lpain(7) = 'PCASECT'
    lchin(7) = chcara(8)
    lpain(8) = 'PCAARPO'
    lchin(8) = chcara(9)
    lpain(9) = 'PCAGNBA'
    lchin(9) = chcara(11)
    lpain(10) = 'PCAGEPO'
    lchin(10) = chcara(5)
    lpain(11) = 'PNBSP_I'
    lchin(11) = chcara(16)
    lpain(12) = 'PFIBRES'
    lchin(12) = chcara(17)
    lpain(13) = 'PCOMPOR'
    lchin(13) = mate(1:8)//'.COMPOR'
    lpain(14) = 'PCAPOUF'
    lchin(14) = chcara(13)
    lpain(15) = 'PCINFDI'
    lchin(15) = chcara(15)
    lpain(16) = 'PCACABL'
    lchin(16) = chcara(10)
    nb = 16
!
    lpaout(1) = 'PMASSINE'
    lchout(1) = chelem
!
    call calcul('S', 'MASS_INER', ligrel, nb, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    mxvale = mxval1
    nbparr = nbpar1
    do iocc = 1, nbocc
        call getvr8('MASS_INER', 'ORIG_INER', iocc=iocc, nbval=0, nbret=nr)
        if (nr .ne. 0) then
            mxvale = mxval2
            nbparr = nbpar2
            goto 20
        end if
    end do
20  continue
!
!     --- CREATION DE LA TABLE ---
    call tbcrsd(resu, 'G')
    call tbajpa(resu, nbparr, noparr, typarr)
!
    AS_ALLOCATE(vr=trav1, size=mxvale)
    do iocc = 1, nbocc
        iorig = 0
        orig(1) = zero
        orig(2) = zero
        orig(3) = zero
        call getvtx('MASS_INER', 'TOUT', iocc=iocc, nbval=0, nbret=nt)
        call getvem(noma, 'GROUP_MA', 'MASS_INER', 'GROUP_MA', iocc, &
                    0, k8b, ng)
        call getvem(noma, 'MAILLE', 'MASS_INER', 'MAILLE', iocc, &
                    0, k8b, nm)

        call getvr8('MASS_INER', 'ORIG_INER', iocc=iocc, nbval=0, nbret=nr)
        if (nr .ne. 0) then
            iorig = 1
            nre = -nr
            call getvr8('MASS_INER', 'ORIG_INER', iocc=iocc, nbval=nre, vect=orig, &
                        nbret=nr)
        end if
        if (nt .ne. 0) then
            call pemica(chelem, mxvale, trav1, 0, [ibid], &
                        orig, iorig, icage)
            valk(1) = noma
            valk(2) = 'TOUT'
            call tbajli(resu, nbparr, noparr, [ibid], trav1, &
                        [c16b], valk, 0)
        end if
        if (ng .ne. 0) then
            nbgrma = -ng
            call wkvect('&&PEMAIN_GROUPM', 'V V K24', nbgrma, jgr)
            call getvem(noma, 'GROUP_MA', 'MASS_INER', 'GROUP_MA', iocc, &
                        nbgrma, zk24(jgr), ng)
            valk2(2) = 'GROUP_MA'
            do ig = 1, nbgrma
                call jeexin(jexnom(mlggma, zk24(jgr+ig-1)), iret)
                if (iret .eq. 0) then
                    call utmess('A', 'UTILITAI3_46', sk=zk24(jgr+ig-1))
                    goto 30
                end if
                call jelira(jexnom(mlggma, zk24(jgr+ig-1)), 'LONUTI', nbma)
                if (nbma .eq. 0) then
                    call utmess('A', 'UTILITAI3_47', sk=zk24(jgr+ig-1))
                    goto 30
                end if
                call jeveuo(jexnom(noma//'.GROUPEMA', zk24(jgr+ig-1)), 'L', jad)
                call pemica(chelem, mxvale, trav1, nbma, zi(jad), &
                            orig, iorig, icage)
                valk2(1) = zk24(jgr+ig-1)
                call tbajli(resu, nbparr, noparr, [ibid], trav1, &
                            [c16b], valk2, 0)
30              continue
            end do
!
!
! --- UNION
            if (nbgrma > 1) then
                call umalma(noma, zk24(jgr), nbgrma, v_allma, nbtot)
                ASSERT(nbtot > 0)
                !
                call pemica(chelem, mxvale, trav1, nbtot, v_allma, orig, iorig, icage)
                valk2(1) = "UNION_GROUP_MA"
                call tbajli(resu, nbparr, noparr, [ibid], trav1, [c16b], valk2, 0)
                !
                AS_DEALLOCATE(vi=v_allma)
            end if
            call jedetr('&&PEMAIN_GROUPM')
        end if
        if (nm .ne. 0) then
            nbmail = -nm
            call wkvect('&&PEMAIN_MAILLE', 'V V K8', nbmail, jma)
            call getvem(noma, 'MAILLE', 'MASS_INER', 'MAILLE', iocc, &
                        nbmail, zk8(jma), nm)
            valk(2) = 'MAILLE'
            call jelira(noma//'.TYPMAIL', 'LONMAX', nbMaiT)
            do im = 1, nbmail
                nume = char8_to_int(zk8(jma+im-1))
                if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                    call utmess('A', 'UTILITAI3_49', sk=zk8(jma+im-1))
                    goto 40
                end if
                call pemica(chelem, mxvale, trav1, 1, [nume], &
                            orig, iorig, icage)
                valk(1) = zk8(jma+im-1)
                call tbajli(resu, nbparr, noparr, [ibid], trav1, &
                            [c16b], valk, 0)
40              continue
            end do
            call jedetr('&&PEMAIN_MAILLE')
        end if
    end do
!
! --- MENAGE
    call detrsd('CHAM_ELEM', '&&PEMAIN.MASS_INER')
    call detrsd('CHAMP_GD', '&&PEMAIN.CH_GEOMER')
    AS_DEALLOCATE(vr=trav1)
!
60  continue
!
    call jedema()
end subroutine
