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
subroutine usupu2(nbpt, nbpair, coef, ang, isupp, &
                  nbinst, temps, puusur, vustub, vusob, &
                  pus, pmoye, pourpu, poupre)

    use DynaGene_module

    implicit none
!     CALCULE LA PUISSANCE D'USURE AU SENS D'ARCHARD
!                    PU  =  FN * VT
!
! OUT : PUUSUR : PUISSANCE USURE
!-----------------------------------------------------------------------
#include "jeveux.h"
#include "nldef.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/impus.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/stapu2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=8) :: noeu
    character(len=16) :: nomk16
    character(len=19) :: trange
    integer(kind=8) :: nbpair, nbinst
    real(kind=8) :: coef(*), ang(*), temps(*)
    real(kind=8) :: vustub(nbpair, nbinst), vusob(nbpair, nbinst)
    real(kind=8) :: pus(*), pmoye, pourpu(*), poupre(*)
!-----------------------------------------------------------------------
    integer(kind=8) :: ichoc, idebut, ifin, ifires, isupp, j
    integer(kind=8) :: nbtot, jwk1, i, ic
    integer(kind=8) :: jwk2, jwk3, jwk4, jwk5, jwk6, lg
    integer(kind=8) :: n1, n2, n3, n4, nbchoc, nbloc
    integer(kind=8) :: nbpas, nbpt, nbval, nt, nbvint, dec, nbnoli, shift, long
    real(kind=8) :: puusur, tdebut, tfin, tmax, tmin
!
    integer(kind=8), pointer :: chindx(:) => null()
    real(kind=8), pointer :: fcho(:) => null()
    real(kind=8), pointer :: dloc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    character(len=8), pointer :: ncho(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: vcho(:) => null()
    real(kind=8), pointer :: vint(:) => null()
    integer(kind=8), pointer :: nltype(:) => null()
    integer(kind=8), pointer :: vindx(:) => null()
    character(len=24), pointer :: nlname(:) => null()
    real(kind=8), pointer :: v_disc(:) => null()
    type(DynaGene) :: dyna_gene
!-----------------------------------------------------------------------

    call jemarq()
    ifires = iunifi('RESULTAT')
    nbpt = 0
!
    call getvid(' ', 'RESU_GENE', scal=trange, nbret=nt)
    if (nt .eq. 0) goto 999
!
    call getvr8(' ', 'PUIS_USURE', scal=puusur, nbret=n1)
    if (n1 .ne. 0) then
        call impus(ifires, 0, puusur)
        goto 999
    end if

    call jeveuo(trange//'.DESC', 'L', vi=desc)
    nbnoli = desc(3)
!
    nomk16 = trange(1:16)
    call jeveuo(nomk16//'.NL.TYPE', 'L', vi=nltype)
    call jeveuo(nomk16//'.NL.VIND', 'L', vi=vindx)
    call jeveuo(nomk16//'.NL.INTI', 'L', vk24=nlname)
    nbvint = vindx(nbnoli+1)-1
!
    AS_ALLOCATE(vi=chindx, size=nbnoli)
!
    nbchoc = 0
    do i = 1, nbnoli
        if (nltype(i) .eq. NL_CHOC) then
            nbchoc = nbchoc+1
            chindx(nbchoc) = i
        end if
    end do
!
    if (nbchoc .eq. 0) call utmess('F', 'PREPOST4_84')
!
    nbtot = nbchoc
!
!
!
    call getvis(' ', 'NB_BLOC', scal=nbloc, nbret=n1)
    if (n1 .eq. 0) nbloc = 1
    call getvr8(' ', 'INST_INIT', scal=tdebut, nbret=n2)
    call getvr8(' ', 'INST_FIN', scal=tfin, nbret=n3)
    call getvtx(' ', 'NOEUD', scal=noeu, nbret=n4)
!
!           --- RECHERCHE DU NOEUD DE CHOC ---

    AS_ALLOCATE(vk8=ncho, size=2*nbtot)

    do ic = 1, nbchoc
        i = chindx(ic)
        ncho(ic) = nlname((i-1)*5+2) (1:8)
        ncho(nbtot+ic) = nlname((i-1)*5+3) (1:8)
    end do

    do ichoc = 1, nbchoc
        if (ncho(ichoc) .eq. noeu) goto 12
    end do
!
    lg = max(1, lxlgut(noeu))
    call utmess('F', 'UTILITAI_87', sk=noeu(1:lg))
    goto 998
!
12  continue
!

    call dyna_gene%init(trange(1:8))

    if (dyna_gene%n_bloc .eq. 0) then
!
        call jeveuo(trange//'.DISC', 'L', vr=disc)
        call jelira(trange//'.DISC', 'LONMAX', nbpt)

        tmax = disc(nbpt)
        tmin = disc(1)
!
        if (n2 .eq. 0) then
            tdebut = tmin
        else
            if (tdebut .lt. tmin) tdebut = tmin
        end if
        if (n3 .eq. 0) then
            tfin = tmax
        else
            if (tfin .gt. tmax) tfin = tmax
        end if
        if (tdebut .ge. tfin) then
            call utmess('F', 'PREPOST4_47')
        end if
!
        do j = 1, nbpt
            if (disc(j) .ge. tdebut) then
                idebut = j
                goto 15
            end if
        end do
15      continue
!
        do j = 1, nbpt
            if (disc(j) .ge. tfin) then
                ifin = j
                goto 17
            end if
        end do
17      continue
!
    else
        nbpt = dyna_gene%length
        if (n2 .eq. 0) then
            idebut = 1
        else
            call dyna_gene%get_values_by_disc(dyna_gene%disc, tdebut, shift, long, vr=disc)
            idebut = 0
            do i = 1, long
                if (tdebut .le. disc(i)) then
                    idebut = shift+i
                    exit
                end if
            end do
            if (idebut .eq. 0) then
                idebut = shift+long
            end if
        end if
        if (n3 .eq. 0) then
            ifin = nbpt
        else
            call dyna_gene%get_values_by_disc(dyna_gene%disc, tfin, shift, long, vr=disc)
            ifin = 0
            do i = 1, long
                if (tfin .le. disc(i)) then
                    ifin = shift+i
                    exit
                end if
            end do
            if (ifin .eq. 0) then
                ifin = shift+long
            end if
        end if
    end if
    nbpas = ifin-idebut+1
    if (nbloc .eq. 0) nbloc = 1
    nbval = nbpas/nbloc

    AS_ALLOCATE(vr=fcho, size=3*nbtot*nbpas)
    AS_ALLOCATE(vr=dloc, size=2*3*nbtot*nbpas)
    AS_ALLOCATE(vr=vcho, size=3*nbtot*nbpas)
    AS_ALLOCATE(vr=v_disc, size=nbpas)
!

    do j = idebut, ifin
        call dyna_gene%get_values_by_index(dyna_gene%disc, j, shift, vr=disc)
        v_disc(j-idebut+1) = disc(j-shift)
    end do

    do ic = 1, nbchoc
        i = chindx(ic)
        do j = idebut, ifin
            call dyna_gene%get_values_by_index(dyna_gene%vint, j, shift, vr=vint)
!
            fcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+1)
            fcho((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+2)
            fcho((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+3)
!
            dloc((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+4)
            dloc((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+5)
            dloc((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+6)
            dec = 3*nbtot*nbpt
            dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+7)
            dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+8)
            dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+9)
!
            vcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+10)
            vcho((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+11)
            vcho((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+12)

        end do
    end do
!
    call wkvect('&&USURPU.WK1', 'V V R', nbpt, jwk1)
    call wkvect('&&USURPU.WK2', 'V V R', nbpt, jwk2)
    call wkvect('&&USURPU.WK3', 'V V R', nbpt, jwk3)
!
    call wkvect('&&USURPU.WK4', 'V V R', nbpt, jwk4)
    call wkvect('&&USURPU.WK5', 'V V R', nbpt, jwk5)
    call wkvect('&&USURPU.WK6', 'V V R', nbpt, jwk6)
!
    call stapu2(nbchoc, nbpas, nbpair, v_disc, fcho, &
                vcho, dloc, coef, ang, zr(jwk1), &
                zr(jwk2), zr(jwk3), zr(jwk4), zr(jwk5), zr(jwk6), &
                1, nbloc, nbval, ichoc, isupp, &
                nbinst, temps, puusur, vustub, vusob, &
                pus, pmoye, pourpu, poupre)
!
    call jedetr('&&USURPU.WK1')
    call jedetr('&&USURPU.WK2')
    call jedetr('&&USURPU.WK3')
    call jedetr('&&USURPU.WK4')
    call jedetr('&&USURPU.WK5')
    call jedetr('&&USURPU.WK6')
!
    AS_DEALLOCATE(vi=chindx)
    AS_DEALLOCATE(vr=dloc)
    AS_DEALLOCATE(vr=fcho)
    AS_DEALLOCATE(vr=vcho)
    AS_DEALLOCATE(vr=v_disc)

    call dyna_gene%free()

998 continue
    AS_DEALLOCATE(vk8=ncho)
999 continue
    call jedema()

end subroutine
