! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine pochoc(trange, nbbloc, tdebut, tfin, offset, &
                  trepos, nbclas, nomres, loptio)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nlget.h"
#include "asterfort/statch.h"
#include "asterfort/statim.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"

    character(len=*) :: trange, nomres
!
!     CALCUL ET IMPRESSION DES STATISTIQUES DE CHOC
!     DEUX OPTIONS PREVUES STATISTIQUES POUR VIBRATIONS USURE
!     ET STATISTIQUES POUR LES IMPACTS SOUS SEISME
!
! ----------------------------------------------------------------------
    character(len=19) :: nomk19
    character(len=16) :: nomk16
    aster_logical :: loptio
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!     --- RECUPERATION DES VECTEURS CONTENANT LES RESULTATS ---
!-----------------------------------------------------------------------
    integer :: idwk1, idwk2, idwk3, idwk4, nbbloc
    integer :: nbnoli, nbclas, nbpt, ifm, info, ic, nbchoc, nbflam, i, j, ifl, nbtot, dec, nbvint
    real(kind=8) :: offset, tdebut, tfin, tmax, tmin, trepos
!-----------------------------------------------------------------------
    integer, pointer :: desc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    integer, pointer :: nltype(:) => null()
    integer, pointer :: vindx(:) => null()
    character(len=24), pointer :: nlname(:) => null()
    real(kind=8), pointer :: vint(:) => null()
!-----------------------------------------------------------------------
    integer, pointer :: chindx(:) => null()
    integer, pointer :: flindx(:) => null()
    real(kind=8), pointer :: vin(:) => null()
    real(kind=8), pointer :: vcho(:) => null()
    real(kind=8), pointer :: fcho(:) => null()
    real(kind=8), pointer :: dloc(:) => null()
    character(len=24), pointer :: inti(:) => null()
    integer, pointer :: icho(:) => null()
    character(len=8), pointer :: ncho(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, info)
!
    nomk16 = '                '
    nomk19 = '                   '
    nomk19(1:8) = trange
    nomk16(1:8) = trange
!
    call jeveuo(nomk19//'.DESC', 'L', vi=desc)
    nbnoli = desc(3)
!
    call jeveuo(nomk19//'.DISC', 'L', vr=disc)
    call jelira(nomk19//'.DISC', 'LONMAX', nbpt)
    write (ifm, *) ' NB DE PAS DE TEMPS :', nbpt
!
    tmax = disc(nbpt)
    tmin = disc(1)
    if (tfin .gt. tmax) tfin = tmax
    if (tdebut .lt. tmin) tdebut = tmin
    if (tdebut .ge. tfin) then
        call utmess('F', 'PREPOST4_47')
    end if
!
    call jeveuo(nomk16//'.NL.TYPE', 'L', vi=nltype)
    call jeveuo(nomk16//'.NL.VIND', 'L', vi=vindx)
    call jeveuo(nomk16//'.NL.INTI', 'L', vk24=nlname)
    call jeveuo(nomk16//'.NL.VINT', 'L', vr=vint)
    nbvint = vindx(nbnoli+1)-1
!
    AS_ALLOCATE(vi=chindx, size=nbnoli)
    AS_ALLOCATE(vi=flindx, size=nbnoli)

    nbchoc = 0
    do i = 1, nbnoli
        if (nltype(i) .eq. NL_CHOC) then
            nbchoc = nbchoc+1
            chindx(nbchoc) = i
        end if
    end do

    nbflam = 0
    do i = 1, nbnoli
        if (nltype(i) .eq. NL_BUCKLING) then
            nbflam = nbflam+1
            flindx(nbflam) = i
        end if
    end do

    nbtot = nbchoc+nbflam

    AS_ALLOCATE(vr=fcho, size=3*nbtot*nbpt)
    AS_ALLOCATE(vr=dloc, size=2*3*nbtot*nbpt)
    AS_ALLOCATE(vr=vcho, size=3*nbtot*nbpt)
    AS_ALLOCATE(vr=vin, size=nbtot*nbpt)
    AS_ALLOCATE(vi=icho, size=nbtot*nbpt)
    AS_ALLOCATE(vk8=ncho, size=2*nbtot)
    AS_ALLOCATE(vk24=inti, size=nbtot)

    do ic = 1, nbchoc
        i = chindx(ic)
        do j = 1, nbpt
            fcho((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+1)
            fcho((j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+2)
            fcho((j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+3)

            dloc((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+4)
            dloc((j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+5)
            dloc((j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+6)
            dec = 3*nbtot*nbpt
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+7)
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+8)
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+9)

            vcho((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+10)
            vcho((j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+11)
            vcho((j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+12)

            vin((j-1)*nbtot+(ic-1)+1) = 0.d0

            icho((j-1)*nbtot+(ic-1)+1) = nint(vint((j-1)*nbvint+vindx(i)-1+13))
        end do
        inti(ic) = nlname((i-1)*5+1) (1:24)
        ncho(ic) = nlname((i-1)*5+2) (1:8)
        ncho(nbtot+ic) = nlname((i-1)*5+3) (1:8)
    end do

    do ifl = 1, nbflam
        i = flindx(ifl)
        ic = nbchoc+ifl
        do j = 1, nbpt
            fcho((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+1)
            fcho((j-1)*3*nbtot+(ic-1)*3+2) = 0.d0
            fcho((j-1)*3*nbtot+(ic-1)*3+3) = 0.d0

            dloc((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+2)
            dloc((j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+3)
            dloc((j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+4)
            dec = 3*nbtot*nbpt
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+5)
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+2) = vint((j-1)*nbvint+vindx(i)-1+6)
            dloc(dec+(j-1)*3*nbtot+(ic-1)*3+3) = vint((j-1)*nbvint+vindx(i)-1+7)

            vcho((j-1)*3*nbtot+(ic-1)*3+1) = vint((j-1)*nbvint+vindx(i)-1+8)
            vcho((j-1)*3*nbtot+(ic-1)*3+2) = 0.d0
            vcho((j-1)*3*nbtot+(ic-1)*3+3) = 0.d0

            vin((j-1)*nbtot+(ic-1)+1) = vint((j-1)*nbvint+vindx(i)-1+9)

            icho((j-1)*nbtot+(ic-1)+1) = 0
        end do
        inti(ic) = nlname((i-1)*5+1) (1:24)
        ncho(ic) = nlname((i-1)*5+2) (1:8)
        ncho(nbtot+ic) = nlname((i-1)*5+3) (1:8)
    end do

    AS_DEALLOCATE(vi=chindx)
    AS_DEALLOCATE(vi=flindx)
    call jelibe(nomk16//'.NL.TYPE')
    call jelibe(nomk16//'.NL.VIND')
    call jelibe(nomk16//'.NL.INTI')
    call jelibe(nomk16//'.NL.VINT')

    call wkvect('&&OP0130.WK1', 'V V R', nbpt, idwk1)
    call wkvect('&&OP0130.WK2', 'V V R', nbpt, idwk2)
    call wkvect('&&OP0130.WK3', 'V V R', nbpt, idwk3)
    call wkvect('&&OP0130.IWK4', 'V V I', nbpt, idwk4)
!
    if (loptio) then
!       --- Wear analysis
        call statch(nbtot, nbpt, disc, dloc, fcho, &
                    vcho, icho, zr(idwk1), zr(idwk2), zr(idwk3), &
                    zi(idwk4), tdebut, tfin, nbbloc, offset, &
                    trepos, ncho, inti, nomres)
    else
!       --- Impact analysis
        call statim(nbtot, nbpt, disc, fcho, vcho, &
                    vin, zr(idwk1), zr(idwk2), zr(idwk3), tdebut, &
                    tfin, nbbloc, offset, trepos, nbclas, &
                    ncho, inti, nomres, nbtot)
    end if
!
    call jedetr('&&OP0130.WK1')
    call jedetr('&&OP0130.WK2')
    call jedetr('&&OP0130.WK3')
    call jedetr('&&OP0130.IWK4')

    AS_DEALLOCATE(vr=dloc)
    AS_DEALLOCATE(vr=fcho)
    AS_DEALLOCATE(vr=vcho)
    AS_DEALLOCATE(vr=vin)
    AS_DEALLOCATE(vi=icho)
    AS_DEALLOCATE(vk8=ncho)
    AS_DEALLOCATE(vk24=inti)
!
    call jedema()
end subroutine
