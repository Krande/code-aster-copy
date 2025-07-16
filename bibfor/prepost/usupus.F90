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
subroutine usupus(puusur, kforn, kvgli, nbpt)
!
    use DynaGene_module
!
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
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/impus.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/reliem.h"
#include "asterfort/statpu.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
!
    character(len=8) :: noeu
    character(len=19) :: trange, kforn, kvgli
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ichoc, idebut, idwk4, ifin, ifires, nbvint
    integer(kind=8) :: ifl, nbtot, nbflam, ic, dec, nbchoc
    integer(kind=8) :: impr, j, jfn, jnomno, nbno
    integer(kind=8) :: jvg, jwk1, jwk2, jwk3, lg
    integer(kind=8) :: n1, n2, n3, nbnoli, nbloc
    integer(kind=8) :: nbpas, nbpt, nbval, nt, shift, long
    character(len=8) :: maillage, modele, base
    character(len=24) :: nomno
    character(len=16) :: motcle(2), typmcl(2)
    real(kind=8) :: puusur, tdebut, tfin, tmax, tmin
!-----------------------------------------------------------------------
    integer(kind=8), pointer :: desc(:) => null()
    real(kind=8), pointer :: disc(:) => null()
    integer(kind=8), pointer :: nltype(:) => null()
    integer(kind=8), pointer :: vindx(:) => null()
    character(len=24), pointer :: nlname(:) => null()
    real(kind=8), pointer :: vint(:) => null()
!-----------------------------------------------------------------------
    integer(kind=8), pointer :: chindx(:) => null()
    integer(kind=8), pointer :: flindx(:) => null()
    real(kind=8), pointer :: vcho(:) => null()
    real(kind=8), pointer :: fcho(:) => null()
    real(kind=8), pointer :: dloc(:) => null()
    character(len=8), pointer :: inti(:) => null()
    integer(kind=8), pointer :: icho(:) => null()
    character(len=8), pointer :: ncho(:) => null()
    real(kind=8), pointer :: v_disc(:) => null()
    type(DynaGene) :: dyna_gene
    blas_int :: b_incx, b_incy, b_n
!
!
!-----------------------------------------------------------------------
    data motcle/'NOEUD', 'GROUP_NO'/
    data typmcl/'NOEUD', 'GROUP_NO'/
!   ------------------------------------------------------------------
!
    call jemarq()
    nomno = '&&USUPUS.MES_NOEUDS'
    ifires = iunifi('RESULTAT')
    nbpt = 0
    impr = 2
!
    call getvr8(' ', 'PUIS_USURE', scal=puusur, nbret=n1)
    if (n1 .ne. 0) then
        call impus(ifires, 0, puusur)
        goto 999
    end if
!
!
    call getvid(' ', 'RESU_GENE', scal=trange, nbret=nt)
    if (nt .ne. 0) then
!
        call dyna_gene%init(trange(1:8))
        call jeveuo(trange//'.DESC', 'L', vi=desc)
        if (desc(1) .eq. 2 .or. desc(1) .eq. 3) then
            nbnoli = desc(3)
            call getvis(' ', 'NB_BLOC', scal=nbloc, nbret=n1)
            if (n1 .eq. 0) nbloc = 1
            call getvr8(' ', 'INST_INIT', scal=tdebut, nbret=n2)
            call getvr8(' ', 'INST_FIN', scal=tfin, nbret=n3)
!
            call dismoi('BASE_MODALE', trange, 'RESU_DYNA', repk=base, arret='F')
            call dismoi('NOM_MODELE', base, 'RESULTAT', repk=modele)
            call dismoi('NOM_MAILLA', base, 'RESULTAT', repk=maillage)
!
            call reliem(modele, maillage, 'NO_NOEUD', ' ', 1, &
                        2, motcle, typmcl, nomno, nbno)
            call jeveuo(nomno, 'L', jnomno)
            noeu = zk8(jnomno)
!
            call jeveuo(trange(1:16)//'.NL.INTI', 'L', vk24=nlname)
!           --- RECHERCHE DU NOEUD DE CHOC ---
            do ichoc = 1, nbnoli
                if (nlname((ichoc-1)*5+2) (1:8) .eq. noeu) goto 12
            end do
            lg = max(1, lxlgut(noeu))
            call utmess('F', 'UTILITAI_87', sk=noeu(1:lg))
12          continue
!
!
            if (dyna_gene%n_bloc .eq. 0) then
                call jeveuo(trange//'.DISC', 'L', vr=disc)
                call jelira(trange//'.DISC', 'LONMAX', nbpt)
                tmax = disc(nbpt)
                tmin = disc(1)
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
                do j = 1, nbpt
                    if (disc(j) .ge. tdebut) then
                        idebut = j
                        goto 15
                    end if
                end do
15              continue
                do j = 1, nbpt
                    if (disc(j) .ge. tfin) then
                        ifin = j
                        goto 17
                    end if
                end do
17              continue
!
            else
                nbpt = dyna_gene%length
                if (n2 .eq. 0) then
                    idebut = 1
                else
                    call dyna_gene%get_values_by_disc(dyna_gene%disc, tdebut, shift, long, &
                                                      vr=disc)
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
!
            end if
!
            nbpas = ifin-idebut+1
            if (nbloc .eq. 0) nbloc = 1
            nbval = nbpas/nbloc
!
            call jeveuo(trange(1:16)//'.NL.TYPE', 'L', vi=nltype)
            call jeveuo(trange(1:16)//'.NL.VIND', 'L', vi=vindx)
            nbvint = vindx(nbnoli+1)-1
!
            AS_ALLOCATE(vi=chindx, size=nbnoli)
            AS_ALLOCATE(vi=flindx, size=nbnoli)
!
            nbchoc = 0
            do i = 1, nbnoli
                if (nltype(i) .eq. NL_CHOC) then
                    nbchoc = nbchoc+1
                    chindx(nbchoc) = i
                end if
            end do
!
            nbflam = 0
            do i = 1, nbnoli
                if (nltype(i) .eq. NL_BUCKLING) then
                    nbflam = nbflam+1
                    flindx(nbflam) = i
                end if
            end do
!
            nbtot = nbchoc+nbflam
!
            AS_ALLOCATE(vr=fcho, size=3*nbtot*nbpas)
            AS_ALLOCATE(vr=dloc, size=2*3*nbtot*nbpas)
            AS_ALLOCATE(vr=vcho, size=3*nbtot*nbpas)
            AS_ALLOCATE(vi=icho, size=nbtot*nbpas)
            AS_ALLOCATE(vk8=ncho, size=2*nbtot)
            AS_ALLOCATE(vk8=inti, size=nbtot)
            AS_ALLOCATE(vr=v_disc, size=nbpas)
!
            do j = idebut, ifin
                call dyna_gene%get_values_by_index(dyna_gene%disc, j, shift, vr=disc)
                v_disc(j-idebut+1) = disc(j-shift)
            end do
!
            do ic = 1, nbchoc
                i = chindx(ic)
                do j = idebut, ifin
!
                    call dyna_gene%get_values_by_index(dyna_gene%vint, j, shift, vr=vint)
!
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+1)
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+2)
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+3)
!
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+4)
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+5)
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+6)
                    dec = 3*nbtot*nbpas
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+1) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+7)
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+2) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+8)
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+3) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+9)
!
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+10)
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+11)
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+12)
!
                    icho((j-idebut)*nbtot+(ic-1)+1) = nint(vint((j-1-shift)*nbvint+vindx(i)-1+13) &
                                                           )
                end do
                inti(ic) = nlname((i-1)*5+1) (1:8)
                ncho(ic) = nlname((i-1)*5+2) (1:8)
                ncho(nbtot+ic) = nlname((i-1)*5+3) (1:8)
            end do
!
            do ifl = 1, nbflam
                i = flindx(ifl)
                ic = nbchoc+ifl
                do j = idebut, ifin
!
                    call dyna_gene%get_values_by_index(dyna_gene%vint, j, shift, vr=vint)
!
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+1)
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+2) = 0.d0
                    fcho((j-idebut)*3*nbtot+(ic-1)*3+3) = 0.d0
!
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+2)
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+2) = vint((j-1-shift)*nbvint+vindx(i)-1+3)
                    dloc((j-idebut)*3*nbtot+(ic-1)*3+3) = vint((j-1-shift)*nbvint+vindx(i)-1+4)
                    dec = 3*nbtot*nbpas
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+1) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+5)
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+2) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+6)
                    dloc(dec+(j-idebut)*3*nbtot+(ic-1)*3+3) = vint( &
                                                              (j-1-shift)*nbvint+vindx(i)-1+7)
!
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+1) = vint((j-1-shift)*nbvint+vindx(i)-1+8)
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+2) = 0.d0
                    vcho((j-idebut)*3*nbtot+(ic-1)*3+3) = 0.d0
!
                    icho((j-idebut)*nbtot+(ic-1)+1) = 0
                end do
                inti(ic) = nlname((i-1)*5+1) (1:8)
                ncho(ic) = nlname((i-1)*5+2) (1:8)
                ncho(nbtot+ic) = nlname((i-1)*5+3) (1:8)
            end do
!
            AS_DEALLOCATE(vi=chindx)
            AS_DEALLOCATE(vi=flindx)
            call jelibe(trange(1:16)//'.NL.TYPE')
            call jelibe(trange(1:16)//'.NL.VIND')
            call jelibe(trange(1:16)//'.NL.INTI')
!
            call wkvect('&&USURPU.WK1', 'V V R', nbpt, jwk1)
            call wkvect('&&USURPU.WK2', 'V V R', nbpt, jwk2)
            call wkvect('&&USURPU.WK3', 'V V R', nbpt, jwk3)
            call wkvect('&&USURPU.IWK4', 'V V I', nbpt, idwk4)
!
            call statpu(nbnoli, nbpas, v_disc, fcho, vcho, &
                        icho, zr(jwk1), zr(jwk2), zr(jwk3), zi(idwk4), &
                        1, nbloc, nbval, ifires, ichoc, &
                        impr, puusur)
!
            call wkvect(kforn, 'V V R', nbpt, jfn)
            call wkvect(kvgli, 'V V R', nbpt, jvg)
            b_n = to_blas_int(nbpt)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jwk1), b_incx, zr(jfn), b_incy)
            do i = 0, nbpt-1
                zr(jvg+i) = sqrt(zr(jwk2+i)**2+zr(jwk3+i)**2)
            end do
!
            call jedetr('&&USURPU.WK1')
            call jedetr('&&USURPU.WK2')
            call jedetr('&&USURPU.WK3')
            call jedetr('&&USURPU.IWK4')
!
!
            AS_DEALLOCATE(vr=dloc)
            AS_DEALLOCATE(vr=fcho)
            AS_DEALLOCATE(vr=vcho)
            AS_DEALLOCATE(vi=icho)
            AS_DEALLOCATE(vk8=ncho)
            AS_DEALLOCATE(vk8=inti)
            AS_DEALLOCATE(vr=v_disc)
!
        else
            call utmess('F', 'PREPOST4_84')
        end if
!
        call dyna_gene%free()
!
    end if
!
!
999 continue
    call jedetr(nomno)
    call jedema()
!
end subroutine
