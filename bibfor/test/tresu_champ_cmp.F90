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

subroutine tresu_champ_cmp(chamgd, typtes, typres, nbref, tbtxt, &
                           refi, refr, refc, epsi, lign1, &
                           lign2, crit, ific, nbcmp, nocmp, &
                           llab, ssigne, ignore, compare)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesred.h"
#include "asterfort/cnsred.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/asmpi_comm_vect.h"
!
    character(len=*), intent(in) :: chamgd
    character(len=8), intent(in) :: typtes
    character(len=*), intent(in) :: typres
    integer(kind=8), intent(in) :: nbref
    character(len=16), intent(in) :: tbtxt(2)
    integer(kind=8), intent(in) :: refi(nbref)
    real(kind=8), intent(in) :: refr(nbref)
    complex(kind=8), intent(in) :: refc(nbref)
    real(kind=8), intent(in) :: epsi
    character(len=200), intent(inout) :: lign1
    character(len=200), intent(inout) :: lign2
    character(len=*), intent(in) :: crit
    integer(kind=8), intent(in) :: ific
    integer(kind=8), intent(in) :: nbcmp
    character(len=8), intent(in) :: nocmp(*)
    aster_logical, intent(in) :: llab
    character(len=*), intent(in) :: ssigne
    aster_logical, intent(in), optional :: ignore
    real(kind=8), intent(in), optional :: compare
! IN  : CHAMGD : NOM DU CHAM_GD
! IN  : TYPTES : TYPE DE TEST A EFFECTUER SUR LE CHAMP
! IN  : REFI   : VALEUR REELLE ENTIERE ATTENDUE
! IN  : REFR   : VALEUR REELLE ATTENDUE
! IN  : REFC   : VALEUR COMPLEXE ATTENDUE
! IN  : CRIT   : 'RELATIF' OU 'ABSOLU'(PRECISION RELATIVE OU ABSOLUE).
! IN  : EPSI   : PRECISION ESPEREE
! IN  : IFIC   : NUMERO LOGIQUE DU FICHIER DE SORTIE
! IN  : TBTXT  : (1)=REFERENCE, (2)=LEGENDE
! IN  : NODDL  : NOM DU DDL A TRAITER
! IN  : LLAB   : FLAG D IMPRESSION DES LABELS
! IN/OUT  : LIGN1  : PREMIERE LIGNE D'IMPRESSION DU RESULTAT
! IN/OUT  : LIGN2  : DEUXIEME LIGNE D'IMPRESSION DU RESULTAT
! OUT : IMPRESSION SUR LISTING
! ----------------------------------------------------------------------
    integer(kind=8) :: vali, neq, i, j, k, iret1, valii, icmp
    integer(kind=8) :: ncmp, vnocmp, jcsd, jcsc, jcsv, jcsl, jcmp, ind
    integer(kind=8) :: nl1, nl11, nl2, nl22, rank, nbno, nbno_list, nbma, nbma_list
    real(kind=8) :: valr, valrr
    complex(kind=8) :: valc
    character(len=1) :: typrez
    character(len=24) :: valk(3)
    character(len=4) :: type
    character(len=8) :: tych, noddl, mesh
    character(len=19) :: cham19, cnsinr, cnsin1
    aster_logical :: skip, l_parallel_mesh
    real(kind=8) :: ordgrd
    mpi_int :: irank
    integer(kind=8), pointer :: v_noex(:) => null()
    integer(kind=8), pointer :: v_maex(:) => null()
    integer(kind=8), pointer :: v_list(:) => null()
!     ------------------------------------------------------------------
!
    skip = .false.
    if (present(ignore)) then
        skip = ignore
    end if
!
    ordgrd = 1.d0
    if (present(compare)) then
        ordgrd = compare
    end if
!
    call jemarq()
!
    cham19 = chamgd
    typrez = typres(1:1)
!
    call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=mesh, arret='F')
    l_parallel_mesh = isParallelMesh(mesh)
!
    if (l_parallel_mesh) then
        call asmpi_info(rank=irank)
        rank = to_aster_int(irank)
        call jeveuo(mesh//'.NOEX', 'L', vi=v_noex)
        call jeveuo(mesh//'.MAEX', 'L', vi=v_maex)
    end if
!
    call wkvect('&&TRESU_CH.CMP', 'V V I', nbcmp, jcmp)
!
!     -- LE CHAMP EXISTE-T-IL ?
!     =========================
    call dismoi('TYPE_CHAMP', cham19, 'CHAMP', repk=tych, arret='C', &
                ier=iret1)

    ASSERT(nbcmp .eq. 1)

    if (tych(1:4) .eq. 'NOEU') then
!   --------------------------------
        cnsin1 = '&&TRESU_CH.CNSIN1'
        cnsinr = '&&TRESU_CH.CNSINR'
        call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nbno)
        call wkvect('&&TRESU_CH.LST', 'V V I', nbno, vi=v_list)
        nbno_list = 0
        if (l_parallel_mesh) then
            do i = 1, nbno
                if (v_noex(i) == rank) then
                    nbno_list = nbno_list+1
                    v_list(nbno_list) = i
                end if
            end do
        else
            do i = 1, nbno
                v_list(i) = i
            end do
            nbno_list = nbno
        end if
        call cnocns(cham19, 'V', cnsin1)
        call cnsred(cnsin1, nbno_list, v_list, 0, ['XXX'], &
                    'V', cnsinr)
        call detrsd('CHAM_NO_S', cnsin1)
        call jedetr('&&TRESU_CH.LST')

        call jeveuo(cnsinr//'.CNSV', 'L', jcsv)
        call jeveuo(cnsinr//'.CNSC', 'L', jcsc)
        call jeveuo(cnsinr//'.CNSL', 'L', jcsl)
        call jeveuo(cnsinr//'.CNSD', 'L', jcsd)
        ncmp = zi(jcsd-1+2)
        do i = 1, nbcmp
            noddl = nocmp(i)
            do j = 1, ncmp
                if (zk8(jcsc-1+j) .eq. noddl) then
                    zi(jcmp-1+i) = j
                    goto 10
                end if
            end do
            call utmess('F', 'CHAMPS_3', sk=noddl)
10          continue
        end do
        call jelira(cnsinr//'.CNSV', 'TYPE', cval=type)
        call jelira(cnsinr//'.CNSV', 'LONMAX', neq)
        neq = neq/ncmp
        if (type(1:1) .ne. typrez) then
            valk(1) = cham19
            valk(2) = type
            valk(3) = typrez
            call utmess('F', 'TEST0_9', nk=3, valk=valk)
            goto 999
        end if

!
    else if (tych(1:2) .eq. 'EL') then
!   -----------------------------------
        cnsin1 = '&&TRESU_CH.CNSIN1'
        cnsinr = '&&TRESU_CH.CNSINR'
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
        call wkvect('&&TRESU_CH.LST', 'V V I', nbma, vi=v_list)
        nbma_list = 0
        if (l_parallel_mesh) then
            do i = 1, nbma
                if (v_maex(i) == rank) then
                    nbma_list = nbma_list+1
                    v_list(nbma_list) = i
                end if
            end do
        else
            do i = 1, nbma
                v_list(i) = i
            end do
            nbma_list = nbma
        end if
        call celces(cham19, 'V', cnsin1)
!       -- tres important : cesred avec nbcmp permet la division neq=neq/ncmp
        call cesred(cnsin1, nbma_list, v_list, 1, nocmp(1), &
                    'V', cnsinr)
        call detrsd('CHAM_ELEM_S', cnsin1)
        call jedetr('&&TRESU_CH.LST')

        call jeveuo(cnsinr//'.CESV', 'L', jcsv)
        call jeveuo(cnsinr//'.CESC', 'L', jcsc)
        call jeveuo(cnsinr//'.CESL', 'L', jcsl)
        call jeveuo(cnsinr//'.CESD', 'L', jcsd)
        ncmp = zi(jcsd-1+2)
        ASSERT(ncmp .eq. 1)
        do i = 1, nbcmp
            noddl = nocmp(i)
            do j = 1, ncmp
                if (zk8(jcsc-1+j) .eq. noddl) then
                    zi(jcmp-1+i) = j
                    goto 20
                end if
            end do
            call utmess('F', 'CHAMPS_3', sk=noddl)
20          continue
        end do
        call jelira(cnsinr//'.CESV', 'TYPE', cval=type)
        call jelira(cnsinr//'.CESV', 'LONMAX', neq)
        neq = neq/ncmp
        if (type(1:1) .ne. typrez) then
            valk(1) = cham19
            valk(2) = type
            valk(3) = typrez
            call utmess('F', 'TEST0_9', nk=3, valk=valk)
            goto 999
        end if
    else
        call utmess('F', 'TEST0_10', sk=cham19)
    end if

    nl1 = lxlgut(lign1)
    lign1(1:nl1+16) = lign1(1:nl1-1)//' NOM_CMP'
    lign1(nl1+17:nl1+17) = '.'

!    ================================================================================
    if (type .eq. 'I') then
        if (typtes .eq. 'SOMM_ABS') then
            vali = 0
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        vali = vali+abs(zi(jcsv-1+ind))
                    end if
                end do
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'I', sci=vali)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'
!
        else if (typtes .eq. 'SOMM') then
            vali = 0
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        vali = vali+zi(jcsv-1+ind)
                    end if
                end do
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'I', sci=vali)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'
!
        else if (typtes .eq. 'MAX') then
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valii = zi(jcsv-1+ind)
                        goto 124
                    end if
                end do
124             continue
                do k = j+1, neq
                    ind = ncmp*(k-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valii = max(valii, zi(jcsv-1+ind))
                    end if
                end do
                if (i .eq. 1) then
                    vali = valii
                    icmp = 1
                else
                    if (valii .gt. vali) then
                        vali = valii
                        icmp = i
                    end if
                end if
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MAX', 'I', sci=vali)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp- &
                                                                 1+icmp))
            lign2(nl2+17:nl2+17) = '.'
!
        else if (typtes .eq. 'MIN') then
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valii = zi(jcsv-1+ind)
                        goto 134
                    end if
                end do
134             continue
                do k = j+1, neq
                    ind = ncmp*(k-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valii = min(valii, zi(jcsv-1+ind))
                        icmp = vnocmp
                    end if
                end do
                if (i .eq. 1) then
                    vali = valii
                    icmp = 1
                else
                    if (valii .lt. vali) then
                        vali = valii
                        icmp = i
                    end if
                end if
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MIN', 'I', sci=vali)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp- &
                                                                 1+icmp))
            lign2(nl2+17:nl2+17) = '.'
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
!

!    ================================================================================
    else if (type .eq. 'R') then
        if (typtes .eq. 'SOMM_ABS') then
            valr = 0.d0
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valr = valr+abs(zr(jcsv-1+ind))
                    end if
                end do
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'R', scr=valr)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'
        else if (typtes .eq. 'SOMM') then
            valr = 0.d0
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valr = valr+zr(jcsv-1+ind)
                    end if
                end do
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'R', scr=valr)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'

        else if (typtes .eq. 'MAX') then
            valrr = -1.d+300
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do k = 1, neq
                    ind = ncmp*(k-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valrr = max(valrr, zr(jcsv-1+ind))
                    end if
                end do
                if (i .eq. 1) then
                    valr = valrr
                    icmp = 1
                else
                    if (valrr .gt. valr) then
                        valr = valrr
                        icmp = i
                    end if
                end if
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MAX', 'R', scr=valr)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp-1+icmp))
            lign2(nl2+17:nl2+17) = '.'

        else if (typtes .eq. 'MIN') then
            valrr = 1.d+300
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do k = 1, neq
                    ind = ncmp*(k-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valrr = min(valrr, zr(jcsv-1+ind))
                    end if
                end do
                if (i .eq. 1) then
                    valr = valrr
                    icmp = 1
                else
                    if (valrr .lt. valr) then
                        valr = valrr
                        icmp = i
                    end if
                end if
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MIN', 'R', scr=valr)
            end if
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp- &
                                                                 1+icmp))
            lign2(nl2+17:nl2+17) = '.'
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
!

!    ================================================================================
    else if (type .eq. 'C') then
        if (l_parallel_mesh) then
            ASSERT(ASTER_FALSE)
        end if
        if (typtes .eq. 'SOMM_ABS') then
            valr = 0.d0
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valr = valr+abs(zc(jcsv-1+ind))
                    end if
                end do
            end do
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'
        else if (typtes .eq. 'SOMM') then
            valc = dcmplx(0.d0, 0.d0)
            do i = 1, nbcmp
                vnocmp = zi(jcmp+i-1)
                do j = 1, neq
                    ind = ncmp*(j-1)+(vnocmp-1)+1
                    if (zl(jcsl-1+ind)) then
                        valc = valc+zc(jcsv-1+ind)
                    end if
                end do
            end do
            nl2 = lxlgut(lign2)
            lign2(1:nl2+16) = lign2(1:nl2-1)//' '//zk8(jcsc-1+zi(jcmp))
            lign2(nl2+17:nl2+17) = '.'
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
    end if

    nl1 = lxlgut(lign1)
    nl11 = lxlgut(lign1(1:nl1-1))
    nl2 = lxlgut(lign2)
    nl22 = lxlgut(lign2(1:nl2-1))

    if (llab) then
        if (nl11 .lt. 80) then
            write (ific, *) lign1(1:nl11)
        else if (nl11 .lt. 160) then
            write (ific, 160) lign1(1:80), lign1(81:nl11)
        else
            write (ific, 120) lign1(1:80), lign1(81:160), lign1(161: &
                                                                nl11)
        end if
!
        if (nl22 .lt. 80) then
            write (ific, *) lign2(1:nl22)
        else if (nl22 .lt. 160) then
            write (ific, 160) lign2(1:80), lign2(81:nl22)
        else
            write (ific, 120) lign2(1:80), lign2(81:160), lign2(161: &
                                                                nl22)
        end if
    end if
!
    call tresu_print_all(tbtxt(1), tbtxt(2), llab, typres, nbref, &
                         crit, epsi, ssigne, refr, valr, &
                         refi, vali, refc, valc, ignore=skip, &
                         compare=ordgrd)
!
    if (tych(1:4) .eq. 'NOEU') then
        call detrsd('CHAM_NO_S', cnsinr)
    else if (tych(1:2) .eq. 'EL') then
        call detrsd('CHAM_ELEM_S', cnsinr)
    end if

999 continue
    call jedetr('&&TRESU_CH.CMP')
!
160 format(1x, a80, a)
120 format(1x, 2(a80), a)
!
    call jedema()
end subroutine
