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

subroutine tresu_champ_all(chamgd, typtes, typres, nbref, tbtxt, &
                           refi, refr, refc, epsi, crit, &
                           llab, ssigne, ignore, compare)
    implicit none
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
#include "asterfort/cesred.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
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
    character(len=*), intent(in) :: crit
    aster_logical, intent(in) :: llab
    character(len=*), intent(in) :: ssigne
    aster_logical, intent(in), optional :: ignore
    real(kind=8), intent(in), optional :: compare
! IN  : CHAMGD : NOM DU CHAM_GD
! IN  : TYPTES : TYPE DE TEST A EFFECTUER SUR LE CHAMP
! IN  : NBREF  : NOMBRE DE VALEURS DE REFERENCE
! IN  : TBTXT  : (1)=REFERENCE, (2)=LEGENDE
! IN  : REFI   : VALEUR REELLE ENTIERE ATTENDUE
! IN  : REFR   : VALEUR REELLE ATTENDUE
! IN  : REFC   : VALEUR COMPLEXE ATTENDUE
! IN  : CRIT   : 'RELATIF' OU 'ABSOLU'(PRECISION RELATIVE OU ABSOLUE).
! IN  : EPSI   : PRECISION ESPEREE
! IN  : LLAB   : FLAG D IMPRESSION DES LABELS
! OUT : IMPRESSION SUR LISTING
! ----------------------------------------------------------------------
    integer(kind=8) :: vali, jvale, neq, i, iret1, iret2, rank, jvale2, neq2
    integer(kind=8) :: nbma, nbma_list
    real(kind=8) :: valr, ordgrd
    complex(kind=8) :: valc
    character(len=1) :: typrez
    character(len=24) :: valk(3)
    character(len=4) :: type
    character(len=8) :: mesh
    character(len=5) :: sufv
    character(len=19) :: cham19, nume_equa, cnsin1, cnsinr
    aster_logical :: skip, l_parallel_mesh, cham_no
    mpi_int :: irank
    integer(kind=8), pointer :: v_noex(:) => null()
    integer(kind=8), pointer :: v_maex(:) => null()
    integer(kind=8), pointer :: v_list(:) => null()
    integer(kind=8), pointer :: v_deeq(:) => null()
!     ------------------------------------------------------------------
    if (present(ignore)) then
        skip = ignore
    else
        skip = .false.
    end if
    if (present(compare)) then
        ordgrd = compare
    else
        ordgrd = 1.d0
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
!     -- LE CHAMP EXISTE-T-IL ?
!     =========================
    sufv = ' '
    call jeexin(cham19//'.VALE', iret1)
    if (iret1 .gt. 0) then
        sufv = '.VALE'
        cham_no = ASTER_TRUE
    else
        cham_no = ASTER_FALSE
        call jeexin(cham19//'.CELV', iret2)
        if (iret2 .gt. 0) then
            sufv = '.CELV'
        end if
    end if
    ASSERT(sufv .ne. ' ')
!
!
    call jelira(cham19//sufv, 'TYPE', cval=type)
    if (type(1:1) .ne. typrez) then
        valk(1) = cham19
        valk(2) = type
        valk(3) = typrez
        call utmess('F', 'TEST0_9', nk=3, valk=valk)
        goto 999
    end if
!
    call jelira(cham19//sufv, 'LONMAX', neq)
    call jeveuo(cham19//sufv, 'L', jvale)
!
! --- Pour un ParallelMesh il faut retailler les champs
!
    if (l_parallel_mesh) then
        call wkvect('&&TRESU_CH.ALL', 'V V '//type(1:1), neq, jvale2)
        neq2 = 0
        if (cham_no) then
            call dismoi('NUME_EQUA', cham19, 'CHAM_NO', repk=nume_equa, arret='F')
            call jeveuo(nume_equa//".DEEQ", 'L', vi=v_deeq)
            do i = 1, neq
                if (v_deeq(2*(i-1)+1) <= 0) then
                    ! On s'arrete en erreur dès que noeud = 0 car on ne sait pas comment le
                    ! sommer
                    ASSERT(v_deeq(2*(i-1)+1) == 0)
                    if (typtes(1:4) == "SOMM") then
                        call utmess('F', 'TEST0_21', sk=typtes)
                    end if
                end if
!
                if (v_deeq(2*(i-1)+1) > 0 .and. v_noex(v_deeq(2*(i-1)+1)) == rank) then
                    neq2 = neq2+1
                    if (type == 'I') then
                        zi(jvale2-1+neq2) = zi(jvale-1+i)
                    elseif (type == 'R') then
                        zr(jvale2-1+neq2) = zr(jvale-1+i)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
                end if
            end do
        else
            cnsin1 = '&&TRESU_CH.CNSIN1'
            cnsinr = '&&TRESU_CH.CNSINR'
            call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
            call wkvect('&&TRESU_CH.LST', 'V V I', nbma, vi=v_list)
            nbma_list = 0
            do i = 1, nbma
                if (v_maex(i) == rank) then
                    nbma_list = nbma_list+1
                    v_list(nbma_list) = i
                end if
            end do
            call celces(cham19, 'V', cnsin1)
            call cesred(cnsin1, nbma_list, v_list, 0, ["XXX"], 'V', cnsinr)
            call jeveuo(cnsinr//'.CESV', 'L', jvale2)
            call jelira(cnsinr//'.CESV', 'LONMAX', neq2)
            call detrsd('CHAM_ELEM_S', cnsin1)
            call jedetr('&&TRESU_CH.LST')
!
        end if
!
        jvale = jvale2
        neq = neq2
    end if
!
!
    if (type .eq. 'I') then
        if (typtes .eq. 'SOMM_ABS') then
            vali = 0
            do i = 1, neq
                vali = vali+abs(zi(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'I', sci=vali)
            end if
        else if (typtes .eq. 'SOMM') then
            vali = 0
            do i = 1, neq
                vali = vali+zi(jvale+i-1)
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'I', sci=vali)
            end if
        else if (typtes .eq. 'MAX') then
            vali = -ismaem()
            do i = 1, neq
                vali = max(vali, zi(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MAX', 'I', sci=vali)
            end if
        else if (typtes .eq. 'MIN') then
            vali = ismaem()
            do i = 1, neq
                vali = min(vali, zi(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MIN', 'I', sci=vali)
            end if
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
!
!
    else if (type .eq. 'R') then
        if (typtes .eq. 'SOMM_ABS') then
            valr = 0.d0
            do i = 1, neq
                valr = valr+abs(zr(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'R', scr=valr)
            end if
        else if (typtes .eq. 'SOMM') then
            valr = 0.d0
            do i = 1, neq
                valr = valr+zr(jvale+i-1)
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_SUM', 'R', scr=valr)
            end if
        else if (typtes .eq. 'MAX') then
            valr = r8miem()
            do i = 1, neq
                valr = max(valr, zr(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MAX', 'R', scr=valr)
            end if
        else if (typtes .eq. 'MIN') then
            valr = r8maem()
            do i = 1, neq
                valr = min(valr, zr(jvale+i-1))
            end do
            if (l_parallel_mesh) then
                call asmpi_comm_vect('MPI_MIN', 'R', scr=valr)
            end if
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
!
!
    else if (type .eq. 'C') then
        if (l_parallel_mesh) then
            ASSERT(ASTER_FALSE)
        end if
        if (typtes .eq. 'SOMM_ABS') then
            valc = dcmplx(0.d0, 0.d0)
            do i = 1, neq
                valc = valc+abs(zc(jvale+i-1))
            end do
        else if (typtes .eq. 'SOMM') then
            valc = dcmplx(0.d0, 0.d0)
            do i = 1, neq
                valc = valc+zc(jvale+i-1)
            end do
        else
            call utmess('F', 'TEST0_8', sk=typtes)
            goto 999
        end if
    end if
!
    call jedetr('&&TRESU_CH.ALL')
!
    call tresu_print_all(tbtxt(1), tbtxt(2), llab, typres, nbref, &
                         crit, epsi, ssigne, refr, valr, &
                         refi, vali, refc, valc, ignore=skip, &
                         compare=ordgrd)
!
999 continue
!
    call jedema()
end subroutine
