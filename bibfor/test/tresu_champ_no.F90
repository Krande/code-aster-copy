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

subroutine tresu_champ_no(cham19, nonoeu, nocmp, nbref, tbtxt, &
                          refi, refr, refc, typres, epsi, &
                          crit, llab, ssigne, ignore, compare)
    implicit none
#include "asterc/indik8.h"
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
#include "jeveux.h"
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
#include "asterf_mpi.h"
#endif
!
    character(len=19), intent(in) :: cham19
    character(len=33), intent(in) :: nonoeu
    character(len=8), intent(in) :: nocmp
    integer(kind=8), intent(in) :: nbref
    character(len=16), intent(in) :: tbtxt(2)
    integer(kind=8), intent(in) :: refi(nbref)
    real(kind=8), intent(in) :: refr(nbref)
    complex(kind=8), intent(in) :: refc(nbref)
    character(len=1), intent(in) :: typres
    real(kind=8), intent(in) :: epsi
    character(len=*), intent(in) :: crit
    aster_logical, intent(in) :: llab
    character(len=*), intent(in) :: ssigne
    aster_logical, intent(in), optional :: ignore
    real(kind=8), intent(in), optional :: compare
!     ENTREES:
!        CHAM19 : NOM DU CHAM_NO DONT ON DESIRE VERIFIER 1 COMPOSANTE
!        NONOEU : NOM DU NOEUD A TESTER
!        NOCMP  : NOM DU DDL A TESTER SUR LE NOEUD NONOEU
!        TBTXT  : (1)=REFERENCE, (2)=LEGENDE
!        NBREF  : NOMBRE DE VALEURS DE REFERENCE
!        REFR   : VALEUR REELLE ATTENDUE SUR LE DDL DU NOEUD
!        REFC   : VALEUR COMPLEXE ATTENDUE SUR LE DDL DU NOEUD
!        CRIT   : 'RELATIF' OU 'ABSOLU'(PRECISION RELATIVE OU ABSOLUE).
!        EPSI   : PRECISION ESPEREE
!        LLAB   : AFFICHAGE DES LABELS
!     SORTIES:
!      LISTING ...
! ----------------------------------------------------------------------
!     FONCTIONS EXTERNES:
!     VARIABLES LOCALES:
    character(len=8) :: nogd, mesh
    character(len=1) :: type
    integer(kind=8) :: gd, iadg, vali, ser, inog, rank
    real(kind=8) :: valr
    complex(kind=8) :: valc
    character(len=8) :: nomma
    character(len=19) :: numeq, valk(3)
    character(len=24) :: nolili
    aster_logical :: skip, l_parallel_mesh, l_ok
    real(kind=8) :: ordgrd
    mpi_int :: irank
    integer(kind=8), pointer :: v_noex(:) => null()
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iancmp, ianueq, iaprno, iavale
    integer(kind=8) :: icmp, idecal, iicmp, ino, ival, ncmp
    integer(kind=8) :: ncmpmx, nec
!-----------------------------------------------------------------------
    call jemarq()
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
    call dismoi('NOM_MAILLA', cham19, 'CHAMP', repk=mesh, arret='F')
    l_parallel_mesh = isParallelMesh(mesh)
!
    call dismoi('NOM_GD', cham19, 'CHAM_NO', repk=nogd)
    call dismoi('NUM_GD', cham19, 'CHAM_NO', repi=gd)
    call dismoi('NOM_MAILLA', cham19, 'CHAM_NO', repk=nomma)
    call dismoi('NUME_EQUA', cham19, 'CHAM_NO', repk=numeq)
!
    call jelira(cham19//'.VALE', 'TYPE', cval=type)
    if (type .ne. typres) then
        valk(1) = cham19
        valk(2) = type
        valk(3) = typres
        call utmess('F', 'CALCULEL6_89', nk=3, valk=valk)
    else if (type .ne. 'R' .and. type .ne. 'C') then
        valk(1) = type
        call utmess('F', 'CALCULEL6_90', sk=valk(1))
    end if
    call jeveuo(cham19//'.VALE', 'L', iavale)
!
    nec = nbec(gd)
!
!     -- ON RECHERCHE LE NUMERO CORRESPONDANT A NOCMP:
    call jelira(jexnum('&CATA.GD.NOMCMP', gd), 'LONMAX', ncmpmx)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd), 'L', iancmp)
    icmp = indik8(zk8(iancmp), nocmp, 1, ncmpmx)
    if (icmp .eq. 0) then
        valk(1) = nocmp
        valk(2) = nogd
        call utmess('F', 'CALCULEL6_91', nk=2, valk=valk)
    end if
!
!        -- RECUPERATION DU NUMERO DU NOEUD:
    l_ok = ASTER_FALSE
    ino = 0
    if (nonoeu(1:8) .ne. ' ') then
        ino = char8_to_int(nonoeu(1:8))
    end if
    inog = ino
    call asmpi_comm_vect('MPI_MAX', 'I', sci=inog)
    if (inog .eq. 0) then
        valk(1) = nonoeu(1:8)
        call utmess('F', 'CALCULEL6_92', sk=valk(1))
    end if
!
    if (ino == 0) then
        go to 100
    end if
!
    if (l_parallel_mesh) then
        call asmpi_info(rank=irank)
        rank = to_aster_int(irank)
        call jeveuo(mesh//'.NOEX', 'L', vi=v_noex)
        if (v_noex(ino) .ne. rank) then
            go to 100
        end if
    end if
!
!        --SI LE CHAMP EST DECRIT PAR 1 "PRNO":
!
    call jenuno(jexnum(numeq//'.LILI', 1), nolili)
    call jeveuo(jexnum(numeq//'.PRNO', 1), 'L', iaprno)
    call jeveuo(numeq//'.NUEQ', 'L', ianueq)
!
!        IVAL : ADRESSE DU DEBUT DU NOEUD INO DANS .NUEQ
!        NCMP : NOMBRE DE COMPOSANTES PRESENTES SUR LE NOEUD
!        IADG : DEBUT DU DESCRIPTEUR GRANDEUR DU NOEUD INO
    ival = zi(iaprno-1+(ino-1)*(nec+2)+1)
    ncmp = zi(iaprno-1+(ino-1)*(nec+2)+2)
    iadg = iaprno-1+(ino-1)*(nec+2)+3
    ASSERT(ncmp .ne. 0)
!
!        -- ON COMPTE LES CMP PRESENTES SUR LE NOEUD AVANT ICMP:
    idecal = 0
    do iicmp = 1, icmp
        if (exisdg(zi(iadg), iicmp)) idecal = idecal+1
    end do
!
    if (exisdg(zi(iadg), icmp)) then
        if (type .eq. 'R') then
            valr = zr(iavale-1+zi(ianueq-1+ival-1+idecal))
        else if (type .eq. 'I') then
            vali = zi(iavale-1+zi(ianueq-1+ival-1+idecal))
        else if (type .eq. 'C') then
            valc = zc(iavale-1+zi(ianueq-1+ival-1+idecal))
        end if
        l_ok = ASTER_TRUE
    end if
!
100 continue
!
!
    if (l_parallel_mesh) then
        ser = 0
        rank = -1
        if (l_ok) then
            ser = 1
            rank = to_aster_int(irank)
        end if
        call asmpi_comm_vect('MPI_MAX', 'I', sci=ser)
        call asmpi_comm_vect('MPI_MAX', 'I', sci=rank)
        if (ser == 1) then
            l_ok = ASTER_TRUE
        else
            l_ok = ASTER_FALSE
        end if
    end if
!
    if (l_ok) then
        if (l_parallel_mesh) then
            if (type == 'R') then
                call asmpi_comm_vect('BCAST', 'R', bcrank=rank, scr=valr)
            elseif (type == 'I') then
                call asmpi_comm_vect('BCAST', 'I', bcrank=rank, sci=vali)
            elseif (type == 'C') then
                call asmpi_comm_vect('BCAST', 'C', bcrank=rank, scc=valc)
            end if
        end if
        call tresu_print_all(tbtxt(1), tbtxt(2), llab, type, nbref, &
                             crit, epsi, ssigne, refr, valr, &
                             refi, vali, refc, valc, ignore=skip, &
                             compare=ordgrd)
    else
        call utmess('F', 'CALCULEL6_93')
    end if
!
    call jedema()
end subroutine
