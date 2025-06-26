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

subroutine tresu_champ_val(cham19, nomail, nonoeu, nupo, nusp, &
                           ivari, nocmp, nbref, tbtxt, refi, &
                           refr, refc, typres, epsi, crit, &
                           llab, ssigne, ignore, compare)
    implicit none
#include "asterf_types.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/tresu_print_all.h"
#include "asterfort/utch19.h"
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
#include "asterf_mpi.h"
#endif
!
    character(len=*), intent(in) :: cham19
    character(len=*), intent(in) :: nomail
    character(len=*), intent(in) :: nonoeu
    integer(kind=8), intent(in) :: nupo
    integer(kind=8), intent(in) :: nusp
    integer(kind=8), intent(in) :: ivari
    character(len=*), intent(in) :: nocmp
    integer(kind=8), intent(in) :: nbref
    character(len=16), intent(in) :: tbtxt(2)
    integer(kind=8), intent(in) :: refi(nbref)
    real(kind=8), intent(in) :: refr(nbref)
    complex(kind=8), intent(in) :: refc(nbref)
    character(len=*), intent(in) :: typres
    real(kind=8), intent(in) :: epsi
    character(len=*), intent(in) :: crit
    aster_logical, intent(in) :: llab
    character(len=*), intent(in) :: ssigne
    aster_logical, intent(in), optional :: ignore
    real(kind=8), intent(in), optional :: compare
! IN  : CHAM19 : NOM DU CHAM_ELEM DONT ON DESIRE VERIFIER 1 COMPOSANTE
! IN  : NOMAIL : NOM DE LA MAILLE A TESTER
! IN  : NONOEU : NOM D'UN NOEUD (POUR LES CHAM_ELEM "AUX NOEUDS").
!                 (SI CE NOM EST BLANC : ON UTILISE NUPO)
! IN  : NUPO   : NUMERO DU POINT A TESTER SUR LA MAILLE NOMAIL
! IN  : NUSP   : NUMERO DU SOUS_POINT A TESTER SUR LA MAILLE NOMAIL
!                (SI NUSP=0 : IL N'Y A PAS DE SOUS-POINT)
! IN  : IVARI   : NUMERO DE LA CMP (POUR VARI_R)
! IN  : NOCMP  : NOM DU DDL A TESTER SUR LE POINT NUPO
! IN  : NBREF  : NOMBRE DE VALEURS DE REFERENCE
! IN  : TBTXT  : (1)=REFERENCE, (2)=LEGENDE
! IN  : REFR   : VALEUR REELLE ATTENDUE SUR LE DDL DU POINT.
! IN  : REFC   : VALEUR COMPLEXE ATTENDUE SUR LE DDL DU POINT.
! IN  : CRIT   : 'RELATIF' OU 'ABSOLU'(PRECISION RELATIVE OU ABSOLUE).
! IN  : EPSI   : PRECISION ESPEREE
! IN  : LLAB   : FLAG D IMPRESSION DE LABELS
! OUT : IMPRESSION SUR LISTING
! ----------------------------------------------------------------------
    integer(kind=8) :: vali, ier, rank
    real(kind=8) :: valr
    complex(kind=8) :: valc
    character(len=8) :: nomma
    aster_logical :: skip, l_parallel_mesh
    real(kind=8) :: ordgrd
    mpi_int :: irank
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
    call dismoi('NOM_MAILLA', cham19, 'CHAM_ELEM', repk=nomma)
    l_parallel_mesh = isParallelMesh(nomma)
!
    call utch19(cham19, nomma, nomail, nonoeu, nupo, &
                nusp, ivari, nocmp, typres, valr, &
                valc, vali, ier)

    if (l_parallel_mesh) then
        rank = -1
        if (ier == 0) then
            call asmpi_info(rank=irank)
            rank = irank
        end if
        call asmpi_comm_vect('MPI_MIN', 'I', sci=ier)
        call asmpi_comm_vect('MPI_MAX', 'I', sci=rank)
    end if

    ASSERT(ier .eq. 0)

    if (l_parallel_mesh) then
        if (typres == 'R') then
            call asmpi_comm_vect('BCAST', 'R', bcrank=rank, scr=valr)
        elseif (typres == 'I') then
            call asmpi_comm_vect('BCAST', 'I', bcrank=rank, sci=vali)
        elseif (typres == 'C') then
            call asmpi_comm_vect('BCAST', 'C', bcrank=rank, scc=valc)
        end if
    end if

    call tresu_print_all(tbtxt(1), tbtxt(2), llab, typres, nbref, &
                         crit, epsi, ssigne, refr, valr, &
                         refi, vali, refc, valc, ignore=skip, &
                         compare=ordgrd)
!
end subroutine
