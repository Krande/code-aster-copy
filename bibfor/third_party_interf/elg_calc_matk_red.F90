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

subroutine elg_calc_matk_red(mat1z, solv1z, mat2z, bas1)
#include "asterf_petsc.h"
    use aster_petsc_module
    implicit none
! aslint: disable=
! person_in_charge: natacha.bereux at edf.fr
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/apetsc.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/elg_calc_matm_red.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
    character(len=*) :: mat1z, mat2z, solv1z
    character(len=1) :: bas1
!--------------------------------------------------------------
! but :
!   calculer la matrice reduite mat2z correspondant a mat1z
!   (solveur/elim_lagr='oui')
! in/jxin  : mat1z : sd_matr_asse avec ses conditions dualisees
!                    a eliminer
! in/jxin  : solv1z : sd_solveur
! in/jxout : mat2z : sd_matr_asse "reduite" (sans lagranges)
! in       : bas1 : 'G'/'V' (pour la création de mat2z)
! remarque : on cree egalement un nume_ddl (sous-terrain) pour
!            mat2z.
!---------------------------------------------------------------
!
#ifdef ASTER_HAVE_PETSC
!
!================================================================
    character(len=19) :: matas1, matas2, solve1
    character(len=1) :: ktyp
    mpi_int :: rang, nbproc
    integer(kind=8) :: iret, ibid, iexi
    real(kind=8) :: rbid(1)
    character(len=24), pointer :: refa(:) => null()
!----------------------------------------------------------------
    call jemarq()
    matas1 = mat1z
    matas2 = mat2z
    solve1 = solv1z
!
    call asmpi_info(rank=rang, size=nbproc)
!    if (nbproc .ne. 1) call utmess('F', 'ELIMLAGR_2')
    call jeveuo(matas1//'.REFA', 'L', vk24=refa)
    if (refa(11) (1:11) .ne. 'MPI_COMPLET') then
        call utmess('F', 'ELIMLAGR_2')
    end if
!
!     -- quelques garde fous :
    call jelira(matas1//'.VALM', 'TYPE', ibid, ktyp)
    if (ktyp .ne. 'R') call utmess('F', 'ELIMLAGR_3')
    call jeexin(matas1//'.CCID', iexi)
    if (iexi .ne. 0) call utmess('F', 'ELIMLAGR_4')
!
!
!   -- mise a jour de matas1.refa(19):
    call jeveuo(matas1//'.REFA', 'E', vk24=refa)
    if (refa(19) .ne. ' ') then
!       Ce n'est peut etre pas tres normal de reduire une matrice qui
!       a deja ete reduite ...
        ASSERT(.false.)
        call detrsd('MATR_ASSE', refa(19))
    end if
    refa(19) = matas2

!
!     1. CALCUL DANS PETSC DES MATRICES NECESSAIRES :
!        Kproj, Tfinal, ...
!     --------------------------------------------------
    call apetsc('ELIM_LAGR', solve1, matas1, rbid, ' ', &
                0_8, 0_8, iret)

    ASSERT(iret .eq. 0)
!
!
!     2. CALCUL DANS L'ESPACE JEVEUX DE LA MATRICE Kproj
!        ET DE SON NUME_DDL => MATAS2 (et NU2)
!     --------------------------------------------------
    call elg_calc_matm_red(matas1, matas2, bas1)
!
!
    call jedema()
#else
    character(len=1) :: kdummy
    call utmess('F', 'ELIMLAGR_1')
    kdummy = mat1z(1:1)//mat2z(1:1)//solv1z(1:1)//bas1(1:1)
#endif
!
end subroutine
