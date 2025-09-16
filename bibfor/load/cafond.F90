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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cafond(load, mesh, model, geomDime, valeType)
!
    implicit none
!
#include "asterf_types.h"
#include "LoadTypes_type.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/char_crea_cart.h"
#include "asterfort/char_read_val.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlim1.h"
#include "asterfort/getelem.h"
#include "asterfort/infniv.h"
#include "asterfort/inical.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mesomm.h"
#include "asterfort/nocart.h"
#include "asterfort/peair1.h"
#include "asterfort/vetyma.h"
#include "asterfort/dismoi.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Keyword = 'EFFE_FOND'
!
! --------------------------------------------------------------------------------------------------
!
!
! In  mesh      : mesh
! In  load      : load
! In  geomDime  : space dimension
! In  model     : model
! In  valeType  : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbout = 1, nbin = 1
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=16), parameter :: keywordfact = 'EFFE_FOND'
    character(len=24), parameter :: listCellHole = '&&CAFOND.LISTHOLE'
    character(len=24), parameter :: listCellSect = '&&CAFOND.LISTSECT'
    character(len=16), parameter :: option = 'CARA_SECT_POUT3'
    character(len=19), parameter :: ligrel = '&&CAFOND.LIGREL'
    integer(kind=8) :: npres, iocc, ii
    integer(kind=8) :: ifm, niv, val_nb, jvalv
    integer(kind=8) :: rang
    aster_logical :: isPMesh
    real(kind=8) :: r8dummy
    real(kind=8) :: hole_area, cara_geom(10), mate_area, coef_mult
    complex(kind=8) :: c16dummy
    character(len=8) :: pres_fonc
    real(kind=8) :: pres_real
    character(len=16) :: k16dummy, answer
    integer(kind=8) :: jvCellHole, jvCellSect
    integer(kind=8) :: nbCellHole, nbCellSect
    character(len=8) :: suffix
    character(len=19) :: map(LOAD_MAP_NBMAX)
    integer(kind=8) :: nbMap, nbCmp(LOAD_MAP_NBMAX)
    integer(kind=8), pointer :: maex(:) => null()
    mpi_int :: mrank
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
    call getfac(keywordfact, npres)
    if (npres .eq. 0) goto 99

! - Warning for COQUE_SOLIDE
    call dismoi('EXI_COQSOL', model, 'MODELE', repk=answer)
    if (answer .eq. 'OUI') then
        call utmess('F', 'SOLIDSHELL1_5')
    end if

    isPMesh = isParallelMesh(mesh)
    if (isPMesh) then
        call jeveuo(mesh//'.MAEX', 'L', vi=maex)
        call asmpi_info(rank=mrank)
        rang = to_aster_int(mrank)
    end if
!
! - Creation and initialization to zero of <CARTE>
!
    call char_crea_cart('MECANIQUE', keywordfact, load, mesh, valeType, &
                        nbMap, map, nbCmp)
    ASSERT(nbMap .eq. 2)
!
! - For computation of geometric caracteristics (elementary)
!
    call inical(nbin, lpain, lchin, nbout, lpaout, lchout)
    lpain(1) = 'PGEOMER'
    lchin(1) = mesh//'.COORDO'
    lpaout(1) = 'PCASECT'
    lchout(1) = '&&CAFOND.PSECT'
    do iocc = 1, npres

! ----- Elements for hole
        suffix = '_INT'
        call getelem(mesh, keywordfact, iocc, 'F', listCellHole, &
                     nbCellHole, suffix=suffix)
        call jeveuo(listCellHole, 'L', jvCellHole)

! ----- Elements for section
        call getelem(mesh, keywordfact, iocc, 'F', listCellSect, &
                     nbCellSect)
        call jeveuo(listCellSect, 'L', jvCellSect)

! ----- Create <LIGREL>
        call exlim1(zi(jvCellSect), nbCellSect, model, 'V', ligrel)

! ----- Get pressure
        call char_read_val(keywordfact, iocc, 'PRES', valeType, val_nb, &
                           pres_real, pres_fonc, c16dummy, k16dummy)
        ASSERT(val_nb .eq. 1)

! ----- Area of hole
        call peair1(mesh, nbCellHole, zi(jvCellHole), hole_area, r8dummy)

! ----- To compute area of material section
        call calcul('S', option, ligrel, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, 'V', &
                    'OUI')
        call mesomm(lchout(1), 10, vr=cara_geom, local=isPMesh)
        call detrsd('LIGREL', ligrel)

! ----- Multiplicative ratio of pressure
        mate_area = cara_geom(1)
        if (isPMesh) then
            call asmpi_comm_vect('MPI_SUM', 'R', scr=mate_area)
        end if
        coef_mult = -hole_area/mate_area
        if (niv .eq. 2) then
            write (ifm, *) 'SURFACE DU TROU    ', hole_area
            write (ifm, *) 'SURFACE DE MATIERE ', mate_area
        end if

! ----- Affectation of values in <CARTE> - Multiplicative ratio of pressure
        call jeveuo(map(1)//'.VALV', 'E', jvalv)
        zr(jvalv-1+1) = coef_mult
        call nocart(map(1), 3, nbCmp(1), mode='NUM', nma=nbCellSect, &
                    limanu=zi(jvCellSect))

! ----- Affectation of values in <CARTE> - Pressure
        call jeveuo(map(2)//'.VALV', 'E', jvalv)
        if (valeType .eq. 'REEL') then
            zr(jvalv-1+1) = pres_real
        else if (valeType .eq. 'FONC') then
            zk8(jvalv-1+1) = pres_fonc
        else
            ASSERT(ASTER_FALSE)
        end if
        call nocart(map(2), 3, nbCmp(2), mode='NUM', nma=nbCellSect, &
                    limanu=zi(jvCellSect))

! ----- Check elements
        call vetyma(mesh, geomDime, keywordfact, listCellSect, nbCellSect)
!
        call jedetr(listCellHole)
        call jedetr(listCellSect)
!
    end do
!
99  continue
!
    call jedema()
end subroutine
