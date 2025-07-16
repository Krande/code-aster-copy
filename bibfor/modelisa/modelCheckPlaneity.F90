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
subroutine modelCheckPlaneity(mesh, model)
!
    use model_module, only: getPlateCell
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/asmpi_barrier_wrap.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_split_comm.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim2.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/nmiret.h"
#include "asterfort/plateGeom_module.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh, model
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_MODELE
!
! Check planeity of plate elements
!
! --------------------------------------------------------------------------------------------------
!
! In  model           : name of the model
! In  mesh           : name of the mesh
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: tabret(0:10)
    character(len=16) :: answer
    aster_logical :: lparallel_mesh
    character(len=1), parameter :: base = "V"
    character(len=19), parameter :: plateLigrel = "&&OP0018.PLATES"
    integer(kind=8), pointer :: cellPlate(:) => null()
    character(len=19) :: modelLigrel
    integer(kind=8) :: nbCellPlate
    character(len=16), parameter :: option = "VERI_PLAN"
    integer(kind=8), parameter :: nbIn = 2, nbOut = 2
    character(len=8) :: lpain(nbIn), lpaout(nbOut)
    character(len=24) :: lchin(nbIn), lchout(nbOut)
    character(len=24) :: chgeom
    character(len=24), parameter :: paraCheck = "&&OP0018.PARACHECK"
    integer(kind=8), parameter :: nbCmp = 2
    character(len=8), parameter :: cmpName(nbCmp) = (/"X1", "X2"/)
    real(kind=8) :: cmpVale(nbCmp)
    character(len=24), parameter :: codret = '&&OP0018.CODRET', indicr = '&&OP0018.INDICR'
    mpi_int :: world, newcom, color, key, ierror
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
    lparallel_mesh = answer .eq. 'OUI'

! - Get list of cells with plate model
    nbCellPlate = 0
    call getPlateCell(model, nbCellPlate, cellPlate)

#ifdef ASTER_HAVE_MPI
! - Split MPI communicator
    call asmpi_comm('GET', world)
    if (nbCellPlate .ne. 0) then
! ----- Plate case
        color = 1
        key = 1
        call asmpi_split_comm(world, color, key, "TMPCOMM", newcom)
        call asmpi_comm('SET', newcom)
    else
! ----- No plate case
        color = 0
        key = 0
        call asmpi_split_comm(world, color, key, "TMPCOMM", newcom)
        call asmpi_comm('SET', newcom)
    end if
#endif

! - Check planeity of quadrangles
    if (nbCellPlate .ne. 0) then
! ----- Create LIGREL for cells with plate model
        call exlim2(cellPlate, nbCellPlate, modelLigrel, base, plateLigrel)

! ----- Create map for parameters
        cmpVale(1) = BASE_TOLE_PLANE
        cmpVale(2) = 1.d0
        call mecact('V', paraCheck, 'LIGREL', plateLigrel, 'NEUT_R', &
                    ncmp=nbCmp, lnomcmp=cmpName, vr=cmpVale)

! ----- Input fields
        call megeom(model, chgeom)
        lchin(1) = chgeom
        lpain(1) = 'PGEOMER'
        lchin(2) = paraCheck
        lpain(2) = 'PCHCKPR'

! ----- Output fields
        lchout(1) = codret
        lpaout(1) = 'PCODRET'
        lchout(2) = indicr
        lpaout(2) = 'PINDICR'

! ----- Compute option
        call calcul('C', option, plateLigrel, nbIn, lchin, &
                    lpain, nbOut, lchout, lpaout, base, &
                    'NON')
        call nmiret(codret, tabret)
        if (tabret(0)) then
            call utmess('A', 'PLATE1_81')
        end if
    end if
!
! - Unsplit MPI communicator
#ifdef ASTER_HAVE_MPI
    call asmpi_barrier_wrap(world, ierror)
    call asmpi_comm('SET', world)
    if (newcom .ne. 0) then
        call asmpi_comm('FREE', newcom)
    end if
#endif
!
    AS_DEALLOCATE(vi=cellPlate)
!
end subroutine
