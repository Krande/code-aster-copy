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
module HHO_postpro_module
!
    use HHO_type
    use HHO_basis_module
    use HHO_utils_module
    use HHO_size_module
    use HHO_eval_module
    use HHO_quadrature_module
    use NonLin_Datastructure_type
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/chpchd.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/infniv.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jevech.h"
#include "asterfort/megeom.h"
#include "asterfort/readVector.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/writeVector.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - post-processing
!
! Post-processing function to vizualize HHO results
!
! --------------------------------------------------------------------------------------------------
    public :: hhoPostMeca, hhoPostDeplMeca, hhoPostTher, hhoPostTherElga
    public :: hhoPostMecaGradVari
    private :: hhoPostDeplMecaOP
!
contains
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostMeca(hhoCell, hhoData, nbnodes)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer(kind=8), intent(in) :: nbnodes
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Evaluate HHO cell unknowns at the nodes
!   In hhoCell         : a HHO Cell
!   In hhoData         : information on HHO methods
!   In nbnodes         : number of nodes
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: max_comp = 81
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: total_dofs, cbs, fbs, jvect, i, ino, ndim, idim, comp_dim
        real(kind=8), dimension(MSIZE_CELL_VEC) :: sol_T
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: sol_T_dim
        real(kind=8), dimension(3, max_comp) :: post_sol
!
        ndim = hhoCell%ndim
        sol_T = 0.d0
        sol_T_dim = 0.d0
        post_sol = 0.d0
!
! --- number of dofs
!
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        ASSERT(ndim*nbnodes .le. max_comp)
        comp_dim = cbs/ndim
!
! --- We get the solution on the cell
!
        call readVector('PDEPLPR', cbs, sol_T, total_dofs-cbs)
!
! --- init cell basis
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- Compute the solution in the cell nodes
!
        do idim = 1, ndim
            sol_T_dim(1:comp_dim) = sol_T(1+(idim-1)*comp_dim:idim*comp_dim)
            do ino = 1, nbnodes
                post_sol(idim, ino) = hhoEvalScalCell( &
                                      hhoBasisCell, hhoData%cell_degree(), &
                                      hhoCell%coorno(1:3, ino), sol_T_dim, comp_dim &
                                      )
            end do
        end do
!
! --- Copy of post_sol in PDEPL_R ('OUT' to fill)
!
        call jevech('PDEPL_R', 'E', jvect)
!
        i = 1
        do ino = 1, nbnodes
            do idim = 1, ndim
                zr(jvect-1+i) = post_sol(idim, ino)
                i = i+1
            end do
        end do
!
    end subroutine
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostMecaGradVari(hhoCell, hhoData, nbnodes)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer(kind=8), intent(in) :: nbnodes
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics - GRAD_VARI
!
!   Evaluate HHO cell unknowns at the nodes
!   In hhoCell         : a HHO Cell
!   In hhoData         : information on HHO methods
!   In nbnodes         : number of nodes
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: mk_total_dofs, mk_cbs, mk_fbs, jvect, i, ino, ndim, idim, comp_dim
        integer(kind=8) :: faces_dofs, gv_cbs, gv_fbs, gv_total_dofs
        real(kind=8), dimension(MSIZE_CELL_VEC) :: sol_U
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: sol_U_dim, sol_V, sol_L
        real(kind=8), dimension(5, 27) :: post_sol
!
        ndim = hhoCell%ndim
        sol_U = 0.d0
        sol_U_dim = 0.d0
        sol_V = 0.d0
        sol_L = 0.d0
        post_sol = 0.d0
!
! --- number of dofs
!
        call hhoMecaDofs(hhoCell, hhoData, mk_cbs, mk_fbs, mk_total_dofs)
        call hhoTherDofs(hhoCell, hhoData, gv_cbs, gv_fbs, gv_total_dofs)
        comp_dim = mk_cbs/ndim
        ASSERT(comp_dim == gv_cbs)
        faces_dofs = mk_total_dofs-mk_cbs+gv_total_dofs-gv_cbs
!
! --- We get the solution on the cell
!
        call readVector('PDEPLPR', mk_cbs, sol_U, faces_dofs)
        call readVector('PDEPLPR', gv_cbs, sol_V, faces_dofs+mk_cbs)
        call readVector('PDEPLPR', gv_cbs, sol_L, faces_dofs+mk_cbs+gv_cbs)
!
! --- init cell basis
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- Compute the solution in the cell nodes
!
        do idim = 1, ndim
            sol_U_dim(1:comp_dim) = sol_U(1+(idim-1)*comp_dim:idim*comp_dim)
            do ino = 1, nbnodes
                post_sol(idim, ino) = hhoEvalScalCell( &
                                      hhoBasisCell, hhoData%cell_degree(), &
                                      hhoCell%coorno(1:3, ino), sol_U_dim, comp_dim &
                                      )
            end do
        end do
        do ino = 1, nbnodes
            post_sol(ndim+1, ino) = hhoEvalScalCell( &
                                    hhoBasisCell, hhoData%cell_degree(), &
                                    hhoCell%coorno(1:3, ino), sol_V, gv_cbs &
                                    )
            post_sol(ndim+2, ino) = hhoEvalScalCell( &
                                    hhoBasisCell, hhoData%cell_degree(), &
                                    hhoCell%coorno(1:3, ino), sol_L, gv_cbs &
                                    )
        end do
!
! --- Copy of post_sol in PDEPL_R ('OUT' to fill)
!
        call jevech('PDEPL_R', 'E', jvect)
!
        i = 1
        do ino = 1, nbnodes
            do idim = 1, ndim+2
                zr(jvect-1+i) = post_sol(idim, ino)
                i = i+1
            end do
        end do
!
    end subroutine
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostDeplMecaOP(model_hho, disp, disp_hho_depl)
!
        implicit none
!
        character(len=24), intent(in) :: model_hho
        character(len=24), intent(in) :: disp
        character(len=24), intent(in) :: disp_hho_depl
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Post-treatment (mechanic)
!
! Compute HHO field at nodes - ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  model_hho        : name of hho model
! In  solalg           : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! In  hhoField         : fields for HHO
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ifm, niv
        integer(kind=8), parameter :: nbin = 3
        integer(kind=8), parameter :: nbout = 1
        character(len=8) :: lpain(nbin), lpaout(nbout)
        character(len=19) :: lchin(nbin), lchout(nbout)
        character(len=19) :: ligrel_model, field_elno, celmod, field_noeu
        character(len=16) :: option
        character(len=1) :: base
        character(len=24) :: chgeom
!
! --------------------------------------------------------------------------------------------------
!
        call infniv(ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'HHO2_5')
        end if
!
! --- Initializations
!
        field_elno = '&&FIELD_ELNO'
        field_noeu = '&&FIELD_NOEU'
        celmod = '&&FIELD_CELMOD'
        base = 'V'
        option = 'HHO_DEPL_MECA'
        call dismoi('NOM_LIGREL', model_hho, 'MODELE', repk=ligrel_model)
!
! --- Init fields
!
        call inical(nbin, lpain, lchin, nbout, lpaout, &
                    lchout)
!
! ---- Geometry field
!
        call megeom(model_hho, chgeom)
!
! ---- Input fields
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom(1:19)
        lpain(2) = 'PDEPLPR'
        lchin(2) = disp(1:19)
        lpain(3) = 'PCHHOBS'
        lchin(3) = model_hho(1:8)//'.HHO.BASE'
!
! ---- Output fields
!
        lpaout(1) = 'PDEPL_R'
        lchout(1) = field_elno
!
! --- Compute
!
        call calcul('S', option, ligrel_model, nbin, lchin, &
                    lpain, nbout, lchout, lpaout, base, &
                    'OUI')
!
! --- Convert to disp_hho_depl
!
        call chpchd(field_elno, 'NOEU', celmod, 'OUI', 'V', &
                    field_noeu)
        call detrsd('CHAMP_GD', disp_hho_depl(1:19))
        call copisd('CHAMP_GD', 'G', field_noeu, disp_hho_depl(1:19))
        call detrsd('CHAMP_GD', field_noeu)
!
    end subroutine
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostDeplMeca(model_hho, result_hho, nume_store)
!
        implicit none
!
        character(len=24), intent(in) :: model_hho
        character(len=8), intent(in) :: result_hho
        integer(kind=8), intent(in) :: nume_store
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Post-treatment (mechanic)
!
! Init HHO field at nodes - ELNO
!
! --------------------------------------------------------------------------------------------------
!
! In  model_hho        : name of hho model
! In  result_hho       : name of hho sd_resultat
! In  nume_store       : index of storage
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ifm, niv, iret
        character(len=24) :: cham_hho, cham
!
! --------------------------------------------------------------------------------------------------
!
        call infniv(ifm, niv)
        if (niv .ge. 2) then
            call utmess('I', 'HHO2_10')
        end if
!
        call jemarq()
!
! ----- Get output fields
!
        call rsexch(' ', result_hho, 'DEPL', nume_store, cham, &
                    iret)
        ASSERT(iret == 0)
        call rsexch(' ', result_hho, 'HHO_DEPL', nume_store, cham_hho, &
                    iret)
        ASSERT(iret == 100)
!
! -----  Compute HHO field at nodes (saved in DEPL)
!
        call hhoPostDeplMecaOP(model_hho, cham, cham_hho)
        call rsnoch(result_hho, 'HHO_DEPL', nume_store)
!
! -----  Compute HHO field at nodes (saved in VITE) - if exists
!
        call rsexch(' ', result_hho, 'VITE', nume_store, cham, &
                    iret)
        if (iret == 0) then
            call rsexch(' ', result_hho, 'HHO_VITE', nume_store, cham_hho, &
                        iret)
            ASSERT(iret == 100)
            call hhoPostDeplMecaOP(model_hho, cham, cham_hho)
            call rsnoch(result_hho, 'HHO_VITE', nume_store)
        end if
!
! -----  Compute HHO field at nodes (saved in ACCE) - if exists
!
        call rsexch(' ', result_hho, 'ACCE', nume_store, cham, &
                    iret)
        if (iret == 0) then
            call rsexch(' ', result_hho, 'HHO_ACCE', nume_store, cham_hho, &
                        iret)
            ASSERT(iret == 100)
            call hhoPostDeplMecaOP(model_hho, cham, cham_hho)
            call rsnoch(result_hho, 'HHO_ACCE', nume_store)
        end if
!
        call jedema()
!
    end subroutine
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostTher(hhoCell, hhoData, nbnodes)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        integer(kind=8), intent(in) :: nbnodes
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Evaluate HHO cell unknowns at the nodes
!   In hhoCell         : a HHO Cell
!   In hhoData         : information on HHO methods
!   In nbnodes         : number of nodes
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), parameter :: max_comp = 27
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: total_dofs, cbs, fbs, ino
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: sol_T
        real(kind=8), dimension(max_comp) :: post_sol
!
        sol_T = 0.d0
        post_sol = 0.d0
!
! --- number of dofs
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        ASSERT(nbnodes .le. max_comp)
!
! --- We get the solution on the cell
!
        call readVector('PTMPCHF', cbs, sol_T, total_dofs-cbs)
!
! --- init cell basis
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- Compute the solution in the cell nodes
!
        do ino = 1, nbnodes
            post_sol(ino) = hhoEvalScalCell( &
                            hhoBasisCell, hhoData%cell_degree(), hhoCell%coorno(1:3, ino), sol_T, &
                            cbs &
                            )
        end do
!
! --- Copy of post_sol in PTEMP_R ('OUT' to fill)
!
        call writeVector("PTEMP_R", nbnodes, post_sol)
!
    end subroutine
!
! ==================================================================================================
! ==================================================================================================
!
    subroutine hhoPostTherElga(hhoCell, hhoData, hhoQuad)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_Quadrature), intent(in) :: hhoQuad
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Evaluate HHO cell unknowns at quadrature points
!   In hhoCell         : a HHO Cell
!   In hhoData         : information on HHO methods
!   In nbnodes         : number of nodes
! --------------------------------------------------------------------------------------------------
!
        type(HHO_basis_cell) :: hhoBasisCell
        integer(kind=8) :: total_dofs, cbs, fbs, ipg
        real(kind=8), dimension(MSIZE_CELL_SCAL) :: sol_T
        real(kind=8), dimension(MAX_QP_CELL) :: post_sol
!
        sol_T = 0.d0
        post_sol = 0.d0
!
! --- number of dofs
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
! --- We get the solution on the cell
!
        call readVector('PTEMPER', cbs, sol_T, total_dofs-cbs)
!
! --- init cell basis
!
        call hhoBasisCell%initialize(hhoCell)
!
! --- Compute the solution in the cell nodes
!
        do ipg = 1, hhoQuad%nbQuadPoints
            post_sol(ipg) = hhoEvalScalCell( &
                            hhoBasisCell, hhoData%cell_degree(), hhoQuad%points(1:3, ipg), sol_T, &
                            cbs &
                            )
        end do
!
! --- Copy of post_sol in PTEMP_R ('OUT' to fill)
!
        call writeVector("PTEMP_R", hhoQuad%nbQuadPoints, post_sol)
!
    end subroutine
!
end module
