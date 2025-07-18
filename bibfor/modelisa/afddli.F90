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
! aslint: disable=W1504
!
subroutine afddli(model, geomDime, gran_cmp_nb, gran_cmp_name, node_nume, node_name, &
                  prnm, repe_type, repe_defi, coef_type, cmp_nb, &
                  cmp_name, cmp_acti, vale_type, vale_real, vale_func, &
                  vale_cplx, cmp_count, list_rela, lxfem, jnoxfl, &
                  jnoxfv, ch_xfem_stat, ch_xfem_lnno, ch_xfem_ltno, connex_inv, &
                  mesh, ch_xfem_heav)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/xddlim.h"
!
    character(len=8), intent(in) :: model
    integer(kind=8), intent(in) :: geomDime
    integer(kind=8), intent(in) :: gran_cmp_nb
    character(len=8), intent(in) :: gran_cmp_name(gran_cmp_nb)
    integer(kind=8), intent(in) :: node_nume
    character(len=8), intent(in) :: node_name
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: prnm(*)
    integer(kind=8), intent(in) :: repe_type
    real(kind=8), intent(in) :: repe_defi(3)
    character(len=4), intent(in) :: coef_type
    integer(kind=8), intent(in) :: cmp_nb
    character(len=16), intent(in) :: cmp_name(cmp_nb)
    integer(kind=8), intent(in) :: cmp_acti(cmp_nb)
    character(len=4), intent(in) :: vale_type
    real(kind=8), intent(in) :: vale_real(cmp_nb)
    character(len=8), intent(in) :: vale_func(cmp_nb)
    complex(kind=8), intent(in) :: vale_cplx(cmp_nb)
    integer(kind=8), intent(inout) :: cmp_count(cmp_nb)
    character(len=19), intent(in) :: list_rela
    aster_logical, intent(in) :: lxfem
    integer(kind=8), intent(in) :: jnoxfl
    integer(kind=8), intent(in) :: jnoxfv
    character(len=19), intent(in) :: connex_inv
    character(len=19), intent(in) :: ch_xfem_stat
    character(len=19), intent(in) :: ch_xfem_lnno
    character(len=19), intent(in) :: ch_xfem_ltno
    character(len=19), intent(in) :: ch_xfem_heav
!
! --------------------------------------------------------------------------------------------------
!
! Loadings - Kinematic
!
! Apply simple kinematic relation
!
! --------------------------------------------------------------------------------------------------
!
! For i=1,cmp_nb
!    coef * component_i = vale_i on node (numnoe,nomnoe) with coef = 1 (real or cplx)
!
! Overload rule: last kinematic relation on the node kept
!
! In  model          : Name of model
! In  gran_cmp_nb    : number of component of <GRANDEUR> (as DEPL_R, TEMP_R, etc)
! In  gran_cmp_name  : names of component of <GRANDEUR> (as DEPL_R, TEMP_R, etc)
! In  node_nume      : number (in mesh) of the node
! In  node_name      : name of the node
! In  prnm(*)        : <GRANDEUR> on node
! In  repe_type      : If 0 -> global reference system
!                      If 2/3 -> local reference system give by repe_defi
! In  repe_defi      : local reference system
! In  coef_type      : type of coefficient (real or complex)
! In  cmp_nb         : number of components
! In  cmp_name       : name of components
! In  cmp_acti       : 1 if component affected, 0 else
! In  vale_type      : affected value type
! In  vale_real      : affected value if real
! In  vale_func      : affected value if function
! In  vale_cplx      : affected value if complex
! In  cmp_count      : how many times components have been affected
! In  list_rela      : JEVEUX object liste_rela for aflrch.F90 subroutine
! In  l_xfem         : .true. if xfem
! In  connex_inv     : inverse connectivity (blank if not xfem)
! In  jnoxfl         : pointer on list of XFEM nodes
! In  jnoxfv         : pointer on list of XFEM nodes
! In  ch_xfem_stat   : status of nodes field (blank if not xfem)
! In  ch_xfem_lnno   : normal level-set field (blank if not xfem)
! In  ch_xfem_ltno   : tangent level-set field (blank if not xfem)
! In  mesh           : name of the mesh
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_cmp, cmp_index
    real(kind=8) :: rbid(3)
    real(kind=8), parameter :: coef_real_unit = 1.d0
    complex(kind=8), parameter :: coef_cplx_unit = dcmplx(1.d0, 0.d0)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Loop on components
    do i_cmp = 1, cmp_nb
!
! ----- Is component exists on this node ?
!
        cmp_index = indik8(gran_cmp_name, cmp_name(i_cmp) (1:8), 1, gran_cmp_nb)
        ASSERT(cmp_index .gt. 0)
        if (.not. exisdg(prnm, cmp_index)) cycle
!
! ----- Apply on component, XFEM case
!
        if (lxfem) then
            ASSERT(coef_type .eq. 'REEL')
            if (zl(jnoxfl-1+2*node_nume)) then
!           ACTUELLEMENT EN XFEM ON NE PEUT IMPOSER DES RELATIONS
!           SUR LES DDLS QU'EN MECA ET EN HM
                if (cmp_name(i_cmp) (1:1) .eq. 'D' .or. cmp_name(i_cmp) .eq. 'PRE1') then
                    call xddlim(model, cmp_name(i_cmp) (1:8), node_name, node_nume, &
                                vale_real(i_cmp), vale_cplx(i_cmp), vale_func(i_cmp), &
                                vale_type, cmp_count(i_cmp), list_rela, geomDime, rbid, &
                                jnoxfv, ch_xfem_stat, ch_xfem_lnno, ch_xfem_ltno, &
                                connex_inv, mesh, ch_xfem_heav)
                    cycle
                end if
            end if
        end if

! ----- Count
        cmp_count(i_cmp) = cmp_count(i_cmp)+1

! ----- Apply on active component
        ASSERT(cmp_acti(i_cmp) .le. 1)
        ASSERT(cmp_acti(i_cmp) .ge. 0)
        if (cmp_acti(i_cmp) .eq. 1) then
            call afrela([coef_real_unit], [coef_cplx_unit], cmp_name(i_cmp) (1:8), node_name, &
                        [repe_type], repe_defi, 1, vale_real(i_cmp), vale_cplx(i_cmp), &
                        vale_func(i_cmp), coef_type, vale_type, 0.d0, &
                        list_rela)
        end if
!
    end do
!
    call jedema()
end subroutine
