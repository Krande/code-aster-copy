! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine thmCompForcNoda(ds_thm)
!
use THM_type
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/thmGetElemDime.h"
#include "asterfort/fnothm.h"
#include "asterfort/thmGetGeneDime.h"
#include "asterfort/thmGetElemInfo.h"
#include "asterfort/thmGetElemRefe.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/thmGetGene.h"
#include "asterfort/thmGetElemIntegration.h"
!
type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! Nodal force (FORC_NODA)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefe, elref2
    integer :: jv_geom, jv_mater, jv_sigm, jv_vect, jv_instm, jv_instp
    integer :: iretm, iretp
    aster_logical :: fnoevo
    integer :: nno, nnos, nnom
    integer :: npi, npi2, npg
    real(kind=8) :: dt
    integer :: dimdep, dimdef, dimcon, dimuel
    integer :: nddls, nddlm
    integer :: nddl_meca, nddl_p1, nddl_p2
    real(kind=8) :: b(21, 120), r(22)
    integer :: jv_poids, jv_poids2
    integer :: jv_func, jv_func2, jv_dfunc, jv_dfunc2, jv_gano
    aster_logical :: l_axi, l_vf, l_steady
    character(len=3) :: inte_type
    integer :: ndim
    integer :: mecani(5), press1(7), press2(7), tempe(5)
!
! --------------------------------------------------------------------------------------------------
!
    fnoevo = .false.
    dt     = 0.d0
!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm, l_axi, l_vf, l_steady, ndim)
!
! - Cannot compute for finite volume
!
    ASSERT(.not.l_vf)
!
! - Get type of integration
!
    call thmGetElemIntegration(l_vf, inte_type)
!
! - Get generalized coordinates
!
    call thmGetGene(ds_thm, l_steady, l_vf  , ndim ,&
                    mecani, press1  , press2, tempe)
!
! - Is transient computation (STAT_NON_LINE or CALC_CHAMP ? )
!
    call tecach('ONO', 'PINSTMR', 'L', iretm, iad=jv_instm)
    call tecach('ONO', 'PINSTPR', 'L', iretp, iad=jv_instp)
    if (iretm .eq. 0 .and. iretp .eq. 0) then
        dt     = zr(jv_instp) - zr(jv_instm)
        fnoevo = .true.
    else
        dt     = 0.d0
        fnoevo = .false.
    endif
!
! - Input/ouput fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mater)
    call jevech('PCONTMR', 'L', jv_sigm)
    call jevech('PVECTUR', 'E', jv_vect)
!
! - Get reference elements
!
    call thmGetElemRefe(l_vf, elrefe, elref2)
!
! - Get informations about element
!
    call thmGetElemInfo(l_vf, elrefe, elref2,&
                        nno, nnos, nnom, &
                        jv_gano, jv_poids, jv_poids2,&
                        jv_func, jv_func2, jv_dfunc, jv_dfunc2,&
                        inte_type, npi   , npi2    , npg)
    ASSERT(npi .le. 27)
    ASSERT(nno .le. 20)
!
! - Get dimensions of generalized vectors
!
    call thmGetGeneDime(ndim  ,&
                        mecani, press1, press2, tempe,&
                        dimdep, dimdef, dimcon)
!
! - Get dimensions about element
!
    call thmGetElemDime(ndim     , nnos   , nnom   ,&
                        mecani   , press1 , press2 , tempe ,&
                        nddls    , nddlm  ,&
                        nddl_meca, nddl_p1, nddl_p2,&
                        dimdep   , dimdef , dimcon , dimuel)
!
! - Compute
!
    call fnothm(ds_thm, zi(jv_mater), ndim     , l_axi    , l_steady , fnoevo ,&
                mecani      , press1   , press2   , tempe    ,&
                nno         , nnos     , npi      , npg    ,&
                zr(jv_geom) , dt       , dimdef   , dimcon   , dimuel ,&
                jv_poids    , jv_poids2,&
                jv_func     , jv_func2 , jv_dfunc , jv_dfunc2,&
                nddls       , nddlm    , nddl_meca, nddl_p1  , nddl_p2,&
                zr(jv_sigm) , b        , r        , zr(jv_vect) )
!
end subroutine
