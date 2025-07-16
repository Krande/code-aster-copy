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
subroutine dilini(option, ivf, ivf2, idfde, &
                  idfde2, jgano, ndim, ipoids, ipoid2, &
                  npi, dimdef, nddls, nddlm, lgpg, &
                  dimcon, typmod, dimuel, nnom, nnos, ds_dil)
!
    use dil_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"

#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/utmess.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"

!
    integer(kind=8) :: ivf, ivf2, idfde, idfde2, jgano, ndim, ipoids, npi
    integer(kind=8) :: ipoid2, dimdef, dimuel, dimcon, nnom, nnos, nddls, nddlm
    integer(kind=8) :: lgpg, iret, itab(7)
    character(len=8) :: typmod(2)
    character(len=16) :: option
    type(dil_modelisation) :: ds_dil
!
! ======================================================================
! --- BUT : INITIALISATION DES GRANDEURS NECESSAIRES POUR LA GESTION ---
! ---       DU CALCUL AVEC REGULARISATION A PARTIR DU MODELE SECOND ----
! ---       GRADIENT A MICRO-DILATATION --------------------------------
! ======================================================================
! =====================================================================
! --- VARIABLES LOCALES ------------------------------------------------
! ======================================================================
    integer(kind=8) :: nno, def1, def2, def3, cont1, cont2, cont3
    character(len=8) :: elrefe, elrf1, elrf2
! ======================================================================
! --- INITIALISATION ---------------------------------------------------
! ======================================================================
    typmod(2) = '        '
    elrf1 = '        '
    elrf2 = '        '
    lgpg = 0
! ======================================================================
! --- TYPE D'ELEMENT ---------------------------------------------------
! ======================================================================
    call elref1(elrefe)
    if (elrefe .eq. 'TR6') then
        elrf1 = 'TR6'
        elrf2 = 'TR3'
    else if (elrefe .eq. 'QU8') then
        elrf1 = 'QU8'
        elrf2 = 'QU4'
    else if (elrefe .eq. 'T10') then
        elrf1 = 'T10'
        elrf2 = 'TE4'
    else if (elrefe .eq. 'P15') then
        elrf1 = 'P15'
        elrf2 = 'PE6'
    else if (elrefe .eq. 'H20') then
        elrf1 = 'H20'
        elrf2 = 'HE8'
    else
        call utmess('F', 'DVP_4', sk=elrefe)
    end if
! ======================================================================
! --- FONCTIONS DE FORME P2 --------------------------------------------
! ======================================================================
    call elrefe_info(elrefe=elrf1, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npi, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
! ======================================================================
! --- FONCTIONS DE FORME P1 --------------------------------------------
! ======================================================================
    call elrefe_info(elrefe=elrf2, fami='RIGI', &
                     jpoids=ipoid2, jvf=ivf2, jdfde=idfde2)
! ======================================================================
! --- RECUPERATION DU TYPE DE LA MODELISATION --------------------------
! ======================================================================
    if (lteatt('DIM_TOPO_MODELI', '2')) then
        typmod(1) = 'D_PLAN'
    else if (lteatt('DIM_TOPO_MODELI', '3')) then
        typmod(1) = '3D'
    end if
! ======================================================================
! --- RECUPERATION DE LA FORMULATION
! ======================================================================
    if (lteatt('FORMULATION', 'DIL')) then
        ds_dil%inco = ASTER_FALSE
    else if (lteatt('FORMULATION', 'DIL_INCO')) then
        ds_dil%inco = ASTER_TRUE
    end if
! ======================================================================
! --- RECUPERATION DU NOMBRE DE VARIABLES INTERNES
! ======================================================================
    if (option(1:9) .ne. 'FORC_NODA') then
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=itab)
        lgpg = max(itab(6), 1)*itab(7)
    end if
!
! ======================================================================
! --- CALCUL DES DIMENSIONS
! ======================================================================
    def1 = 1+2*ndim
    def2 = ndim+1
    def3 = 1
    dimdef = def1+def2+def3
    cont1 = 1+2*ndim
    cont2 = ndim+1
    cont3 = 1
    dimcon = cont1+cont2+cont3
! ======================================================================
    nddls = ndim+2
    nddlm = ndim
! ======================================================================
    nnom = nno-nnos
    dimuel = nnos*nddls+nnom*nddlm
end subroutine
