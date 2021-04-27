! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine op0009()
!
implicit none
!
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cmeGetParameters.h"
#include "asterfort/cmePost.h"
#include "asterfort/cmePrep.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/meamac.h"
#include "asterfort/meamgy.h"
#include "asterfort/meamme.h"
#include "asterfort/medith.h"
#include "asterfort/meimme.h"
#include "asterfort/memaac.h"
#include "asterfort/memame.h"
#include "asterfort/memsth.h"
#include "asterfort/meonme.h"
#include "asterfort/mergth.h"
#include "asterfort/meriac.h"
#include "asterfort/merifs.h"
#include "asterfort/merige.h"
#include "asterfort/merigy.h"
#include "asterfort/merime.h"
#include "asterfort/meriro.h"
#include "asterfort/ntdoch.h"
!
! --------------------------------------------------------------------------------------------------
!
!                       COMMANDE:  CALC_MATR_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: base = 'G'
    character(len=19), parameter :: listLoad = '&&OP0009.LISCHA'
! - Linear case: multi-behaviour for PMF
    character(len=24) :: comporMult
    integer :: nbLoad, modeFourier
    character(len=8) :: model, caraElem, sigm, strx, disp
    character(len=16) :: k8dummy, option
    character(len=19) :: matrElem, rigiMatrElem, massMatrElem
    character(len=24) :: chtime, mate, matr_elem24, mateco
    character(len=8), pointer :: listLoadK8(:) => null()
    character(len=24), pointer :: listLoadK24(:) => null()
    real(kind=8) :: timeCurr, timeIncr
    character(len=24) :: calcElemModel, listElemCalc
    aster_logical, parameter :: hasExteStatVari = ASTER_TRUE
    aster_logical :: onlyDirichlet
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Get results
    call getres(matrElem, k8dummy, k8dummy)

! - Get parameters
    call cmeGetParameters(option       ,&
                          model        , caraElem    ,&
                          mate         , mateco      , comporMult,&
                          listLoadK8   , nbLoad      ,&
                          rigiMatrElem , massMatrElem,&
                          timeCurr     , timeIncr    , modeFourier,&
                          sigm         , strx        , disp       ,&
                          calcElemModel)

! - Preparation
    call cmePrep(option       , model        ,&
                 timeCurr     , timeIncr     , chtime     ,&
                 nbLoad       , listLoadK8   , listLoadK24,&
                 calcElemModel, onlyDirichlet,&
                 listElemCalc)

! - Compute
    if (option .eq. 'RIGI_MECA') then
        call merime(model       , nbLoad         , listLoadK24, mate       , mateco, caraElem,&
                    timeCurr    , comporMult     , matrElem   , modeFourier, base,&
                    listElemCalc, hasExteStatVari, onlyDirichlet)

    else if (option.eq.'RIGI_FLUI_STRU') then
        call merifs(model, nbLoad, listLoadK8, mate, mateco, caraElem,&
                    timeCurr, matrElem, modeFourier)

    else if (option.eq.'RIGI_GEOM') then
        call merige(model, caraElem, sigm      , strx      , matrElem,&
                    base , modeFourier , deplr=disp, mateco=mateco)

    else if (option.eq.'RIGI_ROTA') then
        call meriro(model, caraElem  , nbLoad  , listLoadK8, mate, mateco, &
                    timeCurr , comporMult, matrElem)

    else if (option.eq.'MECA_GYRO') then
        call meamgy(model  , mate, mateco      , caraElem, comporMult, matrElem,&
                    nbLoad, listLoadK8)

    else if (option.eq.'RIGI_GYRO') then
        call merigy(model ,mate,  mateco      , caraElem, comporMult, matrElem,&
                    nbLoad, listLoadK8)

    else if (option .eq. 'MASS_MECA') then
        call memame(option, model   , mate, mateco      , caraElem, timeCurr,&
                    comporMult, matrElem, base, listElemCalc)

    else if (option .eq. 'MASS_FLUI_STRU') then
        call memame(option, model   , mate, mateco      , caraElem, timeCurr,&
                    comporMult, matrElem, base, listElemCalc)

   else if (option .eq. 'MASS_MECA_DIAG') then
        call memame(option    , model   , mate, mateco      , caraElem, timeCurr,&
                    comporMult, matrElem, base, listElemCalc)

    else if (option .eq. 'AMOR_MECA') then
        call meamme(option,&
                    model       , nbLoad      , listLoadK24,&
                    mate        , mateco      , caraElem,&
                    timeCurr    , base        ,&
                    rigiMatrElem, massMatrElem,&
                    matrElem    , listElemCalc)

    else if (option.eq.'IMPE_MECA') then
        call meimme(model, nbLoad, listLoadK8, mate, mateco, matrElem)

    else if (option.eq.'ONDE_FLUI') then
        call meonme(model, nbLoad, listLoadK8, mate, mateco, matrElem)

    else if (option .eq. 'RIGI_MECA_HYST') then
        call meamme(option,&
                    model       , nbLoad      , listLoadK24,&
                    mate        , mateco      , caraElem,&
                    timeCurr    , base        ,&
                    rigiMatrElem, massMatrElem,&
                    matrElem    , listElemCalc)

    else if (option.eq.'RIGI_THER') then
        call ntdoch(listLoad, l_load_user_ = .true._1)
        matr_elem24 = matrElem
        call mergth(model      , listLoad, caraElem, mate, mateco, chtime,&
                    matr_elem24, base,&
                    timeCurr  , nh_ = modeFourier)
        call medith(base, 'CUMU', model, listLoad, matr_elem24)

    else if (option.eq.'MASS_THER') then
        call memsth(model, caraElem, mate, mateco, chtime, matrElem, base, time_curr_ = timeCurr)

    else if (option.eq.'RIGI_ACOU') then
        call meriac(model, nbLoad, listLoadK8, mate, mateco, matrElem, base)

    else if (option.eq.'MASS_ACOU') then
        call memaac(model, mate, mateco, matrElem)

    else if (option.eq.'AMOR_ACOU') then
        call meamac(model, nbLoad, listLoadK8, mate, mateco, matrElem, base)

    else
        ASSERT(ASTER_FALSE)

    endif

! - Post-treatment
    call cmePost(matrElem)

! - Clean
    AS_DEALLOCATE(vk8 = listLoadK8)
    AS_DEALLOCATE(vk24 = listLoadK24)
    call jedetr('&MERITH1           .RELR')
    call jedetr('&MERITH2           .RELR')
    call jedetr('&MERITH3           .RELR')
    call jedetr('&MERITH1           .RERR')
    call jedetr('&MERITH2           .RERR')
    call jedetr('&MERITH3           .RERR')
!
    call jedema()
end subroutine
