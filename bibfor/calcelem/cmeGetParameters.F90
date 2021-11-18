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
subroutine cmeGetParameters(option       ,&
                            model        , caraElem    ,&
                            mate         , mateco      , comporMult,&
                            listLoadK8   , nbLoad      ,&
                            rigiMatrElem , massMatrElem,&
                            timeCurr     , timeIncr    , modeFourier,&
                            sigm         , strx        , disp,&
                            calcElemModel)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/chpver.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/utmess.h"
#include "asterfort/get_load8.h"
!
character(len=16), intent(out) :: option
character(len=8), intent(out) :: model, caraElem
character(len=24), intent(out) :: mate, mateco, comporMult
character(len=8), pointer :: listLoadK8(:)
integer, intent(out) :: nbLoad
character(len=19), intent(out) :: rigiMatrElem, massMatrElem
real(kind=8), intent(out) :: timeCurr, timeIncr
integer, intent(out) :: modeFourier
character(len=8), intent(out) :: sigm, strx, disp
character(len=8), intent(out) :: calcElemModel
!
! --------------------------------------------------------------------------------------------------
!
! CALC_MATR_ELEM
!
! Get parameters
!
! --------------------------------------------------------------------------------------------------
!
! Out option           : option to compute
! Out model            : name of the model
! Out caraElem         : name of elementary characteristics (field)
! Out mate             : name of material characteristics (field)
! Out comporMult       : multi-behaviour for multifibers beams (field)
! Ptr listLoadK8       : pointer to list of loads (K8)
! Out nbLoad           : number of loads
! Out rigiMatrElem     : option for mechanic rigidity (useful for damping)
! Out massMatrElem     : option for mechanic mass (useful for damping)
! Out timeCurr         : current time
! Out timeIncr         : time step
! Out modeFourier      : Fourier mode
! Out sigm             : stress
! Out strx             : fibers information
! Out disp             : displacements
! Out calcElemModel    : value of keyword CALC_ELEM_MODELE
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nocc, ier
    character(len=8) ::  answer
    aster_logical :: l_ther
!
! --------------------------------------------------------------------------------------------------
!
    option        = ' '
    rigiMatrElem  = ' '
    massMatrElem  = ' '
    timeCurr      = 0.d0
    timeIncr      = 0.d0
    modeFourier   = 0
    model         = ' '
    caraElem      = ' '
    mate          = ' '
    mateco        = ' '
    comporMult    = ' '
    listLoadK8    => null()
    nbLoad        = 0
    sigm          = ' '
    strx          = ' '
    disp          = ' '
    calcElemModel = 'OUI'

! - Get parameters
    call getvtx(' ', 'OPTION'   , scal=option, nbret=nocc)
    ASSERT(nocc.eq.1)
    call getvid(' ', 'RIGI_MECA', scal=rigiMatrElem, nbret=nocc)
    if (nocc .eq. 0) then
        rigiMatrElem = ' '
    endif
    call getvid(' ', 'MASS_MECA', scal=massMatrElem, nbret=nocc)
    if (nocc .eq. 0) then
        massMatrElem = ' '
    endif
    call getvr8(' ', 'INST', scal=timeCurr, nbret=nocc)
    if (nocc .eq. 0) then
        timeCurr = 0.d0
    endif
    call getvr8(' ', 'INCR_INST', scal=timeIncr, nbret=nocc)
    if (nocc .eq. 0) then
        timeIncr = 0.d0
    endif
    call getvis(' ', 'MODE_FOURIER', scal=modeFourier, nbret=nocc)
    if (nocc .eq. 0) then
        modeFourier = 0
    endif
    call getvid(' ', 'MODELE', scal=model, nbret=nocc)
    ASSERT(nocc.eq.1)
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=nocc)
    call dismoi('EXI_RDM', model, 'MODELE', repk=answer)
    if ((nocc.eq.0) .and. (answer.eq.'OUI')) then
        call utmess('A', 'MECHANICS1_39')
        caraElem = ' '
    endif
    call getvid(' ', 'CHAM_MATER', scal=mate, nbret=nocc)

    if (nocc .eq. 0) then 
        mate = ' '
        if (option.eq.'RIGI_GEOM')  then
            ! necessaire seulement pour CABLE 
            call dismoi('EXI_CABLE', model, 'MODELE', repk=answer)
            if (answer .eq. 'OUI') then
                call utmess('A', 'MECHANICS1_40')
            endif
        elseif ((option .ne. 'RIGI_ACOU') .and. (option.ne.'RIGI_GEOM') ) then
            ! mater pas besoin pour RIGI_ACOU, 
            ! mater peu besoin pour DIS_/2D_DIS_ pour d'autres options
            call dismoi('BESOIN_MATER', model, 'MODELE', repk=answer)
            if (answer .eq. 'OUI') then
                call utmess('A', 'MECHANICS1_40')
            endif            
        endif
    endif

    l_ther = ASTER_FALSE
    if (option .eq. 'MASS_THER' .or. option.eq. 'RIGI_THER') then
        l_ther = ASTER_TRUE
    endif

    if (mate .ne. ' ') then
        call rcmfmc(mate, mateco, l_ther_ = l_ther)
    else
        mate = ' '
    endif
    call get_load8(model, listLoadK8, nbLoad)

! - Get multi-behaviour for multifibers beams
    comporMult = mate(1:8)//'.COMPOR'
!
    if (option.eq.'RIGI_GEOM') then
        call getvid(' ', 'SIEF_ELGA', scal=sigm, nbret=nocc)
        if (nocc .ne. 0) then
            call chpver('F', sigm, 'ELGA', 'SIEF_R', ier)
        endif
        call getvid(' ', 'STRX_ELGA', scal=strx, nbret=nocc)
        if (nocc .ne. 0) then
            call chpver('F', strx, 'ELGA', 'STRX_R', ier)
        endif
        call getvid(' ', 'DEPL', scal=disp, nbret=nocc)
        if (nocc .ne. 0) then
            call chpver('F', disp, 'NOEU', 'DEPL_R', ier)
        endif
    endif

! - Specific for RIGI_MECA
    if (option .eq. 'RIGI_MECA') then
        call getvtx(' ', 'CALC_ELEM_MODELE', scal=calcElemModel, nbret=nocc)
    endif

!
end subroutine
