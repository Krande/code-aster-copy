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
subroutine cmePrep(optionz      , modelz       ,&
                   timeCurr     , timeIncr     , chtime     ,&
                   nbLoad       , listLoadK8   , listLoadK24,&
                   calcElemModel, onlyDirichlet,&
                   matrElemz    , listElemCalc)
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/mecact.h"
#include "asterfort/getelem.h"
#include "asterfort/jeveuo.h"
#include "asterfort/exlim1.h"
#include "asterfort/jedetr.h"
!
character(len=*), intent(in) :: optionz, modelz
real(kind=8), intent(in) :: timeCurr, timeIncr
character(len=24), intent(out) :: chtime
integer, intent(in) :: nbLoad
character(len=8), pointer :: listLoadK8(:)
character(len=24), pointer :: listLoadK24(:)
character(len=8), intent(in) :: calcElemModel
aster_logical, intent(out) :: onlyDirichlet
character(len=*), intent(in) :: matrElemz
character(len=24), intent(out) :: listElemCalc
!
! --------------------------------------------------------------------------------------------------
!
! CALC_MATR_ELEM
!
! Preparation
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  model            : name of the model
! In  timeCurr         : current time
! In  timeIncr         : time step
! Out chtime           : time parameters (field)
! In  nbLoad           : number of loads
! Ptr listLoadK8       : pointer to list of loads (K8)
! Ptr listLoadK24      : pointer to list of loads (K24)
! In  calcElemModel    : value of keyword CALC_ELEM_MODELE
! Out onlyDirichlet    : flag to compute only Dirichlet [B] matrix
! Out listElemCalc     : list of elements (LIGREL) where matrElem is computed
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywfact = ' '
    integer, parameter :: iocc = 0
    character(len=24), parameter :: listCellJv = '&&CAIMPE.LIST_CELL'
    integer :: nbCell
    integer, pointer :: listCell(:) => null()
    character(len=8) :: model, mesh, matrElem
    integer, parameter :: nbCmp = 6
    character(len=8), parameter :: cmpName(nbCmp) = (/'INST    ','DELTAT  ','THETA   ',&
                                                      'KHI     ','R       ','RHO     '/)
    real(kind=8) :: cmpVale(nbCmp)
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
    option   = optionz
    model    = modelz
    matrElem = matrElemz

! - Preapre field for time
    chtime       = '&&CHTIME'
    cmpVale(1:6) = [timeCurr, timeIncr, 1.d0, 0.d0, 0.d0, 0.d0]
    if ((option.eq.'RIGI_THER') .or. (option.eq.'MASS_THER')) then
        call mecact('V', chtime, 'MODELE', modelz, 'INST_R',&
                    ncmp=nbCmp, lnomcmp=cmpName, vr=cmpVale)
    endif

! - Conversion for names of loads
    listLoadK24 => null()
    if (nbLoad .ne. 0) then
        AS_ALLOCATE(vk24 = listLoadK24, size = nbLoad)
        listLoadK24(1:nbLoad) = listLoadK8(1:nbLoad)
    endif

! - Detect flag to compute only Dirichlet [B] matrix
    onlyDirichlet = ASTER_FALSE
    if (option .eq. 'RIGI_MECA') then
        onlyDirichlet = calcElemModel .eq. 'NON'
    endif

! - Create list of elements (LIGREL) where matrElem is computed

    call dismoi('NOM_MAILLA', model, 'MODELE', repk = mesh)
    if (option .eq. 'RIGI_MECA' .or. option .eq. 'MASS_MECA' .or. &
        option .eq. 'MASS_MECA_DIAG' .or. option .eq. 'AMOR_MECA' .or.&
        option .eq. 'RIGI_MECA_HYST' .or. option .eq. 'MASS_FLUI_STRU') then
        call getelem(mesh, keywFact, iocc, ' ', listCellJv, nbCell)
        if (nbCell .eq. 0) then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk = listElemCalc)
        else
            listElemCalc = matrElem(1:8)//'LIGREL'
            call jeveuo(listCellJv, 'L', vi = listCell)
            call exlim1(listCell, nbCell, model, 'G', listElemCalc)
        endif
        call jedetr(listCellJv)
    else
        call dismoi('NOM_LIGREL', model, 'MODELE', repk = listElemCalc)
    endif
!
end subroutine
