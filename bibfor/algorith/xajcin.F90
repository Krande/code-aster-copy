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

subroutine xajcin(model, option, mxchin, lchin, lpain, &
                  nchin)
!
    implicit none
!
#include "asterfort/assert.h"
!
!
    integer(kind=8), intent(in) :: mxchin
    character(len=*), intent(in) :: model
    character(len=*), intent(in) :: option
    character(len=*), intent(inout) :: lpain(mxchin)
    character(len=*), intent(inout) :: lchin(mxchin)
    integer(kind=8), intent(inout) :: nchin
!
! --------------------------------------------------------------------------------------------------
!
! Add XFEM fields for input fields
!
!  -> OPTIONS : - CHAR_MECA_TEMP_R
!               - CHAR_THER_PARO_F
!               - CHAR_THER_PARO_R
!               - FULL_MECA
!               - RIGI_MECA_*
!               - RAPH_MECA
!               - CHAR_MECA_NEUM
!               - MASS_THER
!
! --------------------------------------------------------------------------------------------------
!
! In  model  : name of model
! In  option : option to select input fields
! In  mxchin : maximum number of input fields
! IO  lpain  : list of parameters
! IO  lchin  : list of fields
! IO  nbin   : number of input fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbadd
!
! --------------------------------------------------------------------------------------------------
!
    if ((option .eq. 'CHAR_MECA_TEMP_R') .or. &
        (option(1:9) .eq. 'FULL_MECA') .or. &
        (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'RIGI_MECA')) then
!
        nbadd = 12
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTTO'
        lchin(nchin+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+2) = 'PCNSETO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+3) = 'PHEAVTO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+4) = 'PLONCHA'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+5) = 'PBASLOR'
        lchin(nchin+5) = model(1:8)//'.BASLOC'
        lpain(nchin+6) = 'PLSN'
        lchin(nchin+6) = model(1:8)//'.LNNO'
        lpain(nchin+7) = 'PLST'
        lchin(nchin+7) = model(1:8)//'.LTNO'
        lpain(nchin+8) = 'PSTANO'
        lchin(nchin+8) = model(1:8)//'.STNO'
        lpain(nchin+9) = 'PPMILTO'
        lchin(nchin+9) = model(1:8)//'.TOPOSE.PMI'
        lpain(nchin+10) = 'PFISNO'
        lchin(nchin+10) = model(1:8)//'.FISSNO'
        lpain(nchin+11) = 'PHEA_NO'
        lchin(nchin+11) = model(1:8)//'.TOPONO.HNO'
        lpain(nchin+12) = 'PHEA_FA'
        lchin(nchin+12) = model(1:8)//'.TOPONO.HFA'
        nchin = nchin+nbadd
!
    elseif (option .eq. 'CHAR_MECA_NEUM') then
!
        nbadd = 19
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTTO'
        lchin(nchin+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+2) = 'PCNSETO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+3) = 'PHEAVTO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+4) = 'PLONCHA'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+5) = 'PLSN'
        lchin(nchin+5) = model(1:8)//'.LNNO'
        lpain(nchin+6) = 'PLST'
        lchin(nchin+6) = model(1:8)//'.LTNO'
        lpain(nchin+7) = 'PSTANO'
        lchin(nchin+7) = model(1:8)//'.STNO'
        lpain(nchin+8) = 'PPMILTO'
        lchin(nchin+8) = model(1:8)//'.TOPOSE.PMI'
        lpain(nchin+9) = 'PFISNO'
        lchin(nchin+9) = model(1:8)//'.FISSNO'
        lpain(nchin+10) = 'PPINTER'
        lchin(nchin+10) = model(1:8)//'.TOPOFAC.OE'
        lpain(nchin+11) = 'PAINTER'
        lchin(nchin+11) = model(1:8)//'.TOPOFAC.AI'
        lpain(nchin+12) = 'PCFACE'
        lchin(nchin+12) = model(1:8)//'.TOPOFAC.CF'
        lpain(nchin+13) = 'PLONGCO'
        lchin(nchin+13) = model(1:8)//'.TOPOFAC.LO'
        lpain(nchin+14) = 'PBASECO'
        lchin(nchin+14) = model(1:8)//'.TOPOFAC.BA'
        lpain(nchin+15) = 'PHEA_NO'
        lchin(nchin+15) = model(1:8)//'.TOPONO.HNO'
        lpain(nchin+16) = 'PHEA_SE'
        lchin(nchin+16) = model(1:8)//'.TOPONO.HSE'
        lpain(nchin+17) = 'PHEA_FA'
        lchin(nchin+17) = model(1:8)//'.TOPONO.HFA'
        lpain(nchin+18) = 'PBASLOR'
        lchin(nchin+18) = model(1:8)//'.BASLOC'
        lpain(nchin+19) = 'PHEAVNO'
        lchin(nchin+19) = model(1:8)//'.HEAVNO'
        nchin = nchin+nbadd
!
    elseif (option .eq. 'REFE_FORC_NODA') then
!
        nbadd = 14
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTTO'
        lchin(nchin+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+2) = 'PCNSETO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+3) = 'PHEAVTO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+4) = 'PLONCHA'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+5) = 'PLSN'
        lchin(nchin+5) = model(1:8)//'.LNNO'
        lpain(nchin+6) = 'PLST'
        lchin(nchin+6) = model(1:8)//'.LTNO'
        lpain(nchin+7) = 'PPMILTO'
        lchin(nchin+7) = model(1:8)//'.TOPOSE.PMI'
        lpain(nchin+8) = 'PPINTER'
        lchin(nchin+8) = model(1:8)//'.TOPOFAC.OE'
        lpain(nchin+9) = 'PAINTER'
        lchin(nchin+9) = model(1:8)//'.TOPOFAC.AI'
        lpain(nchin+10) = 'PCFACE'
        lchin(nchin+10) = model(1:8)//'.TOPOFAC.CF'
        lpain(nchin+11) = 'PBASECO'
        lchin(nchin+11) = model(1:8)//'.TOPOFAC.BA'
        lpain(nchin+12) = 'PHEA_NO'
        lchin(nchin+12) = model(1:8)//'.TOPONO.HNO'
        lpain(nchin+13) = 'PLONFA'
        lchin(nchin+13) = model(1:8)//'.TOPOFAC.LO'
        lpain(nchin+14) = 'PBASLOR'
        lchin(nchin+14) = model(1:8)//'.BASLOC'
        nchin = nchin+nbadd
!
    elseif (option .eq. 'FORC_NODA') then
!
        nbadd = 11
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTTO'
        lchin(nchin+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+2) = 'PCNSETO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+3) = 'PHEAVTO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+4) = 'PLONCHA'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+5) = 'PPMILTO'
        lchin(nchin+5) = model(1:8)//'.TOPOSE.PMI'
        lpain(nchin+6) = 'PBASLOR'
        lchin(nchin+6) = model(1:8)//'.BASLOC'
        lpain(nchin+7) = 'PLSN'
        lchin(nchin+7) = model(1:8)//'.LNNO'
        lpain(nchin+8) = 'PLST'
        lchin(nchin+8) = model(1:8)//'.LTNO'
        lpain(nchin+9) = 'PSTANO'
        lchin(nchin+9) = model(1:8)//'.STNO'
        lpain(nchin+10) = 'PFISNO'
        lchin(nchin+10) = model(1:8)//'.FISSNO'
        lpain(nchin+11) = 'PHEA_NO'
        lchin(nchin+11) = model(1:8)//'.TOPONO.HNO'
        nchin = nchin+nbadd
!
    elseif (option(6:14) .eq. 'THER_PARO') then
!
        nbadd = 9
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTER'
        lchin(nchin+1) = model(1:8)//'.TOPOFAC.OE'
        lpain(nchin+2) = 'PAINTER'
        lchin(nchin+2) = model(1:8)//'.TOPOFAC.AI'
        lpain(nchin+3) = 'PCFACE'
        lchin(nchin+3) = model(1:8)//'.TOPOFAC.CF'
        lpain(nchin+4) = 'PLONGCO'
        lchin(nchin+4) = model(1:8)//'.TOPOFAC.LO'
        lpain(nchin+5) = 'PLST'
        lchin(nchin+5) = model(1:8)//'.LTNO'
        lpain(nchin+6) = 'PSTANO'
        lchin(nchin+6) = model(1:8)//'.STNO'
        lpain(nchin+7) = 'PBASECO'
        lchin(nchin+7) = model(1:8)//'.TOPOFAC.BA'
        lpain(nchin+8) = 'PLSN'
        lchin(nchin+8) = model(1:8)//'.LNNO'
        lpain(nchin+9) = 'PHEA_NO'
        lchin(nchin+9) = model(1:8)//'.TOPONO.HNO'
        nchin = nchin+nbadd
!
    elseif (option .eq. 'RIGI_THER') then
!
        nbadd = 9
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PSTANO'
        lchin(nchin+1) = model(1:8)//'.STNO'
        lpain(nchin+2) = 'PPINTTO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+3) = 'PCNSETO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+4) = 'PHEAVTO'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+5) = 'PLONCHA'
        lchin(nchin+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+6) = 'PBASLOR'
        lchin(nchin+6) = model(1:8)//'.BASLOC'
        lpain(nchin+7) = 'PLSN'
        lchin(nchin+7) = model(1:8)//'.LNNO'
        lpain(nchin+8) = 'PLST'
        lchin(nchin+8) = model(1:8)//'.LTNO'
        lpain(nchin+9) = 'PHEA_NO'
        lchin(nchin+9) = model(1:8)//'.TOPONO.HNO'
        nchin = nchin+nbadd
!
    elseif (option .eq. 'CHAR_THER_EVOL') then
!
        nbadd = 9
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PSTANO'
        lchin(nchin+1) = model(1:8)//'.STNO'
        lpain(nchin+2) = 'PPINTTO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+3) = 'PCNSETO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+4) = 'PHEAVTO'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+5) = 'PLONCHA'
        lchin(nchin+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+6) = 'PBASLOR'
        lchin(nchin+6) = model(1:8)//'.BASLOC'
        lpain(nchin+7) = 'PLSN'
        lchin(nchin+7) = model(1:8)//'.LNNO'
        lpain(nchin+8) = 'PLST'
        lchin(nchin+8) = model(1:8)//'.LTNO'
        lpain(nchin+9) = 'PHEA_NO'
        lchin(nchin+9) = model(1:8)//'.TOPONO.HNO'
        nchin = nchin+nbadd
    elseif (option .eq. 'MASS_THER') then
!
        nbadd = 9
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PSTANO'
        lchin(nchin+1) = model(1:8)//'.STNO'
        lpain(nchin+2) = 'PPINTTO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+3) = 'PCNSETO'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+4) = 'PHEAVTO'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+5) = 'PLONCHA'
        lchin(nchin+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+6) = 'PBASLOR'
        lchin(nchin+6) = model(1:8)//'.BASLOC'
        lpain(nchin+7) = 'PLSN'
        lchin(nchin+7) = model(1:8)//'.LNNO'
        lpain(nchin+8) = 'PLST'
        lchin(nchin+8) = model(1:8)//'.LTNO'
        lpain(nchin+9) = 'PHEA_NO'
        lchin(nchin+9) = model(1:8)//'.TOPONO.HNO'
        nchin = nchin+nbadd
    elseif (option(1:5) .eq. 'MASS_') then
!
        nbadd = 10
        ASSERT(nchin+nbadd .le. mxchin)
        lpain(nchin+1) = 'PPINTTO'
        lchin(nchin+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nchin+2) = 'PHEAVTO'
        lchin(nchin+2) = model(1:8)//'.TOPOSE.HEA'
        lpain(nchin+3) = 'PLONCHA'
        lchin(nchin+3) = model(1:8)//'.TOPOSE.LON'
        lpain(nchin+4) = 'PCNSETO'
        lchin(nchin+4) = model(1:8)//'.TOPOSE.CNS'
        lpain(nchin+5) = 'PBASLOR'
        lchin(nchin+5) = model(1:8)//'.BASLOC'
        lpain(nchin+6) = 'PLSN'
        lchin(nchin+6) = model(1:8)//'.LNNO'
        lpain(nchin+7) = 'PLST'
        lchin(nchin+7) = model(1:8)//'.LTNO'
        lpain(nchin+8) = 'PSTANO'
        lchin(nchin+8) = model(1:8)//'.STNO'
        lpain(nchin+9) = 'PHEA_NO'
        lchin(nchin+9) = model(1:8)//'.TOPONO.HNO'
        lpain(nchin+10) = 'PPMILTO'
        lchin(nchin+10) = model(1:8)//'.TOPOSE.PMI'
        nchin = nchin+nbadd
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
