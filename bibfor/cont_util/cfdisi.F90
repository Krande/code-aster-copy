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
function cfdisi(sdcont_defi_, question_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "Contact_type.h"
!
    integer(kind=8) :: cfdisi
    character(len=*), intent(in) :: sdcont_defi_
    character(len=*), intent(in) :: question_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Get parameter (integer)
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  question         : question to select parameter
! Out cfdisi           : value for selected parameter
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: sdcont_defi, question
    character(len=24) :: sdcont_paraci
    integer(kind=8), pointer :: v_sdcont_paraci(:) => null()
    character(len=24) :: sdcont_ndimco
    integer(kind=8), pointer :: v_sdcont_ndimco(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    sdcont_defi = sdcont_defi_
    question = question_
    cfdisi = 0
!
! - Access to contact datastructure
!
    sdcont_ndimco = sdcont_defi(1:16)//'.NDIMCO'
    sdcont_paraci = sdcont_defi(1:8)//'.PARACI'
!
! - Get parameter
!
    if (question .eq. 'NB_ITER_GEOM') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(1)
    else if (question .eq. 'STOP_SINGULIER') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(2)
    else if (question .eq. 'NB_RESOL') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(3)
    else if (question .eq. 'FORMULATION') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(4)
    else if (question .eq. 'ITER_CONT_MULT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(5)
    else if (question .eq. 'ITER_GEOM_MAXI') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(6)
    else if (question .eq. 'ITER_FROT_MAXI') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(7)
    else if (question .eq. 'ALL_VERIF') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(8)
    else if (question .eq. 'ITER_CONT_MAXI') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(10)
    else if (question .eq. 'ALL_INTERPENETRE') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(11)
    else if (question .eq. 'ITER_GCP_MAXI') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(12)
    else if (question .eq. 'PRE_COND') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(13)
    else if (question .eq. 'ITER_PRE_MAXI') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(14)
    else if (question .eq. 'RECH_LINEAIRE') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(15)
    else if (question .eq. 'AXISYMETRIQUE') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(16)
    else if (question .eq. 'ALGO_CONT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(17)
    else if (question .eq. 'ALGO_UNIL') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(30)
    else if (question .eq. 'ALGO_FROT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(18)
    else if (question .eq. 'LISSAGE') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(19)
    else if (question .eq. 'TYPE_ADAPT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(20)
    else if (question .eq. 'EXIS_PENA') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(22)
    else if (question .eq. 'EXIS_VERIF') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(23)
    else if (question .eq. 'ALL_INTEG_NOEUD') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(24)
    else if (question .eq. 'STOP_INTERP') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(25)
    else if (question .eq. 'EXIS_GLISSIERE') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(26)
    else if (question .eq. 'ALGO_RESO_CONT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(27)
    else if (question .eq. 'ALGO_RESO_FROT') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(28)
    else if (question .eq. 'ALGO_RESO_GEOM') then
        call jeveuo(sdcont_paraci, 'L', vi=v_sdcont_paraci)
        cfdisi = v_sdcont_paraci(9)
    else if (question .eq. 'PROJ_NEWT_ITER') then
        cfdisi = 200
    else if (question .eq. 'FLIP_FLOP_IMAX') then
        cfdisi = 20
    else if (question .eq. 'NDIM') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(1)
    else if (question .eq. 'NZOCO') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(2)
    else if (question .eq. 'NSUCO') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(3)
    else if (question .eq. 'NMACO') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(4)
    else if (question .eq. 'NNOCO') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(5)
    else if (question .eq. 'NTNOE') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(8)
    else if (question .eq. 'NTMAE') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(9)
    else if (question .eq. 'NTNOM') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(10)
    else if (question .eq. 'NTMAM') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(11)
    else if (question .eq. 'NTNOEC') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(12)
    else if (question .eq. 'NTMAEC') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(13)
    else if (question .eq. 'NTNOMC') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(14)
    else if (question .eq. 'NTMAMC') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(15)
    else if (question .eq. 'NTPT') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(16)
    else if (question .eq. 'NTPC') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(17)
    else if (question .eq. 'NTMANO') then
        call jeveuo(sdcont_ndimco, 'L', vi=v_sdcont_ndimco)
        cfdisi = v_sdcont_ndimco(18)
    else
        write (6, *) 'QUESTION: ', question
        ASSERT(ASTER_FALSE)
    end if
!
end function
