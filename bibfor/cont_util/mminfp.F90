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
subroutine mminfp(i_zone, sdcont_defi, question_, vale_i_, vale_r_, &
                  vale_l_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/jeveuo.h"
!
    character(len=24), intent(in) :: sdcont_defi
    character(len=*), intent(in) :: question_
    integer(kind=8), intent(in) :: i_zone
    real(kind=8), optional, intent(out) :: vale_r_
    integer(kind=8), optional, intent(out) :: vale_i_
    aster_logical, optional, intent(out) :: vale_l_
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Utility
!
! Get parameter (all type) - By zone contact
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont_defi      : name of contact definition datastructure (from DEFI_CONTACT)
! In  question         : question to select parameter
! In  i_zone           : index of contact zone
! Out vale_i           : value for selected parameter (integer)
! Out vale_r           : value for selected parameter (real)
! Out vale_l           : value for selected parameter (boolean)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: cont_form, nb_cont_zone
    integer(kind=8) :: zcmcf, zmeth, ztole, zexcl, zdirn, zcmdf
    character(len=24) :: sdcont_caracf, sdcont_caradf, sdcont_dirnor, sdcont_methco
    real(kind=8), pointer :: v_sdcont_dirnor(:) => null()
    integer(kind=8), pointer :: v_sdcont_methco(:) => null()
    real(kind=8), pointer :: v_sdcont_caracf(:) => null()
    real(kind=8), pointer :: v_sdcont_caradf(:) => null()
    character(len=24) :: sdcont_toleco, sdcont_dirapp, sdcont_exclfr
    real(kind=8), pointer :: v_sdcont_toleco(:) => null()
    real(kind=8), pointer :: v_sdcont_dirapp(:) => null()
    real(kind=8), pointer :: v_sdcont_exclfr(:) => null()
    character(len=24) :: sdcont_jeufo1, sdcont_jeufo2
    character(len=8), pointer :: v_sdcont_jeufo1(:) => null()
    character(len=8), pointer :: v_sdcont_jeufo2(:) => null()
    character(len=8) :: jeuf1, jeuf2
    character(len=24) :: question
    real(kind=8) :: vale_r
    integer(kind=8) :: vale_i
    aster_logical :: vale_l
!
! --------------------------------------------------------------------------------------------------
!
    vale_i = 0
    vale_r = 0.d0
    vale_l = .false.
    question = question_
!
! - Access to contact datastructure
!
    sdcont_caracf = sdcont_defi(1:16)//'.CARACF'
    sdcont_caradf = sdcont_defi(1:16)//'.CARADF'
    sdcont_dirapp = sdcont_defi(1:16)//'.DIRAPP'
    sdcont_dirnor = sdcont_defi(1:16)//'.DIRNOR'
    sdcont_methco = sdcont_defi(1:16)//'.METHCO'
    sdcont_toleco = sdcont_defi(1:16)//'.TOLECO'
    sdcont_exclfr = sdcont_defi(1:16)//'.EXCLFR'
    sdcont_jeufo1 = sdcont_defi(1:16)//'.JFO1CO'
    sdcont_jeufo2 = sdcont_defi(1:16)//'.JFO2CO'
    zmeth = cfmmvd('ZMETH')
    ztole = cfmmvd('ZTOLE')
    zcmcf = cfmmvd('ZCMCF')
    zcmdf = cfmmvd('ZCMDF')
    zexcl = cfmmvd('ZEXCL')
    zdirn = cfmmvd('ZDIRN')
!
! - Get parameters
!
    nb_cont_zone = cfdisi(sdcont_defi, 'NZOCO')
    cont_form = cfdisi(sdcont_defi, 'FORMULATION')
    ASSERT(i_zone .gt. 0)
    if (nb_cont_zone .ne. 0) then
        ASSERT(i_zone .le. nb_cont_zone)
    end if
!
    if (question .eq. 'APPARIEMENT') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+1)
    else if (question .eq. 'DIST_POUTRE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_l = v_sdcont_methco(zmeth*(i_zone-1)+2) .eq. 1
    else if (question .eq. 'DIST_COQUE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_l = v_sdcont_methco(zmeth*(i_zone-1)+3) .eq. 1
    else if (question .eq. 'DIST_MAIT') then
        call jeveuo(sdcont_jeufo1, 'L', vk8=v_sdcont_jeufo1)
        jeuf1 = v_sdcont_jeufo1(i_zone)
        if (jeuf1 .eq. ' ') then
            vale_l = .false.
        else
            vale_l = .true.
        end if
    else if (question .eq. 'DIST_ESCL') then
        call jeveuo(sdcont_jeufo2, 'L', vk8=v_sdcont_jeufo2)
        jeuf2 = v_sdcont_jeufo2(i_zone)
        if (jeuf2 .eq. ' ') then
            vale_l = .false.
        else
            vale_l = .true.
        end if
    else if (question .eq. 'NORMALE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+4)
    else if (question .eq. 'MAIT') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+4) .eq. 0) then
            vale_l = .true.
        else
            vale_l = .false.
        end if
    else if (question .eq. 'MAIT_ESCL') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+4) .eq. 1) then
            vale_l = .true.
        else
            vale_l = .false.
        end if
    else if (question .eq. 'ESCL') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+4) .eq. 2) then
            vale_l = .true.
        else
            vale_l = .false.
        end if
    else if (question .eq. 'VECT_MAIT') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+5)
    else if (question .eq. 'VECT_MAIT_DIRX') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+5) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+1)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'VECT_MAIT_DIRY') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+5) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+2)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'VECT_MAIT_DIRZ') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+5) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+3)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'VECT_ESCL') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+6)
    else if (question .eq. 'VECT_ESCL_DIRX') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+6) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+4)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'VECT_ESCL_DIRY') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+6) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+5)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'VECT_ESCL_DIRZ') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+6) .gt. 0) then
            call jeveuo(sdcont_dirnor, 'L', vr=v_sdcont_dirnor)
            vale_r = v_sdcont_dirnor(zdirn*(i_zone-1)+6)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'TYPE_APPA') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+7)
    else if (question .eq. 'TYPE_APPA_FIXE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_l = v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 1
    else if (question .eq. 'TYPE_APPA_DIRX') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 1) then
            call jeveuo(sdcont_dirapp, 'L', vr=v_sdcont_dirapp)
            vale_r = v_sdcont_dirapp(3*(i_zone-1)+1)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'TYPE_APPA_DIRY') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 1) then
            call jeveuo(sdcont_dirapp, 'L', vr=v_sdcont_dirapp)
            vale_r = v_sdcont_dirapp(3*(i_zone-1)+2)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'TYPE_APPA_DIRZ') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+7) .eq. 1) then
            call jeveuo(sdcont_dirapp, 'L', vr=v_sdcont_dirapp)
            vale_r = v_sdcont_dirapp(3*(i_zone-1)+3)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'NBMAE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+8)
    else if (question .eq. 'NBNOE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+9)
    else if (question .eq. 'NBMAM') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+10)
    else if (question .eq. 'NBNOM') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+11)
    else if (question .eq. 'NBMAEC') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+12)
    else if (question .eq. 'NBNOEC') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+13)
    else if (question .eq. 'NBMAMC') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+14)
    else if (question .eq. 'NBNOMC') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+15)
    else if (question .eq. 'JDECME') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+16)
    else if (question .eq. 'JDECMM') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+17)
    else if (question .eq. 'JDECNE') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+18)
    else if (question .eq. 'JDECNM') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+19)
    else if (question .eq. 'NBPT') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+20)
    else if (question .eq. 'NBPC') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+21)
    else if (question .eq. 'VERIF') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+22) .eq. 0) then
            vale_l = .false.
        else if (v_sdcont_methco(zmeth*(i_zone-1)+22) .eq. 1) then
            vale_l = .true.
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'CALCUL') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        if (v_sdcont_methco(zmeth*(i_zone-1)+22) .eq. 0) then
            vale_l = .true.
        else if (v_sdcont_methco(zmeth*(i_zone-1)+22) .eq. 1) then
            vale_l = .false.
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'TYPE_JACOBIEN') then
        call jeveuo(sdcont_methco, 'L', vi=v_sdcont_methco)
        vale_i = v_sdcont_methco(zmeth*(i_zone-1)+23)
    else if (question .eq. 'TOLE_PROJ_EXT') then

        call jeveuo(sdcont_toleco, 'L', vr=v_sdcont_toleco)
        vale_r = v_sdcont_toleco(ztole*(i_zone-1)+1)

    else if (question .eq. 'DIST_APPA') then
        call jeveuo(sdcont_toleco, 'L', vr=v_sdcont_toleco)
        vale_r = v_sdcont_toleco(ztole*(i_zone-1)+2)
    else if (question .eq. 'RESI_APPA') then
        call jeveuo(sdcont_toleco, 'L', vr=v_sdcont_toleco)
        vale_r = v_sdcont_toleco(ztole*(i_zone-1)+2)
    else if (question .eq. 'TOLE_INTERP') then
        call jeveuo(sdcont_toleco, 'L', vr=v_sdcont_toleco)
        vale_r = v_sdcont_toleco(ztole*(i_zone-1)+3)
    else if (question .eq. 'ALGO_CONT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+3)
            vale_i = nint(vale_r)
        end if
    else if (question .eq. 'ALGO_CONT_PENA') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            if (nint(v_sdcont_caracf(zcmcf*(i_zone-1)+3)) .eq. 3) then
                vale_l = .true.
            else
                vale_l = .false.
            end if
        end if
    else if (question .eq. 'ALGO_FROT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+5)
            vale_i = nint(vale_r)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'ALGO_FROT_PENA') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            if (nint(v_sdcont_caracf(zcmcf*(i_zone-1)+5)) .eq. 3) then
                vale_l = .true.
            else
                vale_l = .false.
            end if
        end if
    else if (question .eq. 'GLISSIERE_ZONE') then
        vale_l = .false.
        if (cont_form .eq. 1) then
            call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
            vale_l = nint(v_sdcont_caradf(zcmdf*(i_zone-1)+6)) .eq. 1
        else if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_l = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+9)) .eq. 1
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'INTEGRATION') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_i = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+1))
        else if (cont_form .eq. 5) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_i = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+1))
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'COEF_COULOMB') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+6)
        else if (cont_form .eq. 1) then
            call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
            vale_r = v_sdcont_caradf(zcmdf*(i_zone-1)+4)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'COEF_AUGM_CONT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+2)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'COEF_AUGM_FROT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+4)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'FROTTEMENT_ZONE') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+5)
            vale_i = nint(vale_r)
            vale_l = vale_i .ne. 0
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'SEUIL_INIT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+7)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'SEUIL_AUTO') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_l = (nint(v_sdcont_caracf(zcmcf*(i_zone-1)+13)) .eq. 1)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'PENE_MAXI') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+14)
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'GLIS_MAXI') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+16)
            if (vale_r .ge. 1000) vale_r = -1
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'GRAND_GLIS') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = v_sdcont_caracf(zcmcf*(i_zone-1)+15)
            if (nint(vale_r) .eq. 1) then
                vale_l = .true.
            else
                vale_l = .false.
            end if
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'CONTACT_INIT') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_i = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+8))
        else if (cont_form .eq. 5) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_i = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+8))
        else
            ASSERT(.false.)
        end if
    else if (question .eq. 'CONTACT_INIT_DIST') then
        if (cont_form .eq. 2) then
            call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
            vale_r = (v_sdcont_caracf(zcmcf*(i_zone-1)+17))
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (question .eq. 'COEF_MATR_FROT') then
        call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
        vale_r = v_sdcont_caradf(zcmdf*(i_zone-1)+1)
    else if (question .eq. 'E_N') then
        call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
        vale_r = v_sdcont_caradf(zcmdf*(i_zone-1)+2)
    else if (question .eq. 'E_T') then
        call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
        vale_r = v_sdcont_caradf(zcmdf*(i_zone-1)+3)
    else if (question .eq. 'ALARME_JEU') then
        call jeveuo(sdcont_caradf, 'L', vr=v_sdcont_caradf)
        vale_r = v_sdcont_caradf(zcmdf*(i_zone-1)+5)
    else if (question .eq. 'EXCL_DIR') then
        call jeveuo(sdcont_caracf, 'L', vr=v_sdcont_caracf)
        vale_i = nint(v_sdcont_caracf(zcmcf*(i_zone-1)+12))
    else if (question .eq. 'EXCL_FROT_DIRX') then
        call jeveuo(sdcont_exclfr, 'L', vr=v_sdcont_exclfr)
        vale_r = v_sdcont_exclfr(zexcl*(i_zone-1)+1)
    else if (question .eq. 'EXCL_FROT_DIRY') then
        call jeveuo(sdcont_exclfr, 'L', vr=v_sdcont_exclfr)
        vale_r = v_sdcont_exclfr(zexcl*(i_zone-1)+2)
    else if (question .eq. 'EXCL_FROT_DIRZ') then
        call jeveuo(sdcont_exclfr, 'L', vr=v_sdcont_exclfr)
        vale_r = v_sdcont_exclfr(zexcl*(i_zone-1)+3)
    else
        write (6, *) '   NUM. ZONE    : <', i_zone, '>'
        write (6, *) '   QUESTION     : <', question, '>'
        write (6, *) '   REPONSE  - I : <', vale_i, '>'
        write (6, *) '   REPONSE  - R : <', vale_r, '>'
        write (6, *) '   REPONSE  - L : <', vale_l, '>'
        ASSERT(.false.)
    end if
!
    if (present(vale_i_)) then
        vale_i_ = vale_i
    end if
!
    if (present(vale_r_)) then
        vale_r_ = vale_r
    end if
!
    if (present(vale_l_)) then
        vale_l_ = vale_l
    end if
end subroutine
