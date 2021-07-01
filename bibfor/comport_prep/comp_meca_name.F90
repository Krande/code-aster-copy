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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine comp_meca_name(nbVari     , nbVariMeca,&
                          l_excl     , vari_excl,&
                          l_kit_meta , l_mfront_offi, l_prot_comp,&
                          rela_comp  , defo_comp    , kit_comp     ,&
                          type_cpla  , post_iter    , regu_visc    ,&
                          libr_name  , subr_name    , model_mfront , model_dim,&
                          infoVari)
!
implicit none
!
#include "asterf_types.h"
#include "asterc/lcinfo.h"
#include "asterc/lcvari.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/comp_mfront_vname.h"
#include "asterfort/comp_meca_code.h"
!
integer, intent(in) :: nbVari, nbVariMeca
aster_logical, intent(in) :: l_excl
character(len=16), intent(in) :: vari_excl
aster_logical, intent(in) :: l_kit_meta, l_mfront_offi, l_prot_comp
character(len=16), intent(in) :: rela_comp, defo_comp, kit_comp(4)
character(len=16), intent(in) :: type_cpla, post_iter, regu_visc
character(len=255), intent(in) :: libr_name, subr_name
character(len=16), intent(in) :: model_mfront
integer, intent(in) :: model_dim
character(len=16), pointer :: infoVari(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  nbVari           : number of internal variables
! In  nbVariMeca       : number of internal variables for mechanic
! In  l_excl           : .true. if exception case (no names for internal variables)
! In  vari_excl        : name of internal variables if l_excl
! In  l_kit_meta       : .true. if metallurgy
! In  l_mfront_offi    : .true. if MFront official
! In  rela_comp        : RELATION comportment
! In  defo_comp        : DEFORMATION comportment
! In  kit_comp         : KIT comportment
! In  type_cpla        : plane stress method
! In  post_iter        : type of post_treatment
! In  regu_visc        : keyword for viscuous regularization
! In  libr_name        : name of library
! In  subr_name        : name of comportement in library
! In  model_mfront     : type of modelisation MFront
! In  model_dim        : dimension of modelisation (2D or 3D)
! Ptr infoVari         : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_vari_meta, nb_vari_rela, idummy, idummy2, nbVariOther
    character(len=6) :: phas_name(10)
    character(len=8) :: rela_name(30)
    integer :: iVari, i_vari_meta, i_vari_rela, iVariMeca
    character(len=16) :: comp_code_py, rela_code_py, meta_code_py
!
! --------------------------------------------------------------------------------------------------
!
    if (l_excl) then
        infoVari(1:nbVari) = vari_excl
    else
        call comp_meca_code(rela_comp_  = rela_comp , defo_comp_  = defo_comp ,&
                            type_cpla_  = type_cpla , kit_comp_   = kit_comp,&
                            post_iter_  = post_iter , regu_visc_  = regu_visc ,&
                            l_implex_   = .false._1,&
                            comp_code_py_ = comp_code_py, rela_code_py_ = rela_code_py,&
                            meta_code_py_ = meta_code_py)

! ----- Name of internal state variables
        if (l_prot_comp) then
            nbVariOther = nbVari - nbVariMeca
            do iVariMeca = 1, nbVariMeca
                infoVari(iVariMeca) = 'NoName'
            end do
            if (nbVariOther .ne. 0) then
                call lcvari(comp_code_py, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            endif

        else if (l_mfront_offi) then
            nbVariOther = nbVari - nbVariMeca
            call comp_mfront_vname(nbVariMeca, &
                                   libr_name , subr_name, model_mfront, model_dim,&
                                   infoVari)
            if (nbVariOther .ne. 0) then
                call lcvari(comp_code_py, nbVariOther, infoVari(nbVariMeca+1:nbVari))
            endif

        else
            if (l_kit_meta) then
                call lcinfo(meta_code_py, idummy, nb_vari_meta, idummy2)
                call lcinfo(rela_code_py, idummy, nb_vari_rela, idummy2)
                ASSERT(nb_vari_meta .le. 10)
                ASSERT(nb_vari_rela .le. 30)
                call lcvari(meta_code_py, nb_vari_meta, phas_name)
                call lcvari(rela_code_py, nb_vari_rela, rela_name)
                iVari = 0
                do i_vari_meta = 1, nb_vari_meta
                    do i_vari_rela = 1, nb_vari_rela
                        iVari = iVari + 1
                        infoVari(iVari) = phas_name(i_vari_meta)//'##'//rela_name(i_vari_rela)
                    enddo
                enddo
                do i_vari_rela = 1, nb_vari_rela
                    iVari = iVari + 1
                    infoVari(iVari) = rela_name(i_vari_rela)
                enddo
                iVari = iVari + 1
                infoVari(iVari) = 'INDIPLAS'
                if (defo_comp .eq. 'SIMO_MIEHE') then
                    iVari = iVari + 1
                    infoVari(iVari) = 'TRAC_EPSE'
                endif
                if (defo_comp .eq. 'GDEF_LOG') then
                    iVari = iVari + 1
                    infoVari(iVari) = 'TXX'
                    iVari = iVari + 1
                    infoVari(iVari) = 'TYY'
                    iVari = iVari + 1
                    infoVari(iVari) = 'TZZ'
                    iVari = iVari + 1
                    infoVari(iVari) = 'TXY'
                    iVari = iVari + 1
                    infoVari(iVari) = 'TXZ'
                    iVari = iVari + 1
                    infoVari(iVari) = 'TYZ'
                endif
                ASSERT(iVari .eq. nbVari)
                call lcdiscard(rela_code_py)
                call lcdiscard(meta_code_py)
            else
                call lcvari(comp_code_py, nbVari, infoVari)
            endif
        endif
        call lcdiscard(comp_code_py)
    endif
!
end subroutine
