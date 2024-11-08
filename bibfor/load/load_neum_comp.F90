! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine load_neum_comp(stop, iLoad, loadName, loadNume, loadApply, &
                          ligrel_calc, nb_in_maxi, nb_in_prep, lpain, lchin, &
                          base, resuElem, vectElem, iden_direct, name_inputz)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/copisd.h"
#include "asterfort/exisd.h"
#include "asterfort/reajre.h"
#include "asterfort/gcnco2.h"
#include "asterfort/corich.h"
#include "asterfort/load_neum_spec.h"
!
    character(len=1), intent(in) :: stop
    integer, intent(in) :: iLoad
    character(len=8), intent(in) :: loadName
    integer, intent(in) :: loadNume
    character(len=4), intent(in) :: loadApply
    character(len=19), intent(in) :: ligrel_calc
    integer, intent(in) :: nb_in_maxi, nb_in_prep
    character(len=*), intent(inout) :: lpain(nb_in_maxi), lchin(nb_in_maxi)
    character(len=19), intent(inout) :: resuElem
    character(len=19), intent(in) :: vectElem
    character(len=1), intent(in) :: base
    character(len=*), optional, intent(in) :: iden_direct, name_inputz
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation
!
! Elementary (on one load) - Vector
!
! --------------------------------------------------------------------------------------------------
!
! In  stop           : CALCUL subroutine comportement
! In  iLoad          : index of current load
! In  loadName       : name of current load
! In  loadNume       : identification of load type
! In  loadApply      : type of application for load
!                        'Dead' - Dead loads (not dependent on displacements)
!                        'Pilo' - Loads for continuation (not dependent on displacements)
!                        'Suiv' - Undead loads (dependent on displacements)
! In  ligrel_calc    : LIGREL to compute
! In  nb_in_maxi     : maximum number of input fields
! In  nb_in_prep     : number of input fields before specific ones
! IO  lpain          : list of input parameters
! IO  lchin          : list of input fields
! IO  resuElem       : name of resuElem
! In  vectElem       : name of vectElem
! In  base           : JEVEUX base to create vectElem
! In  iden_direct    : direct identification of type
! In  name_inputz    : direct name of input field
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbTypeNeum = 19
    character(len=8), parameter :: lpaout = 'PVECTUR'
    integer :: iexist, iTypeNeum, nb_in_add
    character(len=16) :: load_option
    character(len=24) :: load_ligrel
    integer :: nbout, nbin
    character(len=8) :: newnom, mesh, answer
!
! --------------------------------------------------------------------------------------------------
!
    do iTypeNeum = 1, nbTypeNeum

! ----- Get information about load
        if (present(iden_direct)) then
            call load_neum_spec(loadName, loadNume, loadApply, ligrel_calc, iTypeNeum, &
                                nbTypeNeum, nb_in_maxi, nb_in_prep, lchin, lpain, &
                                nb_in_add, load_ligrel, load_option, iden_direct=iden_direct, &
                                name_inputz=name_inputz)
        else
            call load_neum_spec(loadName, loadNume, loadApply, ligrel_calc, iTypeNeum, &
                                nbTypeNeum, nb_in_maxi, nb_in_prep, lchin, lpain, &
                                nb_in_add, load_ligrel, load_option)
        end if
!
        if (load_option .ne. 'No_Load') then
! --------- Generate new RESU_ELEM name
            newnom = resuElem(10:16)
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=iLoad)

! --------- Total number of fields
            nbin = nb_in_prep+nb_in_add
            nbout = 1

! --------- Computation (or not)
            if (load_option .eq. 'Copy_Load') then
                call copisd('CHAMP_GD', base, load_ligrel(1:8), resuElem)
            else
                call calcul(stop, load_option, load_ligrel, nbin, lchin, &
                            lpain, nbout, resuElem, lpaout, base, &
                            'OUI')
            end if

! --------- Copying output field
            call exisd('CHAMP_GD', resuElem, iexist)
            if (load_option .ne. 'Copy_Load') then
                call dismoi('NOM_MAILLA', load_ligrel, 'LIGREL', repk=mesh)
                call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
            else
                answer = 'NON'
            end if
            if (answer .eq. 'NON') then
                ASSERT((iexist .gt. 0) .or. (stop .eq. 'C'))
            end if
            call reajre(vectElem, resuElem, base)
        end if
    end do

end subroutine
