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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cfmxr0(mesh, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/cfnumn.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Create CONT_NOEU for storing contact results
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cmp = 30
    integer(kind=8), parameter :: nb_per = 4
    integer(kind=8) :: zresu, zperc
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: i_zone, i_cmp, i_node_slav, ibid
    integer(kind=8) :: nb_node_slav, node_slav_indx(1), node_slav_nume(1)
    integer(kind=8) :: nb_node_mesh, i_node, node_nume
    integer(kind=8) :: nb_cont_zone
    character(len=19) :: cnsinr, cnsper, cnoinr
    integer(kind=8) :: jcnslr
    integer(kind=8) :: jcnslp
    integer(kind=8) :: jdecne
    aster_logical :: l_cont_cont, l_cont_disc
    real(kind=8), pointer :: cnsvp(:) => null()
    real(kind=8), pointer :: cnsvr(:) => null()
    character(len=8), parameter, dimension(nb_cmp) :: list_cmp = (/ &
                                                      'CONT    ', 'JEU     ', 'RN      ', &
                                                      'RNX     ', 'RNY     ', 'RNZ     ', &
                                                      'GLIX    ', 'GLIY    ', 'GLI     ', &
                                                      'RTAX    ', 'RTAY    ', 'RTAZ    ', &
                                                      'RTGX    ', 'RTGY    ', 'RTGZ    ', &
                                                      'RX      ', 'RY      ', 'RZ      ', &
                                                      'R       ', 'HN      ', 'I       ', &
                                                      'IX      ', 'IY      ', 'IZ      ', &
                                                      'PT_X    ', 'PT_Y    ', 'PT_Z    ', &
                                                      'PROJ_X  ', 'PROJ_Y  ', 'PROJ_Z  '/)
    character(len=8), parameter, dimension(nb_per) :: list_cmp_per = (/ &
                                                    'V1      ', 'V2      ', 'V3      ', 'V4      '/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'CONTACT5_11')
    end if
!
! - Parameters
!
    l_cont_cont = cfdisl(ds_contact%sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_disc = cfdisl(ds_contact%sdcont_defi, 'FORMUL_DISCRETE')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node_mesh)
!
! - Get fields name
!
    cnsinr = ds_contact%fields_cont_node
    cnoinr = ds_contact%field_cont_node
    cnsper = ds_contact%field_cont_perc
!
! - Sizes
!
    zresu = cfmmvd('ZRESU')
    zperc = cfmmvd('ZPERC')
    ASSERT(zresu .eq. nb_cmp)
    ASSERT(zperc .eq. nb_per)
!
! - Create CONT_NOEU
!
    call cnscre(mesh, 'INFC_R', zresu, list_cmp, 'V', cnsinr)
    call jeveuo(cnsinr(1:19)//'.CNSV', 'E', vr=cnsvr)
    call jeveuo(cnsinr(1:19)//'.CNSL', 'E', jcnslr)
!
! - Init CONT_NOEU
!
    if (l_cont_cont .or. l_cont_disc) then
        do i_zone = 1, nb_cont_zone
            jdecne = mminfi(ds_contact%sdcont_defi, 'JDECNE', i_zone)
            nb_node_slav = mminfi(ds_contact%sdcont_defi, 'NBNOE', i_zone)
            do i_node_slav = 1, nb_node_slav
                node_slav_indx(1) = i_node_slav+jdecne
                call cfnumn(ds_contact%sdcont_defi, 1, node_slav_indx(1), node_slav_nume(1))
                do i_cmp = 1, zresu
                    cnsvr(zresu*(node_slav_nume(1)-1)+i_cmp) = 0.d0
                    zl(jcnslr-1+zresu*(node_slav_nume(1)-1)+i_cmp) = .true.
                end do
            end do
        end do
    else
        do i_node = 1, nb_node_mesh
            node_nume = i_node
            do i_cmp = 1, zresu
                cnsvr(zresu*(node_nume-1)+i_cmp) = 0.d0
                zl(jcnslr-1+zresu*(node_nume-1)+i_cmp) = .true.
            end do
        end do
    end if
!
! - Create CONT_NOEU_PERC
!
    if (l_cont_cont .or. l_cont_disc) then
        call cnscre(mesh, 'VARI_R', zperc, list_cmp_per, 'V', cnsper)
        call jeveuo(cnsper(1:19)//'.CNSV', 'E', vr=cnsvp)
        call jeveuo(cnsper(1:19)//'.CNSL', 'E', jcnslp)
    end if
!
! - Init CONT_NOEU_PERC
!
    if (l_cont_cont .or. l_cont_disc) then
        do i_zone = 1, nb_cont_zone
            jdecne = mminfi(ds_contact%sdcont_defi, 'JDECNE', i_zone)
            nb_node_slav = mminfi(ds_contact%sdcont_defi, 'NBNOE', i_zone)
            do i_node_slav = 1, nb_node_slav
                node_slav_indx(1) = i_node_slav+jdecne
                call cfnumn(ds_contact%sdcont_defi, 1, node_slav_indx(1), node_slav_nume(1))
                do i_cmp = 1, zperc
                    cnsvp(zperc*(node_slav_nume(1)-1)+i_cmp) = 0.d0
                    zl(jcnslp-1+zperc*(node_slav_nume(1)-1)+i_cmp) = .false.
                end do
            end do
        end do
    end if
!
! - Transform CONT_NOEU in CHAM_NO
!
    call cnscno(cnsinr, ' ', 'NON', 'V', cnoinr, 'F', ibid)
!
    call jedema()
end subroutine
