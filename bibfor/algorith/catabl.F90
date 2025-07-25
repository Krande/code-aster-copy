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
subroutine catabl(table_new, table_old, time, nume_store, nb_obje, &
                  obje_name, obje_sdname)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/memaxm.h"
#include "asterfort/tbacce.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: table_new
    character(len=8), intent(in) :: table_old
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: nume_store
    integer(kind=8), intent(in) :: nb_obje
    character(len=16), intent(in) :: obje_name(nb_obje)
    character(len=24), intent(in) :: obje_sdname(nb_obje)
!
! --------------------------------------------------------------------------------------------------
!
! Command CALCUL
!
! Management of result (TABLE_CONTAINER)
!
! --------------------------------------------------------------------------------------------------
!
! In  table_new        : name of created table
! In  table_old        : name of old table
! In  time             : time
! In  nume_store       : index of current step time
! In  nb_obje          : number of new objects to add
! In  obje_name        : name of new objects to add
! In  obje_sdname      : datastructure name of new objects to add
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbpara = 6
    character(len=19), parameter :: nompar(nbpara) = (/ &
                                    'NOM_OBJET ', 'TYPE_OBJET', &
                                    'NOM_SD    ', 'NUME_ORDRE', &
                                    'INST      ', 'VALE_I    '/)
    character(len=19), parameter :: typpar(nbpara) = (/ &
                                    'K16', 'K16', 'K24', 'I  ', 'R8 ', 'I  '/)
    integer(kind=8) :: prepar(nbpara)
!
    integer(kind=8), parameter :: l_nb_obje = 9
    character(len=16), parameter :: l_obje_name(l_nb_obje) = (/ &
                                    'MATR_TANG_ELEM  ', 'SIEF_ELGA       ', 'VARI_ELGA       ', &
                                    'FORC_INTE_ELEM  ', 'FORC_DIRI_ELEM  ', 'FORC_NODA_ELEM  ', &
                                    'CODE_RETOUR_INTE', 'FORC_VARC_ELEM_M', 'FORC_VARC_ELEM_P'/)
    character(len=16), parameter :: l_obje_type(l_nb_obje) = (/ &
                                    'MATR_ELEM_DEPL_R', 'CHAM_ELEM       ', 'CHAM_ELEM       ', &
                                    'VECT_ELEM_DEPL_R', 'VECT_ELEM_DEPL_R', 'VECT_ELEM_DEPL_R', &
                                    'ENTIER          ', 'VECT_ELEM_DEPL_R', 'VECT_ELEM_DEPL_R'/)
!
    character(len=19) :: nomtab
    aster_logical :: l_new_table, l_repl_object
    integer(kind=8) :: i_repl_object
    integer(kind=8) :: jnobj, jnosd, jnuor, jtobj, jrins, jlins
    integer(kind=8) :: nboldp, nblign, t_nume_store
    integer(kind=8) :: ipara, ilign, i_l_obj, i_obj, ibid, iret, nbcol
    character(len=24) :: vk(3)
    character(len=16) :: k16bid, t_obje_name, obje_type
    real(kind=8) :: r8bid, v_iret(1)
    complex(kind=8) :: c16bid
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    nomtab = table_new
    nboldp = 0
    nblign = 0
    l_new_table = .false.
    prepar(1:nbpara) = 0
!
! - New table or not ?
!
    if (table_old .eq. ' ') then
        l_new_table = .true.
    else
        l_new_table = .false.
    end if
!
! - Create new table
!
    if (l_new_table) then
        call detrsd('TABLE_CONTAINER', table_new)
        call tbcrsd(table_new, 'G')
        call tbajpa(table_new, nbpara, nompar, typpar)
    end if
!
! - Check old table
!
    if (.not. l_new_table) then
        call jeveuo(nomtab//'.TBNP', 'L', vi=tbnp)
        call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
        nboldp = tbnp(1)
        if (nboldp .ne. nbpara) then
            call utmess('F', 'CALCUL1_1')
        end if
        nblign = tbnp(2)
        do ipara = 1, nbpara
            if (tblp(1+(ipara-1)*4) .eq. nompar(ipara)) then
                prepar(ipara) = ipara
            end if
        end do
        do ipara = 1, nbpara
            if (prepar(ipara) .eq. 0) then
                call utmess('F', 'CALCUL1_2')
            end if
        end do
    end if
!
! - Memory pointer on old table
!
    if (.not. l_new_table) then
        call jeveuo(tblp(1+(prepar(5)-1)*4+3), 'L', jlins)
        call jeveuo(tblp(1+(prepar(1)-1)*4+2), 'L', jnobj)
        call jeveuo(tblp(1+(prepar(2)-1)*4+2), 'L', jtobj)
        call jeveuo(tblp(1+(prepar(3)-1)*4+2), 'E', jnosd)
        call jeveuo(tblp(1+(prepar(4)-1)*4+2), 'E', jnuor)
        call jeveuo(tblp(1+(prepar(5)-1)*4+2), 'E', jrins)
    end if
!
! - Loop on objects to add new one or replace old one
!
    do i_obj = 1, nb_obje
!
! ----- Find the type of object
!
        obje_type = ' '
        do i_l_obj = 1, l_nb_obje
            if (l_obje_name(i_l_obj) .eq. obje_name(i_obj)) then
                obje_type = l_obje_type(i_l_obj)
            end if
        end do
        ASSERT(obje_type .ne. ' ')
!
! ----- Find right line in table
!
        l_repl_object = .false.
        i_repl_object = 0
        if (l_new_table) then
            l_repl_object = .false.
            i_repl_object = 0
        else
! --------- Loop on lines in table
            do ilign = 1, nblign
                if (zi(jlins+ilign-1) .eq. 1) then
! ----------------- Current object name
                    call tbacce(nomtab, ilign, 'NOM_OBJET', 'L', ibid, &
                                r8bid, c16bid, t_obje_name)
                    call tbacce(nomtab, ilign, 'NUME_ORDRE', 'L', t_nume_store, &
                                r8bid, c16bid, k16bid)
! ----------------- New object or replace old one ?
                    if (nume_store .eq. t_nume_store .and. t_obje_name .eq. obje_name(i_obj)) then
                        l_repl_object = .true.
                        i_repl_object = ilign
                        goto 50
                    end if
                end if
            end do
        end if
50      continue
!
! ----- Add object (new line) or replace old one ?
!
        if (l_repl_object) then
            ASSERT(i_repl_object .ne. 0)
            call utmess('I', 'CALCUL1_4', sk=obje_name(i_obj), si=t_nume_store)
            call jedetr(zk24(jnosd+i_repl_object-1))
            zk24(jnosd+i_repl_object-1) = obje_sdname(i_obj)
            zi(jnuor+i_repl_object-1) = nume_store
            zr(jrins+i_repl_object-1) = time
        else
            ASSERT(i_repl_object .eq. 0)
            nbcol = nbpara
            vk(1) = obje_name(i_obj)
            vk(2) = obje_type
            vk(3) = obje_sdname(i_obj)
            iret = 0
            if (obje_name(i_obj) .eq. 'CODE_RETOUR_INTE') then
                call memaxm('MAX', obje_sdname(i_obj), 'IRET', 1, ['IRET'], v_iret, 0, [0])
                vk(2) = ' '
                vk(3) = ' '
                iret = nint(v_iret(1))
            else
                nbcol = nbpara-1
            end if
            call tbajli(nomtab, nbcol, nompar, [nume_store, iret], [time], &
                        [c16bid], vk, 0)
        end if
    end do
!
    call jedema()
!
end subroutine
