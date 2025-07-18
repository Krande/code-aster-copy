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

subroutine nunueq(mesh, nume_equa, nb_equa, igds, sd_iden_relaz)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jelira.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/exisdg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: nume_equa
    integer(kind=8), intent(in) :: nb_equa
    integer(kind=8), intent(in) :: igds
    character(len=*), optional, intent(in) :: sd_iden_relaz
!
! --------------------------------------------------------------------------------------------------
!
! Numbering
!
! Set NUEQ object
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh           : name of mesh
! In  nume_equa      : name of NUME_EQUA
! In  nb_equa        : number of equations
! In  igds           : index of GRANDEUR used to numbering
! In  sd_iden_rela   : name of object for identity relations between dof
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: length_prno, iexi
    character(len=8) :: node_name_term, cmp_name_term
    character(len=8) :: cmp_name
    integer(kind=8) :: i_equ, i_dof, i_ligr, i_node, i_rela, i_term, i_equ_old, i_in_rela, i_equ_sav
    integer(kind=8) :: jprno, iadg, istart
    integer(kind=8) :: nb_node, ncmpmx, nec, nb_dof, nb_ligr, nb_term, nb_iden_rela
    integer(kind=8) :: i_cmp_glob, i_cmp_loca, i_cmp
    logical :: l_in_rela, l_new_equa
    character(len=24) :: sd_iden_rela
    character(len=24) :: prno, nueq
    integer(kind=8), pointer :: v_nueq(:) => null()
    integer(kind=8), pointer :: v_rela_dof(:) => null()
    integer(kind=8), pointer :: v_sdiden_info(:) => null()
    integer(kind=8), pointer :: v_sdiden_dime(:) => null()
    integer(kind=8), pointer :: v_sdiden_iset(:) => null()
    integer(kind=8), pointer :: v_sdiden_nueq(:) => null()
    integer(kind=8), pointer :: v_rela_aux(:) => null()
    character(len=8), pointer :: v_sdiden_term(:) => null()
    character(len=8), pointer :: p_cata_cmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    nueq = nume_equa(1:19)//'.NUEQ'
    prno = nume_equa(1:19)//'.PRNO'
    call jeveuo(nueq, 'E', vi=v_nueq)
!
! - Information about GRANDEUR
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', igds), 'L', vk8=p_cata_cmp)
    call jelira(jexnum('&CATA.GD.NOMCMP', igds), 'LONMAX', ncmpmx)
    nec = nbec(igds)
    ASSERT(ncmpmx .ne. 0)
    ASSERT(nec .ne. 0)
!
!
! - Informations about identity relation
!
    nb_iden_rela = 0
    sd_iden_rela = ' '
    if (present(sd_iden_relaz)) then
        sd_iden_rela = sd_iden_relaz
        if (sd_iden_rela .ne. ' ') then
            call jeveuo(sd_iden_rela(1:19)//'.COLL', 'L', vk8=v_sdiden_term)
            call jeveuo(sd_iden_rela(1:19)//'.INFO', 'L', vi=v_sdiden_info)
            call jeveuo(sd_iden_rela(1:19)//'.DIME', 'L', vi=v_sdiden_dime)
            nb_iden_rela = v_sdiden_info(1)
            call jeexin(sd_iden_rela(1:19)//'.ISET', iexi)
            if (iexi .eq. 0) then
                call wkvect(sd_iden_rela(1:19)//'.ISET', 'V V I', 1, vi=v_sdiden_iset)
                call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', nb_node)
                AS_ALLOCATE(size=nb_node*ncmpmx, vi=v_rela_aux)
                istart = 0
                do i_rela = 1, nb_iden_rela
                    nb_term = v_sdiden_dime(i_rela)
                    do i_term = 1, nb_term
                        node_name_term = v_sdiden_term(2*(i_term-1)+1+istart)
                        cmp_name_term = v_sdiden_term(2*(i_term-1)+2+istart)
                        do i_cmp_glob = 1, ncmpmx
                            if (p_cata_cmp(i_cmp_glob) .eq. cmp_name_term) then
                                i_cmp = i_cmp_glob
                            end if
                        end do
                        i_node = char8_to_int(node_name_term)
                        v_rela_aux((i_node-1)*ncmpmx+i_cmp) = i_rela
                    end do
                    istart = istart+nb_term*2
                end do
            else
                call jeveuo(sd_iden_rela(1:19)//'.ISET', 'L', vi=v_sdiden_iset)
                if (v_sdiden_iset(1) .eq. 1) then
                    call jeveuo(sd_iden_rela(1:19)//'.NUEQ', 'L', vi=v_sdiden_nueq)
                    v_nueq(:) = v_sdiden_nueq(:)
                    go to 100
                end if
            end if
        end if
    end if
!
    if (nb_iden_rela .ne. 0) then
        AS_ALLOCATE(vi=v_rela_dof, size=nb_iden_rela)
    end if
!
! - Trivial bijection
!
    if (nb_iden_rela .eq. 0) then
        nb_dof = nb_equa
        do i_dof = 1, nb_dof
            i_equ = i_dof
            v_nueq(i_dof) = i_equ
        end do
    end if
!
! - If identity relations between dof exists
!
    i_equ = 0
    if (nb_iden_rela .ne. 0) then
        call jelira(prno, 'NMAXOC', nb_ligr)
        do i_ligr = 1, nb_ligr
            call jelira(jexnum(prno, i_ligr), 'LONMAX', length_prno)
            if (length_prno .gt. 0) then
                call jeveuo(jexnum(prno, i_ligr), 'L', jprno)
                nb_node = length_prno/(nec+2)
                do i_node = 1, nb_node
                    i_dof = zi(jprno-1+(i_node-1)*(nec+2)+1)-1
                    iadg = jprno-1+(i_node-1)*(nec+2)+3
                    i_cmp_loca = 0
                    do i_cmp_glob = 1, ncmpmx
                        if (exisdg(zi(iadg), i_cmp_glob)) then
                            cmp_name = p_cata_cmp(i_cmp_glob)
                            i_cmp_loca = i_cmp_loca+1
                            i_dof = i_dof+1
                            l_new_equa = .false.
                            if (i_ligr .gt. 1) then
!
! ----------------------------- For non-physical nodes: not identity relation possible !
!
                                l_new_equa = .true.
                                l_in_rela = .false.
                            else
!
! ----------------------------- Find this dof in identity relation
!
                                i_in_rela = v_rela_aux((i_node-1)*ncmpmx+i_cmp_glob)
                                if (i_in_rela .gt. 0) then
                                    l_in_rela = .true.
                                else
                                    l_in_rela = .false.
                                end if
!
! ----------------------------- This dof in identity relation
!
                                if (l_in_rela) then
                                    ASSERT(i_in_rela .gt. 0)
                                    i_equ_old = v_rela_dof(i_in_rela)
                                    if (i_equ_old .eq. 0) then
                                        l_new_equa = .true.
                                    else
                                        l_new_equa = .false.
                                        i_equ_sav = i_equ
                                    end if
                                else
                                    l_new_equa = .true.
                                end if
                            end if
!
! ------------------------- Set equation number
!
                            if (l_new_equa) then
                                i_equ = i_equ+1
                                if (l_in_rela) then
                                    v_rela_dof(i_in_rela) = i_equ
                                end if
                            else
                                i_equ = i_equ_old
                            end if
                            v_nueq(i_dof) = i_equ
                            if (.not. l_new_equa) then
                                i_equ = i_equ_sav
                            end if
                        end if
                    end do
                end do
            end if
        end do
        ASSERT(i_equ .lt. nb_equa)
    end if
!
    AS_DEALLOCATE(vi=v_rela_dof)
!
!
! - Save .NUEQ
!
    if (sd_iden_rela .ne. ' ') then
        call jelira(nueq, 'LONMAX', nb_dof)
        call wkvect(sd_iden_rela(1:19)//'.NUEQ', 'V V I', nb_dof, vi=v_sdiden_nueq)
        v_sdiden_nueq(:) = v_nueq(:)
        v_sdiden_iset(1) = 1
        AS_DEALLOCATE(vi=v_rela_aux)
    end if
100 continue
    call jedema()
!
end subroutine
