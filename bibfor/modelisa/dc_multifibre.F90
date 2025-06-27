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

subroutine dc_multifibre(nbocci, sdcomp)
!
! --------------------------------------------------------------------------------------------------
!
!                   DEFI_COMPOR / MULTIFIBRE
!
! --------------------------------------------------------------------------------------------------
!
! person_in_charge: jean-luc.flejou at edf.fr
!
!
    implicit none
!
    integer(kind=8), intent(in) :: nbocci
    character(len=8), intent(in) :: sdcomp
!
#include "jeveux.h"
#include "MultiFiber_type.h"
#include "asterf_types.h"
#include "asterc/lccree.h"
#include "asterc/lctest.h"
#include "asterc/lcdiscard.h"
#include "asterfort/assert.h"
#include "asterfort/comp_nbvari_std.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbDeborst, imi, jcprk, iocc, irett
    integer(kind=8) :: ibid, nbg, nbgrfib, jnmgrfib, ig, ig1, jnbfig, iaff
    integer(kind=8) :: nbVariMax, nbVari, icp, numeLaw
    integer(kind=8) :: valmi(5)
!
    aster_logical :: l_kit, l_cristal
!
    character(len=8)  :: materi, sdgf, mator
    character(len=16) :: rela_comp, defo_comp, type_cpla, rela_comp_py, kit_comp(4)
    character(len=16) :: regu_visc, post_iter, post_incr
    character(len=24) :: vnbfig, vnmfig, kgroup
    character(len=80) :: valmk(5)
!
    character(len=16), parameter :: keywordfact = 'MULTIFIBRE'
!
! --------------------------------------------------------------------------------------------------
    call jemarq()
    !
    ! on récupère les renseignements dans la sd_group_fibre
    !   noms de tous les groupes, nb maxi de groupes, nb de fibres par groupe
    nbVariMax = 0
    !
    call getvid(' ', 'GEOM_FIBRE', scal=sdgf, nbret=ibid)
    vnbfig = sdgf//'.NB_FIBRE_GROUPE'
    vnmfig = sdgf//'.NOMS_GROUPES'
    call jeveuo(vnbfig, 'L', jnbfig)
    call jelira(vnbfig, 'LONMAX', nbgrfib)
    ! le +1 c'est pour le matériau de torsion : MATER_SECT fin de routine
    call wkvect(sdcomp//'.CPRK', 'G V K24', MULTI_FIBER_SIZEK*nbgrfib+1, jcprk)
    call wkvect('&&OP0059.NOMS_GROUPES', 'V V K24', nbgrfib, jnmgrfib)
    call wkvect('&&OP0059.VERIF_AFFECT', 'V V I', nbgrfib, iaff)
    do ig = 1, nbgrfib
        zi(iaff-1+ig) = 0
    end do
    nbDeborst = 0
    !
    do iocc = 1, nbocci
        call getvtx(keywordfact, 'GROUP_FIBRE', iocc=iocc, nbval=0, nbret=nbg)
        nbg = -nbg
        if (nbg .gt. nbgrfib) then
            valmi(1) = iocc
            valmi(2) = nbg
            valmi(3) = nbgrfib
            valmk(1) = keywordfact//' / GROUP_FIBRE'
            call utmess('F', 'MODELISA8_19', ni=3, vali=valmi, nk=1, valk=valmk)
        end if
        call getvtx(keywordfact, 'GROUP_FIBRE', iocc=iocc, nbval=nbg, vect=zk24(jnmgrfib))
        call getvid(keywordfact, 'MATER', iocc=iocc, scal=materi)
        !
        ! Get name of RELATION
        call compGetRelation(keywordfact, iocc, rela_comp)
        call comp_meca_l(rela_comp, 'KIT', l_kit)
        call comp_meca_l(rela_comp, 'CRISTAL', l_cristal)
        ASSERT(.not. l_kit)
        ASSERT(.not. l_cristal)
        !
        ! Select unidimensional algorithm
        call lccree(1, rela_comp, rela_comp_py)
        call lctest(rela_comp_py, 'MODELISATION', '1D', irett)
        type_cpla = 'ANALYTIQUE'
        if (irett .eq. 0) then
            type_cpla = 'DEBORST'
            nbDeborst = nbDeborst+1
        end if
        call lcdiscard(rela_comp_py)
        !
        ! Other parameters
        defo_comp = 'VIDE'
        kit_comp(:) = 'VIDE'
        post_iter = 'VIDE'
        regu_visc = 'VIDE'
        post_incr = 'VIDE'
        ! Get number of internal state variables
        call comp_nbvari_std(rela_comp, defo_comp, type_cpla, &
                             kit_comp, post_iter, &
                             regu_visc, post_incr, &
                             nbVari, numeLaw)
        !
        do ig = 1, nbg
            ! Numéro correspondant au nom
            call jenonu(jexnom(vnmfig, zk24(jnmgrfib+ig-1)), ig1)
            if (ig1 .eq. 0) then
                call utmess('F', 'MODELISA8_8', sk=zk24(jnmgrfib+ig-1))
            end if
            icp = jcprk-1+(ig1-1)*MULTI_FIBER_SIZEK
            ! Information si déjà affecté ==> Surcharge
            if (zi(iaff-1+ig1) .eq. 1) then
                valmk(1) = zk24(icp+MULTI_FIBER_NAME)
                valmk(2) = zk24(icp+MULTI_FIBER_MATER)
                valmk(3) = zk24(icp+MULTI_FIBER_RELA)
                valmk(4) = materi
                valmk(5) = rela_comp
                call utmess('I', 'COMPOR5_19', nk=5, valk=valmk)
            end if
            zk24(icp+MULTI_FIBER_NAME) = zk24(jnmgrfib+ig-1)
            zk24(icp+MULTI_FIBER_MATER) = materi
            zk24(icp+MULTI_FIBER_RELA) = rela_comp
            zk24(icp+MULTI_FIBER_ALGO) = type_cpla
            zk24(icp+MULTI_FIBER_DEFO) = defo_comp
            write (zk24(icp+MULTI_FIBER_NBFI), '(I24)') zi(jnbfig-1+ig1)
            write (zk24(icp+MULTI_FIBER_NBVARI), '(I24)') nbVari
            zi(iaff-1+ig1) = 1
        end do
        ! on met à jour le nombre de variables internes maxi
        nbVariMax = max(nbVariMax, nbVari)
    end do
    !
    ! vérification que tout est affecté au moins une fois. Les groupes non affectes 'VIDE'
    do ig = 1, nbgrfib
        if (zi(iaff-1+ig) .eq. 0) then
            call jenuno(jexnum(vnmfig, ig), kgroup)
            icp = jcprk-1+(ig-1)*MULTI_FIBER_SIZEK
            zk24(icp+MULTI_FIBER_NAME) = kgroup
            zk24(icp+MULTI_FIBER_MATER) = 'VIDE'
        end if
    end do
    if (nbDeborst .ge. 1) then
        valmk(1) = 'double DEBORST'
        call utmess('I', 'COMPOR5_20', sk=valmk(1))
    end if
    ! on récupère le nom du matériau pour la torsion, mis à la fin
    call getvid(' ', 'MATER_SECT', scal=mator, nbret=ibid)
    icp = jcprk-1+nbgrfib*MULTI_FIBER_SIZEK
    zk24(icp+1) = mator
    call wkvect(sdcomp//'.CPRI', 'G V I', MULTI_FIBER_SIZEI, imi)
    ! type 3 = MULTIFIBRE
    zi(imi-1+MULTI_FIBER_TYPE) = 3
    zi(imi-1+MULTI_FIBER_NBVARMAX) = nbVariMax
    zi(imi-1+MULTI_FIBER_NBGRFIBR) = nbgrfib
    !
    call jedema()
end subroutine
