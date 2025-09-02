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

subroutine mstget(matrix, keywordfactz, nbocc, ddlsta)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/select_dof.h"
#include "asterfort/getnode.h"
#include "asterfort/pteddl.h"
#include "asterfort/rgndas.h"
#include "asterfort/typddl.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
    character(len=*), intent(in) :: keywordfactz
    character(len=*), intent(in) :: matrix
    integer(kind=8), intent(in) :: nbocc
    integer(kind=8), intent(inout) :: ddlsta(*)

!     ------------------------------------------------------------------
!     OPERATEUR : MODE_STATIQUE
!     RECUPERATION DES DDL SUR LESQUELS IL FAUT CALCULER DES MODES STATS
!     ------------------------------------------------------------------
! IN  : MATRIX : NOM DE LA MATRICE ASSEMBLEE DU SYSTEME
! IN  : KEYWORDFACT : MOT FACTEUR  'MODE_STAT', 'FORCE_NODALE', 'PSEUDO_MODE'
! IN  : NBOCC  : NOMBRE DE MOT CLE FACTEUR
! OUT : DDLSTA : TABLEAU DES DDL
!                DDLSTA(I) = 0  PAS DE MODE STATIQUE POUR LE DDL I
!                DDLSTA(I) = 1  MODE STATIQUE POUR LE DDL I
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: neq
    character(len=8) :: mesh
    character(len=14) :: nume_ddl
    character(len=16) :: keywordfact
    character(len=24) :: list_node
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iocc, ic, ieq
    integer(kind=8) :: iii, imode, ind, jcmp
    integer(kind=8) :: jind1, jind2, lacb, lact, lblo, lcmp
    integer(kind=8) :: llag, na, nac, nba
    integer(kind=8) :: nbb, nbl, nbliai, ncmp, nd, nt
    integer(kind=8) :: nb_node, nsc, ntc
    integer(kind=8) :: nb_mode_appui, idx_gd, icmp, nbcmp, nb_cmp_appui, nume_cmp(6)
    character(len=8) :: nom_appui
    character(len=19) :: nume_equa
    integer(kind=8), pointer :: list_equa(:) => null()
    integer(kind=8), pointer :: p_list_node(:) => null()
    integer(kind=8), pointer :: p_deeq(:) => null()
    integer(kind=8), pointer :: list_cmp(:) => null()
    character(len=8), pointer :: p_cata_nomcmp(:) => null()
    character(len=16), pointer :: nom_appui_cmp(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    keywordfact = keywordfactz
    call dismoi('NOM_MAILLA', matrix, 'MATR_ASSE', repk=mesh)
    call dismoi('NOM_NUME_DDL', matrix, 'MATR_ASSE', repk=nume_ddl)
    call dismoi('NB_EQUA', matrix, 'MATR_ASSE', repi=neq)
    call wkvect('&&MSTGET.LISTE.LAGRAN', 'V V I', neq, llag)
    call wkvect('&&MSTGET.LISTE.BLOQUE', 'V V I', neq, lblo)
    call wkvect('&&MSTGET.LISTE.ACTIF', 'V V I', neq, lact)
    call wkvect('&&MSTGET.LISTE.ACTBLO', 'V V I', neq, lacb)
    call typddl('LAGR', nume_ddl, neq, zi(llag), nba, &
                nbb, nbl, nbliai)
    call typddl('BLOQ', nume_ddl, neq, zi(lblo), nba, &
                nbb, nbl, nbliai)
    call typddl('ACTI', nume_ddl, neq, zi(lact), nba, &
                nbb, nbl, nbliai)
    call typddl('ACBL', nume_ddl, neq, zi(lacb), nba, &
                nbb, nbl, nbliai)
    list_node = '&&MSTGET.LIST_NODE'
!
    if (keywordfact .eq. 'MODE_STAT') then
        jind1 = lblo
        jind2 = lact
    else if (keywordfact .eq. 'FORCE_NODALE') then
        jind1 = lact
        jind2 = lblo
    else if (keywordfact .eq. 'PSEUDO_MODE') then
        jind1 = lacb
        jind2 = llag
        nb_mode_appui = 0
        do iocc = 1, nbocc
            call getvtx(keywordfact, 'NOM_APPUI', iocc=iocc, nbval=0, nbret=na)
            if (na .ne. 0) then
                nb_mode_appui = nb_mode_appui + 6
            end if
            ! FIXME : vérifier que le nom appui n'apparait pas 2 fois le même
        end do
        if ( nb_mode_appui .gt. 0) then
            call dismoi('NUME_EQUA', nume_ddl, 'NUME_DDL', repk=nume_equa)
            call dismoi('NUM_GD_SI', nume_ddl, 'NUME_DDL', repi=idx_gd)
            call jeveuo(nume_equa(1:19)//'.DEEQ', 'L', vi=p_deeq)
            call jeveuo(jexnum('&CATA.GD.NOMCMP', idx_gd), 'L', vk8=p_cata_nomcmp)
            call jelira(jexnum('&CATA.GD.NOMCMP', idx_gd), 'LONMAX', nbcmp)
            call wkvect('&&MSTGET.NOM.APPUI_CMP', 'V V K16', nb_mode_appui, vk16=nom_appui_cmp)
            nb_mode_appui = 0
        end if
    else if (keywordfact .eq. 'MODE_INTERF') then
        jind1 = lblo
        jind2 = lact
    else
        ASSERT(.false.)
    end if
!
    do iocc = 1, nbocc
        nom_appui = ' '
        if (keywordfact .eq. 'PSEUDO_MODE') then
            call getvtx(keywordfact, 'AXE', iocc=iocc, nbval=0, nbret=na)
            call getvtx(keywordfact, 'DIRECTION', iocc=iocc, nbval=0, nbret=nd)
            if ((na+nd) .ne. 0) goto 10
            call getvtx(keywordfact, 'NOM_APPUI', iocc=iocc, nbval=0, nbret=na)
            if (na .ne. 0) then
                call getvtx(keywordfact, 'NOM_APPUI', iocc=iocc, nbval=-na, scal=nom_appui)
            end if
        end if
!
! ----- Get nodes
!
        AS_ALLOCATE(vi=list_equa, size=neq)
        call getnode(mesh, keywordfact, iocc, ' ', list_node, &
                     nb_node, elem_excl=.true._1)
        if (nb_node .eq. 0) then
            call utmess('F', 'MODESTAT1_1')
        end if
        call jeveuo(list_node, 'L', vi=p_list_node)
!
! ----- Select dof
!
        call select_dof(listEqua_=list_equa, &
                        numeDofZ_=nume_ddl, &
                        nbNodeToSelect_=nb_node, listNodeToSelect_=p_list_node)
!
! ----- Get all components
!
        call getvtx(keywordfact, 'TOUT_CMP', iocc=iocc, nbval=0, nbret=ntc)
        if (ntc .ne. 0) then
            call wkvect('&&MSTGET.LISTE.CMP', 'V V I', neq, lcmp)
            do ieq = 1, neq
                zi(lcmp+ieq-1) = 1
            end do
        end if
!
! ----- Get list of components
!
        call getvtx(keywordfact, 'AVEC_CMP', iocc=iocc, nbval=0, nbret=nac)
        if (nac .ne. 0) then
            ncmp = -nac
            call wkvect('&&MSTGET.NOM.CMP', 'V V K8', ncmp, jcmp)
            call getvtx(keywordfact, 'AVEC_CMP', iocc=iocc, nbval=ncmp, vect=zk8(jcmp))
            call wkvect('&&MSTGET.LISTE.CMP', 'V V I', neq*ncmp, lcmp)
            call pteddl('NUME_DDL', nume_ddl, ncmp, zk8(jcmp), neq, &
                        tabl_equa=zi(lcmp))
            do ic = 2, ncmp
                ind = (ic-1)*neq
                do ieq = 1, neq
                    zi(lcmp+ieq-1) = max(zi(lcmp+ind+ieq-1), zi(lcmp+ieq-1))
                end do
            end do
        end if
!
! ----- Exclude list of components
!
        call getvtx(keywordfact, 'SANS_CMP', iocc=iocc, nbval=0, nbret=nsc)
        if (nsc .ne. 0) then
            ncmp = -nsc
            call wkvect('&&MSTGET.NOM.CMP', 'V V K8', ncmp, jcmp)
            call getvtx(keywordfact, 'SANS_CMP', iocc=iocc, nbval=ncmp, vect=zk8(jcmp))
            ncmp = ncmp+1
            zk8(jcmp+ncmp-1) = 'LAGR'
            call wkvect('&&MSTGET.LISTE.CMP', 'V V I', neq*ncmp, lcmp)
            call pteddl('NUME_DDL', nume_ddl, ncmp, zk8(jcmp), neq, &
                        tabl_equa=zi(lcmp))
            do ic = 2, ncmp
                ind = (ic-1)*neq
                do ieq = 0, neq-1
                    zi(lcmp+ieq) = max(zi(lcmp+ind+ieq), zi(lcmp+ieq))
                end do
            end do
            do ieq = 1, neq
                zi(lcmp+ieq-1) = 1-zi(lcmp+ieq-1)
            end do
        end if
!
! ----- If TOUT = 'OUI', deselect nodes
!
        call getvtx(keywordfact, 'TOUT', iocc=iocc, nbval=0, nbret=nt)
        if (nt .ne. 0) then
            do ieq = 1, neq
                if (zi(jind1+ieq-1) .ne. 1) then
                    zi(lcmp+ieq-1) = 0
                end if
            end do
        end if
!
! ----- Update list_equa
!
        do ieq = 1, neq
            list_equa(ieq) = list_equa(ieq)*zi(lcmp+ieq-1)
        end do
!
! ----- pour un appui, le nombre de modes correspond au nombre de composantes actives sur les noeuds
!
        if (nom_appui .ne. ' ') then
            AS_ALLOCATE(vi=list_cmp, size=nbcmp)
            do ieq = 1, neq
                if (list_equa(ieq) .eq. 1) then
                    icmp = p_deeq(2*(ieq-1)+2)
                    ASSERT(icmp .gt. 0)
                    list_cmp(icmp) = 1
                end if
            end do
            ic = 0
            do icmp = 1, nbcmp
                if (list_cmp(icmp) .eq. 1 ) then
                    ic = ic + 1
                    ASSERT(ic .le. 6)
                    nume_cmp(ic) = icmp
                    nom_appui_cmp(nb_mode_appui+ic) = nom_appui//p_cata_nomcmp(icmp)
                end if
            end do
            nb_cmp_appui = ic
            AS_DEALLOCATE(vi=list_cmp)
        end if
!
! ----- Checking:
! -----   MODE_STAT: dof must been blocked
! -----   FORCE_NODALE: dof must been free
! -----   PSEUDO_MODE: dof must been physical type
! -----   MODE_INTERF: dof must been blocked
!
        do ieq = 1, neq
            imode = list_equa(ieq)
            iii = zi(jind2+ieq-1)*imode
            if (iii .ne. 0) then
                call rgndas(nume_ddl, ieq, l_print=.true.)
                if (keywordfact .eq. 'MODE_STAT') then
                    call utmess('E', 'MODESTAT1_2')
                else if (keywordfact .eq. 'FORCE_NODALE') then
                    call utmess('E', 'MODESTAT1_3')
                else if (keywordfact .eq. 'PSEUDO_MODE') then
                    call utmess('E', 'MODESTAT1_4')
                else if (keywordfact .eq. 'MODE_INTERF') then
                    call utmess('E', 'MODESTAT1_5')
                else
                    ASSERT(.false.)
                end if
                imode = 0
            end if
            if (imode .gt. 0) then
                ! FIXME : emmettre une alarme ou une erreur si le ddl est surchargé
                ASSERT(ddlsta(ieq) .eq. 0)
                ddlsta(ieq) = imode
                if (nom_appui .ne. ' ') then
                    icmp = p_deeq(2*(ieq-1)+2)
                    do ic = 1, nb_cmp_appui
                        if (nume_cmp(ic) .eq. icmp) then
                            ddlsta(ieq) = -(nb_mode_appui+ic)*ddlsta(ieq)
                        end if
                    end do
                end if
            end if
        end do
        if (nom_appui .ne. ' ') then
            nb_mode_appui = nb_mode_appui + nb_cmp_appui
        end if
!
!        --- NETTOYAGE ---
!
        call jedetr(list_node)
        AS_DEALLOCATE(vi=list_equa)
        call jedetr('&&MSTGET.LISTE.CMP')
        call jedetr('&&MSTGET.NOM.CMP')
!
10      continue
    end do
    if (keywordfact .eq. 'PSEUDO_MODE' .and. nb_mode_appui .gt. 0) then
        call juveca('&&MSTGET.NOM.APPUI_CMP', nb_mode_appui)
    end if
    call jedetr('&&MSTGET.LISTE.LAGRAN')
    call jedetr('&&MSTGET.LISTE.BLOQUE')
    call jedetr('&&MSTGET.LISTE.ACTIF')
    call jedetr('&&MSTGET.LISTE.ACTBLO')
!
    call jedema()
end subroutine
