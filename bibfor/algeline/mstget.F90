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
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
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
    integer(kind=8), pointer :: list_equa(:) => null()
    integer(kind=8), pointer :: p_list_node(:) => null()
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
    else if (keywordfact .eq. 'MODE_INTERF') then
        jind1 = lblo
        jind2 = lact
    else
        ASSERT(.false.)
    end if
!
    do iocc = 1, nbocc
        if (keywordfact .eq. 'PSEUDO_MODE') then
            call getvtx(keywordfact, 'AXE', iocc=iocc, nbval=0, nbret=na)
            call getvtx(keywordfact, 'DIRECTION', iocc=iocc, nbval=0, nbret=nd)
            if ((na+nd) .ne. 0) goto 10
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
! ----- Checking:
! -----   MODE_STAT: dof must been blocked
! -----   FORCE_NODALE: dof must been free
! -----   PSEUDO_MODE: dof must been physical type
! -----   MODE_INTERF: dof must been blocked
!
        do ieq = 1, neq
            imode = list_equa(ieq)*zi(lcmp+ieq-1)
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
            ddlsta(ieq) = max(ddlsta(ieq), imode)
        end do
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
    call jedetr('&&MSTGET.LISTE.LAGRAN')
    call jedetr('&&MSTGET.LISTE.BLOQUE')
    call jedetr('&&MSTGET.LISTE.ACTIF')
    call jedetr('&&MSTGET.LISTE.ACTBLO')
!
    call jedema()
end subroutine
