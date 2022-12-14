! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine cpysol(nomat, numddl, rsolu, debglo, vecpet, nbval)
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_recv_r.h"
#include "asterc/asmpi_send_r.h"
#include "asterc/asmpi_sendrecv_r.h"
#include "asterc/loisem.h"
#include "asterf_config.h"
#include "asterf_petsc.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterf_debug.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/codlet.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/mrconl.h"
#include "asterfort/wkvect.h"
#include "asterfort/crnustd.h"
#include "asterfort/create_graph_comm.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
#ifdef ASTER_HAVE_PETSC
    PetscInt :: debglo, nbval
#else
    integer(kind=4) :: debglo, nbval
#endif
    real(kind=8) :: rsolu(*), vecpet(*)
    character(len=14) :: numddl
    character(len=19) :: nomat
#ifdef ASTER_HAVE_MPI
#include "mpif.h"
!
    integer :: rang, nbproc, numpro, jjointr, jjointe, lmat
    integer :: lgenvo, lgrecep, jvaleue, jvaleur, iaux, jaux, jnulg
    integer :: jprddl, jnequ, nloc, nlili, ili, iret, ijoin
    integer :: numglo, jdeeq, jrefn, jmlogl, nuno1, nucmp1, numloc
    integer :: iret1, iret2, jjoine, nbnoee, idprn1, idprn2, nec, dime
    integer :: nunoel, l, jjoinr, jnujoi1, jnujoi2, nbnoer, nddll, ntot
    integer :: numnoe, step, nb_comm
    aster_logical :: ldebug
    integer, pointer :: v_nuls(:) => null()
    integer, pointer :: v_deeg(:) => null()
    integer, pointer :: v_comm(:) => null()
    integer, pointer :: v_tag(:) => null()
    integer, save :: nstep = 0
!
    mpi_int :: n4r, n4e, iaux4, tag4, numpr4
    mpi_int :: mrank, msize, iermpi, mpicou
!
    character(len=4) :: chnbjo
    character(len=8) :: k8bid, noma
    character(len=19) :: nomlig, comm_name, tag_name
    character(len=24) :: nojoinr, nojoine, nonulg, join
!----------------------------------------------------------------------
    integer :: zzprno
!
!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE
!
    zzprno(ili, nunoel, l) = zi(idprn1 - 1 + zi(idprn2 + ili - 1) + (nunoel - 1)*(nec + 2) + l - 1)
!
    call jemarq()
!
    step = -1
    ldebug = ASTER_FALSE .and. step == nstep
    nstep = nstep + 1
!
    call asmpi_comm('GET', mpicou)
!
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    ASSERT(nbproc <= MT_DOMMAX)
    DEBUG_MPI('cpysol', rang, nbproc)
!
    comm_name = '&&CPYSOL.COMM'
    tag_name = '&&CPYSOL.TAG'
    call create_graph_comm(numddl, "NUME_DDL", nb_comm, comm_name, tag_name)
    call jeveuo(comm_name, 'L', vi=v_comm)
    call jeveuo(tag_name, 'L', vi=v_tag)
!
    call jeveuo(numddl//'.NUME.NULG', 'L', jnulg)
    call jeveuo(numddl//'.NUME.PDDL', 'L', jprddl)
    call jeveuo(numddl//'.NUME.NEQU', 'L', jnequ)
    nloc = zi(jnequ)
!
    do iaux = 0, nloc
        if (zi(jprddl + iaux) .eq. rang) then
            numglo = zi(jnulg + iaux)
            rsolu(iaux + 1) = vecpet(numglo - debglo + 1)
        end if
    end do
!
    do iaux = 1, nb_comm
        numpro = v_comm(iaux)
        call codlet(numpro, 'G', chnbjo)
        nojoinr = numddl//'.NUMER'//chnbjo
        nojoine = numddl//'.NUMEE'//chnbjo
        call jeexin(nojoine, iret1)
        call jeexin(nojoinr, iret2)
        if ((iret1 + iret2) .ne. 0) then
            call jeveuo(nojoinr, 'L', jjointr)
            call jelira(nojoinr, 'LONMAX', lgrecep, k8bid)
            call jeveuo(nojoine, 'L', jjointe)
            call jelira(nojoine, 'LONMAX', lgenvo, k8bid)
            ASSERT(lgenvo .gt. 0 .and. lgrecep .gt. 0)
!
            call wkvect('&&CPYSOL.TMP1E', 'V V R', lgenvo, jvaleue)
            call wkvect('&&CPYSOL.TMP1R', 'V V R', lgrecep, jvaleur)
            do jaux = 0, lgenvo - 1
                numloc = zi(jjointe + jaux)
                ASSERT(zi(jprddl + numloc - 1) .eq. rang)
                zr(jvaleue + jaux) = rsolu(numloc)
            end do
!
            n4e = lgenvo
            n4r = lgrecep
            tag4 = to_mpi_int(v_tag(iaux))
            numpr4 = numpro
            call asmpi_sendrecv_r(zr(jvaleue), n4e, numpr4, tag4, &
                                  zr(jvaleur), n4r, numpr4, tag4, mpicou)

            do jaux = 0, lgrecep - 1
                numloc = zi(jjointr + jaux)
                ASSERT(zi(jprddl + numloc - 1) .eq. numpro)
                rsolu(numloc) = zr(jvaleur + jaux)
            end do
            call jedetr('&&CPYSOL.TMP1E')
            call jedetr('&&CPYSOL.TMP1R')
        end if
    end do
!
    call jedetr(comm_name)
    call jedetr(tag_name)

!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(numddl//'.NUME.PRNO', 'E', idprn1)
    call jeveuo(jexatr(numddl//'.NUME.PRNO', 'LONCUM'), 'L', idprn2)
    call jelira(jexnum(numddl//'.NUME.PRNO', 1), 'LONMAX', ntot, k8bid)

!   RECUPERATION DU NOM DU MAILLAGE DANS LE BUT D'OBTENIR LE JOINT
    call jeveuo(numddl//'.NUME.REFN', 'L', jrefn)
    noma = zk24(jrefn)

    call jeveuo(noma//'.DIME', 'L', dime)

!   !!! VERIFIER QU'IL N'Y A PAS DE MACRO-ELTS
!   CALCUL DU NOMBRE D'ENTIERS CODES A PARTIR DE LONMAX
    nec = ntot/zi(dime) - 2
    call jelira(numddl//'.NUME.PRNO', 'NMAXOC', nlili, k8bid)
    do ili = 2, nlili
        call jenuno(jexnum(numddl//'.NUME.LILI', ili), nomlig)
        call create_graph_comm(nomlig, "LIGREL", nb_comm, comm_name, tag_name)
        call jeveuo(comm_name, 'L', vi=v_comm)
        call jeveuo(tag_name, 'L', vi=v_tag)
        do ijoin = 1, nb_comm
            numpro = v_comm(ijoin)
            numpr4 = to_mpi_int(numpro)
            tag4 = to_mpi_int(v_tag(ijoin))
            call codlet(numpro, 'G', chnbjo)
            nojoine = nomlig//'.E'//chnbjo
            nojoinr = nomlig//'.R'//chnbjo

            call jeexin(nojoine, iret1)
            if (iret1 .ne. 0) then
                call jeveuo(nojoine, 'L', jjoine)
                call jelira(nojoine, 'LONMAX', nbnoee, k8bid)
                call wkvect('&&CRNUGL.NUM_DDL_GLOB_E', 'V V R', nbnoee, jnujoi1)
                do jaux = 1, nbnoee
                    numnoe = -zi(jjoine + jaux - 1)
                    nddll = zzprno(ili, numnoe, 1)
                    zr(jnujoi1 + jaux - 1) = rsolu(nddll)
                end do
                n4e = to_mpi_int(nbnoee)
            end if

            call jeexin(nojoinr, iret2)
            if (iret2 .ne. 0) then
                call jeveuo(nojoinr, 'L', jjoinr)
                call jelira(nojoinr, 'LONMAX', nbnoer, k8bid)
                call wkvect('&&CRNUGL.NUM_DDL_GLOB_R', 'V V R', nbnoer, jnujoi2)
                n4r = to_mpi_int(nbnoer)
            end if

            if (rang .lt. numpro) then
                if (iret1 .ne. 0) then
                    call asmpi_send_r(zr(jnujoi1), n4e, numpr4, tag4, mpicou)
                end if
                if (iret2 .ne. 0) then
                    call asmpi_recv_r(zr(jnujoi2), n4r, numpr4, tag4, mpicou)
                end if
            else if (rang .gt. numpro) then
                if (iret2 .ne. 0) then
                    call asmpi_recv_r(zr(jnujoi2), n4r, numpr4, tag4, mpicou)
                end if
                if (iret1 .ne. 0) then
                    call asmpi_send_r(zr(jnujoi1), n4e, numpr4, tag4, mpicou)
                end if
            end if

            if (iret2 .ne. 0) then
                do jaux = 1, nbnoer
                    numnoe = -zi(jjoinr + jaux - 1)
                    nddll = zzprno(ili, numnoe, 1)
                    rsolu(nddll) = zr(jnujoi2 + jaux - 1)
                end do
            end if
            call jedetr('&&CRNUGL.NUM_DDL_GLOB_E')
            call jedetr('&&CRNUGL.NUM_DDL_GLOB_R')
        end do
        call jedetr(comm_name)
        call jedetr(tag_name)
    end do
!
! -- REMISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION
    call jeveuo(nomat//'.&INT', 'L', lmat)
    call mrconl('MULT', lmat, 0, 'R', rsolu, 1)
!
!
! -- debug
    if (ldebug) then
        print*, "DEBUG IN CPYSOL"
        call jeexin(numddl//'.NUME.NULS', iret)
        if(iret == 0) then
            call crnustd(numddl)
        end if
        call jeveuo(numddl//'.NUME.NULS', 'L', vi=v_nuls)
        call jeveuo(numddl//'.NUME.DEEG', 'L', vi=v_deeg)
        do iaux = 1, nloc
            if (zi(jprddl + iaux - 1) .eq. rang) then
                nuno1  = v_deeg(2*(iaux-1) + 1)
                nucmp1 = v_deeg(2*(iaux-1) + 2)
                write (30 + rang, *) nuno1, nucmp1, v_nuls(iaux), rsolu(iaux)
                !,zi(jnulg - 1 + iaux)
            end if
        end do
        flush (30 + rang)
    end if

    call jedema()
#endif
!
end subroutine
