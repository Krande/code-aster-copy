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

subroutine crnggn(chamnz)
    implicit none
#include "asterc/asmpi_allgather_i.h"
#include "asterc/asmpi_comm.h"
#include "asterf_config.h"
#include "asterf_types.h"
#include "asterf.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
    character(len=*) :: chamnz

#ifdef ASTER_HAVE_MPI
#include "mpif.h"
!
    integer(kind=8) :: ili, idprn1, idprn2, ntot, lonmax, nbno_prno
    integer(kind=8) :: nbddll, i_proc, ino, iret, nbcmp
    integer(kind=8) :: numero_noeud, numero_cmp, rang, nbproc, jrefn
    integer(kind=8) :: nec, numloc, nbddl_lag
    integer(kind=8) :: pos, i_ddl, jnbddl, gd, nec_max
    mpi_int :: mrank, msize, mpicou
    mpi_int, parameter :: one4 = to_mpi_int(1)
    integer(kind=8), pointer :: v_noext(:) => null()
    integer(kind=8), pointer :: v_deeq(:) => null()
    integer(kind=8), pointer :: v_nequ(:) => null()
    integer(kind=8), pointer :: v_delg(:) => null()
    integer(kind=8), pointer :: v_nugll(:) => null()
    integer(kind=8), pointer :: v_posdd(:) => null()
    integer(kind=8), pointer :: v_mult(:) => null()
    integer(kind=8), pointer :: v_mults(:) => null()
    integer(kind=8), pointer :: v_owner(:) => null()
    integer(kind=8), pointer :: v_mult1(:) => null()
    integer(kind=8), pointer :: v_mult2(:) => null()
    integer(kind=8), pointer :: v_mdlag(:) => null()
!
    character(len=8) :: k8bid, mesh
    character(len=19) :: nomlig, nume_equa, chamno
    character(len=24) :: owner, mult1, mult2
!
!----------------------------------------------------------------------
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
#define zzprno(ili,nunoel,l)  zi(idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
    call jemarq()
!
    chamno = chamnz
!
    call asmpi_comm('GET', mpicou)
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
!   RECUPERATION DU NOM DU MAILLAGE DANS LE BUT D'OBTENIR LE JOINT
    call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=nume_equa)
    call jeveuo(nume_equa//'.REFN', 'L', jrefn)
    mesh = zk24(jrefn) (1:8)
!
!   SI LE MAILLAGE N'EST PAS PARALLELE, ON SORT
    if (.not. isParallelMesh(mesh)) goto 999
!
    call dismoi('NUM_GD_SI', nume_equa, 'NUME_EQUA', repi=gd)
    nec = nbec(gd)
    nec_max = nec
    call asmpi_comm_vect('MPI_MAX', 'I', sci=nec_max)
    ASSERT(nec == nec_max)
!
    call jeveuo(mesh//'.NOEX', 'L', vi=v_noext)
!
    call jeveuo(nume_equa//'.DEEQ', 'L', vi=v_deeq)
    call jeveuo(nume_equa//'.NEQU', 'L', vi=v_nequ)
    call jeveuo(nume_equa//'.DELG', 'L', vi=v_delg)
    nbddll = v_nequ(1)

!   Creation de la numerotation globale
    call wkvect(nume_equa//'.NULG', 'G V I', max(nbddll, 1), vi=v_nugll)
    call wkvect(nume_equa//'.PDDL', 'G V I', max(nbddll, 1), vi=v_posdd)
    call wkvect('&&CRNUGL.MULT_DDL', 'V V I', max(nbddll, 1), vi=v_mult)
    call wkvect('&&CRNUGL.MULT_DDL2', 'V V I', max(nbddll, 1), vi=v_mults)
!
! --- Il ne faut pas changer la valeur d'initialisation car on s'en sert pour detecter
!     qui est propriétaire d'un noeud (-1 si pas propriétaire)
    v_nugll(1:nbddll) = -1
    v_posdd(1:nbddll) = -1
    v_mult(1:nbddll) = -1
    v_mults(1:nbddll) = -1
    numloc = 0
! --- On numérote les ddls physiques si le noeud est propriétaire
    do i_ddl = 1, nbddll
        numero_noeud = v_deeq((i_ddl-1)*2+1)
        numero_cmp = v_deeq((i_ddl-1)*2+2)
        if (numero_noeud .gt. 0 .and. numero_cmp .gt. 0) then
            if (v_noext(numero_noeud) == rang) then
                v_nugll(i_ddl) = numloc
                v_posdd(i_ddl) = rang
                numloc = numloc+1
            end if
        end if
    end do
!
!   RECHERCHE DES ADRESSES DU .PRNO DE .NUME
    call jeveuo(nume_equa//'.PRNO', 'E', idprn1)
    call jeveuo(jexatr(nume_equa//'.PRNO', 'LONCUM'), 'L', idprn2)
    call jelira(nume_equa//'.PRNO', 'NMAXOC', ntot, k8bid)
!
! --- On numérote les lagranges des noeuds tardifs
    nbddl_lag = 0
    do ili = 2, ntot
        call jeexin(jexnum(nume_equa//'.PRNO', ili), iret)
        if (iret .ne. 0) then
            call jelira(jexnum(nume_equa//'.PRNO', ili), 'LONMAX', lonmax)
            nbno_prno = lonmax/(nec+2)
            call jenuno(jexnum(nume_equa//'.LILI', ili), nomlig)
            owner = nomlig//'.PNOE'
            mult1 = nomlig//'.MULT'
            mult2 = nomlig//'.MUL2'
            call jeveuo(owner, 'L', vi=v_owner)
            call jeveuo(mult1, 'L', vi=v_mult1)
            call jeveuo(mult2, 'L', vi=v_mult2)
            do ino = 1, nbno_prno
                ! Le proc est proprio du noeud
                if (v_owner(ino) == rang) then
                    i_ddl = zzprno(ili, ino, 1)
                    nbcmp = zzprno(ili, ino, 2)
                    ASSERT(nbcmp .eq. 1)
                    ASSERT(v_nugll(i_ddl) == -1)
                    ASSERT(v_posdd(i_ddl) == -1)
                    ASSERT(v_mult(i_ddl) == -1)
                    ASSERT(v_mults(i_ddl) == -1)
                    v_nugll(i_ddl) = numloc
                    v_posdd(i_ddl) = rang
                    v_mult(i_ddl) = v_mult2(ino)
                    v_mults(i_ddl) = v_mult1(ino)
                    nbddl_lag = nbddl_lag+1
                    numloc = numloc+1
                end if
            end do
        end if
    end do
!
    call wkvect('&&CRNULG.NBDDLL', 'V V I', nbproc, jnbddl)
!
!   ON CHERCHE LE NOMBRE TOTAL DE DDL AINSI QUE LE DECALAGE
!   DE NUMEROTATION A APPLIQUER POUR CHAQUE PROC
!
!   Chacun envoie le nb de ddl qu'il possède
    call asmpi_allgather_i([numloc], one4, zi(jnbddl), one4, mpicou)

    do i_proc = 2, nbproc
        zi(jnbddl-1+i_proc) = zi(jnbddl-1+i_proc)+zi(jnbddl-1+i_proc-1)
    end do
!   Nombre total de degré de liberté
    v_nequ(2) = zi(jnbddl-1+nbproc)
    do i_proc = nbproc, 2, -1
        zi(jnbddl-1+i_proc) = zi(jnbddl-1+i_proc-1)
    end do
    zi(jnbddl-1+1) = 0
!
    if (nbddl_lag .ne. 0) then
        call wkvect(nume_equa//'.MDLA', 'G V I', 3*nbddl_lag, vi=v_mdlag)
    end if

    pos = 0
!   Decalage de la numerotation
    do i_ddl = 1, nbddll
        if (v_delg(i_ddl) .ne. 0 .and. v_posdd(i_ddl) == rang) then
            v_mdlag(3*pos+1) = i_ddl
            v_mdlag(3*pos+2) = v_mult(i_ddl)
            v_mdlag(3*pos+3) = v_mults(i_ddl)
            pos = pos+1
        end if
        if (v_nugll(i_ddl) .ne. -1) then
            v_nugll(i_ddl) = v_nugll(i_ddl)+zi(jnbddl-1+rang+1)
        end if
    end do
    ASSERT(nbddl_lag .eq. pos)
!
! --- Pour debuggage en hpc
    if (ASTER_FALSE) then
        print *, "DEBUG IN CRNLGN"
        do i_ddl = 1, nbddll
! numero ddl local, numéro noeud local,  num composante du noeud,
!            num ddl global, num proc proprio
            write (120+rang, *) i_ddl, v_deeq((i_ddl-1)*2+1), &
                v_deeq((i_ddl-1)*2+2), v_nugll(i_ddl), v_posdd(i_ddl)
        end do
        flush (120+rang)
    end if
!
    call jedetr("&&CRNUGL.MULT_DDL")
    call jedetr("&&CRNUGL.MULT_DDL2")
    call jedetr('&&CRNULG.NBDDLL')
!
999 continue
!
    call jedema()
#else
    character(len=14) :: k14
    k14 = chamnz
#endif
!
end subroutine
