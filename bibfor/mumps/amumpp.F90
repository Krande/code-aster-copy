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

subroutine amumpp(option, nbsol, kxmps, ldist, type,&
                  impr, ifmump, eli2lg, rsolu, csolu,&
                  vcine, prepos, lpreco, lmhpc)
!
!
    use iso_c_binding, only: c_ptr, c_f_pointer
    implicit none
!-----------------------------------------------------------------------
! BUT : ROUTINE DE PRE/POST-TRAITEMENT DE LA SOLUTION ET DU
!       SECOND MEMBRE POUR AMUMPS/C/D/Z
!
! IN  OPTION :   IN   : OPTION D'UTILISATION.
! IN  NBSOL  :   IN   : NBRE DE SYSTEMES A RESOUDRE
! IN  KXMPS  :   IN   : INDICE DE L'INSTANCE MUMPS DANS DMPS
! IN  LDIST  :  LOG   : LOGICAL MUMPS DISTRIBUE OR NOT
! IN  TYPE   :   K1   : TYPE DU POINTEUR R OU C
! IN  IMPR   :  K14   : FLAG POUR IMPRESSION MATRICE
! IN  IFMUMP :   IN   : UNITE LOGIQUE POUR IMPRESSION FICHIER
! IN  ELI2LG :  LOG   : LOGICAL POUR NE LAISSER QU'1 LAGRANGE ACTIF
! I/O RSOLU  :    R   : EN ENTREE : SECONDS MEMBRES REELS
!                     : EN SORTIE : SOLUTIONS
! I/O CSOLU  :    C   : EN ENTREE : SECONDS MEMBRES COMPLEXES
!                     : EN SORTIE : SOLUTIONS
! IN  VCINE  :  K19   : NOM DU CHAM_NO DE CHARGEMENT CINEMATIQUE
! IN  PREPOS :  LOG   : SI .TRUE. ON FAIT LES PRE ET POSTTRAITEMENTS DE
!           MISE A L'ECHELLE DU RHS ET DE LA SOLUTION (MRCONL) ET DE LA
!           PRISE EN COMPTE DES AFFE_CHAR_CINE (CSMBGG).
!           SI .FALSE. ON NE LES FAIT PAS (PAR EXEMPLE EN MODAL).
! IN  LPRECO :  LOG   : MUMPS EST-IL UTILISE COMME PRECONDITIONNEUR ?
!-----------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
!
#include "asterf_types.h"
#include "asterf.h"
#include "jeveux.h"
#include "asterc/r4maem.h"
#include "asterc/r4miem.h"
#include "asterc/r8maem.h"
#include "asterc/r8miem.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/csmbgg.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jgetptc.h"
#include "asterfort/mcconl.h"
#include "asterfort/mrconl.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nudlg2.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
    integer :: option, nbsol, kxmps, ifmump
    aster_logical :: ldist, eli2lg, prepos, lpreco, lmhpc
    character(len=1) :: type
    character(len=14) :: impr
    character(len=19) :: vcine
    real(kind=8) :: rsolu(*)
    complex(kind=8) :: csolu(*)
!
#ifdef ASTER_HAVE_MUMPS
#include "asterf_mumps.h"
    type(smumps_struc), pointer :: smpsk => null()
    type(cmumps_struc), pointer :: cmpsk => null()
    type(dmumps_struc), pointer :: dmpsk => null()
    type(zmumps_struc), pointer :: zmpsk => null()
    integer :: n, nnbsol, rang, lmat, i, ierd, idvalc, k, ifm, niv
    integer :: jj, nbeql, nuglo, jrhs, j
    integer :: jrefn, jmlogl, jdeeq, nuno, nucmp
    character(len=1) :: rouc
    character(len=4) :: etam
    character(len=8) :: mesh
    character(len=14) :: nonu
    character(len=19) :: nomat, nosolv
    character(len=24) :: vcival, nonulg
    aster_logical :: ltypr, lverif
    aster_logical, parameter :: l_debug = ASTER_FALSE
    real(kind=8) :: rr4max, raux, rmin, rmax, rtest, valr(2)
    complex(kind=8) :: cbid, caux
    integer, pointer :: delg(:) => null()
    integer, pointer :: dlg2(:) => null()
    integer, pointer :: nequ(:) => null()
    integer, pointer :: pddl(:) => null()
    integer, pointer :: nulg(:) => null()
    real(kind=8), pointer :: rsolu2(:) => null()
    complex(kind=8), pointer :: csolu2(:) => null()
    type(c_ptr) :: pteur_c
    cbid = dcmplx(0.d0, 0.d0)
!
!-----------------------------------------------------------------------
    call jemarq()
    call infdbg('SOLVEUR',ifm, niv)
! pour verifier la norme de la solution produite
    lverif=.true.
    lverif=.false.
!
!     ------------------------------------------------
!     INITS
!     ------------------------------------------------
    rr4max=r4maem()
    if (type .eq. 'S') then
        smpsk=>smps(kxmps)
        rang=smpsk%myid
        n=smpsk%n
        smpsk%nrhs=to_mumps_int(nbsol)
        smpsk%lrhs=to_mumps_int(n)
        ltypr=.true.
        rmax=r4maem()*0.5
        rmin=r4miem()*2.0
    else if (type.eq.'C') then
        cmpsk=>cmps(kxmps)
        rang=cmpsk%myid
        n=cmpsk%n
        cmpsk%nrhs=to_mumps_int(nbsol)
        cmpsk%lrhs=to_mumps_int(n)
        ltypr=.false.
        rmax=r4maem()*0.5
        rmin=r4miem()*2.0
    else if (type.eq.'D') then
        dmpsk=>dmps(kxmps)
        rang=dmpsk%myid
        n=dmpsk%n
        dmpsk%nrhs=to_mumps_int(nbsol)
        dmpsk%lrhs=to_mumps_int(n)
        ltypr=.true.
        rmax=r8maem()*0.5
        rmin=r8miem()*2.0
    else if (type.eq.'Z') then
        zmpsk=>zmps(kxmps)
        rang=zmpsk%myid
        n=zmpsk%n
        zmpsk%nrhs=to_mumps_int(nbsol)
        zmpsk%lrhs=to_mumps_int(n)
        ltypr=.false.
        rmax=r8maem()*0.5
        rmin=r8miem()*2.0
    else
        ASSERT(.false.)
    endif
    nnbsol=n*nbsol
    nomat=nomats(kxmps)
    nosolv=nosols(kxmps)
    nonu=nonus(kxmps)
    etam=etams(kxmps)
!
    vcival=vcine//'.VALE'
    call mtdscr(nomat)
    call jeveuo(nomat//'.&INT', 'E', lmat)
!
    if (lmhpc) then
        ASSERT(nbsol.eq.1)
        call jeveuo(nonu//'.NUME.PDDL', 'L', vi=pddl)
        call jeveuo(nonu//'.NUME.NULG', 'L', vi=nulg)
        call jeveuo(nonu//'.NUME.NEQU', 'L', vi=nequ)
        nbeql=nequ(1)
    else
        nbeql=n
    endif
!   Adresses needed to get the solution wrt nodes and dof numbers (see below)
    if (l_debug) then
        call jeveuo(nonu//'.NUME.REFN', 'L', jrefn)
        call jeveuo(nonu//'.NUME.DEEQ', 'L', jdeeq)
        mesh = zk24(jrefn)(1:8)
        if (lmhpc) then
            nonulg = mesh//'.NULOGL'
            call jeveuo(nonulg, 'L', jmlogl)
        endif
    endif
!
!
!
    if (option .eq. 0) then
!
        if (rang .eq. 0) then
            if (type .eq. 'S') then
                allocate(smpsk%rhs(nnbsol))
            else if (type.eq.'C') then
                allocate(cmpsk%rhs(nnbsol))
            else if (type.eq.'D') then
                allocate(dmpsk%rhs(nnbsol))
            else if (type.eq.'Z') then
                allocate(zmpsk%rhs(nnbsol))
            else
                ASSERT(.false.)
            endif
        endif
!
!       ------------------------------------------------
!        PRETRAITEMENTS ASTER DU/DES SECONDS MEMBRES :
!       ------------------------------------------------
!
!        --- PAS DE PRETRAITEMENT SI NON DEMANDE
        if (.not.lpreco .and. prepos) then
!
            if (rang .eq. 0 .or. lmhpc) then
!           --- MISE A L'ECHELLE DES LAGRANGES DANS LE SECOND MEMBRE
!           --- RANG 0 UNIQUEMENT
                if (ltypr) then
                    call mrconl('MULT', lmat, nbeql, 'R', rsolu,&
                                nbsol)
                else
                    call mcconl('MULT', lmat, nbeql, 'C', csolu,&
                                nbsol)
                endif
            endif
!
!           --- PRISE EN COMPTE DES CHARGES CINEMATIQUES :
            call jeexin(vcival, ierd)
            if (ierd .ne. 0) then
!              --- ON RAZ RSOLU SUR LES RANGS > 0 POUR NE CUMULER QUE
!              --- LA CONTRIBUTION DES CHARGES CINEMATIQUES EN DISTRIBUE
                if (ldist .and. rang .ne. 0) then
                    do i = 1, nnbsol
                        if (ltypr) then
                            rsolu(i)=0.d0
                        else
                            csolu(i)=dcmplx(0.d0,0.d0)
                        endif
                    enddo
                endif
                call jeveuo(vcival, 'L', idvalc)
                call jelira(vcival, 'TYPE', cval=rouc)
                if (ltypr) then
                    ASSERT(rouc.eq.'R')
                    do i = 1, nbsol
                        call csmbgg(lmat, rsolu(nbeql*(i-1)+1), zr(idvalc), [cbid], [cbid],&
                                    'R')
                    enddo
                else
                    ASSERT(rouc.eq.'C')
                    do i = 1, nbsol
                        call csmbgg(lmat, [0.d0], [0.d0], csolu(nbeql*(i-1)+1), zc(idvalc),&
                                    'C')
                    enddo
                endif
!
!         --- REDUCTION DU SECOND MEMBRE AU PROC MAITRE EN DISTRIBUE
!         --- POUR ETRE COHERENT AVEC LA MATRICE QUI CONTIENT DES N
!         --- SUR LA DIAGONALE
                if (ldist) then
                    if (ltypr) then
                        call asmpi_comm_vect('REDUCE', 'R', nbval=nnbsol, vr=rsolu)
                    else
                        call asmpi_comm_vect('REDUCE', 'C', nbval=nnbsol, vc=csolu)
                    endif
                endif
!
            endif
!
        endif
!
        if (lmhpc) then
            ASSERT(nbsol.eq.1)
            if (ltypr) then
                call wkvect('&&AMUMPP.RHS', 'V V R', nnbsol, jrhs)
            else
                call wkvect('&&AMUMPP.RHS', 'V V C', nnbsol, jrhs)
            endif
            if (.not.lpreco) then
                do j = 1, nbeql
                    if ( pddl(j).eq.rang ) then
                        nuglo = nulg(j)
                        if(ltypr) then
                            if(l_debug) then
                                nuno  = zi(jdeeq+2*(j-1))
                                if( nuno.ne.0 ) nuno = zi(jmlogl + nuno - 1) + 1
                                nucmp = zi(jdeeq+2*(j-1) + 1)
!                    num??ro noeud global, num comp du noeud, rhs
                                write(101+rang,*) nuno, nucmp, nuglo, rsolu(j)
                            end if
                            zr(jrhs+nuglo) = rsolu(j)
                        else
                            zc(jrhs+nuglo) = csolu(j)
                        endif
                    endif
                    if(l_debug) then
                        nuno  = zi(jdeeq+2*(j-1))
                        if( nuno.ne.0 ) nuno = zi(jmlogl + nuno - 1) + 1
                        nucmp = zi(jdeeq+2*(j-1) + 1)
!                    num??ro noeud global, num comp du noeud, rhs
                        write(201+rang,*) nuno, nucmp, rsolu(j), pddl(j), j ,nulg(j)
                    end if
                enddo
                if(l_debug) flush(101+rang)
                if(l_debug) flush(201+rang)
                if (ltypr) then
                    call asmpi_comm_vect('REDUCE', 'R', nbval=nnbsol, vr=zr(jrhs))
                    call jgetptc(jrhs, pteur_c, vr=zr(1))
                    call c_f_pointer(pteur_c, rsolu2, [nnbsol])
                else
                    call asmpi_comm_vect('REDUCE', 'C', nbval=nnbsol, vc=zc(jrhs))
                    call jgetptc(jrhs, pteur_c, vc=zc(1))
                    call c_f_pointer(pteur_c, csolu2, [nnbsol])
                endif
            else
!               if mumps is used as a preconditionner in HPC mode (asterxx),
!               each process knows the RHS and the numberings are the same
                do j = 1, nnbsol
                    if(ltypr) then
                        zr(jrhs+j-1) = rsolu(j)
                    else
                        zc(jrhs+j-1) = csolu(j)
                    endif
                enddo
                if (ltypr) then
                    call jgetptc(jrhs, pteur_c, vr=zr(1))
                    call c_f_pointer(pteur_c, rsolu2, [nnbsol])
                else
                    call jgetptc(jrhs, pteur_c, vc=zc(1))
                    call c_f_pointer(pteur_c, csolu2, [nnbsol])
                endif
            endif
        else
            if (ltypr) then
                if(l_debug) then
                    do j = 1, nnbsol
                        nuno  = zi(jdeeq+2*(j-1))
                        nucmp = zi(jdeeq+2*(j-1) + 1)
        !                num??ro noeud, num comp du noeud, solution
                        write(49,*) nuno, nucmp, j, rsolu(j)
                    end do
                    flush(49)
                end if
                call jgetptc(1, pteur_c, vr=rsolu(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call jgetptc(1, pteur_c, vc=csolu(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            endif
        endif
!
!        --- COPIE DE RSOLU DANS %RHS:
        if (rang .eq. 0) then
            if (type .eq. 'S') then
                do i = 1, nnbsol
                    raux=rsolu2(i)
                    rtest=abs(raux)
                    if (rtest .lt. rmin) then
                        raux=0.d0
                    else if (rtest.gt.rmax) then
                        raux=rmax*sign(1.d0,raux)
                    endif
                    smpsk%rhs(i)=real(raux, kind=4)
                enddo
            else if (type.eq.'C') then
                do i = 1, nnbsol
                    caux=csolu2(i)
                    rtest=abs(caux)
                    if (rtest .lt. rmin) then
                        caux=dcmplx(0.d0,0.d0)
                    else if (rtest.gt.rmax) then
                        caux=dcmplx(rmax*sign(1.d0,dble(caux)),0.d0)
                        caux=rmax*dcmplx(1.d0*sign(1.d0,dble(caux)), 1.d0*sign(1.d0,imag(caux)))
                    endif
                    cmpsk%rhs(i)=cmplx(caux, kind=4)
                enddo
            else if (type.eq.'D') then
                do i = 1, nnbsol
                    raux=rsolu2(i)
                    rtest=abs(raux)
                    if (rtest .lt. rmin) then
                        raux=0.d0
                    else if (rtest.gt.rmax) then
                        valr(1)=rtest
                        valr(2)=rmax
                        call utmess('F', 'FACTOR_79', si=i, nr=2, valr=valr)
                    endif
                    dmpsk%rhs(i)=raux
                enddo
            else if (type.eq.'Z') then
                do i = 1, nnbsol
                    caux=csolu2(i)
                    rtest=abs(caux)
                    if (rtest .lt. rmin) then
                        caux=dcmplx(0.d0,0.d0)
                    else if (rtest.gt.rmax) then
                        valr(1)=rtest
                        valr(2)=rmax
                        call utmess('F', 'FACTOR_79', si=i, nr=2, valr=valr)
                    endif
                    zmpsk%rhs(i)=caux
                enddo
            else
                ASSERT(.false.)
            endif
        endif
        if (lmhpc) call jedetr('&&AMUMPP.RHS')
!
!         -- IMPRESSION DU/DES SECONDS MEMBRES (SI DEMANDE) :
        if (impr(1:3) .eq. 'OUI') then
            if (rang .eq. 0) then
                if (type .eq. 'S') then
                    do k = 1, nnbsol
                        write(ifmump,*) k,smpsk%rhs(k)
                    enddo
                else if (type.eq.'C') then
                    do k = 1, nnbsol
                        write(ifmump,*) k,cmpsk%rhs(k)
                    enddo
                else if (type.eq.'D') then
                    do k = 1, nnbsol
                        write(ifmump,*) k,dmpsk%rhs(k)
                    enddo
                else if (type.eq.'Z') then
                    do k = 1, nnbsol
                        write(ifmump,*) k,zmpsk%rhs(k)
                    enddo
                else
                    ASSERT(.false.)
                endif
                write(ifmump,*) 'MUMPS FIN RHS'
            endif
            if (impr(1:11) .eq. 'OUI_NOSOLVE') then
                call utmess('F', 'FACTOR_71', si=ifmump)
            endif
        endif
!
!
    else if (option.eq.2) then
!
        if (lmhpc) then
            if (ltypr) then
                call wkvect('&&AMUMPP.RHS', 'V V R', nnbsol, jrhs)
                call jgetptc(jrhs, pteur_c, vr=zr(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call wkvect('&&AMUMPP.RHS', 'V V C', nnbsol, jrhs)
                call jgetptc(jrhs, pteur_c, vc=zc(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            endif
        else
            if (ltypr) then
                call jgetptc(1, pteur_c, vr=rsolu(1))
                call c_f_pointer(pteur_c, rsolu2, [nnbsol])
            else
                call jgetptc(1, pteur_c, vc=csolu(1))
                call c_f_pointer(pteur_c, csolu2, [nnbsol])
            endif
        endif
!
!       ------------------------------------------------
!        POST-TRAITEMENTS ASTER DE LA SOLUTION :
!       ------------------------------------------------
        if (rang .eq. 0 .or. lmhpc) then
            if (rang .eq. 0) then
                if (type .eq. 'S') then
                    do i = 1, nnbsol
                        rsolu2(i)=smpsk%rhs(i)
                    enddo
                    deallocate(smpsk%rhs)
                else if (type.eq.'C') then
                    do i = 1, nnbsol
                        csolu2(i)=cmpsk%rhs(i)
                    enddo
                    deallocate(cmpsk%rhs)
                else if (type.eq.'D') then
                    call dcopy(nnbsol, dmpsk%rhs, 1, rsolu2, 1)
                    deallocate(dmpsk%rhs)
                else if (type.eq.'Z') then
                    call zcopy(nnbsol, zmpsk%rhs, 1, csolu2, 1)
                    deallocate(zmpsk%rhs)
                else
                    ASSERT(.false.)
                endif
            endif
!
            if (lmhpc) then
                if (ltypr) then
                    call asmpi_comm_vect('BCAST', 'R', nbval=nnbsol, bcrank=0, vr=rsolu2)
                    if(l_debug) then
                        do j = 1, nbeql
                            nuno  = zi(jdeeq+2*(j-1))
                            if( nuno.ne.0 ) nuno = zi(jmlogl + nuno - 1) + 1
                            nucmp = zi(jdeeq+2*(j-1) + 1)
!                    numero ddl local, num??ro noeud local, num??ro noeud global, num comp du noeud,
!                                num ddl global, num proc proprio, solution
                            ! write(51+rang,*) j, zi(jdeeq+2*(j-1)), nuno, nucmp,  &
                            !                      nulg(j), pddl(j), rsolu(j)
                            write(51+rang,*) nuno, nucmp, rsolu(j)
                        end do
                        flush(51+rang)
                    end if
                    if (.not.lpreco) then
                        do j = 1, nbeql
                            rsolu(j) = rsolu2(nulg(j)+1)
                        enddo
                    else
!               if mumps is used as a preconditionner in HPC mode (asterxx),
!               each process knows the SOL and the numberings are the same
                        do j = 1, nnbsol
                            rsolu(j) = rsolu2(j)
                        enddo
                    endif
                else
                    call asmpi_comm_vect('BCAST', 'C', nbval=nnbsol, bcrank=0, vc=csolu2)
                    if (.not.lpreco) then
                        do j = 1, nbeql
                            csolu(j) = csolu2(nulg(j)+1)
                        enddo
                    else
                        do j = 1, nnbsol
                            csolu(j) = csolu2(j)
                        enddo
                    endif
                endif
            endif
!
            if (eli2lg) then
!           -- PRISE EN COMPTE DES LAGRANGES "2" :
!           -- EN SORTIE DE RESOLUTION AVEC ELIM_LAGR='LAGR2' ON A :
!           -- LAGR1 = LAGR1 + LAGR2, ON DOIT RECTIFIER CELA :
!           -- LAGR1 = LAGR1/2 PUIS LAGR2 = LAGR1
!           -- VALIDE QUE SUR PROC 0, MAIS C'EST OK CAR ON
!           -- BROADCAST LA SOLUTION APRES
                call jeveuo(nonu//'.NUME.DELG', 'L', vi=delg)
                call nudlg2(nonu)
                call jeveuo(nonu//'.NUME.DLG2', 'L', vi=dlg2)
                if (ltypr) then
                    do i = 1, nbsol
                        do k = 1, nbeql
                            if (delg(k) .eq. -1) then
                                rsolu((i-1)*nbeql+k)= 0.5d0 * rsolu((i-1)*nbeql+k)
                                jj = dlg2(k)
                                rsolu((i-1)*nbeql+jj) = rsolu((i-1)*nbeql+k)
                            endif
                        enddo
                    enddo
                else
                    do i = 1, nbsol
                        do k = 1, nbeql
                            if (delg(k) .eq. -1) then
                                csolu((i-1)*nbeql+k) = 0.5d0*csolu((i-1)*nbeql+k)
                                jj = dlg2(k)
                                csolu((i-1)*nbeql+jj) = csolu((i-1)*nbeql+k)
                            endif
                        enddo
                    enddo
                endif
            endif
!
!         --- MISE A L'ECHELLE DES LAGRANGES DANS LA SOLUTION :
!         ON NE LE FAIT PAS SI NON DEMANDE
            if (.not.lpreco .and. prepos) then
                if (ltypr) then
                    call mrconl('MULT', lmat, nbeql, 'R', rsolu,&
                                nbsol)
                else
                    call mcconl('MULT', lmat, nbeql, 'C', csolu,&
                                nbsol)
                endif
            endif
        endif
!
!       -- BROADCAST DE SOLU A TOUS LES PROC
        if (.not.lmhpc) then
            if (ltypr) then
                call asmpi_comm_vect('BCAST', 'R', nbval=nnbsol, bcrank=0, vr=rsolu)
            else
                call asmpi_comm_vect('BCAST', 'C', nbval=nnbsol, bcrank=0, vc=csolu)
            endif
        endif
!
!       -- IMPRESSION DU/DES SOLUTIONS (SI DEMANDE) :
        if ((lverif).and.(rang .eq. 0)) then
          raux=0.d0
          if (ltypr) then
            do k = 1, nnbsol
              raux=raux+abs(rsolu2(k))
            enddo
          else
            do k = 1, nnbsol
              raux=raux+abs(csolu2(k))
            enddo
          endif
          write(ifm,*) 'NORME L1 MUMPS SOLUTION=',raux
        endif
        if (impr(1:9) .eq. 'OUI_SOLVE') then
            if (rang .eq. 0) then
                if (ltypr) then
                    do k = 1, nnbsol
                        write(ifmump,*) k,rsolu2(k)
                    enddo
                else
                    do k = 1, nnbsol
                        write(ifmump,*) k,csolu2(k)
                    enddo
                endif
                write(ifmump,*) 'MUMPS FIN SOLUTION'
            endif
        endif

        if(l_debug .and. .not. lmhpc) then
            do j = 1, nnbsol
                nuno  = zi(jdeeq+2*(j-1))
                nucmp = zi(jdeeq+2*(j-1) + 1)
!                num??ro noeud, num comp du noeud, solution
                write(50,*) nuno, nucmp, rsolu(j)
            end do
            flush(50)
        end if
        if (lmhpc) call jedetr('&&AMUMPP.RHS')
    else
!       ------------------------------------------------
!        MAUVAISE OPTION
!       ------------------------------------------------
        ASSERT(.false.)
    endif
    call jedema()
#endif
end subroutine
