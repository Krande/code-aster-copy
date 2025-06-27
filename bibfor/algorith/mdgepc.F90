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
!
subroutine mdgepc(neq, nbmode, bmodal, xgene, u, &
                  kacce, kprof, inst, instt, indice, &
                  taille, kcham)
!     CONVERTIT EN BASE PHYSIQUE MODE VECTEUR
!     REM U EST INITIALISE A 0.
!     OPTIONNELLEMENT: DIFFERENTES OPTIONS D'ACCELERATION
!         - PARALLELISME MPI (SUR DDL) x OPENMP (SUR DDLxMODE)
!         - SANS BLAS, AVEC BLAS2, AVEC BLAS3
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
!        NOM        MODE                    ROLE
!  ________________ ____ ______________________________________________
! IN    NEQ            <--   NB D'EQUATIONS DU SYSTEME ASSEMBLE
! IN    NBMODE         <--   NB DE MODES NORMAUX CONSIDERES
! IN    BMODAL         <--   BASE MODALE CONSIDEREE
! IN    XGENE          <--   VECTEUR DES COORDONNEES GENERALISEES
! OUT   U              <--   VECTEUR DES COORDONNEES PHYSIQUES
!
!    OPTIONNELS
! IN    KACCE          <--   PILOTAGE DE L'ACCELERATION: 'CAS0' A 'CAS3'
!                      CAS0: SANS PARALLÉLISME ET SANS BLAS (PAR DEFAUT)
!                      CAS1: AVEC MPI (DDL) X OPENMP (DDL X MODAL)
!                            JOINDRE INDICE/TAILLE ET, SI POSSIBLE, INST.
!                            SINON PREPARATION MPI RECALCULEE A CHAQUE PASSAGE
!                      CAS2: AVEC MPI (DDL) X OPENMP (DDL X MODAL) + BLAS2
!                            JOINDRE INDICE/TAILLE/INST
!                      CAS3: AVEC MPI (DDL) X OPENMP (DDL X MODAL) + BLAS3
!                            JOINDRE INDICE/TAILLE/INST/INSTT/KCHAM
!
!
! IN    KPROF          <--   PROFILING: 'OUI', 'NON', 'OUIA' (OUI + AFFICHAGE)
!                            'NON' PAR DEFAUT
!
! IN    INST           <--   INSTANT DE CALCUL
!                            SI -NB_MAX_INST, C'EST LE DERNIER, ON FAIT LE MENAGE
!                            OBLIGATOIRE SI KACCE='CAS2' OU 'CAS3'
!
! IN    INSTT          <--   NBRE TOTAL D'INSTANT DE CALCUL
!                            OBLIGATOIRE SI KACCE='CAS3'
!
! IN    KCHAM          <--   NOM DU CHAMNO RESULTAT, OBLIGATOIRE SI KACCE='CAS3'
!
! INOUT  INDICE        <--   INDICE DE DECALAGE DEPENDANT DU RANG MPI
!                            OBLIGATOIRE SI KACCE='CAS1' OU 'CAS2' OU 'CAS3'
! INOUT  TAILLE        <--   TAILLE DU BLOC DEPENDANT DU RANG MPI
!                            OBLIGATOIRE SI KACCE='CAS1' OU 'CAS2' OU 'CAS3'
! .________________.____.______________________________________________.
    implicit none
! aslint: disable=C1513
!
#include "jeveux.h"
#include "blas/zcopy.h"
#include "blas/zgemv.h"
#include "blas/zgemm.h"
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vecinc.h"
#include "asterfort/wkvect.h"
! PARAMETRES
    integer(kind=8), intent(in) :: neq, nbmode
    real(kind=8), intent(in) :: bmodal(neq, nbmode)
    complex(kind=8), intent(in) :: xgene(nbmode)
    complex(kind=8), intent(out) :: u(neq)
    integer(kind=8), intent(in), optional :: inst, instt
    integer(kind=8), intent(inout), optional :: indice, taille
    character(len=4), intent(in), optional :: kacce, kprof
    character(len=24), intent(in), optional :: kcham
! VARIABLES LOCALES
! Taille du bloc d'instants si KACCE='CAS3'
    integer(kind=8) :: nbinst
!
    integer(kind=8) :: i, j, rang, nbproc, ja, jb, jc, nbloc, iaux, nrest, m
    integer(kind=8) :: ietdeb, ietrat, ietmax, ietfin, ifm, niv, instl, insttl
    integer(kind=8) :: nbloci, nresti, jchamno, jkchamno, borne1, borne2, iloc
    real(kind=8) :: tempsl, tempst
    complex(kind=8) :: zero, un
    character(len=4) :: kaccl, kprol
    character(len=24) :: kgemv, kgemmb, kgemmc, kchamno, k24j
    mpi_int :: mpicou, mpicow, mrang, mnbproc
    aster_logical :: ltest
    blas_int :: b_incx, b_incy, b_n
    blas_int :: b_k, b_lda, b_ldb, b_ldc, b_m
    data tempst/0.d0/
    save tempst
!-----------------------------------------------------------------------
    call jemarq()
!---------------------------------------------------
! Préparation des données
!---------------------------------------------------
! Init. globales
    zero = dcmplx(0.d0, 0.d0)
    un = dcmplx(1.d0, 0.d0)
    call infniv(ifm, niv)
    if (.not. present(kacce)) then
        kaccl = 'CAS0'
    else
        ltest = ( &
                (kacce .eq. 'CAS0') .or. (kacce .eq. 'CAS1') .or. (kacce .eq. 'CAS2') .or. &
                (kacce .eq. 'CAS3') &
                )
        ASSERT(ltest)
        kaccl = kacce
    end if
    if (.not. present(kprof)) then
        kprol = 'NON'
    else
        ltest = ( &
                (kprof(1:3) .eq. 'OUI') .or. (kprof(1:3) .eq. 'NON') .or. &
                (kprof(1:4) .eq. 'OUIA') &
                )
        ASSERT(ltest)
        kprol = kprof
    end if
    if (present(inst)) then
        instl = inst
    else
        instl = 1
    end if
    if (present(instt)) then
        insttl = instt
    else
        insttl = 1
    end if
! Affichage/profiling
    if (kprol(1:3) .eq. 'OUI') call system_clock(ietdeb, ietrat, ietmax)
!
! Init. par option
    ASSERT((neq .ge. 1) .and. (nbmode .ge. 1))
    rang = -999
    nbproc = -999
    nbloc = -999
    m = neq
    iaux = 1
    if ((kaccl .eq. 'CAS1') .or. (kaccl .eq. 'CAS2') .or. (kaccl .eq. 'CAS3')) then
        ASSERT(present(indice) .and. present(taille))
        if ((kaccl .eq. 'CAS2') .or. (kaccl .eq. 'CAS3')) then
            ASSERT(present(inst))
        end if
        if (kaccl .eq. 'CAS3') then
            ASSERT(present(kcham) .and. present(instt))
        end if
        if (abs(instl) .eq. 1) then
! On prend abs(instl) pour ne pas oublier le cas particulier: un seul pas de temps
! Premier passage, on calcul indice (iaux) et bloc (m)
            call asmpi_comm('GET_WORLD', mpicow)
            call asmpi_comm('GET', mpicou)
            if (mpicow .ne. mpicou) then
                ASSERT(.False.)
            end if
            call asmpi_info(mpicou, mrang, mnbproc)
            rang = to_aster_int(mrang)
            nbproc = to_aster_int(mnbproc)
            nbloc = neq/nbproc
            ASSERT(nbloc .ge. 0)
            if (nbloc .ge. 1) then
! Parallelisme MPI en ddl
                nrest = neq-(nbloc*nbproc)
                if (rang .eq. 0) then
                    m = nbloc+nrest
                    iaux = 1
                else
                    m = nbloc
                    iaux = nbloc+nrest+(rang-1)*nbloc+1
                end if
            else
! Pas de parallelisme MPI en ddl, pas assez de ddl
                m = neq
                iaux = 1
            end if
! Pour transmettre routine appelante et ainsi mutualiser pour les autres passages
            taille = m
            indice = iaux
        else
! Tous les autres passages, on les récupère pour éviter de les recalculer
            m = taille
            iaux = indice
        end if
        ASSERT((m .ge. 1) .and. (iaux .ge. 1))
    end if
!
!---------------------------------------------------
! Calcul proprement dit, différent suivant l'option
!---------------------------------------------------
!
!    write(ifm,*)'<mdgepc> kacce/m/iaux/neq/nbmode/instl/insttl=', &
!                     kaccl,m,iaux,neq,nbmode,instl, insttl
    if (kaccl(1:4) .eq. 'CAS0') then
! Sans parallélisme et sans BLAS
!
        do i = 1, neq
            u(i) = zero
            do j = 1, nbmode
                u(i) = u(i)+bmodal(i, j)*xgene(j)
            end do
        end do
!
    else if (kaccl(1:4) .eq. 'CAS1') then
! Avec MPI (ddl) x OpenMP (ddl x modal) + sans BLAS
!
        call vecinc(neq, zero, u)
        borne1 = iaux
        borne2 = iaux+m-1
        !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
        do i = borne1, borne2
            do j = 1, nbmode
                u(i) = u(i)+bmodal(i, j)*xgene(j)
            end do
        end do
        !$OMP END PARALLEL DO
! COM MPI que si au moins 1 ddl par processus
        if (m .ne. neq) call asmpi_comm_vect('MPI_SUM', 'C', nbval=neq, vc=u)
!
    else if (kaccl(1:4) .eq. 'CAS2') then
! Avec MPI (ddl) x OpenMP (ddl x modal) + avec BLAS2
!
        kgemv = '&&MDGEPH.DGEMV'
! Premier pas de temps global: on créé les buffers
        if (abs(instl) .eq. 1) then
            call wkvect(kgemv, 'V V C', m*nbmode, ja)
! A la main car pas de BLAS1 pour faire une copie de reels à complexes
            !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
            do j = 1, nbmode
                do i = 1, m
                    zc(ja+(j-1)*m+i-1) = bmodal(iaux+i-1, j)*un
                end do
            end do
            !$OMP END PARALLEL DO
        else
            call jeveuo(kgemv, 'L', ja)
        end if
        call vecinc(neq, zero, u)
        b_lda = to_blas_int(m)
        b_m = to_blas_int(m)
        b_n = to_blas_int(nbmode)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zgemv('N', b_m, b_n, un, zc(ja), &
                   b_lda, xgene, b_incx, zero, u(iaux), &
                   b_incy)
! COM MPI que si au moins 1 ddl par processus
        if (m .ne. neq) call asmpi_comm_vect('MPI_SUM', 'C', nbval=neq, vc=u)
!
! dernier pas de temps global: on détruit les buffers
        if (instl .lt. 0) call jedetr(kgemv)
!
    else if (kaccl(1:4) .eq. 'CAS3') then
! Avec MPI (ddl) x OpenMP (ddl x modal) + avec BLAS3
!
        kgemv = '&&MDGEPH.DGEMV'
        kgemmb = '&&MDGEPH.DGEMMB'
        kgemmc = '&&MDGEPH.DGEMMC'
        kchamno = '&&MDGEPH.CHAMNO'
! Taille de bloc d'instant
        nbinst = min(32, insttl)
! Nbre de blocs d'instants jusqu'au pas de temps courant
        nbloci = abs(instl)/nbinst
! Nbre d'instants dans le dernier bloc, si=0, calcul effectif
        nresti = abs(instl)-(nbloci*nbinst)
! Indice local dans le dernier bloc d'instants
        if (nresti .eq. 0) then
            iloc = nbinst
        else
            iloc = nresti
        end if
! Cas particulier de fin de transitoire
        if (instl .lt. 0) nresti = 0
!
! Premier pas de temps global: on créé les buffers
        if (abs(instl) .eq. 1) then
! On rempli le buffer correspondant à bmodal: matrice A
            call wkvect(kgemv, 'V V C', m*nbmode, ja)
! A la main car pas de BLAS1 pour faire une copie de reels à complexes
            !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
            do j = 1, nbmode
                do i = 1, m
                    zc(ja+(j-1)*m+i-1) = bmodal(iaux+i-1, j)*un
                end do
            end do
            !$OMP END PARALLEL DO
            call wkvect(kgemmb, 'V V C', nbmode*nbinst, jb)
            call wkvect(kgemmc, 'V V C', m*nbinst, jc)
            call wkvect(kchamno, 'V V K24', nbinst, jkchamno)
        else
            call jeveuo(kgemv, 'L', ja)
            call jeveuo(kgemmb, 'E', jb)
            call jeveuo(kchamno, 'E', jkchamno)
        end if
        zk24(jkchamno-1+iloc) = kcham
!
! On rempli le buffer correspondant à xgene: matrice B
        b_n = to_blas_int(nbmode)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call zcopy(b_n, xgene, b_incx, zc(jb+(iloc-1)*nbmode), b_incy)
!
! Calcul effectif
        if (nresti .eq. 0) then
            call jeveuo(kgemmc, 'E', jc)
            call vecinc(m*iloc, zero, zc(jc))
            b_ldc = to_blas_int(m)
            b_ldb = to_blas_int(nbmode)
            b_lda = to_blas_int(m)
            b_m = to_blas_int(m)
            b_n = to_blas_int(iloc)
            b_k = to_blas_int(nbmode)
            call zgemm('N', 'N', b_m, b_n, b_k, &
                       un, zc(ja), b_lda, zc(jb), b_ldb, &
                       zero, zc(jc), b_ldc)
            do j = 1, iloc
                k24j = zk24(jkchamno-1+j)
                call jeveuo(k24j, 'E', jchamno)
                call vecinc(neq, zero, zc(jchamno))
                b_n = to_blas_int(m)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call zcopy(b_n, zc(jc+(j-1)*m), b_incx, zc(jchamno-1+iaux), b_incy)
! COM MPI que si au moins 1 ddl par processus
                if (m .ne. neq) call asmpi_comm_vect('MPI_SUM', 'C', nbval=neq, vc=zc(jchamno))
                call jelibe(k24j)
            end do
        end if
!
! dernier pas de temps global: on détruit les buffers
        if (instl .lt. 0) then
            call jedetr(kgemv)
            call jedetr(kgemmb)
            call jedetr(kgemmc)
            call jedetr(kchamno)
        end if
!
    else
! Mauvaise option
        ASSERT(.false.)
!
    end if
!
! Affichage/profiling
    if (kprol(1:3) .eq. 'OUI') then
        call system_clock(ietfin)
        tempsl = real(ietfin-ietdeb)/real(ietrat)
        tempst = tempst+tempsl
        if (kprol(1:4) .eq. 'OUIA') then
            write (ifm, *) '<mdgepc> kacce/m/iaux/neq/nbmode/instl=', &
                kaccl, m, iaux, neq, nbmode, instl
            write (ifm, *) '<mdgepc> temps_local/temps_total=', tempsl, tempst
        end if
    end if
    call jedema()
end subroutine
