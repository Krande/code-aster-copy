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
subroutine mdgep3(neq, nbexci, psidel, temps, nomfon, tab, kacce, kprof, &
                  inst, indice, taille)
!    MULTI-APPUIS :
!    CONVERSION LES DDL GENERALISES EN BASE PHYSIQUE : CONTRIBUTION
!    DES DEPLACEMENTS DIFFERENTIELS DES ANCRAGES
!     OPTIONNELLEMENT: ACCELERATION VIA PARALLELISME MPI (SUR DDL)
!-----------------------------------------------------------------------
! IN  : NEQ    : NB D'EQUATIONS DU SYSTEME ASSEMBLE
! IN  : NBEXCI : NOMBRE D'ACCELERO DIFFERENTS
! IN  : PSIDEL : VALEUR DU VECTEUR PSI*DELTA
! IN  : TEMPS  : INSTANT DE CALCUL DES DEPL_IMPO
! IN  : NOMFON : NOM DE LA FONCTION DEPL_IMPO
! OUT : TAB    : VALEUR DE PSIDEL*VALE_NOMFOM(TEMPS)
!
!    OPTIONNELS
! IN    KACCE          <--   PILOTAGE DE L'ACCELERATION: 'CAS0' A 'CAS2'
!                      CAS0: SANS PARALLÉLISME ET SANS BLAS (PAR DEFAUT)
!                      CAS1/CAS2: MPI (DDL) x OPENMP (SUR DDL)
!                            JOINDRE INDICE/TAILLE ET, SI POSSIBLE, INST.
!                            SINON PREPARATION MPI RECALCULEE A CHAQUE PASSAGE
!                      CAS3: INTERDIT CAR RESTE A FAIRE, CF. MDGEPH/MDGEPC
!
! IN    KPROF          <--   PROFILING: 'OUI', 'NON', 'OUIA' (OUI + AFFICHAGE)
!                            'NON' PAR DEFAUT
!
! IN    INST           <--   INSTANT DE CALCUL
!
! INOUT  INDICE        <--   INDICE DE DECALAGE DEPENDANT DU RANG MPI
!                            OBLIGATOIRE SI KACCE DIFFERENT DE CAS0
! INOUT  TAILLE        <--   TAILLE DU BLOC DEPENDANT DU RANG MPI
!                            OBLIGATOIRE SI KACCE DIFFERENT DE CAS0
! .________________.____.______________________________________________.
    implicit none
#include "asterc/asmpi_comm.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/vecini.h"
! PARAMETRES
    integer(kind=8), intent(in) :: neq, nbexci
    real(kind=8), intent(in) :: psidel(neq, nbexci), temps
    character(len=8), intent(in) :: nomfon(2*nbexci)
    real(kind=8), intent(out) :: tab(neq)
    character(len=4), intent(in), optional :: kacce, kprof
    integer(kind=8), intent(in), optional :: inst
    integer(kind=8), intent(inout), optional :: indice, taille
! VARIABLES LOCALES
    integer(kind=8) :: ieq, ier, iex, cntr, ietdeb, ietrat, ietmax, ietfin
    integer(kind=8) :: rang, nbproc, nbloc, nrest, ifm, niv, instl, borne1, borne2
    integer(kind=8) :: m, iaux
    real(kind=8) :: coef, tempsl, tempst, zero
    character(len=4) :: kaccl, kprol
    character(len=8) :: nompar, k8bid
    mpi_int :: mpicou, mpicow, mrang, mnbproc
    data tempst/0.d0/
    save tempst
!-----------------------------------------------------------------------
    call jemarq()
! Init. globales
    call infniv(ifm, niv)
    k8bid = '        '
    nompar = 'INST'
    cntr = 0
    zero = 0.d0
!
    if (.not. present(kacce)) then
        kaccl = 'CAS0'
    else
        ASSERT((kacce .eq. 'CAS0') .or. (kacce .eq. 'CAS1') .or. (kacce .eq. 'CAS2'))
        if (kacce .eq. 'CAS0') then
            kaccl = 'CAS0'
        else
            kaccl = 'CAS1'
        end if
    end if
    if (.not. present(kprof)) then
        kprol = 'NON'
    else
        ASSERT((kprof(1:3) .eq. 'OUI') .or. (kprof(1:3) .eq. 'NON') .or. (kprof(1:4) .eq. 'OUIA'))
        kprol = kprof
    end if
    if (present(inst)) then
        instl = inst
    else
        instl = 1
    end if
! Affichage/profiling
    if (kprol(1:3) .eq. 'OUI') call system_clock(ietdeb, ietrat, ietmax)
!
! Init. par option
    ASSERT((nbexci .ge. 1) .and. (neq .ge. 1))
    rang = -999
    nbproc = -999
    nbloc = -999
    m = neq
    iaux = 1
    if (kaccl .eq. 'CAS1') then
        ASSERT(present(indice) .and. present(taille))
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
! Parallelisme MPI sur accelero
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
            iaux = indice
            m = taille
        end if
        ASSERT((m .ge. 1) .and. (iaux .ge. 1))
    end if
!
!---------------------------------------------------
! Calcul proprement dit, différent suivant l'option
!---------------------------------------------------
!
!    write(ifm,*)'<mdhep3> kacce/m/iaux/neq/nbexci/instl=', &
!                     kaccl,m,iaux,neq,nbexci,instl

    call vecini(neq, zero, tab)
    if (kaccl(1:4) .eq. 'CAS0') then
! Sans parallélisme et sans BLAS
        do iex = 1, nbexci
            if (nomfon(iex) .eq. k8bid) goto 10
            cntr = cntr+1
            call fointe('F ', nomfon(iex), 1, [nompar], [temps], coef, ier)
            do ieq = 1, neq
                tab(ieq) = tab(ieq)+psidel(ieq, iex)*coef
            end do
10          continue
        end do
!
    else if (kaccl(1:4) .eq. 'CAS1') then
! Avec MPI (ddl) et sans OpenMP
        do iex = 1, nbexci
            if (nomfon(iex) .eq. k8bid) goto 15
            cntr = cntr+1
            call fointe('F ', nomfon(iex), 1, [nompar], [temps], coef, ier)
            borne1 = iaux
            borne2 = iaux+m-1
            do ieq = borne1, borne2
                tab(ieq) = tab(ieq)+psidel(ieq, iex)*coef
            end do
15          continue
        end do
! COM MPI que si au moins 1 ddl par processus
        if (m .ne. neq) call asmpi_comm_vect('MPI_SUM', 'R', nbval=neq, vr=tab)
    else
! Mauvaise option
        ASSERT(.false.)
!
    end if
!
    if (cntr .eq. 0) call utmess('A', 'ALGORITH13_44')
!
! Affichage/profiling
    if (kprol(1:3) .eq. 'OUI') then
        call system_clock(ietfin)
        tempsl = real(ietfin-ietdeb)/real(ietrat)
        tempst = tempst+tempsl
        if (kprol(1:4) .eq. 'OUIA') then
            write (ifm, *) '<mdgep3> kacce/iaux/m/instl=', kaccl, iaux, m, instl
            write (ifm, *) '<mdgep3> temps_local/temps_total=', tempsl, tempst
        end if
    end if
    call jedema()
end subroutine
