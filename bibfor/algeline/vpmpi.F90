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

subroutine vpmpi(option, eigsol, icom1_, icom2_, lcomod_, &
                 mpicou_, mpicow_, nbvecg_, nfreqg_, rangl_, &
                 omemax_, omemin_, vpinf_, vpmax_)
!
! ROUTINE ORGANISANT LE PARALLELISME MULTI-NIVEAUX DANS MODE_ITER_SIMULT (+ APPELS DS VPPARA).
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "asterc/asmpi_comm.h"
#include "asterc/asmpi_split_comm.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/utmess.h"
#include "asterfort/vpleci.h"
!
!
! --- INPUT
!
    integer(kind=8), intent(in) :: option
    character(len=19), optional, intent(in) :: eigsol
!
! --- OUTPUT
!
    integer(kind=8), optional, intent(out) :: icom1_, icom2_, nbvecg_, nfreqg_
!
! --- INPUT/OUTPUT
!
    mpi_int, optional, intent(inout) :: mpicou_, mpicow_
    integer(kind=8), optional, intent(inout) :: rangl_
    aster_logical, optional, intent(inout) :: lcomod_
    real(kind=8), optional, intent(inout) :: omemax_, omemin_, vpinf_, vpmax_
!
! --- VARIABLES LOCALES
!
    mpi_int :: mrang, mnbproc
    integer(kind=8) :: l1, l2, l3, nbvect, nbproc, nfreq, rang, typeco, vali(5)
    real(kind=8) :: rbid
    character(len=24) :: k24bid, valk(5)

    integer(kind=8) :: icom1, icom2, nbvecg, nfreqg
    mpi_int :: mpicou, mpicow
    integer(kind=8) :: rangl
    aster_logical :: lcomod
    real(kind=8) :: omemax, omemin, vpinf, vpmax
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!

!   traitement des arguments optionnels en entrée
    if (option .ne. 2) then
        ASSERT(.not. present(eigsol))
    end if

    if (present(mpicou_)) mpicou = mpicou_
    if (present(mpicow_)) mpicow = mpicow_
    if (present(rangl_)) rangl = rangl_
    if (present(lcomod_)) lcomod = lcomod_
    if (present(omemax_)) omemax = omemax_
    if (present(omemin_)) omemin = omemin_
    if (present(vpinf_)) vpinf = vpinf_
    if (present(vpmax_)) vpmax = vpmax_

    select case (option)
    case (1)
! ---  STEP 1: RECUPERATION DES PARAMETRES MPI + TESTS
! --- INPUT: option
! --- OUTPUT: mpicou, mpicow, icom1, icom2, rangl, lcomod
        icom1 = -999
        icom2 = -999
        call asmpi_comm('GET_WORLD', mpicow)
        call asmpi_comm('GET', mpicou)
! --  ON EST CENSE FONCTIONNER EN COMM_WORLD
        if (mpicow .ne. mpicou) then
            ASSERT(.false.)
        end if
        call asmpi_info(mpicow, mrang, mnbproc)
        rang = to_aster_int(mrang)
        nbproc = to_aster_int(mnbproc)
!
        call getvis('PARALLELISME_MACRO', 'TYPE_COM', iocc=1, scal=typeco, nbret=l1)
        call getvis('PARALLELISME_MACRO', 'IPARA1_COM', iocc=1, scal=icom1, nbret=l2)
        call getvis('PARALLELISME_MACRO', 'IPARA2_COM', iocc=1, scal=icom2, nbret=l3)
        valk(1) = 'TYPE_COM'
        valk(2) = 'IPARA1_COM'
        valk(3) = 'IPARA2_COM'
        valk(4) = 'RANG'
        valk(5) = 'NBPROC'
        vali(1) = typeco
        vali(2) = icom1
        vali(3) = icom2
        vali(4) = rang
        vali(5) = nbproc
        if (l1*l2*l3 .ne. 1) call utmess('F', 'APPELMPI_6', nk=3, valk=valk, ni=3, &
                                         vali=vali)
!
        if (( &
            ((typeco .ne. 1) .and. (typeco .ne. -999)) .or. &
            ((icom1 .ne. -999) .and. ((icom1 .lt. 1) .or. (icom1 .gt. nbproc))) .or. &
       ((icom2 .ne. -999) .and. ((icom2 .lt. 1) .or. (icom2 .gt. nbproc))) .or. (icom1 .gt. icom2) &
            .or. (nbproc .lt. 1) .or. (rang .lt. 0) &
            )) call utmess('F', 'APPELMPI_8', nk=5, valk=valk, ni=5, &
                           vali=vali)
!
!
        if ((typeco .eq. 1) .and. (nbproc .gt. 1)) then
            lcomod = .true.
! --  DECOMPOSE LE COM GLOBAL MPICOW EN COM LOCAL MPICOU
! --  PLUS AFFECTATION DE CE NOUVEAU COM AFIN DE NE PAS PERTURBER LA FACTO DE LA DEMI-BANDE
            call asmpi_split_comm(mpicow, to_mpi_int(icom1), to_mpi_int(0), 'ipara1', mpicou)
            if (mpicow .eq. mpicou) then
                ASSERT(.false.)
            end if
            call asmpi_barrier()
            call asmpi_comm('SET', mpicou)
! --  RANG DANS LE SOUS-COMM MPICOU LIE A CHAQUE OCCURENCE MUMPS: RANGL
            call asmpi_info(comm=mpicou, rank=mrang)
            rangl = to_aster_int(mrang)
        else
            rangl = -9999
            mpicou = -9999
            lcomod = .false.
        end if
!
    case (2)
! --- STEP 2: REDIMENSIONNEMENT DES BUFFERS DE COMMUNICATION
! --- INPUT: option, eigsol, lcomod, mpicou, mpicow, rangl
! --- OUTPUT: nbvecg, nfreqg
        nbvecg = -9999
        nfreqg = -9999
        if (lcomod) then
            ASSERT(present(eigsol))
            call vpleci(eigsol, 'I', 1, k24bid, rbid, &
                        nfreq)
            call vpleci(eigsol, 'I', 2, k24bid, rbid, &
                        nbvect)
! --  ON REMET LE COM WORLD POUR COMMUNIQUER NBVECT/NBFREQ
            call asmpi_comm('SET', mpicow)
            call asmpi_barrier()
! --  EST-ON LE PROCESSUS MAITRE DU COM LOCAL: RANGL=0 ?
! --  SI OUI, ON ENVOI LES BONNES VALEURS DE NBVECT/NFREQ SUR LE COM GLOBAL MPICOW,
! --  SINON ON RENVOI ZERO POUR NE PAS COMPTER PLUSIEURS FOIS L'INFO.
            if (rangl .eq. 0) then
                nbvecg = nbvect
                nfreqg = nfreq
            else
                nbvecg = 0
                nfreqg = 0
            end if
            call asmpi_comm_vect('MPI_SUM', 'I', sci=nbvecg)
            call asmpi_comm_vect('MPI_SUM', 'I', sci=nfreqg)
! --  ON REMET LE COM LOCAL POUR LES FACTO ET SOLVES A SUIVRE
            call asmpi_barrier()
            call asmpi_comm('SET', mpicou)
        end if
!
    case (3)
! --- STEP 3: POUR MEMOIRE, ETAPE GEREE DANS VPPARA
        ASSERT(.false.)
!
    case (4)
! --- STEP 4:
! --- EN CAS DE TEST DE STURM LOCAL A CHAQUE SOUS-BANDE, REMISE A JOUR DES BORNES VIA LE COM WORLD.
! --- PUIS ON REMET LE COMCOU POUR NE PAS GENER LES FACTOS EVENTUELLES DE VPCNTL.
! --- INPUT: lcomod, mpicow, mpicou
! --- INPUT/OUTPUT: omemin, omemax, vpinf, vpmax
        if (lcomod) then
            call asmpi_comm('SET', mpicow)
            call asmpi_barrier()
            call asmpi_comm_vect('MPI_MIN', 'R', scr=omemin)
            call asmpi_comm_vect('MPI_MIN', 'R', scr=vpinf)
            call asmpi_comm_vect('MPI_MAX', 'R', scr=omemax)
            call asmpi_comm_vect('MPI_MAX', 'R', scr=vpmax)
            call asmpi_barrier()
            call asmpi_comm('SET', mpicou)
        end if
!
    case (5)
! --- STEP 5:
! --- AVANT DE QUITTER L'OP. ON REMET LE COM WORLD (AU CAS OU).
! --- DESTRUCTION DES SOUS-COMMUNICATEURS EVENTUELLEMENT ASSOCIES A UNE OCCURENCE MUMPS
! --- UNE OCCURENCE MUMPS (APRES CELLE DE LADITE OCCURENCE).
! --- INPUT: lcomod, mpicow, mpicou
        if (lcomod) then
            call asmpi_comm('SET', mpicow)
            call asmpi_barrier()
            call asmpi_comm('FREE', mpicou)
        end if
!
    case default
        ASSERT(.false.)
    end select

!   traitement des arguments optionnels en sortie

    if (present(mpicou_)) mpicou_ = mpicou
    if (present(mpicow_)) mpicow_ = mpicow
    if (present(rangl_)) rangl_ = rangl
    if (present(lcomod_)) lcomod_ = lcomod
    if (present(omemax_)) omemax_ = omemax
    if (present(omemin_)) omemin_ = omemin
    if (present(vpinf_)) vpinf_ = vpinf
    if (present(vpmax_)) vpmax_ = vpmax

    if (present(icom1_)) icom1_ = icom1
    if (present(icom2_)) icom2_ = icom2
    if (present(nbvecg_)) nbvecg_ = nbvecg
    if (present(nfreqg_)) nfreqg_ = nfreqg
!
!     FIN DE VPMPI
!
end subroutine
