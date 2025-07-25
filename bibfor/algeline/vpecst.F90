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

subroutine vpecst(ifm, typres, omgmin, omgmax, nbfre1, &
                  nbfre2, nbfreq, nblagr, typep, typcon, &
                  dimc1, zimc1)
    implicit none
#include "asterfort/assert.h"
#include "asterfort/freqom.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ifm, nbfre1, nbfre2, nbfreq, nblagr
    real(kind=8) :: omgmin, omgmax, dimc1
    complex(kind=8) :: zimc1
    character(len=1) :: typep
    character(len=8) :: typcon
    character(len=16) :: typres
!   PRINTING OF THE NUMBER OF THE EIGENVALUES THAT BELONG TO A CHOOSEN
!   PATTERN
!     ------------------------------------------------------------------
! IN  OMGMIN : R8 : PULSATION MIN
! IN  OMGMAX : R8 : PULSATION MAX
! IN  NBFRE1 : IS : NB DE PULSATION INFERIEURE A OMGMIN
! IN  NBFRE2 : IS : NB DE PULSATION INFERIEURE A OMGMAX
! IN/OUT NBFREQ : IS : OUTPUT: NBRE DE FREQ DANS LA BANDE OMGMIN OMGMAX
!                    INPUT: NUMERO DE LA BANDE (SEULEMENT SI TYPEP='D')
! IN  NBLAGR : IS : NB DE DDLS DE LAGRANGE
! IN  TYPRES : TX : TYPE DE CALCUL (DYNAMIQUE OU FLAMBEMENT)
! IN  TYPEP  : K1 : TYPE D EIGENVALUE-PROBLEM
!    'R' GEP REEL POUR MODE_ITER_SIMULT/BANDE OU INFO_MODE SUR UNE
!        SEULE BANDE OU MODE_ITER_INV
!    'F' IDEM 'R' MAIS SANS CALCUL NBFREQ, IL EST DONNE EN INPUT.
!        UTILISE POUR OPTION 'SIMULT/'BANDE' AVEC TABLE ISSUE
!        D'INFO_MODE
!    'D' GEP REEL POUR INFO_MODE SUR LA BANDE NUMERO NFREQ
!    'C' GEP COMPLEXE OU QEP POUR INFO_MODE
!    'S' GEP REEL POUR TEST DE STURM POSTTRAITEMENT MODE_ITER_SIMULT
! IN  TYPCON : K8 : TYPE DE CONTOUR (LICITE SI TYPEP='C')
! IN  DIMC1  : R8 : DIMENSION CHARACTERISTIQUE REEL N 1 DU CONTOUR
!                   (LICITE SI TYPEP='C')
! IN  ZIMC1  : C16: DIMENSION CHARACTERISTIQUE COMPLEXE N 1 DU CONTOUR
!                   (LICITE SI TYPEP='C')
!     ------------------------------------------------------------------
!     REMARQUE:  NBFRE1 ET NBFRE2  SONT CALCULES PAR VPSTUR
!     ------------------------------------------------------------------
    real(kind=8) :: fmin, fmax, valr(3)
    integer(kind=8) :: ibande, vali(2)
!     ------------------------------------------------------------------
!
!   --- ON RECUPERE LE NUMERO DE LA BANDE SI NECESSAIRE
    if (typep .eq. 'D') ibande = nbfreq
!
!   --- WE ARE ON THE REEL PLANE (STURM TEST)
    if ((typep .eq. 'R') .or. (typep .eq. 'S') .or. (typep .eq. 'D') .or. (typep .eq. 'F')) then
!
!   --- TEST DE STURM
        if (typep .ne. 'F') then
            if (typres .eq. 'DYNAMIQUE') then
                nbfreq = abs(nbfre2-nbfre1)
            else
                if ((omgmin*omgmax) .ge. 0.d0) then
                    nbfreq = abs(nbfre2-nbfre1)
                else
                    nbfreq = abs(nbfre2+nbfre1-2*nblagr)
                end if
            end if
        end if
        if (nbfreq .gt. 9999) then
            call utmess('A', 'ALGELINE3_64', si=nbfreq)
        end if
!
!   --- AFFICHAGE SI MODE_ITER_SIMULT+BANDE OU MODE_ITER_INV+
!   --- SEPARE/AJUSTE OU INFO_MODE+MIN/MAX
        if ((typep .eq. 'R') .or. (typep .eq. 'F')) then
            if (typres .eq. 'DYNAMIQUE') then
                fmin = freqom(omgmin)
                fmax = freqom(omgmax)
                write (ifm, 950)
                if (typep .eq. 'R') then
                    call utmess('I', 'ALGELINE6_33')
                else if (typep .eq. 'F') then
                    call utmess('I', 'ALGELINE6_34')
                end if
                if (nbfreq .eq. 0) then
                    valr(1) = fmin
                    valr(2) = fmax
                    vali(1) = 1
                    call utmess('I', 'ALGELINE6_35', si=vali(1), nr=2, valr=valr)
                else if (fmin .eq. 0.d0) then
                    valr(1) = fmax
                    vali(1) = nbfreq
                    call utmess('I', 'ALGELINE6_36', si=vali(1), sr=valr(1))
                else
                    valr(1) = fmin
                    valr(2) = fmax
                    vali(1) = 1
                    vali(2) = nbfreq
                    call utmess('I', 'ALGELINE6_37', ni=2, vali=vali, nr=2, &
                                valr=valr)
                end if
            else
                if (typep .eq. 'R') then
                    call utmess('I', 'ALGELINE6_28')
                else if (typep .eq. 'F') then
                    call utmess('I', 'ALGELINE6_29')
                end if
                if (nbfreq .eq. 0) then
                    valr(1) = -omgmax
                    valr(2) = -omgmin
                    vali(1) = 1
                    call utmess('I', 'ALGELINE6_30', si=vali(1), nr=2, valr=valr)
                else if (omgmin .eq. 0.d0) then
                    valr(1) = -omgmax
                    vali(1) = nbfreq
                    call utmess('I', 'ALGELINE6_31', si=vali(1), sr=valr(1))
                else
                    valr(1) = -omgmax
                    valr(2) = -omgmin
                    vali(1) = 1
                    vali(2) = nbfreq
                    call utmess('I', 'ALGELINE6_32', ni=2, vali=vali, nr=2, &
                                valr=valr)
                end if
            end if
!
!   --- AFFICHAGE SI INFO_MODE+LISTE
        else if (typep .eq. 'D') then
            if (typres .eq. 'DYNAMIQUE') then
                fmin = freqom(omgmin)
                fmax = freqom(omgmax)
                if (ibande .eq. 1) then
                    write (ifm, 950)
                    call utmess('I', 'ALGELINE6_33')
                end if
                valr(1) = fmin
                valr(2) = fmax
                vali(1) = ibande
                if (nbfreq .eq. 0) then
                    call utmess('I', 'ALGELINE6_35', si=vali(1), nr=2, valr=valr)
                else
                    vali(2) = nbfreq
                    call utmess('I', 'ALGELINE6_37', ni=2, vali=vali, nr=2, &
                                valr=valr)
                end if
            else
                if (ibande .eq. 1) then
                    write (ifm, 950)
                    call utmess('I', 'ALGELINE6_28')
                end if
                valr(1) = -omgmax
                valr(2) = -omgmin
                vali(1) = ibande
                if (nbfreq .eq. 0) then
                    call utmess('I', 'ALGELINE6_30', si=vali(1), nr=2, valr=valr)
                else
                    vali(2) = nbfreq
                    call utmess('I', 'ALGELINE6_32', ni=2, vali=vali, nr=2, &
                                valr=valr)
                end if
            end if
!
        else if (typep .eq. 'S') then
!   --- AFFICHAGE DEDIE DS UNE DES ROUTINES APPELLANTES
        end if
!
!   --- WE ARE ON THE COMPLEX PLANE (APM TEST)
    else if (typep .eq. 'C') then
        nbfreq = nbfre2
        if (nbfreq .gt. 9999) then
            call utmess('A', 'ALGELINE3_64')
            call utmess('I', 'ALGELINE7_19', si=nbfreq)
        end if
        write (ifm, 950)
        call utmess('I', 'ALGELINE6_38')
        if (nbfreq .eq. 0) then
            if (typcon(1:6) .eq. 'CERCLE') then
                valr(1) = dble(zimc1)
                valr(2) = dimag(zimc1)
                valr(3) = dimc1
                call utmess('I', 'ALGELINE6_39', nr=3, valr=valr)
            end if
        else
            if (typcon(1:6) .eq. 'CERCLE') then
                valr(1) = dble(zimc1)
                valr(2) = dimag(zimc1)
                valr(3) = dimc1
                vali(1) = nbfreq
                call utmess('I', 'ALGELINE6_40', si=vali(1), nr=3, valr=valr)
            end if
        end if
!
!   --- ILLEGAL OPTION
    else
        ASSERT(.false.)
    end if
    if (typep .ne. 'S') write (ifm, 950)
!
950 format(72('-'),/)
end subroutine
