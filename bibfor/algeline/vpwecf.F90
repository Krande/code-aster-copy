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

subroutine vpwecf(option, typres, nfreq, mxfreq, resufi, &
                  resufr, resufk, lamor, ktyp, lns)
!     ECRITURE DES FREQUENCES RELATIVEMENT A LA METHODE UTILISEE
!     IMPRESSION D'OFFICE SUR "MESSAGE"
!-----------------------------------------------------------------------
    implicit none
!
! PARAMETRES D'APPEL
#include "asterf_types.h"
#include "asterc/isnnem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nfreq, mxfreq, resufi(mxfreq, *), lamor
    real(kind=8) :: resufr(mxfreq, *)
    character(len=*) :: option, resufk(mxfreq, *), typres
    character(len=1) :: ktyp
    aster_logical :: lns
!
! VARIABLES LOCALES
    integer(kind=8) :: ifm, ifreq, indf, niv, vali(4), ideb, ifin, ipas, counter
    real(kind=8) :: am, undf, erc, errmoy
    real(kind=8) :: valr(5)
    character(len=27) :: valk(4)
!     ------------------------------------------------------------------
    call infniv(ifm, niv)
    undf = r8vide()
    indf = isnnem()
    errmoy = 0.d0
    if (nfreq .eq. 0) then
        ASSERT(.false.)
    end if

    if (typres .eq. 'DYNAMIQUE') then
        ideb = 1
        ifin = nfreq
        ipas = 1
    else
        ideb = nfreq
        ifin = 1
        ipas = -1
    end if
    counter = 0
    if (resufk(nfreq, 2) .eq. 'BATHE_WILSON') then
        if (typres .eq. 'DYNAMIQUE') then
            call utmess('I', 'ALGELINE6_59')
        else
            call utmess('I', 'ALGELINE6_60')
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            am = resufr(ifreq, 4)
            errmoy = errmoy+abs(am)
            valr(1) = am
            vali(2) = resufi(ifreq, 3)
            vali(3) = resufi(ifreq, 5)
            if (typres .eq. 'DYNAMIQUE') then
                valr(2) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(2) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            call utmess('I', 'ALGELINE6_61', ni=3, vali=vali, nr=2, &
                        valr=valr)
        end do
        valr(1) = errmoy/nfreq
        call utmess('I', 'ALGELINE6_58', sr=valr(1))
!
    else if (resufk(nfreq, 2) .eq. 'LANCZOS') then
        if (lamor .eq. 0) then
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_62')
            else
                call utmess('I', 'ALGELINE6_63')
            end if
        else
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_64')
            else
                call utmess('I', 'ALGELINE6_65')
            end if
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            vali(2) = resufi(ifreq, 2)
            if (lamor .eq. 0) then
                am = resufr(ifreq, 4)
            else
                am = resufr(ifreq, 3)
            end if
            errmoy = errmoy+abs(am)
            valr(1) = am
            if (typres .eq. 'DYNAMIQUE') then
                valr(2) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(2) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            call utmess('I', 'ALGELINE6_66', ni=2, vali=vali, nr=2, &
                        valr=valr)
        end do
        if (lamor .eq. 0) then
            valr(1) = errmoy/nfreq
            call utmess('I', 'ALGELINE6_58', sr=valr(1))
        end if
!
    else if (resufk(nfreq, 2) .eq. 'SORENSEN') then
        if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_67')
            else
                call utmess('I', 'ALGELINE6_68')
            end if
        else
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_70')
            else
                call utmess('I', 'ALGELINE6_71')
            end if
        end if

        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
                am = resufr(ifreq, 4)
                errmoy = errmoy+abs(am)
            else
                am = resufr(ifreq, 3)
                erc = resufr(ifreq, 4)
                errmoy = errmoy+abs(erc)
            end if
            valr(1) = am
            if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
                if (typres .eq. 'DYNAMIQUE') then
                    valr(2) = resufr(ifreq, 1)
                    vali(1) = resufi(ifreq, 1)
                else
                    valr(2) = -resufr(ifreq, 2)
                    vali(1) = counter
                end if
                call utmess('I', 'ALGELINE6_69', si=vali(1), nr=2, valr=valr)
            else
                if (typres .eq. 'DYNAMIQUE') then
                    valr(2) = resufr(ifreq, 1)
                    vali(1) = resufi(ifreq, 1)
                else
                    valr(2) = -resufr(ifreq, 2)
                    vali(1) = counter
                end if
                valr(3) = erc
                call utmess('I', 'ALGELINE6_72', si=vali(1), nr=3, valr=valr)
            end if
        end do
        valr(1) = errmoy/nfreq
        call utmess('I', 'ALGELINE6_58', sr=valr(1))
!
    else if (resufk(nfreq, 2) (1:2) .eq. 'QZ') then
        valk(1) = resufk(nfreq, 2) (1:16)
        if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_73', sk=valk(1))
            else
                call utmess('I', 'ALGELINE6_74', sk=valk(1))
            end if
        else
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_75', sk=valk(1))
            else
                call utmess('I', 'ALGELINE6_76', sk=valk(1))
            end if
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
                am = resufr(ifreq, 4)
                errmoy = errmoy+abs(am)
            else
                am = resufr(ifreq, 3)
                erc = resufr(ifreq, 4)
                errmoy = errmoy+abs(erc)
            end if
            valr(1) = am
            if ((lamor .eq. 0) .and. (ktyp .eq. 'R') .and. (.not. lns)) then
                if (typres .eq. 'DYNAMIQUE') then
                    valr(2) = resufr(ifreq, 1)
                    vali(1) = resufi(ifreq, 1)
                else
                    valr(2) = -resufr(ifreq, 2)
                    vali(1) = counter
                end if
                call utmess('I', 'ALGELINE6_69', si=vali(1), nr=2, valr=valr)
            else
                if (typres .eq. 'DYNAMIQUE') then
                    valr(2) = resufr(ifreq, 1)
                    vali(1) = resufi(ifreq, 1)
                else
                    valr(2) = -resufr(ifreq, 2)
                    vali(1) = counter
                end if
                valr(3) = erc
                call utmess('I', 'ALGELINE6_72', si=vali(1), nr=3, valr=valr)
            end if
        end do
        valr(1) = errmoy/nfreq
        call utmess('I', 'ALGELINE6_58', sr=valr(1))
!
    elseif ((resufk(nfreq, 2) .eq. 'INVERSE_R' .or. resufk(nfreq, 2) &
             .eq. 'INVERSE_C') .and. (option(1:6) .eq. 'PROCHE')) then
        if (typres .eq. 'DYNAMIQUE') then
            call utmess('I', 'ALGELINE6_77')
        else
            call utmess('I', 'ALGELINE6_78')
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if (typres .eq. 'DYNAMIQUE') then
                valr(1) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(1) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            valr(2) = resufr(ifreq, 3)
            vali(2) = resufi(ifreq, 4)
            valr(3) = resufr(ifreq, 15)
            valr(4) = resufr(ifreq, 4)
!
            call utmess('I', 'ALGELINE6_79', ni=2, vali=vali, nr=4, &
                        valr=valr)
            resufr(ifreq, 14) = undf
            resufr(ifreq, 15) = undf
            resufi(ifreq, 2) = indf
            resufi(ifreq, 3) = indf
            resufi(ifreq, 4) = indf
            resufi(ifreq, 8) = resufi(ifreq, 4)
        end do
        write (ifm, 7777)
!
    elseif (resufk(nfreq, 2) .eq. 'INVERSE_R' .and. option(1:6) .eq. &
            'AJUSTE') then
        if (typres .eq. 'DYNAMIQUE') then
            call utmess('I', 'ALGELINE6_80')
        else
            call utmess('I', 'ALGELINE6_81')
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if (typres .eq. 'DYNAMIQUE') then
                valr(1) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(1) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            valr(2) = resufr(ifreq, 3)
            vali(2) = resufi(ifreq, 2)
            vali(3) = resufi(ifreq, 3)
            valr(3) = resufr(ifreq, 14)
            vali(4) = resufi(ifreq, 4)
            valr(4) = resufr(ifreq, 15)
            valr(5) = resufr(ifreq, 4)
            call utmess('I', 'ALGELINE6_82', ni=4, vali=vali, nr=5, &
                        valr=valr)
!
            resufr(ifreq, 14) = undf
            resufr(ifreq, 15) = undf
            resufi(ifreq, 2) = indf
            resufi(ifreq, 3) = indf
            resufi(ifreq, 4) = indf
            resufi(ifreq, 7) = resufi(ifreq, 4)
        end do
        write (ifm, 7777)
!
    elseif (resufk(nfreq, 2) .eq. 'INVERSE_R' .and. option(1:6) .eq. &
            'SEPARE') then
        if (typres .eq. 'DYNAMIQUE') then
            call utmess('I', 'ALGELINE6_83')
        else
            call utmess('I', 'ALGELINE6_84')
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if (typres .eq. 'DYNAMIQUE') then
                valr(1) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(1) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            vali(1) = resufi(ifreq, 1)
            valr(2) = resufr(ifreq, 3)
            vali(2) = resufi(ifreq, 2)
            vali(3) = resufi(ifreq, 4)
            valr(3) = resufr(ifreq, 15)
            valr(4) = resufr(ifreq, 4)
            call utmess('I', 'ALGELINE6_85', ni=3, vali=vali, nr=4, &
                        valr=valr)
!
            resufr(ifreq, 14) = undf
            resufr(ifreq, 15) = undf
            resufi(ifreq, 2) = indf
            resufi(ifreq, 3) = indf
            resufi(ifreq, 4) = indf
            resufi(ifreq, 6) = resufi(ifreq, 4)
        end do
        write (ifm, 7777)
!
    elseif (resufk(nfreq, 2) .eq. 'INVERSE_C' &
            .and. (option(1:6) .eq. 'AJUSTE' .or. option(1:6) .eq. 'SEPARE')) then
        if (typres .eq. 'DYNAMIQUE') then
            call utmess('I', 'ALGELINE6_86')
        else
            call utmess('I', 'ALGELINE6_87')
        end if
        do ifreq = ideb, ifin, ipas
            counter = counter+1
            if (typres .eq. 'DYNAMIQUE') then
                valr(1) = resufr(ifreq, 1)
                vali(1) = resufi(ifreq, 1)
            else
                valr(1) = -resufr(ifreq, 2)
                vali(1) = counter
            end if
            vali(1) = resufi(ifreq, 1)
            valr(2) = resufr(ifreq, 3)
            vali(2) = resufi(ifreq, 2)
            valr(3) = resufr(ifreq, 14)
            vali(3) = resufi(ifreq, 4)
            valr(4) = resufr(ifreq, 15)
            valr(5) = resufr(ifreq, 4)
            call utmess('I', 'ALGELINE6_88', ni=3, vali=vali, nr=5, &
                        valr=valr)
!
            resufr(ifreq, 14) = undf
            resufr(ifreq, 15) = undf
            resufi(ifreq, 2) = indf
            resufi(ifreq, 3) = indf
            resufi(ifreq, 4) = indf
            resufi(ifreq, 8) = resufi(ifreq, 4)
        end do
        write (ifm, 7777)
!
    end if
!
7777 format(/)
!
!     ------------------------------------------------------------------
end subroutine
