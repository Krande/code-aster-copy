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
subroutine affbar(tmp, tmpf, fcx, nommai, isec, &
                  car, val, exp, ncar, kioc, &
                  ier)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8maem.h"
#include "asterc/r8pi.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: isec, ncar, ier
    real(kind=8) :: val(*)
    character(len=6) :: kioc
    character(len=8) :: fcx, nommai, car(*), exp(*)
    character(len=24) :: tmp, tmpf
!     VERIFICATION DE LA BONNE AFFECTATION DES SECTIONS DE BARRE :
!        - INTERDICTION D ECRASER UNE GEOMETRIE DE SECTION PAR UNE AUTRE
!     AFFECTATION DES CARACTERISTIQUES GENERALES ET GEOMETRIQUES
!                         AUX MAILLES DE TYPE BARRE DANS L OBJET TAMPON
!
!     L OBJET TAMPON (TMP) CONTIENT (8*NBBARRE) VALEURS
!     EXP  : A    HY   HZ   EPY  EPZ  EP   R   TSEC
!     TSEC : TYPE  GEOMETRIQUE DE SECTION : 0 = GENERALE
!                                           1 = RECTANGLE
!                                           2 = CERCLE
!-----------------------------------------------------------------------
    real(kind=8) :: tst, pi, zero
    real(kind=8) :: hy, hz, epy, epz, hyi, hzi, e, re, ri
    character(len=8) :: resu
    character(len=24) :: valk(2)
    character(len=16) :: concep, cmd
    aster_logical :: secple
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iisec, j, jdge, jdgef, num
!-----------------------------------------------------------------------
    call jemarq()
    call getres(resu, concep, cmd)
    tst = r8maem()
    pi = r8pi()
    zero = 0.d0
    secple = .true.
!
    call jenonu(jexnom(tmp, nommai), num)
!
! --- TESTS D ECRASEMENT DE SECTION
    if (num .ne. 0) then
        call jeveuo(jexnom(tmp, nommai), 'E', jdge)
        iisec = nint(zr(jdge+ncar-1))
        if (iisec .ne. isec) then
            valk(1) = kioc
            valk(2) = nommai
            call utmess('A', 'MODELISA_69', nk=2, valk=valk)
            ier = ier+1
            goto 999
        end if
    else
        call jecroc(jexnom(tmp, nommai))
        call jeveuo(jexnom(tmp, nommai), 'E', jdge)
        do i = 1, ncar
            zr(jdge+i-1) = tst
        end do
    end if
!
!     --- NOM DE LA FONCTION DU CX
    call jenonu(jexnom(tmpf, nommai), num)
    if (num .eq. 0) then
        call jecroc(jexnom(tmpf, nommai))
    end if
    call jeveuo(jexnom(tmpf, nommai), 'E', jdgef)
    zk8(jdgef) = fcx
!
    if (isec .eq. 0) then
        do j = 1, ncar
            if (car(j) .eq. 'A       ') zr(jdge) = val(j)
        end do
    else if (isec .eq. 1) then
        do j = 1, ncar
            if (car(j) .eq. 'HY      ') then
                zr(jdge+1) = val(j)
            else if (car(j) .eq. 'HZ      ') then
                zr(jdge+2) = val(j)
            else if (car(j) .eq. 'EPY     ') then
                zr(jdge+3) = val(j)
                secple = .false.
            else if (car(j) .eq. 'EPZ     ') then
                zr(jdge+4) = val(j)
                secple = .false.
            else if (car(j) .eq. 'H       ') then
                zr(jdge+1) = val(j)
                zr(jdge+2) = val(j)
            else if (car(j) .eq. 'EP      ') then
                zr(jdge+3) = val(j)
                zr(jdge+4) = val(j)
                secple = .false.
            end if
        end do
    else if (isec .eq. 2) then
        do j = 1, ncar
            if (car(j) .eq. 'R       ') then
                zr(jdge+5) = val(j)
            else if (car(j) .eq. 'EP      ') then
                zr(jdge+6) = val(j)
                secple = .false.
            end if
        end do
    end if
    zr(jdge+7) = isec
!
! --- COMPLETUDE DES DONNES GENERALES
    if (isec .eq. 0) then
        if (zr(jdge) .eq. tst) then
            valk(1) = nommai
            valk(2) = exp(1)
            call utmess('A', 'MODELISA_70', nk=2, valk=valk)
            ier = ier+1
        end if
        if (zr(jdge) .le. zero) then
            valk(1) = nommai
            valk(2) = exp(1)
            call utmess('A', 'MODELISA_71', nk=2, valk=valk)
            ier = ier+1
        end if
!
! --- COMPLETUDE DES DONNES GEOMETRIQUES RECTANGLE
    else if (isec .eq. 1) then
        do j = 1, 2
            if (zr(jdge+j) .eq. tst) then
                valk(1) = nommai
                valk(2) = exp(1+j)
                call utmess('A', 'MODELISA_72', nk=2, valk=valk)
                ier = ier+1
            end if
            if (zr(jdge+j) .le. zero) then
                valk(1) = nommai
                valk(2) = exp(1+j)
                call utmess('A', 'MODELISA_73', nk=2, valk=valk)
                ier = ier+1
            end if
        end do
        if (.not. secple) then
            do j = 3, 4
                if (zr(jdge+j) .eq. tst) then
                    valk(1) = nommai
                    valk(2) = exp(1+j)
                    call utmess('A', 'MODELISA_72', nk=2, valk=valk)
                    ier = ier+1
                end if
                if (zr(jdge+j) .le. zero) then
                    valk(1) = nommai
                    valk(2) = exp(1+j)
                    call utmess('A', 'MODELISA_73', nk=2, valk=valk)
                    ier = ier+1
                end if
            end do
        end if
!
! --- COMPLETUDE DES DONNES GEOMETRIQUES CERCLE
    else if (isec .eq. 2) then
        if (zr(jdge+5) .eq. tst) then
            valk(1) = nommai
            valk(2) = exp(5)
            call utmess('A', 'MODELISA_74', nk=2, valk=valk)
            ier = ier+1
        end if
        if (zr(jdge+5) .le. zero) then
            valk(1) = nommai
            valk(2) = exp(5)
            call utmess('A', 'MODELISA_75', nk=2, valk=valk)
            ier = ier+1
        end if
        if (.not. secple) then
            if (zr(jdge+6) .le. zero) then
                valk(1) = nommai
                valk(2) = exp(6)
                call utmess('A', 'MODELISA_76', nk=2, valk=valk)
                ier = ier+1
            end if
        end if
    end if
!
!
    if (ier .ne. 0) goto 999
!
! --- AFFECTATION DES VALEURS PAR DEFAUT POUR LES DONNEES RECTANGLE
    if (isec .eq. 1) then
        hy = zr(jdge+1)
        hz = zr(jdge+2)
        if (secple) then
            zr(jdge) = hy*hz
        else
            epy = zr(jdge+3)
            epz = zr(jdge+4)
            hyi = hy-2.d0*epy
            hzi = hz-2.d0*epz
            zr(jdge) = hy*hz-hyi*hzi
        end if
        if (zr(jdge+3) .eq. tst) zr(jdge+3) = zr(jdge+1)
        if (zr(jdge+4) .eq. tst) zr(jdge+4) = zr(jdge+2)
!
! --- AFFECTATION DES VALEURS PAR DEFAUT POUR LES DONNEES CERCLE
    else if (isec .eq. 2) then
        re = zr(jdge+5)
        if (secple) then
            zr(jdge) = pi*re*re
        else
            e = zr(jdge+6)
            ri = re-e
            zr(jdge) = pi*(re*re-ri*ri)
        end if
        if (zr(jdge+6) .eq. tst) zr(jdge+6) = zr(jdge+5)
    end if
    do i = 1, ncar
        if (zr(jdge+i-1) .eq. tst) zr(jdge+i-1) = zero
    end do
!
999 continue
    call jedema()
end subroutine
