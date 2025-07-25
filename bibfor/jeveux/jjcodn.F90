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
function jjcodn(icre, nomrep, nomec, irep, crep, &
                nmax, nuti)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "asterfort/jxhcod.h"
#include "asterfort/utmess.h"
    character(len=*) :: nomrep, nomec, crep(*)
    integer(kind=8) :: icre, irep(*), nmax, nuti
! ----------------------------------------------------------------------
    integer(kind=8) :: ilorep, ideno, ilnom, iluti, idehc
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idehco, idenom, iin, in, j, jin
    integer(kind=8) :: jjcodn, k, ll, lnom, lorep, ne
!-----------------------------------------------------------------------
    parameter(ilorep=1, ideno=2, ilnom=3, iluti=5, idehc=6)
    integer(kind=8) :: iret
    character(len=32) :: cle, nom, valk(2)
    aster_logical :: rinser
! DEB ------------------------------------------------------------------
    jjcodn = 0
    rinser = .false.
    iret = 0
    lorep = irep(ilorep)
    idenom = irep(ideno)
    lnom = irep(ilnom)
    idehco = irep(idehc)
    ll = min(lnom, len(nomec))
    nom = nomec(1:ll)
    i = jxhcod(nom, lorep)
    ne = 1
    valk(1) = nom
    valk(2) = nomrep
!
5   continue
    if (irep(idehco+i) .eq. 0 .and. .not. rinser) then
        if (icre .eq. 3) then
            if (nuti .ge. nmax) then
                call utmess('F', 'JEVEUX1_33', sk=valk(2))
            else
                j = nuti+1
                do k = 1, ll
                    crep(idenom+lnom*(j-1)+k) = nomec(k:k)
                end do
                nuti = nuti+1
                irep(iluti) = nuti
                irep(idehco+i) = j
                iret = j
            end if
        else
            if (icre .eq. 0) then
                iret = 0
            else
                call utmess('F', 'JEVEUX1_34', nk=2, valk=valk)
            end if
        end if
    else
        j = irep(idehco+i)
        do k = 1, ll
            cle(k:k) = crep(idenom+lnom*(abs(j)-1)+k)
        end do
        do k = ll+1, 32
            cle(k:k) = ' '
        end do
        if (cle .eq. nom) then
            if (icre .eq. 3) then
                call utmess('F', 'JEVEUX1_35', nk=2, valk=valk)
            else if (icre .eq. 0) then
                iret = j
            else if (icre .eq. -3) then
                irep(idehco+i) = -j
                crep(idenom+lnom*(j-1)+1) = '?'
                iret = -j
            end if
        else
            if (j .lt. 0 .and. .not. rinser) then
                if (icre .eq. 3) then
                    rinser = .true.
                    jin = j
                    iin = i
                end if
            end if
            if (ne .eq. 1) in = jxhcod(nom, lorep-2)
            ne = ne+1
            i = 1+mod(i+in, lorep)
            if (ne .le. lorep) then
                j = irep(idehco+i)
                if (j .eq. 0 .and. rinser) goto 10
                goto 5
            else
                if (icre .eq. 3) then
                    call utmess('F', 'JEVEUX1_36', nk=2, valk=valk)
                else if (icre .eq. 0) then
                    iret = 0
                end if
            end if
        end if
    end if
10  continue
    if (rinser) then
        irep(idehco+iin) = -jin
        do k = 1, ll
            crep(idenom+lnom*(-jin-1)+k) = nomec(k:k)
        end do
        iret = -jin
    end if
    jjcodn = iret
! FIN ------------------------------------------------------------------
end function
