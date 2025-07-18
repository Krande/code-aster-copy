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
subroutine te0403(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/fcent.h"
#include "asterfort/fointe.h"
#include "asterfort/fpesa.h"
#include "asterfort/fpres.h"
#include "asterfort/fsurf.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecael.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
    character(len=16) :: option, nomte
!
!
    integer(kind=8) :: nb1
    real(kind=8) :: vecl(51)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectt(9, 2, 3), vectpt(9, 3, 3)
    real(kind=8) :: valpar(4), pr
    character(len=8) :: nompar(4), nomail
    character(len=24) :: valk
! DEB ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadzi, iazk24, ib, ier, itemps, j
    integer(kind=8) :: jgeom, jpres, jvecg, lzi, lzr, nb2
!-----------------------------------------------------------------------
    call jevech('PGEOMER', 'L', jgeom)
!
    call jevech('PVECTUR', 'E', jvecg)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    call vectan(nb1, nb2, zr(jgeom), zr(lzr), vecta, &
                vectn, vectt)
!
    do ib = 1, nb2
        do i = 1, 2
            do j = 1, 3
                vectpt(ib, i, j) = vectt(ib, i, j)
            end do
        end do
        vectpt(ib, 3, 1) = vectn(ib, 1)
        vectpt(ib, 3, 2) = vectn(ib, 2)
        vectpt(ib, 3, 3) = vectn(ib, 3)
    end do
!
!
    if (option .eq. 'CHAR_MECA_FRCO3D' .or. option .eq. 'CHAR_MECA_FFCO3D') then
!------------------------------------------------------
!      PAS DE CHANGEMENT DE SIGNE POUR LES FORCES REPARTIES
!------------------------------------------------------
        call fsurf(option, nomte, zr(jgeom), nb1, vecl, &
                   vectpt)
!
    else if (option .eq. 'CHAR_MECA_PESA_R') then
        call fpesa(nomte, zr(jgeom), nb1, vecl)
!
    else if (option .eq. 'CHAR_MECA_ROTA_R') then
        call fcent(nomte, zr(jgeom), nb1, vecl)
!
    else if (option .eq. 'CHAR_MECA_PRES_R') then
!------------------------------------------------------
!      CHANGEMENT DE SIGNE POUR LES PRESSIONS DANS FPRES
!------------------------------------------------------
        call fpres(nomte, zr(jgeom), nb1, vecl, vectpt)
!
    else if (option .eq. 'CHAR_MECA_PRES_F') then
        call jevech('PPRESSF', 'L', jpres)
        if (zk8(jpres) .eq. '&FOZERO') goto 999
        call jevech('PINSTR', 'L', itemps)
        valpar(4) = zr(itemps)
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        do j = 0, nb1-1
            valpar(1) = zr(jgeom+3*j)
            valpar(2) = zr(jgeom+3*j+1)
            valpar(3) = zr(jgeom+3*j+2)
            call fointe('FM', zk8(jpres), 4, nompar, valpar, &
                        pr, ier)
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                valk = nomail
                call utmess('F', 'ELEMENTS4_92', sk=valk)
            end if
        end do
        goto 999
    end if
!
    call trnflg(nb2, vectpt, vecl, zr(jvecg))
!
999 continue
end subroutine
