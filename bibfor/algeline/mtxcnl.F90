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
subroutine mtxcnl(cumul, typcst, const, typmat, lmat, &
                  typres, lres, neq)
    implicit none
#include "jeveux.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: lmat, lres
    character(len=*) :: cumul, typcst
    character(len=1) :: typmat, typres
    real(kind=8) :: const(2)
!     MIXAGE DES .CONL    ---> LIRE AVERTISSEMENT CI DESSOUS
!     ------------------------------------------------------------------
!     CECI EST UNE ROUTINE INTERNE VOUS N'AVEZ PAS LE DROIT DE L'APPELER
!     DIRECTEMENT, NI MEME DE CRITIQUER.
!     NEANMOINS SI VOUS VOULEZ LA REECRIRE A VOTRE AISE.
!     ------------------------------------------------------------------
!
    character(len=24) :: valk(3)
!
    real(kind=8) :: un, zero, rcum
    complex(kind=8) :: cun, c8cst
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ival, neq
!-----------------------------------------------------------------------
    zero = 0.d0
    un = 1.d0
    cun = dcmplx(un, 0.d0)
!
    if (cumul .eq. 'CUMU') then
        rcum = un
    else
        rcum = zero
        if (typres .eq. 'R') then
            do ival = 0, neq-1
                zr(lres+ival) = un
            end do
        else if (typres .eq. 'C') then
            do ival = 0, neq-1
                zc(lres+ival) = cun
            end do
        end if
    end if
!
! --- MATRICE REELLE EN RESULTAT
!
    if (typres .eq. 'R') then
        if (typmat .eq. 'R') then
            if (typcst(1:1) .eq. 'R') then
                do ival = 0, neq-1
                    if (zr(lmat+ival) .ne. un) then
                        zr(lres+ival) = rcum*zr(lres+ival)+const(1)*zr(lmat+ival)
                    end if
                end do
            else
                valk(1) = typres
                valk(2) = typmat
                valk(3) = typcst(1:1)
                call utmess('F', 'ALGELINE4_25', nk=3, valk=valk)
            end if
        else if (typmat .eq. 'C') then
            if (typcst(1:1) .eq. 'R') then
                do ival = 0, neq-1
                    if (zc(lmat+ival) .ne. cun) then
                        zr(lres+ival) = rcum*zr(lres+ival)+const(1)*dble(zc(lmat+ival))
                    end if
                end do
            else if (typcst(1:1) .eq. 'C') then
                c8cst = dcmplx(const(1), const(2))
                do ival = 0, neq-1
                    if (zc(lmat+ival) .ne. cun) then
                        zr(lres+ival) = rcum*zr(lres+ival)+dble(c8cst*zc(lmat+ival))
                    end if
                end do
            else
                valk(1) = typres
                valk(2) = typmat
                valk(3) = typcst(1:1)
                call utmess('F', 'ALGELINE4_25', nk=3, valk=valk)
            end if
        else
            valk(1) = typres
            valk(2) = typmat
            call utmess('F', 'ALGELINE4_27', nk=2, valk=valk)
        end if
!
! --- MATRICE COMPLEXE EN RESULTAT
!
    else if (typres .eq. 'C') then
        if (typmat .eq. 'C') then
            if (typcst(1:1) .eq. 'R') then
                do ival = 0, neq-1
                    if (zc(lmat+ival) .ne. cun) then
                        zc(lres+ival) = rcum*zc(lres+ival)+const(1)*zc(lmat+ival)
                    end if
                end do
            else if (typcst(1:1) .eq. 'C') then
                c8cst = dcmplx(const(1), const(2))
                do ival = 0, neq-1
                    if (zc(lmat+ival) .ne. cun) then
                        zc(lres+ival) = rcum*zc(lres+ival)+c8cst*zc(lmat+ival)
                    end if
                end do
            else
                valk(1) = typres
                valk(2) = typmat
                valk(3) = typcst(1:1)
                call utmess('F', 'ALGELINE4_25', nk=3, valk=valk)
            end if
        else if (typmat .eq. 'R') then
            if (typcst(1:1) .eq. 'R') then
                do ival = 0, neq-1
                    if (zr(lmat+ival) .ne. un) then
                        zc(lres+ival) = rcum*zc(lres+ival)+const(1)*zr(lmat+ival)
                    end if
                end do
            else if (typcst(1:1) .eq. 'C') then
                c8cst = dcmplx(const(1), const(2))
                do ival = 0, neq-1
                    if (zr(lmat+ival) .ne. un) then
                        zc(lres+ival) = rcum*zc(lres+ival)+c8cst*zr(lmat+ival)
                    end if
                end do
            else
                valk(1) = typres
                valk(2) = typmat
                valk(3) = typcst(1:1)
                call utmess('F', 'ALGELINE4_25', nk=3, valk=valk)
            end if
        else
            valk(1) = typres
            valk(2) = typmat
            call utmess('F', 'ALGELINE4_27', nk=2, valk=valk)
        end if
    else
        valk(1) = typres
        call utmess('F', 'ALGELINE4_31', sk=valk(1))
!
    end if
!
end subroutine
