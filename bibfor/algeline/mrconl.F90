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
subroutine mrconl(oper, lmat, neq2, typev, rvect, &
                  nvect)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: lmat, neq2, nvect
    character(len=1) :: typev
    character(len=4) :: oper
    real(kind=8) :: rvect(*)
!     TENIR COMPTE DU CONDITIONNEMENT DES LAGRANGE SUR DES VECTEURS
!     ------------------------------------------------------------------
! IN  : OPER    : 'MULT' : ON MULTIPLIE RVECT PAR RCOEF
!                 'DIVI' : ON DIVISE RVECT PAR RCOEF
! IN  : LMAT    : ADRESSE DU DESCRIPTEUR DE LA MATRICE
! IN  : NEQ2    : NOMBRE D'EQUATIONS DU VECTEUR RVECT(NEQ2,NVECT)
!                 SI NEQ2.LE.0 ALORS NEQ2 = LMAT(2)
! IN  : NVECT   : NOMBRE DE VECTEURS DU VECTEUR RVECT(NEQ2,NVECT)
! IN  : TYPEV   : TYPE DES COEFFICIENTS DU VECTEUR
!               = 'R'  : A COEFFICIENTS REELS
!               = 'C'  : A COEFFICIENTS COMPLEXES
!               = ' '  : LES COEFFICIENTS DU VECTEURS SONT DU MEME TYPE
!                        QUE CEUX DE LA MATRICE.
! VAR :  RVECT  : VECTEUR A MODIFIER
!               REMARQUE : RVECT DOIT ETRE DU MEME TYPE QUE LA MATRICE
!                          SOIT REAL*8
!                          SOIT COMPLEX*16
!     ------------------------------------------------------------------
!
!
    character(len=1) :: ftype(2), type, typecn
    character(len=24) :: conl
    complex(kind=8) :: c8cst
    integer(kind=8) :: ieq, ii, ind, iret, ive, jconl, neq
!     ------------------------------------------------------------------
    data ftype/'R', 'C'/
!     ------------------------------------------------------------------
!
    call jemarq()
    ASSERT(oper .eq. 'MULT' .or. oper .eq. 'DIVI')
    if (typev .eq. ' ') then
        type = ftype(zi(lmat+3))
    else
        type = typev
    end if
    neq = neq2
    if (neq2 .le. 0) neq = zi(lmat+2)
!
!
    conl = zk24(zi(lmat+1)) (1:19)//'.CONL'
    call jeexin(conl, iret)
    if (iret .ne. 0) then
        call jelira(conl, 'TYPE', cval=typecn)
        call jeveuo(conl, 'L', jconl)
        jconl = jconl-1
!
!
        if (type .eq. 'R' .and. typecn .eq. 'R') then
            do ive = 1, nvect
                ind = neq*(ive-1)
                if (oper .eq. 'MULT') then
                    do ieq = 1, neq
                        rvect(ind+ieq) = rvect(ind+ieq)*zr(jconl+ieq)
                    end do
                else
                    do ieq = 1, neq
                        rvect(ind+ieq) = rvect(ind+ieq)/zr(jconl+ieq)
                    end do
                end if
            end do
!
!
        else if (type .eq. 'C') then
            if (typecn .eq. 'R') then
                do ive = 1, nvect
                    ind = neq*(ive-1)
                    if (oper .eq. 'MULT') then
                        do ieq = 1, neq
                            ii = ind+2*ieq
                            rvect(ii-1) = rvect(ii-1)*zr(jconl+ieq)
                            rvect(ii) = rvect(ii)*zr(jconl+ieq)
                        end do
                    else
                        do ieq = 1, neq
                            ii = ind+2*ieq
                            rvect(ii-1) = rvect(ii-1)/zr(jconl+ieq)
                            rvect(ii) = rvect(ii)/zr(jconl+ieq)
                        end do
                    end if
                end do
!
            else if (typecn .eq. 'C') then
                do ive = 1, nvect
                    ind = neq*(ive-1)
                    if (oper .eq. 'MULT') then
                        do ieq = 1, neq
                            ii = ind+2*ieq
                            c8cst = dcmplx(rvect(ii-1), rvect(ii))
                            c8cst = c8cst*zc(jconl+ieq)
                            rvect(ii-1) = dble(c8cst)
                            rvect(ii) = dimag(c8cst)
                        end do
                    else
                        do ieq = 1, neq
                            ii = ind+2*ieq
                            c8cst = dcmplx(rvect(ii-1), rvect(ii))
                            c8cst = c8cst/zc(jconl+ieq)
                            rvect(ii-1) = dble(c8cst)
                            rvect(ii) = dimag(c8cst)
                        end do
                    end if
                end do
            end if
        end if
    end if
    call jedema()
end subroutine
