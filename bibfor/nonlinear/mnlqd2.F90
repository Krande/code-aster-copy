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
subroutine mnlqd2(ind, imat, neq, ninc, nd, &
                  nchoc, h, hf, parcho, xcdl, &
                  adime, xvect, xtemp)
    implicit none
!
!
!     MODE_NON_LINE -- MATRICE JACOBIENNE (Q(V,E_I))
!     -    -                -            -   -
! ----------------------------------------------------------------------
!
! CALCUL PARTIELLE DE  LA MATRICE JACOBIENNE POUR UN CERTAIN
!                                                       VECTEUR SOLUTION
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlaft.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer(kind=8) :: ind, imat(2), neq, ninc, nd, nchoc, h, hf
    character(len=14) :: parcho, xcdl, adime, xvect, xtemp
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: jeu, alpha, coef
    integer(kind=8) :: iq2, itemp1, itemp2, itemp3, itemp4
    integer(kind=8) :: iddl, i, nddl
    integer(kind=8) :: icdl, iadim, itemp, k, ivec, nt, ih, puismax
    integer(kind=8) :: neqs, deb, hind, ddl, nddlx, nddly
    aster_logical :: stp
    integer(kind=8), pointer :: vneqs(:) => null()
    real(kind=8), pointer :: jeumax(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    character(len=8), pointer :: type(:) => null()
    integer(kind=8), pointer :: vnddl(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
    puismax = int(dlog(4.d0*dble(hf)+1.d0)/dlog(2.d0)+1.d0)
    nt = 2**puismax
    call wkvect('&&mnlqd2.q2', 'V V R', ninc-1, iq2)
    call wkvect('&&mnlqd2.temp1', 'V V R', neq, itemp1)
    call wkvect('&&mnlqd2.temp2', 'V V R', neq, itemp2)
    call wkvect('&&mnlqd2.temp3', 'V V R', 2*hf+1, itemp3)
    call wkvect('&&mnlqd2.temp4', 'V V R', 2*hf+1, itemp4)
    stp = .true.
!
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', vr=jeumax)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(xcdl, 'L', icdl)
    call jeveuo(adime, 'L', iadim)
    call jeveuo(xvect, 'L', ivec)
    call jeveuo(xtemp, 'E', itemp)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(itemp), b_incx)
! ----------------------------------------------------------------------
! --- INCONNUE DU SYSTEME DYNAMIQUE i.e. ND+1:ND*(2*H+1)
! ----------------------------------------------------------------------
    if (ind .le. nd*(2*h+1) .and. ind .gt. nd) then
        ih = int((ind-1)/nd)
        iddl = ind-nd*int((ind-1)/nd)
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                if (i .eq. iddl) then
                    zr(itemp1-1+k) = 1.d0
                end if
            end if
        end do
        call mrmult('ZERO', imat(2), zr(itemp1), zr(itemp2), 1, &
                    .false._1)
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                if (ih .le. h) then
                    coef = dble(ih)*dble(ih)
                else
                    coef = dble(ih-h)*dble(ih-h)
                end if
                zr(iq2-1+ih*nd+i) = -zr(ivec-1+ninc-2)*coef*zr(itemp2-1+k)/zr(iadim+1)
            end if
        end do
        if (ih .le. h) then
            zr(iq2-1+(h+ih)*nd+iddl) = dble(ih)*zr(ivec-1+ninc-3)
        else
            zr(iq2-1+(ih-h)*nd+iddl) = -dble(ih)*zr(ivec-1+ninc-3)
        end if
    end if
! ----------------------------------------------------------------------
! --- EQUATIONS SUPPLEMENTAIRES
! ----------------------------------------------------------------------
    neqs = 0
    deb = nd*(2*h+1)
    do i = 1, nchoc
        alpha = raid(i)/zr(iadim-1+1)
        jeu = vjeu(i)/jeumax(1)
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
            if ((ind .le. nd*(2*h+1)) .or. ((ind .gt. deb) .and. (ind .le. (deb+(2*hf+1))))) then
! ---     (F/ALPHA-XG))
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp4), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec-1+deb+1), b_incx, zr(itemp4), b_incy)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 1.d0/alpha, zr(itemp4), b_incx)
                b_n = to_blas_int(h+1)
                b_incx = to_blas_int(nd)
                b_incy = to_blas_int(1)
                call daxpy(b_n, -1.d0/jeu, zr(ivec-1+nddl), b_incx, zr(itemp4), &
                           b_incy)
                b_n = to_blas_int(h)
                b_incx = to_blas_int(nd)
                b_incy = to_blas_int(1)
                call daxpy(b_n, -1.d0/jeu, zr(ivec-1+nd*(h+1)+nddl), b_incx, zr(itemp4-1+hf+2), &
                           b_incy)
            end if
            if (ind .le. nd*(2*h+1)) then
                hind = int((ind-1)/nd)
!            WRITE(6,*) 'HIND',HIND
                ddl = ind-nd*hind
! ---     -(F/ALPHA-XG)*(F/ALPHA-XG))
                if (ddl .eq. nddl) then
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                    if (hind .le. h) then
                        zr(itemp3-1+hind+1) = -1.d0/jeu
                    else
                        zr(itemp3-1+hf+1+hind-h) = -1.d0/jeu
                    end if
!              WRITE(6,*) 'TEMP3',TEMP3(1:2*HF+1)
!              WRITE(6,*) 'TEMP4',TEMP4(1:2*HF+1)
                    call mnlaft(zr(itemp4), zr(itemp3), hf, nt, zr(iq2-1+deb+(2*hf+1)+1))
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, -1.d0, zr(iq2-1+deb+(2*hf+1)+1), b_incx)
!              WRITE(6,*) 'Q1',IND,DEB,Q1(DEB+(2*HF+1)+1:DEB+2*(2*HF+1))
                end if
            else if ((ind .gt. deb) .and. (ind .le. (deb+(2*hf+1)))) then
! ---     -(F/ALPHA-XG)*(F/ALPHA-XG))
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                zr(itemp3-1+ind-deb) = 1.d0/alpha
                call mnlaft(zr(itemp4), zr(itemp3), hf, nt, zr(iq2-1+deb+(2*hf+1)+1))
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, -1.d0, zr(iq2-1+deb+(2*hf+1)+1), b_incx)
            else if ((ind .gt. (deb+2*hf+1) .and. ind .le. (deb+4*hf+2))) &
                then
! ---     -F*Z
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp4), b_incx)
                zr(itemp3-1+ind-deb-(2*hf+1)) = -1.d0
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec-1+deb+1), b_incx, zr(itemp4), b_incy)
                call mnlaft(zr(itemp4), zr(itemp3), hf, nt, zr(iq2-1+deb+1))
            end if
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
            if (ind .le. nd*(2*h+1)) then
                hind = int((ind-1)/nd)
                ddl = ind-nd*hind
                if ((ddl .eq. nddlx) .or. (ddl .eq. nddly)) then
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 0.d0, zr(itemp4), b_incx)
                    if (hind .le. h) then
                        zr(itemp4-1+hind+1) = 1.d0/jeu
                    else
                        zr(itemp4-1+hf+1+hind-h) = 1.d0/jeu
                    end if
! ---         FX*R - FN*([UX]/JEU)
! ---         FY*R - FN*([UY]/JEU)
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(ivec+deb+3*(2*hf+1)), b_incx, zr(itemp3), b_incy)
                    if (ddl .eq. nddlx) then
                        call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+1))
                        b_n = to_blas_int(2*hf+1)
                        b_incx = to_blas_int(1)
                        call dscal(b_n, -1.d0, zr(iq2-1+deb+1), b_incx)
                    else if (ddl .eq. nddly) then
                        call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+(2*hf+1)+1))
                        b_n = to_blas_int(2*hf+1)
                        b_incx = to_blas_int(1)
                        call dscal(b_n, -1.d0, zr(iq2-1+deb+(2*hf+1)+1), b_incx)
                    end if
! ---         R*R - ([UX]/JEU)^2 - ([UY]/JEU)^2
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                    b_n = to_blas_int(h+1)
                    b_incx = to_blas_int(nd)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(ivec-1+ddl), b_incx, zr(itemp3), b_incy)
                    b_n = to_blas_int(h)
                    b_incx = to_blas_int(nd)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(ivec-1+nd*(h+1)+ddl), b_incx, zr(itemp3-1+hf+2), b_incy)
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 1.d0/jeu, zr(itemp3), b_incx)
                    call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+2*(2*hf+1)+1))
                    b_n = to_blas_int(2*hf+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, -1.d0, zr(iq2-1+deb+2*(2*hf+1)+1), b_incx)
                end if
            else if (ind .gt. deb+2*(2*hf+1) .and. ind .le. deb+3*(2*hf+1)) then
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp4), b_incx)
                zr(itemp4-1+ind-deb-2*(2*hf+1)) = 1.d0
! ---       FX*[R] - FN*(UX/JEU)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec+deb), b_incx, zr(itemp3), b_incy)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+1))
! ---       FY*[R] - FN*(UY/JEU)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec+deb+(2*hf+1)), b_incx, zr(itemp3), b_incy)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+(2*hf+1)+1))
! ---       R*[R] - (UX/JEU)^2 - (UY/JEU)^2
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec+deb+2*(2*hf+1)), b_incx, zr(itemp3), b_incy)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+2*(2*hf+1)+1))
            else if (ind .gt. deb+3*(2*hf+1) .and. ind .le. deb+4*(2*hf+1)) then
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp4), b_incx)
                zr(itemp4-1+ind-deb-3*(2*hf+1)) = 1.d0
! ---       (FN/ALPHA - R)*[FN]
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, -1.d0, zr(ivec+deb+2*(2*hf+1)), b_incx, zr(itemp3), &
                           b_incy)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, 1.d0/alpha, zr(ivec+deb+3*(2*hf+1)), b_incx, zr(itemp3), &
                           b_incy)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+3*(2*hf+1)+1))
            end if
        else if (type(i) (1:4) .eq. 'PLAN') then
            nddl = vnddl(6*(i-1)+1)
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itemp3), b_incx)
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itemp4), b_incx)
            if (ind .gt. deb .and. ind .le. deb+(2*hf+1)) then
! ---       (F/ALPHA - XG)*[F]
                b_n = to_blas_int(h+1)
                b_incx = to_blas_int(nd)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec-1+nddl), b_incx, zr(itemp3-1+1:h+1), b_incy)
                b_n = to_blas_int(h)
                b_incx = to_blas_int(nd)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(ivec-1+nd*(h+1)+nddl), b_incx, zr(itemp3-1+hf+2:hf+h+1), &
                           b_incy)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, -1.d0, zr(itemp3), b_incx)
                b_n = to_blas_int(2*hf+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call daxpy(b_n, 1.d0/alpha, zr(ivec+deb), b_incx, zr(itemp3), &
                           b_incy)
                zr(itemp4-1+ind-deb) = 1.d0
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq2-1+deb+1))
            end if
        end if
        neqs = neqs+vneqs(i)
        deb = deb+vneqs(i)*(2*hf+1)
    end do
!
! ----------------------------------------------------------------------
! --- GAMMA1 i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+1
! --- GAMMA2 i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+2
! ----------------------------------------------------------------------
    if (ind .eq. ninc) then
        zr(iq2-1+ninc-3) = -1.d0*zr(ivec-1+ninc-1)
        zr(iq2-1+ninc-2) = -1.d0*zr(ivec-1+ninc)
    end if
! ----------------------------------------------------------------------
! --- EQUATION DE PHASE i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+3
! ----------------------------------------------------------------------
    do k = 1, h
        if (ind .eq. (nd*(h+k)+1)) then
            zr(iq2-1+ninc-1) = k*zr(ivec-1+ninc)
        end if
    end do
!
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(iq2), b_incx, zr(itemp), b_incy)
!
    call jedetr('&&mnlqd2.q2')
    call jedetr('&&mnlqd2.temp1')
    call jedetr('&&mnlqd2.temp2')
    call jedetr('&&mnlqd2.temp3')
    call jedetr('&&mnlqd2.temp4')
!
    call jedema()
!
end subroutine
