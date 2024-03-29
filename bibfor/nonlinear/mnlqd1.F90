! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine mnlqd1(ind, imat, neq, ninc, nd, &
                  nchoc, h, hf, parcho, xcdl, &
                  adime, xvect, xtemp)
    implicit none
!
!
!     MODE_NON_LINE -- MATRICE JACOBIENNE (Q(E_I,V))
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
    integer :: ind, imat(2), neq, ninc, nd, nchoc, h, hf
    character(len=14) :: parcho, xcdl, adime, xvect, xtemp
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: jeu, alpha, coef
    integer :: iq1, itemp1, itemp2, itemp3, itemp4
    integer :: j, i, nddl
    integer :: icdl, iadim, itemp, k, ivec, nt, puismax
    integer :: neqs, deb, hind, ddl, nddlx, nddly
    aster_logical :: stp
    character(len=8), pointer :: type(:) => null()
    integer, pointer :: vnddl(:) => null()
    integer, pointer :: vneqs(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: jeumax(:) => null()
    real(kind=8), pointer :: raid(:) => null()
!
    call jemarq()
!
    puismax = int(dlog(4.d0*dble(hf)+1.d0)/dlog(2.d0)+1.d0)
    nt = 2**puismax
    call wkvect('&&mnlqd1.q1', 'V V R', ninc-1, iq1)
    call wkvect('&&mnlqd1.temp1', 'V V R', neq*2*h, itemp1)
    call wkvect('&&mnlqd1.temp2', 'V V R', neq*2*h, itemp2)
    call wkvect('&&mnlqd1.temp3', 'V V R', 2*hf+1, itemp3)
    call wkvect('&&mnlqd1.temp4', 'V V R', 2*hf+1, itemp4)
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
    call dscal(ninc-1, 0.d0, zr(itemp), 1)
! ----------------------------------------------------------------------
! --- INCONNUE DU SYSTEME DYNAMIQUE i.e. ND+1:ND*(2*H+1)
! ----------------------------------------------------------------------
    if (ind .eq. (ninc-2)) then
        do j = 1, 2*h
            i = 0
            do k = 1, neq
                if (zi(icdl-1+k) .eq. 0) then
                    i = i+1
                    zr(itemp1-1+(j-1)*neq+k) = zr(ivec-1+j*nd+i)
                end if
            end do
        end do
        call mrmult('ZERO', imat(2), zr(itemp1), zr(itemp2), 2*h, &
                    .false._1)
        do j = 1, 2*h
            i = 0
            do k = 1, neq
                if (zi(icdl-1+k) .eq. 0) then
                    i = i+1
                    if (j .le. h) then
                        coef = dble(j)*dble(j)
                    else
                        coef = dble(j-h)*dble(j-h)
                    end if
                    zr(iq1-1+j*nd+i) = -coef*zr(itemp2-1+(j-1)*neq+k)/zr(iadim+1)
                end if
            end do
        end do
    else if (ind .eq. (ninc-3)) then
        call dcopy(nd*h, zr(ivec+nd*(h+1)), 1, zr(iq1-1+nd+1), 1)
        call dscal(nd*h, -1.d0, zr(iq1-1+nd+1), 1)
        call dcopy(nd*h, zr(ivec+nd), 1, zr(iq1-1+nd*(h+1)+1), 1)
        do k = 1, h
            call dscal(nd, dble(k), zr(iq1-1+k*nd+1), 1)
            call dscal(nd, dble(k), zr(iq1-1+(h+k)*nd+1), 1)
        end do
    end if
! ----------------------------------------------------------------------
! --- EQUATIONS SUPPLEMENTAIRES
! ----------------------------------------------------------------------
    neqs = 0
    deb = nd*(2*h+1)
    hind = int((ind-1)/nd)
    ddl = ind-nd*hind
    do i = 1, nchoc
        alpha = raid(i)/zr(iadim-1+1)
        jeu = vjeu(i)/jeumax(1)
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
!          WRITE(6,*) 'NDDL',NDDL
            if (ind .le. deb+(2*hf+1)) then
! ---     (F/ALPHA-XG))
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(2*hf+1, zr(ivec-1+deb+1), 1, zr(itemp4), 1)
                call dscal(2*hf+1, 1.d0/alpha, zr(itemp4), 1)
                call daxpy(h+1, -1.d0/jeu, zr(ivec-1+nddl), nd, zr(itemp4), &
                           1)
                call daxpy(h, -1.d0/jeu, zr(ivec-1+nd*(h+1)+nddl), nd, zr(itemp4-1+hf+2), &
                           1)
            end if
            if (ind .le. nd*(2*h+1)) then
! ---     -(F/ALPHA-[XG])*(F/ALPHA-XG))
                if (ddl .eq. nddl) then
                    call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                    if (hind .le. h) then
                        zr(itemp3-1+hind+1) = -1.d0/jeu
                    else
                        zr(itemp3-1+hf+1+hind-h) = -1.d0/jeu
                    end if
!              WRITE(6,*) 'TEMP3',TEMP3(1:2*HF+1)
!              WRITE(6,*) 'TEMP4',TEMP4(1:2*HF+1)
                    call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+(2*hf+1)+1))
                    call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+(2*hf+1)+1), 1)
                end if
            else if ((ind .gt. deb) .and. (ind .le. (deb+(2*hf+1)))) then
! ---     -(F/ALPHA-XG)*(F/ALPHA-XG))
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                zr(itemp3-1+ind-deb) = 1.d0/alpha
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+(2*hf+1)+1))
                call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+(2*hf+1)+1), 1)
! ---     -F*Z
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                zr(itemp3-1+ind-deb) = 1.d0
                call dcopy(2*hf+1, zr(ivec-1+deb+(2*hf+1)+1), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+1))
                call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+1), 1)
            end if
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
            if (ind .le. nd*(2*h+1)) then
                if ((ddl .eq. nddlx) .or. (ddl .eq. nddly)) then
! ---         R*R - ([UX]/JEU)^2 - ([UY]/JEU)^2
                    call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                    call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                    if (hind .le. h) then
                        zr(itemp3-1+hind+1) = 1.d0/jeu
                    else
                        zr(itemp3-1+hf+1+hind-h) = 1.d0/jeu
                    end if
                    call dcopy(h+1, zr(ivec-1+ddl), nd, zr(itemp4), 1)
                    call dcopy(h, zr(ivec-1+nd*(h+1)+ddl), nd, zr(itemp4-1+hf+2), 1)
                    call dscal(2*hf+1, 1.d0/jeu, zr(itemp4), 1)
                    call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+2*(2*hf+1)+1))
                    call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+2*(2*hf+1)+1), 1)
                end if
            else if (ind .gt. deb .and. ind .le. deb+(2*hf+1)) then
! ---       [FX]*R - FN*(UX/JEU)
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                zr(itemp3-1+ind-deb) = 1.d0
                call dcopy(2*hf+1, zr(ivec+deb+2*(2*hf+1)), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+1))
            else if (ind .gt. deb+(2*hf+1) .and. ind .le. deb+2*(2*hf+1)) &
                then
! ---       [FY]*R - FN*(UY/JEU)
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                zr(itemp3-1+ind-deb-(2*hf+1)) = 1.d0
                call dcopy(2*hf+1, zr(ivec+deb+2*(2*hf+1)), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+(2*hf+1)+1))
            else if (ind .gt. deb+2*(2*hf+1) .and. ind .le. deb+3*(2*hf+1)) &
                then
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                zr(itemp3-1+ind-deb-2*(2*hf+1)) = 1.d0
! ---       [R]*R - (UX/JEU)^2 - (UY/JEU)^2
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(2*hf+1, zr(ivec+deb+2*(2*hf+1)), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+2*(2*hf+1)+1))
! ---       (FN/ALPHA - [R])*FN
                zr(itemp3-1+ind-deb-2*(2*hf+1)) = -1.d0
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(2*hf+1, zr(ivec+deb+3*(2*hf+1)), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+3*(2*hf+1)+1))
            else if (ind .gt. deb+3*(2*hf+1) .and. ind .le. deb+4*(2*hf+1)) &
                then
                call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
                zr(itemp3-1+ind-deb-3*(2*hf+1)) = 1.d0
! ---       FX*R - [FN]*(UX/JEU)
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(h+1, zr(ivec-1+nddlx), nd, zr(itemp4), 1)
                call dcopy(h, zr(ivec-1+nd*(h+1)+nddlx), nd, zr(itemp4-1+hf+2), 1)
                call dscal(2*hf+1, 1.d0/jeu, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+1))
                call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+1), 1)
! ---       FY*R - [FN]*(UY/JEU)
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(h+1, zr(ivec-1+nddly), nd, zr(itemp4), 1)
                call dcopy(h, zr(ivec-1+nd*(h+1)+nddly), nd, zr(itemp4-1+hf+2), 1)
                call dscal(2*hf+1, 1.d0/jeu, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+(2*hf+1)+1))
                call dscal(2*hf+1, -1.d0, zr(iq1-1+deb+(2*hf+1)+1), 1)
! ---       ([FN/ALPHA] - R)*FN
                zr(itemp3-1+ind-deb-3*(2*hf+1)) = 1.d0/alpha
                call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
                call dcopy(2*hf+1, zr(ivec+deb+3*(2*hf+1)), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+3*(2*hf+1)+1))
            end if
        else if (type(i) (1:4) .eq. 'PLAN') then
            nddl = vnddl(6*(i-1)+1)
            call dscal(2*hf+1, 0.d0, zr(itemp3), 1)
            call dscal(2*hf+1, 0.d0, zr(itemp4), 1)
            if (ind .le. nd*(2*h+1)) then
! ---       (F/ALPHA - [XG])*F
                if (ddl .eq. nddl) then
                    if (hind .le. h) then
                        zr(itemp3-1+hind+1) = -1.d0/jeu
                    else
                        zr(itemp3-1+hf+1+hind-h) = -1.d0/jeu
                    end if
                    call dcopy(2*hf+1, zr(ivec+deb), 1, zr(itemp4), 1)
                    call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+1))
                end if
            else if (ind .gt. deb .and. ind .le. deb+(2*hf+1)) then
! ---       ([F]/ALPHA - XG)*F
                zr(itemp3-1+ind-deb) = 1.d0/alpha
                call dcopy(2*hf+1, zr(ivec+deb), 1, zr(itemp4), 1)
                call mnlaft(zr(itemp3), zr(itemp4), hf, nt, zr(iq1-1+deb+1))
            end if
        end if
        neqs = neqs+vneqs(i)
        deb = deb+vneqs(i)*(2*hf+1)
    end do
! ----------------------------------------------------------------------
! --- GAMMA1 i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+1
! ----------------------------------------------------------------------
    if (ind .eq. (ninc-1)) then
        zr(iq1-1+ninc-3) = -1.d0*zr(ivec-1+ninc)
    end if
! ----------------------------------------------------------------------
! --- GAMMA2 i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+2
! ----------------------------------------------------------------------
    if (ind .eq. ninc) then
        zr(iq1-1+ninc-2) = -1.d0*zr(ivec-1+ninc)
    end if
! ----------------------------------------------------------------------
! --- EQUATION DE PHASE i.e. ND*(2*H+1)+2*NCHOC(2*HF+1)+3
! ----------------------------------------------------------------------
    if (ind .eq. ninc) then
        do k = 1, h
            zr(iq1-1+ninc-1) = zr(iq1-1+ninc-1)+k*zr(ivec-1+(h+k)*nd+1)
        end do
    end if
!
    call dcopy(ninc-1, zr(iq1), 1, zr(itemp), 1)
!
    call jedetr('&&mnlqd1.q1')
    call jedetr('&&mnlqd1.temp1')
    call jedetr('&&mnlqd1.temp2')
    call jedetr('&&mnlqd1.temp3')
    call jedetr('&&mnlqd1.temp4')
!
    call jedema()
!
end subroutine
