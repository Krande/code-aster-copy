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
subroutine mnleng(imat, xcdl, parcho, xus, ninc, &
                  nd, nchoc, h, nbpt, xeng)
    implicit none
!
!
!     MODE_NON_LINE CALCUL DE L'ENERGIE MECANIQUE
!     -    -   -    -           -       -
! ----------------------------------------------------------------------
!
! EFFECTUE LE PRODUIT DE DEUX SIGNAUX FREQUENTIELS X ET Y PAR LA METHODE
! AFT  : IFFT -> FFT -> IFFT
! LES COEFFICIENTS SONT RANGES AINSI : Z = [Z0 ZC1...ZCH ZS1...ZSH]
! X ET Y PEUVENT CONTENIR N VECTEURS, PAR EX : X = [Z1 Z2 ...ZN]
! ----------------------------------------------------------------------
! IN  IMAT   : I(2)          : DESCRIPTEUR DES MATRICES :
!                               - IMAT(1) => MATRICE DE RAIDEUR
!                               - IMAT(2) => MATRICE DE MASSE
! IN  XCDL   : K14           : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14           : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  XUS    : K14           : BRANCHE SOLUTION
! IN  IND    : I             : INDICE DISCRETISATION
! IN  OMEGA  : R8            : PULSATION OMEGA
! IN  NINC   : I             : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I             : NOMBRE DE DDLS ACTIFS
! IN  NCHOC  : I             : NOMBRE DE CONTACTEURS
! IN  H      : I             : NOMBRE D'HARMONIQUES
! IN  NBPT   : I             : NOMBRE DE POINT DE DISCRETISATION DE LA
!                                                                BRANCHE
! OUT XENG   : K14           : ENERGIE MECANIQUE
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer(kind=8) :: imat(2), ind, ninc, nd, h, nbpt, nchoc
    character(len=14) :: xcdl, parcho, xus, xeng
    real(kind=8) :: e
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: pi
    real(kind=8) :: omega, alpha, jeu, rayon, origx, origy, ratio
    integer(kind=8) :: ix, iy, idy, imdy, iky
    integer(kind=8) :: ius, ieng, k, icdl, neq, i
    integer(kind=8) :: ireg, nddl, nddlx, nddly
    real(kind=8), pointer :: dye(:) => null()
    real(kind=8), pointer :: kye(:) => null()
    real(kind=8), pointer :: mdye(:) => null()
    real(kind=8), pointer :: ye(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    integer(kind=8), pointer :: vnddl(:) => null()
    character(len=8), pointer :: type(:) => null()
    real(kind=8), pointer :: orig(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
!
! ----------------------------------------------------------------------
! --- RECUPERATION POINTEUR ET TAILLE DE LA MATRICE
! ----------------------------------------------------------------------
    call jeveuo(xus, 'L', ius)
    call jeveuo(xeng, 'E', ieng)
    b_n = to_blas_int(nbpt-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(ieng), b_incx)
    call jeveuo(xcdl, 'L', icdl)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.REG', 'L', ireg)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.ORIG', 'L', vr=orig)
    neq = zi(imat(1)+2)
! ----------------------------------------------------------------------
! --- DECLARATION VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call wkvect('&&MNLENG.X', 'V V R', nd*(2*h+1), ix)
    call wkvect('&&MNLENG.Y', 'V V R', nd, iy)
    call wkvect('&&MNLENG.DY', 'V V R', nd, idy)
    call wkvect('&&MNLENG.MDY', 'V V R', nd, imdy)
    call wkvect('&&MNLENG.KY', 'V V R', nd, iky)
!
    AS_ALLOCATE(vr=ye, size=neq)
    AS_ALLOCATE(vr=dye, size=neq)
    AS_ALLOCATE(vr=kye, size=neq)
    AS_ALLOCATE(vr=mdye, size=neq)
    do ind = 1, nbpt-1
        b_n = to_blas_int(nd*(2*h+1))
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(ix), b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(iy), b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(idy), b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(imdy), b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(iky), b_incx)
!
        omega = zr(ius-1+ind*ninc)
        b_n = to_blas_int(nd*(2*h+1))
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ius+(ind-1)*ninc), b_incx, zr(ix), b_incy)
! ----------------------------------------------------------------------
! --- PASSAGE EN TEMPOREL (t=T/4)
! ----------------------------------------------------------------------
! ---   PI
        pi = r8pi()
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ix), b_incx, zr(iy), b_incy)
        ratio = 4.d0
        do k = 1, h
! ---     COS
            b_n = to_blas_int(nd)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, dcos(2*k*pi/ratio), zr(ix-1+nd*k+1), b_incx, zr(iy), &
                       b_incy)
            b_n = to_blas_int(nd)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, k*omega*dcos(2*k*pi/ratio), zr(ix-1+nd*(h+k)+1), b_incx, zr(idy), &
                       b_incy)
! ---     SIN
            b_n = to_blas_int(nd)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, dsin(2*k*pi/ratio), zr(ix-1+nd*(h+k)+1), b_incx, zr(iy), &
                       b_incy)
            b_n = to_blas_int(nd)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -k*omega*dsin(2*k*pi/ratio), zr(ix-1+nd*k+1), b_incx, zr(idy), &
                       b_incy)
        end do
! ----------------------------------------------------------------------
! --- CALCUL DE K*Y ET M*DY
! ----------------------------------------------------------------------
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, ye, b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, dye, b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, kye, b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, mdye, b_incx)
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                ye(k) = zr(iy-1+i)
                dye(k) = zr(idy-1+i)
            end if
        end do
        call mrmult('ZERO', imat(1), ye, kye, 1, &
                    .false._1)
        call mrmult('ZERO', imat(2), dye, mdye, 1, &
                    .false._1)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(iky), b_incx)
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, 0.d0, zr(imdy), b_incx)
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(iky-1+i) = kye(k)
                zr(imdy-1+i) = mdye(k)
            end if
        end do
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        e = ddot(b_n, zr(iy), b_incx, zr(iky), b_incy)/2
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        e = e+ddot(b_n, zr(idy), b_incx, zr(imdy), b_incy)/2
        do k = 1, nchoc
            alpha = raid(k)
            jeu = vjeu(k)
            if (type(k) (1:4) .eq. 'PLAN') then
                nddl = vnddl(6*(k-1)+1)
                if (zr(iy-1+nddl) .gt. jeu) then
                    e = e+0.5*alpha*(zr(iy-1+nddl)-jeu)**2
                end if
            else if (type(k) (1:7) .eq. 'BI_PLAN') then
                nddl = vnddl(6*(k-1)+1)
                if (zr(iy-1+nddl) .gt. jeu) then
                    e = e+0.5*alpha*(zr(iy-1+nddl)-jeu)**2
                else if (zr(iy-1+nddl) .lt. (-1.d0*jeu)) then
                    e = e+0.5*alpha*(zr(iy-1+nddl)+jeu)**2
                end if
            else if (type(k) (1:6) .eq. 'CERCLE') then
                nddlx = vnddl(6*(k-1)+1)
                nddly = vnddl(6*(k-1)+2)
                origx = orig(3*(k-1)+1)
                origy = orig(3*(k-1)+2)
                rayon = sqrt((zr(iy-1+nddlx)-origx)**2+(zr(iy-1+nddly)-origy)**2)
                if (rayon .gt. jeu) then
                    e = e+alpha*(rayon-jeu)**2
                end if
            end if
        end do
        zr(ieng-1+ind) = e
    end do
!
    call jedetr('&&MNLENG.X')
    call jedetr('&&MNLENG.Y')
    call jedetr('&&MNLENG.DY')
    call jedetr('&&MNLENG.MDY')
    call jedetr('&&MNLENG.KY')
!
    AS_DEALLOCATE(vr=ye)
    AS_DEALLOCATE(vr=dye)
    AS_DEALLOCATE(vr=kye)
    AS_DEALLOCATE(vr=mdye)
!
    call jedema()
!
end subroutine
