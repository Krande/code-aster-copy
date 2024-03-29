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
    integer :: imat(2), ind, ninc, nd, h, nbpt, nchoc
    character(len=14) :: xcdl, parcho, xus, xeng
    real(kind=8) :: e
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: pi
    real(kind=8) :: omega, alpha, jeu, rayon, origx, origy, ratio
    integer :: ix, iy, idy, imdy, iky
    integer :: ius, ieng, k, icdl, neq, i
    integer :: ireg, nddl, nddlx, nddly
    real(kind=8), pointer :: dye(:) => null()
    real(kind=8), pointer :: kye(:) => null()
    real(kind=8), pointer :: mdye(:) => null()
    real(kind=8), pointer :: ye(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    integer, pointer :: vnddl(:) => null()
    character(len=8), pointer :: type(:) => null()
    real(kind=8), pointer :: orig(:) => null()
!
    call jemarq()
!
! ----------------------------------------------------------------------
! --- RECUPERATION POINTEUR ET TAILLE DE LA MATRICE
! ----------------------------------------------------------------------
    call jeveuo(xus, 'L', ius)
    call jeveuo(xeng, 'E', ieng)
    call dscal(nbpt-1, 0.d0, zr(ieng), 1)
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
        call dscal(nd*(2*h+1), 0.d0, zr(ix), 1)
        call dscal(nd, 0.d0, zr(iy), 1)
        call dscal(nd, 0.d0, zr(idy), 1)
        call dscal(nd, 0.d0, zr(imdy), 1)
        call dscal(nd, 0.d0, zr(iky), 1)
!
        omega = zr(ius-1+ind*ninc)
        call dcopy(nd*(2*h+1), zr(ius+(ind-1)*ninc), 1, zr(ix), 1)
! ----------------------------------------------------------------------
! --- PASSAGE EN TEMPOREL (t=T/4)
! ----------------------------------------------------------------------
! ---   PI
        pi = r8pi()
        call dcopy(nd, zr(ix), 1, zr(iy), 1)
        ratio = 4.d0
        do k = 1, h
! ---     COS
            call daxpy(nd, dcos(2*k*pi/ratio), zr(ix-1+nd*k+1), 1, zr(iy), &
                       1)
            call daxpy(nd, k*omega*dcos(2*k*pi/ratio), zr(ix-1+nd*(h+k)+1), 1, zr(idy), &
                       1)
! ---     SIN
            call daxpy(nd, dsin(2*k*pi/ratio), zr(ix-1+nd*(h+k)+1), 1, zr(iy), &
                       1)
            call daxpy(nd, -k*omega*dsin(2*k*pi/ratio), zr(ix-1+nd*k+1), 1, zr(idy), &
                       1)
        end do
! ----------------------------------------------------------------------
! --- CALCUL DE K*Y ET M*DY
! ----------------------------------------------------------------------
        call dscal(nd, 0.d0, ye, 1)
        call dscal(nd, 0.d0, dye, 1)
        call dscal(nd, 0.d0, kye, 1)
        call dscal(nd, 0.d0, mdye, 1)
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
        call dscal(nd, 0.d0, zr(iky), 1)
        call dscal(nd, 0.d0, zr(imdy), 1)
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(iky-1+i) = kye(k)
                zr(imdy-1+i) = mdye(k)
            end if
        end do
        e = ddot(nd, zr(iy), 1, zr(iky), 1)/2
        e = e+ddot(nd, zr(idy), 1, zr(imdy), 1)/2
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
