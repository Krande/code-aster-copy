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
subroutine mnline(imat, xcdl, parcho, adime, xvect, &
                  ninc, nd, nchoc, h, hf, &
                  xline)
    implicit none
!
!
!     MODE_NON_LINE PARTIE LINEAIRE
!     -    -               ---
! ----------------------------------------------------------------------
!
! REGROUPE LES TERMES LINEAIRES DU PROBLEME A RESOUDRE
! ----------------------------------------------------------------------
! IN  IMAT   : I(2) : DESCRIPTEUR DES MATRICES :
!                        - IMAT(1) => MATRICE DE RAIDEUR
!                        - IMAT(2) => MATRICE DE MASSE
! IN  XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  XVECT  : K14  : NOM DU VECTEUR SOLUTION
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I    : NOMBRE DE DEGRES DE LIBERTE
! IN  NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN  H      : I    : NOMBRE D'HARMONIQUES de X
! IN  HF     : I    : NOMBRE D'HARMONIQUES POUR F
! OUT XLINE  : I    : NOM DU VECTEUR DES TERMES LINEAIRES
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf
    character(len=14) :: xcdl, parcho, adime, xvect, xline
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: alpha, eta, jeu
    integer(kind=8) :: neq, ivect, icdl, iline, ivect1, ivect2, j, i, k
    integer(kind=8) :: nddl, iadim
    integer(kind=8) :: neqs, ncmp, nddlx, nddly
    real(kind=8), pointer :: orig(:) => null()
    real(kind=8), pointer :: reg(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: jeumax(:) => null()
    integer(kind=8), pointer :: vnddl(:) => null()
    character(len=8), pointer :: type(:) => null()
    integer(kind=8), pointer :: vneqs(:) => null()
    integer(kind=8), pointer :: vncmp(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
! ----------------------------------------------------------------------
! --- RECUPERATION DU NOM DE LA MATRICE ET TAILLE DE LA MATRICE
! ----------------------------------------------------------------------
    neq = zi(imat(1)+2)
! ----------------------------------------------------------------------
! --- RECUPERATION POINTEUR DE XVECT, L(XVECT) ET XCDL
! ----------------------------------------------------------------------
    call jeveuo(xvect, 'L', ivect)
    call jeveuo(xcdl, 'L', icdl)
    call jeveuo(adime, 'L', iadim)
    call jeveuo(xline, 'E', iline)
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.REG', 'L', vr=reg)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', vr=jeumax)
    call jeveuo(parcho//'.NCMP', 'L', vi=vncmp)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(parcho//'.ORIG', 'L', vr=orig)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(iline), b_incx)
! ----------------------------------------------------------------------
! --- CREATION D'UN VECTEUR TEMPORAIRE
! ----------------------------------------------------------------------
    call wkvect('&&MNLINE.VECT1', 'V V R', neq*(2*h+1), ivect1)
    call wkvect('&&MNLINE.VECT2', 'V V R', neq*(2*h+1), ivect2)
! --- COPIE DE XVECT DANS UN VECTEUR AVEC DDLS NON-ACTIFS
    do j = 1, 2*h+1
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(ivect1-1+(j-1)*neq+k) = zr(ivect-1+(j-1)*nd+i)
                if (abs(zr(ivect-1+(j-1)*nd+i)-1.d0) .lt. 1.d-16) then
                end if
            end if
        end do
    end do
! ----------------------------------------------------------------------
! --- EQUATION DE LA DYNAMIQUE (TRAITEMENT DU DEPLACEMENT)
! --- XLINE(1:ND*(2*H+1))=K*XVECT(1:ND*(2*H+1))
! ----------------------------------------------------------------------
    call mrmult('ZERO', imat(1), zr(ivect1), zr(ivect2), 2*h+1, &
                .false._1)
! --- COPIE DE XVECT1 DANS XVECT EN SUPPRIMANT LES DDLS NON ACTIFS
    do j = 1, 2*h+1
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(iline-1+(j-1)*nd+i) = zr(ivect2-1+(j-1)*neq+k)/zr(iadim)
            end if
        end do
    end do
! ----------------------------------------------------------------------
! --- EQUATION DE LA DYNAMIQUE (TRAITEMENT DE LA FORCE NON-LINEAIRE)
! --- XLINE(1:ND*(2*H+1))=XLINE(1:ND*(2*H+1))+F
! ----------------------------------------------------------------------
    neqs = 0
    do i = 1, nchoc
        eta = reg(i)
        jeu = vjeu(i)/jeumax(1)
        ncmp = vncmp(i)
        do j = 1, ncmp
            nddl = vnddl(6*(i-1)+j)
! ---       CSTE & COS
            b_n = to_blas_int(h+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(nd)
            call daxpy(b_n, jeu, zr(ivect+nd*(2*h+1)+(neqs+j-1)*(2*hf+1)), b_incx, &
                       zr(iline-1+nddl), b_incy)
! ---       SIN
            b_n = to_blas_int(h)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(nd)
            call daxpy(b_n, jeu, zr(ivect+nd*(2*h+1)+(neqs+j-1)*(2*hf+1)+hf+1), b_incx, &
                       zr(iline-1+nd*(h+1)+nddl), b_incy)
        end do
        neqs = neqs+vneqs(i)
    end do
! ----------------------------------------------------------------------
! --- EQUATION SUPPLEMENTAIRE POUR DEFINIR LA FORCE NON-LINEAIRE
! --- XLINE(ND*(2*H+1)+1:ND*(2*H+1)+2*NCHOC*(2*HF+1))
! ----------------------------------------------------------------------
    neqs = 0
    do i = 1, nchoc
        alpha = raid(i)/zr(iadim-1+1)
        eta = reg(i)
        jeu = vjeu(i)/jeumax(1)
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
! ---     F -ETA*XG
! ---       -ETA*XG (CSTE & COS)
            b_n = to_blas_int(h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -eta/jeu, zr(ivect-1+nddl), b_incx, &
                       zr(iline-1+nd*(2*h+1)+neqs*(2*hf+1)+1), b_incy)
! ---       -ETA*XG (SIN)
            b_n = to_blas_int(h)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -eta/jeu, zr(ivect-1+nd*(h+1)+nddl), b_incx, &
                       zr(iline-1+nd*(2*h+1)+neqs*(2*hf+1)+hf+2), b_incy)
! ---      + F
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, zr(ivect-1+nd*(2*h+1)+neqs*(2*hf+1)+1), b_incx, &
                       zr(iline-1+nd*(2*h+1)+neqs*(2*hf+1)+1), b_incy)
! ---     Z
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, zr(ivect-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1), b_incx, &
                       zr(iline-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1), b_incy)
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
! ---     + ORIG1*[FN]
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, orig(1+3*(i-1))/jeu, zr(ivect+nd*(2*h+1)+(neqs+3)*(2*hf+1)), b_incx, &
                       zr(iline-1+nd*(2*h+1)+neqs*(2*hf+1)+1), b_incy)
! ---     + ORIG2*[FN]
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, orig(1+3*(i-1)+1)/jeu, zr(ivect+nd*(2*h+1)+(neqs+3)*(2*hf+1)), &
                       b_incx, zr(iline-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1), b_incy)
! ---     + 2*ORIG1*UX + 2*ORIG2*UY (CSTE & COS)
            b_n = to_blas_int(h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 2*orig(1+3*(i-1))/jeu**2, zr(ivect-1+nddlx), b_incx, &
                       zr(iline-1+nd*(2*h+1)+(neqs+2)*(2*hf+1)+1), b_incy)
            b_n = to_blas_int(h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 2*orig(1+3*(i-1)+1)/jeu**2, zr(ivect-1+nddly), b_incx, &
                       zr(iline-1+nd*(2*h+1)+(neqs+2)*(2*hf+1)+1), b_incy)
! ---     + 2*ORIG1*UX + 2*ORIG2*UY (SIN)
            b_n = to_blas_int(h)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 2*orig(1+3*(i-1))/jeu**2, zr(ivect-1+nd*(h+1)+nddlx), b_incx, &
                       zr(iline-1+nd*(2*h+1)+(neqs+2)*(2*hf+1)+hf+2), b_incy)
            b_n = to_blas_int(h)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 2*orig(1+3*(i-1)+1)/jeu**2, zr(ivect-1+nd*(h+1)+nddly), b_incx, &
                       zr(iline-1+nd*(2*h+1)+(neqs+2)*(2*hf+1)+hf+2), b_incy)
! ---     FN
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect+nd*(2*h+1)+(neqs+3)*(2*hf+1)), b_incx, &
                       zr(iline+nd*(2*h+1)+(neqs+3)*(2*hf+1)), b_incy)
        else if (type(i) (1:4) .eq. 'PLAN') then
! ---     F
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)), b_incx, &
                       zr(iline+nd*(2*h+1)+neqs*(2*hf+1)), b_incy)
        end if
        neqs = neqs+vneqs(i)
    end do
! ----------------------------------------------------------------------
! --- AUTRES EQUATIONS
! ----------------------------------------------------------------------
! --- GAMMA1
    zr(iline-1+ninc-3) = zr(ivect-1+ninc-3)
! --- GAMMA2
    zr(iline-1+ninc-2) = zr(ivect-1+ninc-2)
! --- EQUATION DE PHASE
    zr(iline-1+ninc-1) = 0.d0
! ----------------------------------------------------------------------
! --- DESTRUCTION DES VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call jedetr('&&MNLINE.VECT1')
    call jedetr('&&MNLINE.VECT2')
!
    call jedema()
!
end subroutine
