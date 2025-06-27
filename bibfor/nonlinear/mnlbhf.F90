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
subroutine mnlbhf(xvect, parcho, adime, ninc, nd, &
                  nchoc, h, hf, err)
    implicit none
!
!
!       MODE NON LINEAIRE - VERIFICATION BONNE VALEUR HF
!       -         -        -             --           -
! ----------------------------------------------------------------------
!
! INITIALISATION DU POINT DE DEPART DE LA MAN
! ----------------------------------------------------------------------
! IN   XVECT  : K14  : NOM DU VECTEUR A INITIALISER
! IN   PARCHO : K14  : SD PARAMETRE DES CONTACTEURS
! IN   ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN   NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN   ND     : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN   NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN   H      : I    : NOMBRE D'HARMONIQUES POUR X
! IN   HF     : I    : NOMBRE D'HARMONIQUES POUR F
! OUT  ERR    : R8   : ERREUR ENTRE LES VALEURS DE REF. DE F ET CALC.
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlbil.h"
#include "asterfort/mnlcir.h"
#include "asterfort/mnluil.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer(kind=8) :: ninc, nd, nchoc, h, hf
    character(len=14) :: parcho, adime, xvect
    real(kind=8) :: err
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    character(len=14) :: xdep1, xdep2, xtemp
    real(kind=8) :: omega, alpha, eta, jeu
    real(kind=8) :: nrm
    integer(kind=8) :: ivect
    integer(kind=8) :: iadim, i, j
    integer(kind=8) :: neqs, nddlx, nddly
    integer(kind=8) :: nt, nddl, idep1, idep2, itemp
    real(kind=8), pointer :: tep2(:) => null()
    integer(kind=8), pointer :: vneqs(:) => null()
    real(kind=8), pointer :: jeumax(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: orig(:) => null()
    real(kind=8), pointer :: reg(:) => null()
    character(len=8), pointer :: type(:) => null()
    integer(kind=8), pointer :: vnddl(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
! ----------------------------------------------------------------------
! --- QUELQUES VALEURS UTILES
! ----------------------------------------------------------------------
    call jeveuo(adime, 'L', iadim)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.REG', 'L', vr=reg)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', vr=jeumax)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.ORIG', 'L', vr=orig)
    call jeveuo(xvect, 'L', ivect)
    omega = zr(ivect-1+ninc)
! ----------------------------------------------------------------------
! --- REMPLISSAGE DES FORCES DE CHOCS ET DES VECTEURS AUXILIAIRES
! ----------------------------------------------------------------------
    xdep1 = '&&MNLBHF.DEP1'
    xdep2 = '&&MNLBHF.DEP2'
    xtemp = '&&MNLBHF.TEMP'
    call wkvect(xdep1, 'V V R', 2*h+1, idep1)
    call wkvect(xdep2, 'V V R', 2*h+1, idep2)
    AS_ALLOCATE(vr=tep2, size=2*h+1)
    call wkvect(xtemp, 'V V R', ninc, itemp)
    nt = int(2**int(dlog(2.d0*dble(3*hf)+1.d0)/dlog(2.d0)+1.d0))
    neqs = 0
    err = 0.d0
    do i = 1, nchoc
! ---   ON RECUPERE LES PARAMETRES DE CHOCS
        alpha = raid(i)/zr(iadim)
        eta = reg(i)
        jeu = vjeu(i)/jeumax(1)
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0/jeu, zr(ivect-1+nddl), b_incx, zr(idep1), &
                       b_incy)
            call mnlbil(zr(idep1), omega, alpha, eta, h, &
                        hf, nt, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)))
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep2), b_incx)
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itemp), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect-1+nddlx), b_incx, zr(idep1), b_incy)
            zr(idep1) = zr(idep1)-orig(3*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/jeu, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect-1+nddly), b_incx, zr(idep2), b_incy)
            zr(idep2) = zr(idep2)-orig(3*(i-1)+2)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/jeu, zr(idep2), b_incx)
!
            call mnlcir(xdep1, xdep2, omega, alpha, eta, &
                        h, hf, nt, xtemp)
!
            b_n = to_blas_int(4*(2*hf+1))
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)), b_incx, zr(itemp), &
                       b_incy)
!
            b_n = to_blas_int(4*(2*hf+1))
            b_incx = to_blas_int(1)
            nrm = dnrm2(b_n, zr(itemp), b_incx)
            if (nrm .gt. 0.d0) then
                nrm = 0.d0
                do j = 1, 2
                    b_n = to_blas_int(2*h+1)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, 0.d0, tep2, b_incx)
                    b_n = to_blas_int(h+1)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(itemp+(j-1)*(2*hf+1)), b_incx, tep2, b_incy)
                    b_n = to_blas_int(h)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call dcopy(b_n, zr(itemp+(j-1)*(2*hf+1)+hf+1), b_incx, tep2, b_incy)
                    b_n = to_blas_int(2*h+1)
                    b_incx = to_blas_int(1)
                    nrm = nrm+dnrm2(b_n, tep2, b_incx)
                end do
                err = err+nrm/2.d0
            end if
        else if (type(i) (1:4) .eq. 'PLAN') then
            nddl = vnddl(6*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itemp), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0/jeu, zr(ivect-1+nddl), b_incx, zr(idep1), &
                       b_incy)
            call mnluil(zr(idep1), omega, alpha, eta, h, &
                        hf, nt, zr(itemp))
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)), b_incx, zr(itemp), &
                       b_incy)
            b_n = to_blas_int(2*hf+1)
            b_incx = to_blas_int(1)
            nrm = dnrm2(b_n, zr(itemp), b_incx)
            if (nrm .gt. 0.d0) then
                b_n = to_blas_int(2*h+1)
                b_incx = to_blas_int(1)
                call dscal(b_n, 0.d0, tep2, b_incx)
                b_n = to_blas_int(h+1)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(itemp), b_incx, tep2, b_incy)
                b_n = to_blas_int(h)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, zr(itemp+hf+1), b_incx, tep2, b_incy)
                b_n = to_blas_int(2*h+1)
                b_incx = to_blas_int(1)
                err = err+dnrm2(b_n, tep2, b_incx)
            end if
        end if
        neqs = neqs+vneqs(i)
    end do
!
    err = err/dble(nchoc)
!
    call jedetr(xtemp)
    AS_DEALLOCATE(vr=tep2)
    call jedetr(xdep1)
    call jedetr(xdep2)
!
    call jedema()
!
end subroutine
