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
! person_in_charge: jean-luc.flejou at edf.fr
!
subroutine dielas(DD, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    DD : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: DD
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, jdc, irep, neq, ii, jj, ifono, icontp, icontm
    real(kind=8) :: r8bid, klv(78), klc(12, 12), fl(12), dulth(12)
    character(len=8) :: k8bid
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PCONTMR', 'L', icontm)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   absolu vers local ?
!   irep = 1 . la matrice est dans le repère global ==> passer en local
    if (irep .eq. 1) then
        if (DD%ndim .eq. 3) then
            call utpsgl(DD%nno, DD%nc, DD%pgl, zr(jdc), klv)
        else if (DD%ndim .eq. 2) then
            call ut2mgl(DD%nno, DD%nc, DD%pgl, zr(jdc), klv)
        end if
    else
        b_n = to_blas_int(DD%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
!   calcul de la matrice tangente
    if (DD%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (DD%ndim .eq. 3) then
            call utpslg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
        else if (DD%ndim .eq. 2) then
            call ut2mlg(DD%nno, DD%nc, DD%pgl, klv, zr(imat))
        end if
    end if
    neq = DD%nno*DD%nc
    !
    if (DD%lVect .or. DD%lSigm) then
! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, DD%nbt, klc, neq)
! calcul de fl = klc.dul (incrément d'effort)
!   Correction thermique : DeltaDilatation = alpha*(T+ - T-)*xl
        fl(:) = 0.0
        if (abs(DD%DeltaDilatation) .gt. r8prem()) then
            dulth(:) = 0.0
            dulth(1) = -DD%DeltaDilatation
            dulth(DD%nc+1) = -dulth(1)
            do ii = 1, DD%nc
                do jj = 1, DD%nc
                    fl(ii) = fl(ii)-klc(ii, jj)*dulth(jj)
                    fl(ii+DD%nc) = fl(ii+DD%nc)-klc(ii+DD%nc, jj+DD%nc)*dulth(jj+DD%nc)
                end do
            end do
        end if
        call pmavec('CUMU', neq, klc, DD%dul, fl)
    end if
! calcul des efforts généralisés
    if (DD%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (DD%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (DD%nno .eq. 2) then
            do ii = 1, DD%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+DD%nc) = fl(ii+DD%nc)+zr(icontm-1+ii+DD%nc)
            end do
        end if
    end if
! calcul des forces nodales
    if (DD%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (DD%nno .eq. 1) then
            do ii = 1, neq
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (DD%nno .eq. 2) then
            do ii = 1, DD%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+DD%nc) = fl(ii+DD%nc)+zr(icontm-1+ii+DD%nc)
            end do
        end if
! forces nodales aux noeuds 1 et 2 (repère global)
        if (DD%nc .ne. 2) then
            call utpvlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
        else
            call ut2vlg(DD%nno, DD%nc, DD%pgl, fl, zr(ifono))
        end if
    end if
end subroutine
