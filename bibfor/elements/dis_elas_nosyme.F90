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
subroutine dis_elas_nosyme(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut3mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, jdc, irep, neq, ii, ifono, icontp, icontm
    real(kind=8) :: r8bid, dfg(12), dfl(12), fgm(12)
    character(len=8) :: k8bid
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PCONTMR', 'L', icontm)
!   La matrice de raideur est donnée dans quel repère ? 1:global, 2:local
    call infdis('REPK', irep, r8bid, k8bid)
!
    neq = for_discret%nno*for_discret%nc
!   Traitement du cas : repère GLOBAL
    if (irep .eq. 1) then
! Matrice tangente, dans le repère global
        if (for_discret%lMatr) then
            call jevech('PMATUNS', 'E', imat)
            b_n = to_blas_int(for_discret%nbt)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jdc), b_incx, zr(imat), b_incy)
        end if
! calcul de df = K.du (incrément d'effort, dans le repère global)
        if (for_discret%lVect .or. for_discret%lSigm) then
            call pmavec('ZERO', neq, zr(jdc), for_discret%dug, dfg)
        end if
! calcul des efforts généralisés. Ils doivent être dans le repère local.
        if (for_discret%lSigm) then
            call jevech('PCONTPR', 'E', icontp)
! Passage des incréments d'efforts du repère global vers local
            if (for_discret%ndim .eq. 3) then
                call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, dfg, dfl)
            else
                call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, dfg, dfl)
            end if
! Attention aux signes des efforts sur le 1er noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
            if (for_discret%nno .eq. 1) then
                do ii = 1, for_discret%nc
                    zr(icontp-1+ii) = dfl(ii)+zr(icontm-1+ii)
                end do
            else if (for_discret%nno .eq. 2) then
                do ii = 1, for_discret%nc
                    zr(icontp-1+ii) = -dfl(ii)+zr(icontm-1+ii)
                    zr(icontp-1+ii+for_discret%nc) = dfl(ii+for_discret%nc)+zr(icontm-1+ii+for_d&
                                                     &iscret%nc)
                end do
            end if
        end if
! calcul des forces nodales. Elles doivent être dans le repère global.
        if (for_discret%lVect) then
            call jevech('PVECTUR', 'E', ifono)
! Passage des forces nodales (t-) du repère local vers global
            if (for_discret%ndim .eq. 3) then
                call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, zr(icontm), fgm)
            else
                call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, zr(icontm), fgm)
            end if
! Attention aux signes des efforts sur le 1er noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
            if (for_discret%nno .eq. 1) then
                do ii = 1, for_discret%nc
                    zr(ifono-1+ii) = dfg(ii)+fgm(ii)
                end do
            else if (for_discret%nno .eq. 2) then
                do ii = 1, for_discret%nc
                    zr(ifono-1+ii) = dfg(ii)-fgm(ii)
                    zr(ifono-1+ii+for_discret%nc) = dfg(ii+for_discret%nc)+fgm(ii+for_discret%nc)
                end do
            end if
        end if
!
!   Traitement du cas : repère LOCAL
    else
! Matrice tangente, dans le repère global
        if (for_discret%lMatr) then
            call jevech('PMATUNS', 'E', imat)
            if (for_discret%ndim .eq. 3) then
                call ut3mlg(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), zr(imat))
            else
                ASSERT(.false.)
            end if
        end if
! calcul de df = K.du (incrément d'effort, dans le repère local)
        if (for_discret%lVect .or. for_discret%lSigm) then
            call pmavec('ZERO', neq, zr(jdc), for_discret%dul, dfl)
        end if
! calcul des efforts généralisés. Ils doivent être dans le repère local.
        if (for_discret%lSigm) then
            call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le 1er noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
            if (for_discret%nno .eq. 1) then
                do ii = 1, for_discret%nc
                    zr(icontp-1+ii) = dfl(ii)+zr(icontm-1+ii)
                end do
            else if (for_discret%nno .eq. 2) then
                do ii = 1, for_discret%nc
                    zr(icontp-1+ii) = -dfl(ii)+zr(icontm-1+ii)
                    zr(icontp-1+ii+for_discret%nc) = dfl(ii+for_discret%nc)+zr(icontm-1+ii+for_d&
                                                     &iscret%nc)
                end do
            end if
        end if
! calcul des forces nodales. Elles sont dans le repère global.
        if (for_discret%lVect) then
            call jevech('PVECTUR', 'E', ifono)
! Passage des forces nodales (t-) du repère local vers global
            if (for_discret%ndim .eq. 3) then
                call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, zr(icontm), fgm)
                call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, dfl, dfg)
            else
                call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, zr(icontm), fgm)
                call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, dfl, dfg)
            end if
! Attention aux signes des efforts sur le 1er noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
            if (for_discret%nno .eq. 1) then
                do ii = 1, for_discret%nc
                    zr(ifono-1+ii) = dfg(ii)+fgm(ii)
                end do
            else if (for_discret%nno .eq. 2) then
                do ii = 1, for_discret%nc
                    zr(ifono-1+ii) = dfg(ii)-fgm(ii)
                    zr(ifono-1+ii+for_discret%nc) = dfg(ii+for_discret%nc)+fgm(ii+for_discret%nc)
                end do
            end if
        end if
    end if
!
end subroutine
