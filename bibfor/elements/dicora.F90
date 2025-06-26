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
subroutine dicora(for_discret, iret)
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
#include "jeveux.h"
#include "asterfort/dicorn.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iiter, imat, ivarim, neq, iterat, ii, ifono, icontp, ivarip, icontm, irmetg
    real(kind=8) :: ulp(12), klv(78), klv2(78), varipc(7), klc(144), fl(12)
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PITERAT', 'L', iiter)
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
!
    neq = for_discret%nno*for_discret%nc
    ulp(1:12) = for_discret%ulm(1:12)+for_discret%dul(1:12)
!   relation de comportement de la cornière
    irmetg = 0
    if (for_discret%lMatrPred) irmetg = 1
    iterat = zi(iiter)
    call dicorn(irmetg, for_discret%nbt, neq, iterat, zi(imat), &
                for_discret%ulm, for_discret%dul, ulp, zr(icontm), zr(ivarim), &
                klv, klv2, varipc)
!   actualisation de la matrice tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv2, zr(imat))
    end if
!
    if (for_discret%lVect .or. for_discret%lSigm) then
!       demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
!       calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    end if

!   calcul des efforts généralisés et des forces nodales
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)+zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)+ &
                                                 zr(icontm-1+ii+for_discret%nc)
            end do
        end if
    end if
    ! calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fl(ii) = fl(ii)+zr(icontm-1+ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                fl(ii) = fl(ii)-zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+ &
                                        zr(icontm-1+ii+for_discret%nc)
            end do
        end if
!       forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
!   mise à jour des variables internes
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, 7
            zr(ivarip+ii-1) = varipc(ii)
            zr(ivarip+ii+6) = varipc(ii)
        end do
    end if
end subroutine
