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

subroutine discret_sief(for_discret, klv, dul, sim, ilogic, sip, fono, force)
!
    use te0047_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8)         :: ilogic
    real(kind=8)    :: klv(*), dul(*), sim(*)
    real(kind=8)    :: sip(*), fono(*), force(*)
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES EFFORTS GÉNÉRALISÉS (REPÈRE LOCAL)
!     ET DES FORCES NODALES (REPÈRE GLOBAL). COMME ON TRAITE DES
!     ÉLÉMENTS DISCRETS, CES QUANTITÉS SONT ÉGALES, AU REPÈRE PRÈS.
!
! --------------------------------------------------------------------------------------------------
!
! IN
!       klv    : matrice de "raideur tangente"
!       dul    : incrément de déplacement local
!       sim    : efforts généralisés a l'instant précédent
!       ilogic :
!       duly   :
!
! OUT
!       sip    : efforts generalises actualises
!       fono   : forces nodales
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8)      :: ii, neq
    real(kind=8) :: klc(144), fl(12)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT((ilogic .eq. 0) .or. (ilogic .eq. 2))
!
    neq = for_discret%nno*for_discret%nc
!   demi-matrice klv transformée en matrice pleine klc
    call vecma(klv, for_discret%nbt, klc, neq)
!   calcul de fl = klc.dul (incrément d'effort)
    call pmavec('ZERO', neq, klc, dul, fl)
!   efforts généralisés aux noeuds 1 et 2 (repère local)
!   on change le signe des efforts sur le premier noeud pour les MECA_DIS_TR_L et MECA_DIS_T_L
    if (for_discret%nno .eq. 1) then
        do ii = 1, neq
            sip(ii) = fl(ii)+sim(ii)
            fl(ii) = fl(ii)+sim(ii)
        end do
    else if (for_discret%nno .eq. 2) then
        do ii = 1, for_discret%nc
            sip(ii) = -fl(ii)+sim(ii)
            sip(ii+for_discret%nc) = fl(ii+for_discret%nc)+sim(ii+for_discret%nc)
            fl(ii) = fl(ii)-sim(ii)
            fl(ii+for_discret%nc) = fl(ii+for_discret%nc)+sim(ii+for_discret%nc)
        end do
    end if
!
    if (ilogic .eq. 2) then
        if (for_discret%nno .eq. 1) then
            sip(1) = force(1)
            sip(2) = force(2)
            fl(1) = force(1)
            fl(2) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = force(3)
                sip(3) = force(3)
            end if
        else if (for_discret%nno .eq. 2) then
            sip(1) = force(1)
            sip(2) = force(2)
            sip(1+for_discret%nc) = force(1)
            sip(2+for_discret%nc) = force(2)
            fl(1) = -force(1)
            fl(2) = -force(2)
            fl(1+for_discret%nc) = force(1)
            fl(2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                sip(3) = force(3)
                sip(3+for_discret%nc) = force(3)
                fl(3) = -force(3)
                fl(3+for_discret%nc) = force(3)
            end if
        end if
        if (abs(force(1)) .lt. r8prem()) then
            fl(1:neq) = 0.0
            sip(1:neq) = 0.0
        end if
    end if
!
!   forces nodales aux noeuds 1 et 2 (repère global)
    if (for_discret%ndim .eq. 3) then
        call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, fono)
    else
        call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, fono)
    end if
end subroutine
