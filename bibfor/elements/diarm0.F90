! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine diarm0(for_discret, iret)
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
#include "asterfort/diarme.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
type(te0047_dscr), intent(in) :: for_discret
integer, intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jdc, irep, imat, ivarim, ifono, icontp, ivarip, icontm, neq, ii
    real(kind=8) :: r8bid, klv(78), force(3), klc(144), fl(12), duly, ulp(12), varip
    character(len=8) :: k8bid
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    ! param??tres en entr??e
    call jevech('PCADISK', 'L', jdc)
    call infdis('REPK', irep, r8bid, k8bid)
    call dcopy(for_discret%nbt, zr(jdc), 1, klv, 1)
    if (irep .eq. 1) then
        call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
    endif
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PCONTMR', 'L', icontm)
!   relation de comportement de l'armement
    force(:) = 0.0d0
    ulp(1:12) = for_discret%ulm(1:12) + for_discret%dul(1:12)
!
    call diarme(for_discret%nbt, neq, zi(imat), for_discret%ulm, for_discret%dul,&
                ulp, zr(icontm), zr(ivarim), klv, varip,&
                force(1), duly)
!   actualisation de la matrice tangente
    if ( for_discret%lMatr ) then
        call jevech('PMATUUR', 'E', imat)
        call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
    endif
    neq = for_discret%nno*for_discret%nc
    !
    if ( for_discret%lVect .or. for_discret%lSigm) then
!       demi-matrice klv transform??e en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
!       calcul de fl = klc.dul (incr??ment d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    endif
    ! calcul des efforts g??n??ralis??s, des forces nodales et des variables internes
    if ( for_discret%lSigm ) then
        call jevech('PCONTPR', 'E', icontp)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii) + zr(icontm-1+ii)
            enddo
        else if (for_discret%nno.eq.2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii)                = -fl(ii) + zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) =  fl(ii+for_discret%nc) + &
                                                  zr(icontm-1+ii+for_discret%nc)
            enddo
        endif
!       modif pour les armements
        zr(icontp-1+2) = zr(icontm-1+2) + force(1)*duly
        zr(icontp-1+8) = zr(icontm-1+8) + force(1)*duly
    endif
    ! calcul des forces nodales
    if ( for_discret%lVect ) then
        call jevech('PVECTUR', 'E', ifono)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fl(ii) = fl(ii) + zr(icontm-1+ii)
            enddo
        else if (for_discret%nno.eq.2) then
            do ii = 1, for_discret%nc
                fl(ii)                = fl(ii) - zr(icontm-1+ii)
                fl(ii+for_discret%nc) = fl(ii+for_discret%nc) + &
                                        zr(icontm-1+ii+for_discret%nc)
            enddo
        endif
!       modif pour les armements
        fl(2) = -zr(icontm-1+2) - force(1)*duly
        fl(8) =  zr(icontm-1+8) + force(1)*duly
!       forces nodales aux noeuds 1 et 2 (rep??re global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        endif
    endif
    ! mise ?? jour des variables internes
    if ( for_discret%lVari ) then
        call jevech('PVARIPR', 'E', ivarip)
        zr(ivarip)   = varip
        zr(ivarip+1) = varip
    endif
end subroutine
