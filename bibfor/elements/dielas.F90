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
subroutine dielas(for_discret, iret)
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
type(te0047_dscr), intent(in) :: for_discret
integer, intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer :: imat, jdc, irep, neq, ii, ifono, icontp, icontm
    real(kind=8) :: r8bid, klv(78), klc(144), fl(12)
    character(len=8) :: k8bid
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
        if (for_discret%ndim .eq. 3) then
            call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        else if (for_discret%ndim.eq.2) then
            call ut2mgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        endif
    else
        call dcopy(for_discret%nbt, zr(jdc), 1, klv, 1)
    endif
!   calcul de la matrice tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim.eq.2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        endif
    endif
    neq = for_discret%nno*for_discret%nc
    !
    if ( for_discret%lVect .or. for_discret%lSigm ) then
        ! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
        ! calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    endif
    ! calcul des efforts généralisés
    if ( for_discret%lSigm ) then
        call jevech('PCONTPR', 'E', icontp)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii) + zr(icontm-1+ii)
            enddo
        else if (for_discret%nno.eq.2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii)                = -fl(ii)    + zr(icontm-1+ii)
                zr(icontp-1+ii+for_discret%nc) =  fl(ii+for_discret%nc) + &
                                                  zr(icontm-1+ii+for_discret%nc)
            enddo
        endif
    endif
    ! calcul des forces nodales
    if ( for_discret%lVect ) then
        call jevech('PVECTUR', 'E', ifono)
        ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                fl(ii)          = fl(ii) + zr(icontm-1+ii)
            enddo
        else if (for_discret%nno.eq.2) then
            do ii = 1, for_discret%nc
                fl(ii)                         =  fl(ii)    - zr(icontm-1+ii)
                fl(ii+for_discret%nc)          =  fl(ii+for_discret%nc) +  &
                                                  zr(icontm-1+ii+for_discret%nc)
            enddo
        endif
        ! forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        endif
    endif
end subroutine
