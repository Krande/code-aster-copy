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
subroutine dicho0(for_discret, iret)
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
#include "asterc/r8prem.h"
#include "asterfort/dichoc.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
!
type(te0047_dscr), intent(in) :: for_discret
integer, intent(out)          :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jdc, irep, imat, ivarim, ii, jinst, ivitp, idepen, iviten, neq, igeom, ivarip
    integer :: iretlc, ifono
    integer :: icontm, icontp
!
    real(kind=8)    :: klv(78), varmo(8), varpl(8)
    real(kind=8)    :: force(3)
    real(kind=8)    :: klc(144), fl(12), dvl(12), dpe(12), dve(12), ulp(12)
!
    real(kind=8)        :: r8bid, duly
    character(len=8)    :: k8bid
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTMR', 'L', icontm)
!
    call infdis('REPK', irep, r8bid, k8bid)
!   absolu vers local ? ---
!   irep = 1 = matrice en repère global ==> passer en local ---
    if (irep .eq. 1) then
        if (for_discret%ndim .eq. 3) then
            call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        else if (for_discret%ndim.eq.2) then
            call ut2mgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        endif
    else
        call dcopy(for_discret%nbt, zr(jdc), 1, klv, 1)
    endif
    call jevech('PMATERC', 'L', imat)
    call jevech('PVARIMR', 'L', ivarim)
!
    do ii = 1, 8
        varmo(ii) = zr(ivarim+ii-1)
        varpl(ii) = 0.
    enddo
!
    call jevech('PINSTPR', 'L', jinst)
!
    call tecach('ONO', 'PVITPLU', 'L', iretlc, iad=ivitp)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        else if (for_discret%ndim.eq.2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(ivitp), dvl)
        endif
    else
        dvl(:)    = 0.0
    endif
!
    call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        else if (for_discret%ndim.eq.2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(idepen), dpe)
        endif
    else
        dpe(:) = 0.0d0
    endif
!
    call tecach('ONO', 'PVITENT', 'L', iretlc, iad=iviten)
    if (iretlc .eq. 0) then
        if (for_discret%ndim .eq. 3) then
            call utpvgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        else if (for_discret%ndim.eq.2) then
            call ut2vgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(iviten), dve)
        endif
    else
        dve(:) = 0.d0
    endif
!
    neq = for_discret%nno*for_discret%nc
    ulp(:) = for_discret%ulm(:) + for_discret%dul(:)
!   relation de comportement de choc
    call dichoc(for_discret%nbt, neq, for_discret%nno, for_discret%nc, zi(imat),&
                for_discret%dul, ulp, zr(igeom), for_discret%pgl, klv,&
                duly, dvl, dpe, dve, force,&
                varmo, varpl, for_discret%ndim)
!   actualisation de la matrice tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim.eq.2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        endif
    endif
    !
    if ( for_discret%lVect .or. for_discret%lSigm ) then
        ! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
        ! calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    endif
!   calcul des efforts généralisés
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
        if (for_discret%nno .eq. 1) then
            zr(icontp-1+1) = force(1)
            zr(icontp-1+2) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3) = force(3)
            endif
        else if (for_discret%nno.eq.2) then
            zr(icontp-1+1)                = force(1)
            zr(icontp-1+1+for_discret%nc) = force(1)
            zr(icontp-1+2)                = force(2)
            zr(icontp-1+2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                zr(icontp-1+3)               = force(3)
                zr(icontp-1+3+for_discret%nc) = force(3)
            endif
        endif
        if (abs(force(1)) .lt. r8prem()) then
            do ii = 1, neq
                zr(icontp-1+ii) = 0.0d0
            enddo
        endif
    endif
    ! calcul des efforts généralisés et des forces nodales
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
        if (for_discret%nno .eq. 1) then
            fl(1) = force(1)
            fl(2) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3) = force(3)
            endif
        else if (for_discret%nno.eq.2) then
            fl(1)                = -force(1)
            fl(1+for_discret%nc) = force(1)
            fl(2)                = -force(2)
            fl(2+for_discret%nc) = force(2)
            if (for_discret%ndim .eq. 3) then
                fl(3)                = -force(3)
                fl(3+for_discret%nc) =  force(3)
            endif
        endif
        if (abs(force(1)) .lt. r8prem()) then
            do ii = 1, neq
                fl(ii) = 0.0d0
            enddo
        endif
!       forces nodales aux noeuds 1 et 2 (repère global)
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        endif
    endif
    ! mise à jour des variables internes
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        do ii = 1, 8
            zr(ivarip+ii-1) = varpl(ii)
            if (for_discret%nno .eq. 2) zr(ivarip+ii+7) = varpl(ii)
        enddo
    endif
end subroutine
