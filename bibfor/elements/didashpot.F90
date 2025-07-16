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
subroutine didashpot(for_discret, iret)
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
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
#include "blas/dcopy.h"
#include "asterfort/Behaviour_type.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imat, jdc, irep, neq, ii, ifono, icontp, icontm
    real(kind=8) :: r8bid, klv(78), klc(144), fl(12)
    character(len=16), pointer :: compor(:) => null()
    character(len=8) :: k8bid
!
    integer(kind=8) :: iadzi, iazk24
    character(len=24) :: messak(5)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
! Seulement en 3D et Translation
    if (for_discret%nomte(1:11) .ne. 'MECA_DIS_T_') then
        call jevech('PCOMPOR', 'L', vk16=compor)
        messak(1) = for_discret%nomte
        messak(2) = compor(INCRELAS)
        messak(3) = compor(RELA_NAME)
        call tecael(iadzi, iazk24)
        messak(4) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_23', nk=4, valk=messak)
    end if
    !
! paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PCONTMR', 'L', icontm)
    !
    call infdis('REPK', irep, r8bid, k8bid)
! absolu vers local ?
! irep = 1 . la matrice est dans le repère global ==> passer en local
    if (irep .eq. 1) then
        if (for_discret%ndim .eq. 3) then
            call utpsgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        else if (for_discret%ndim .eq. 2) then
            call ut2mgl(for_discret%nno, for_discret%nc, for_discret%pgl, zr(jdc), klv)
        end if
    else
        b_n = to_blas_int(for_discret%nbt)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)
    end if
! calcul de la matrice tangente
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        end if
    end if
    neq = for_discret%nno*for_discret%nc
    !
    if (for_discret%lVect .or. for_discret%lSigm) then
! demi-matrice klv transformée en matrice pleine klc
        call vecma(klv, for_discret%nbt, klc, neq)
! calcul de fl = klc.dul (incrément d'effort)
        call pmavec('ZERO', neq, klc, for_discret%dul, fl)
    end if
! calcul des efforts généralisés et des forces nodales
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nno .eq. 1) then
            do ii = 1, neq
                zr(icontp-1+ii) = fl(ii)
            end do
        else if (for_discret%nno .eq. 2) then
            do ii = 1, for_discret%nc
                zr(icontp-1+ii) = -fl(ii)
                zr(icontp-1+ii+for_discret%nc) = fl(ii+for_discret%nc)
            end do
        end if
    end if
! calcul des forces nodales
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_TR_L et MECA_DIS_T_L
        if (for_discret%nc .ne. 2) then
            call utpvlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        else
            call ut2vlg(for_discret%nno, for_discret%nc, for_discret%pgl, fl, zr(ifono))
        end if
    end if
end subroutine
