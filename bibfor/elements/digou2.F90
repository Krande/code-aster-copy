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
subroutine digou2(for_discret, iret)
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
#include "asterfort/digouj.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpslg.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdc, irep, imat, ivarim, ifono, icontp, ivarip, icontm, neq
    real(kind=8) :: r8bid, klv(78), klv2(78)
    character(len=8) :: k8bid
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!   paramètres en entrée
    call jevech('PCADISK', 'L', jdc)
    call jevech('PMATERC', 'L', imat)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
!   matrice de rigidité en repère local
    call infdis('REPK', irep, r8bid, k8bid)
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
!
    ifono = 1
    icontp = 1
    ivarip = 1
    if (for_discret%lVect) then
        call jevech('PVECTUR', 'E', ifono)
    end if
    if (for_discret%lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (for_discret%lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
!   relation de comportement : élastique partout
!   sauf suivant Y local : élasto-plastique VMIS_ISOT_TRAC
    neq = for_discret%nno*for_discret%nc
    call digouj(for_discret%option, for_discret%rela_comp, for_discret%nno, for_discret%nbt, neq, &
                for_discret%nc, zi(imat), for_discret%dul, zr(icontm), zr(ivarim), &
                for_discret%pgl, klv, klv2, zr(ivarip), zr(ifono), &
                zr(icontp), for_discret%nomte)
!
    if (for_discret%lMatr) then
        call jevech('PMATUUR', 'E', imat)
        if (for_discret%ndim .eq. 3) then
            call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        else if (for_discret%ndim .eq. 2) then
            call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klv, zr(imat))
        end if
    end if
!
end subroutine
