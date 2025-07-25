!    -----------------------------------------------------------------
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
!    -----------------------------------------------------------------
!
subroutine check_aniso(propname, objname)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dortvp.h"
#include "asterfort/indk16.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: propname, objname
!
!    -----------------------------------------------------------------------------------------------
!
! DEFI_MATERIAU
!
! Compute eigenvalues for Hooke matrix (check stability)
!
!    -----------------------------------------------------------------------------------------------
!
! In  mate             : name of output datastructure
! In  v_mate_func      : pointer to flags for function
!
!    -----------------------------------------------------------------------------------------------
!
    character(len=2) :: m2blan, k8bid
    character(len=16) :: nomrc
    character(len=19) :: noobrc
    real(kind=8) :: dorth(6, 6)
    real(kind=8) :: nu12, nu21, nu13, nu31, nu23, nu32
    integer(kind=8) :: iel, ien, iet, igln, iglt, igtn
    integer(kind=8) :: inuln, inult, inutn
    integer(kind=8) ::  nbr, ndim
    real(kind=8) :: c1, delta, deux, e1, e2, e3, g12
    real(kind=8) :: g13, g23, un, undemi, zero
    real(kind=8), pointer :: valr(:) => null()
    character(len=16), pointer :: valk(:) => null()
!
!    -----------------------------------------------------------------------------------------------
    nomrc = propname
    noobrc = objname
!
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    deux = 2.0d0
!
    m2blan = ' '
    k8bid = ' '
!
    dorth(:, :) = zero
!
    e1 = zero
    e2 = zero
    e3 = zero
    g12 = zero
    g23 = zero
    g13 = zero
    nu12 = zero
    nu23 = zero
    nu13 = zero
!
    call jemarq()
!
!   on ne traite que les cas isotrope-transverse et orthotrope
    ASSERT(nomrc .eq. 'ELAS_ISTR' .or. nomrc .eq. 'ELAS_ORTH')
!
!   recuperation du nom des composantes et des valeurs
!   definissant le materiau :
    call jeveuo(noobrc//'.VALR', 'L', vr=valr)
    call jeveuo(noobrc//'.VALK', 'L', vk16=valk)
!
!   longueur du tableau des composantes :
    call jelira(noobrc//'.VALR', 'LONUTI', nbr)
!
!   recuperation des indices des composantes relatives
!   a l'orthotropie et a l'isotropie transverse dans
!   le tableau du nom des composantes :
    iel = indk16(valk, 'E_L', 1, nbr)
    iet = indk16(valk, 'E_T', 1, nbr)
    ien = indk16(valk, 'E_N', 1, nbr)
!
    iglt = indk16(valk, 'G_LT', 1, nbr)
    igtn = indk16(valk, 'G_TN', 1, nbr)
    igln = indk16(valk, 'G_LN', 1, nbr)
!
    inult = indk16(valk, 'NU_LT', 1, nbr)
    inutn = indk16(valk, 'NU_TN', 1, nbr)
    inuln = indk16(valk, 'NU_LN', 1, nbr)
!
!   recuperation des composantes relatives a l'orthotropie
!   et a l'isotropie transverse :
    if (iel .ne. 0) e1 = valr(iel)
    if (iet .ne. 0) e2 = valr(iet)
    if (ien .ne. 0) e3 = valr(ien)
!
    if (iglt .ne. 0) g12 = valr(iglt)
    if (igtn .ne. 0) g23 = valr(igtn)
    if (igln .ne. 0) g13 = valr(igln)
!
    if (inult .ne. 0) nu12 = valr(inult)
    if (inutn .ne. 0) nu23 = valr(inutn)
    if (inuln .ne. 0) nu13 = valr(inuln)
!
!   traitement du cas de l'isotropie transverse :
    if (nomrc .eq. 'ELAS_ISTR') then
!
!       si g13 = 0 , on peut supposer que l'on est en 2d
!       on ne traite que le cas deformations planes ou
!       axisymetrique car le cas contraintes planes revient
!       a l'elasticite isotrope :
!       (comme on n'a accès qu'au matériau, les différents cas
!       2d, 3d... ne peuvent être discriminés qu'à partir de ce que
!       l'utilisateur a spécifié)
        if (igln .eq. 0) then
            ndim = 2
            if (ien .eq. 0 .or. e3 .le. r8prem()) then
                goto 999
            end if
!
            c1 = e1/(un+nu12)
            delta = un-nu12-deux*nu13*nu13*e3/e1
!
            dorth(1, 1) = c1*(un-nu13*nu13*e3/e1)/delta
            dorth(1, 2) = c1*((un-nu13*nu13*e3/e1)/delta-un)
            dorth(1, 3) = e3*nu13/delta
            dorth(2, 1) = dorth(1, 2)
            dorth(2, 2) = dorth(1, 1)
            dorth(2, 3) = dorth(1, 3)
            dorth(3, 1) = dorth(1, 3)
            dorth(3, 2) = dorth(2, 3)
            dorth(3, 3) = e3*(un-nu12)/delta
            dorth(4, 4) = undemi*c1
!
!           calcul des valeurs propres de la matrice dorth :
            call dortvp(ndim, nomrc, dorth, 'DP')
!
!       traitement du cas 3d :
        else
            ndim = 3
            if (ien .eq. 0 .or. e3 .le. r8prem()) then
                goto 999
            end if
!
            c1 = e1/(un+nu12)
            delta = un-nu12-deux*nu13*nu13*e3/e1
!
            dorth(1, 1) = c1*(un-nu13*nu13*e3/e1)/delta
            dorth(1, 2) = c1*((un-nu13*nu13*e3/e1)/delta-un)
            dorth(1, 3) = e3*nu13/delta
            dorth(2, 1) = dorth(1, 2)
            dorth(2, 2) = dorth(1, 1)
            dorth(2, 3) = dorth(1, 3)
            dorth(3, 1) = dorth(1, 3)
            dorth(3, 2) = dorth(2, 3)
            dorth(3, 3) = e3*(un-nu12)/delta
            dorth(4, 4) = undemi*c1
            dorth(5, 5) = g13
            dorth(6, 6) = dorth(5, 5)
!
!           calcul des valeurs propres de la matrice dorth :
            call dortvp(ndim, nomrc, dorth, m2blan)
!
        end if
!
!   traitement du cas de l'orthotropie :
    else if (nomrc .eq. 'ELAS_ORTH') then
!
!       si g13 = 0 , on peut supposer que l'on est en 2d :
        if (igln .eq. 0) then
            ndim = 2
            if (iet .eq. 0 .or. e2 .le. r8prem() .or. e3 .le. r8prem()) then
                goto 999
            end if
!
            if (ien .eq. 0) then
                call utmess('A', 'MATERIAL2_10')
            else
!               traitement des cas des deformations planes et de l'axisymetrie :
                nu21 = e2*nu12/e1
                nu31 = e3*nu13/e1
                nu32 = e3*nu23/e2
                delta = un-nu23*nu32-nu31*nu13-nu21*nu12-deux*nu23*nu31*nu12
!
                dorth(1, 1) = (un-nu23*nu32)*e1/delta
                dorth(1, 2) = (nu21+nu31*nu23)*e1/delta
                dorth(1, 3) = (nu31+nu21*nu32)*e1/delta
                dorth(2, 2) = (un-nu31*nu13)*e2/delta
                dorth(2, 3) = (nu32+nu31*nu12)*e2/delta
                dorth(3, 3) = (un-nu21*nu12)*e3/delta
                dorth(2, 1) = dorth(1, 2)
                dorth(3, 1) = dorth(1, 3)
                dorth(3, 2) = dorth(2, 3)
                !
                dorth(4, 4) = g12
!
!               calcul des valeurs propres de la matrice dorth :
                call dortvp(ndim, nomrc, dorth, 'DP')
            end if
!
!           traitement du cas des contraintes planes :
!
            dorth(:, :) = zero
!
            nu21 = e2*nu12/e1
            delta = un-nu12*nu21
!
            dorth(1, 1) = e1/delta
            dorth(1, 2) = nu12*e2/delta
            dorth(2, 2) = e2/delta
            dorth(2, 1) = dorth(1, 2)
!
            dorth(4, 4) = g12
!
!           calcul des valeurs propres de la matrice dorth :
            call dortvp(ndim, nomrc, dorth, 'CP')
!
!       traitement du cas 3d :
        else
            ndim = 3
            if (iet .eq. 0 .or. e2 .le. r8prem() .or. e3 .le. r8prem()) then
                goto 999
            end if
            if (ien .eq. 0) then
                ndim = 2
!
                dorth(:, :) = zero
!
                nu21 = e2*nu12/e1
                delta = un-nu12*nu21
!
                dorth(1, 1) = e1/delta
                dorth(1, 2) = nu12*e2/delta
                dorth(2, 2) = e2/delta
                dorth(2, 1) = dorth(1, 2)
!
                dorth(4, 4) = g12
!
!               calcul des valeurs propres de la matrice dorth :
                call dortvp(ndim, nomrc, dorth, 'CP')
            end if
!
            nu21 = e2*nu12/e1
            nu31 = e3*nu13/e1
            nu32 = e3*nu23/e2
            delta = un-nu23*nu32-nu31*nu13-nu21*nu12-deux*nu23*nu31*nu12
!
            dorth(1, 1) = (un-nu23*nu32)*e1/delta
            dorth(1, 2) = (nu21+nu31*nu23)*e1/delta
            dorth(1, 3) = (nu31+nu21*nu32)*e1/delta
            dorth(2, 2) = (un-nu31*nu13)*e2/delta
            dorth(2, 3) = (nu32+nu31*nu12)*e2/delta
            dorth(3, 3) = (un-nu21*nu12)*e3/delta
            dorth(2, 1) = dorth(1, 2)
            dorth(3, 1) = dorth(1, 3)
            dorth(3, 2) = dorth(2, 3)
!
            dorth(4, 4) = g12
            dorth(5, 5) = g13
            dorth(6, 6) = g23
!
!           calcul des valeurs propres de la matrice dorth :
            call dortvp(ndim, nomrc, dorth, m2blan)
        end if
    else
        ASSERT(.false.)
    end if
999 continue
    call jedema
!
end subroutine
