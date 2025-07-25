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

subroutine trigd(dg1, deb1, dg2, deb2, cumul, ino, nno)

    use calcul_module, only: ca_iachin_, ca_iachlo_, ca_ianueq_, ca_ilchlo_, &
                             ca_itypgd_, ca_lprno_, ca_ncmpmx_, ca_nec_

    implicit none

! person_in_charge: jacques.pellet at edf.fr

#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/exisdg.h"

    integer(kind=8) :: dg1(*), dg2(*), deb1, deb2, ino, nno
    aster_logical :: cumul
!-----------------------------------------------------------------------
!     entrees:
!       dg1    : dg de la grandeur de chin
!       deb1   : adresse du debut de la grandeur dans chin.vale
!       dg2    : dg de la grandeur de chloc
!       deb2   : adresse du debut de la grandeur dans chloc.vale
!       cumul  : / .f. : on affecte dans deb2 la grandeur lue
!                / .t. : on cumule dans deb2 la grandeur lue
!                cumul= .t. => on veut faire "moyenne" chno -> elem
!       ino    : numero du noeud associe a deb1
!                (ino n'est utilise que si cumul)
!       nno    : nombre de neouds de la maille
!                (nno n'est utilise que si cumul)
!
!     sorties:
!       recopie de chin.vale dans chloc pour les cmps voulus.
!       deb2   : est mis a jour pour le noeud suivant (si exchno).
!                et pour l'element suivant si excart ou exchno.
!-----------------------------------------------------------------------
    aster_logical :: change
    integer(kind=8) :: cmp, ind1, nec2, nsav, ksav, nmax
    parameter(nec2=120)
    parameter(nsav=5)
    parameter(nmax=13)
    integer(kind=8) :: ind2(nsav), necold(nsav)
    integer(kind=8) :: dg1old(nec2, nsav), dg2old(nec2, nsav), poscmp(nec2*nmax, nsav)
    integer(kind=8) :: ieq, i, k
    save dg1old, dg2old, poscmp, ind2, necold
    data necold/nsav*0/

!----------------------------------------------------------------------
!   1. on regarde si on ne trouve pas un poscmp(*,ksav) qui convient.
!        (calcul de ksav)
!   -----------------------------------------------------------------
    do k = 1, nsav
        if (ca_nec_ .ne. necold(k)) cycle
        change = .false.
        do i = 1, ca_nec_
            if (dg1(i) .ne. dg1old(i, k)) change = .true.
            if (dg2(i) .ne. dg2old(i, k)) change = .true.
        end do
        if (change) cycle
!       -- on a trouve un ksav convenable :
        ksav = k
        goto 80
    end do

!   -- on place le "assert" ici pour le faire moins souvent
    ASSERT(ca_nec_ .le. nec2)

!   2.1 on decale les tableaux vers le "fond" :
!         et on ajoute les nouvelles valeurs en ksav=1
!   ------------------------------------------------
    do k = nsav-1, 1, -1
        do i = 1, necold(k)
            dg1old(i, k+1) = dg1old(i, k)
            dg2old(i, k+1) = dg2old(i, k)
        end do
        do i = 1, ind2(k)
            poscmp(i, k+1) = poscmp(i, k)
        end do
        ind2(k+1) = ind2(k)
        necold(k+1) = necold(k)
    end do

    ksav = 1
    necold(ksav) = ca_nec_
    do i = 1, ca_nec_
        dg1old(i, ksav) = dg1(i)
        dg2old(i, ksav) = dg2(i)
    end do

!   2.2 remplissage de poscmp(ksav):
!   -------------------------------
    ind1 = 0
    ind2(ksav) = 0
    do cmp = 1, ca_ncmpmx_
        if (exisdg(dg1, cmp)) ind1 = ind1+1
        if (exisdg(dg2, cmp)) then
            ind2(ksav) = ind2(ksav)+1
            if (exisdg(dg1, cmp)) then
                poscmp(ind2(ksav), ksav) = ind1
            else
                poscmp(ind2(ksav), ksav) = 0
            end if
        end if
    end do

80  continue
!   3. recopie des valeurs dans le champ_local :
!   --------------------------------------------
    ASSERT(ind2(ksav) <= nec2*nmax)
    do cmp = 1, ind2(ksav)
        if (poscmp(cmp, ksav) .gt. 0) then
            ieq = deb1-1+poscmp(cmp, ksav)
            if (ca_lprno_ .eq. 1) ieq = zi(ca_ianueq_-1+ieq)

            if (.not. cumul) then
                if (ca_itypgd_ .eq. 1) then
                    zr(ca_iachlo_-1+deb2-1+cmp) = zr(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 2) then
                    zc(ca_iachlo_-1+deb2-1+cmp) = zc(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 3) then
                    zi(ca_iachlo_-1+deb2-1+cmp) = zi(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 4) then
                    zk8(ca_iachlo_-1+deb2-1+cmp) = zk8(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 5) then
                    zk16(ca_iachlo_-1+deb2-1+cmp) = zk16(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 6) then
                    zk24(ca_iachlo_-1+deb2-1+cmp) = zk24(ca_iachin_-1+ieq)
                else
                    ASSERT(.false.)
                end if

            else
                if (ca_itypgd_ .eq. 1) then
                    zr(ca_iachlo_-1+deb2-1+cmp) = zr(ca_iachlo_-1+deb2-1+cmp)+ &
                                                  zr(ca_iachin_-1+ieq)
                else if (ca_itypgd_ .eq. 2) then
                    zc(ca_iachlo_-1+deb2-1+cmp) = zc(ca_iachlo_-1+deb2-1+cmp)+ &
                                                  zc(ca_iachin_-1+ieq)
                else
                    ASSERT(.false.)
                end if
            end if

            if (cumul) then
                if (ino .eq. 1) then
                    zl(ca_ilchlo_-1+deb2-1+cmp) = .true.
                else
                end if
            else
                zl(ca_ilchlo_-1+deb2-1+cmp) = .true.
            end if

        else
            zl(ca_ilchlo_-1+deb2-1+cmp) = .false.
        end if
    end do

    if (cumul) then
        if (ino .eq. nno) deb2 = deb2+ind2(ksav)
    else
        deb2 = deb2+ind2(ksav)
    end if
end subroutine
