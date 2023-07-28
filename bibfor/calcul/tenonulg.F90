! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
!
subroutine tenonulg(nulg)
!
    use calcul_module, only: ca_ialiel_, ca_illiel_, ca_igr_, ca_iel_, &
                             ca_ilmaco_, ca_iamaco_, ca_ianulg_
!
    implicit none
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
!
    integer, intent(out) :: nulg(*)
!----------------------------------------------------------------------
!
! Sorties:
!     nulg : numéros globaux des noeuds (en std et hpc)
!            tableau de dimension au moins le nombre de noeuds de la maille
!
!   Remarques :
!     la maille ne doit pas être tardive
!----------------------------------------------------------------------
!
    integer :: ima, ino, nno, nunolc, nunogl
!----------------------------------------------------------------------
!
!
!   -- recuperation du numero de la maille et du nombre de noeuds :
!   ---------------------------------------------------------------
!
    ima = zi(ca_ialiel_-1+zi(ca_illiel_+ca_igr_-1)+ca_iel_-1)
    if (ima .gt. 0) then
        nno = zi(ca_ilmaco_-1+ima+1)-zi(ca_ilmaco_-1+ima)
    else
        ASSERT(ASTER_FALSE)
    end if
!
!
!   -- recuperation des numeros globaux des noeuds :
!   -------------------------------------------------
    do ino = 1, nno
        nunolc = zi(ca_iamaco_-1+zi(ca_ilmaco_+ima-1)+ino-1)
        if (ca_ianulg_ .ne. -1) then
            nunogl = zi(ca_ianulg_-1+nunolc)
        else
            nunogl = nunolc
        end if
        nulg(ino) = nunogl
    end do
!
!
end subroutine
