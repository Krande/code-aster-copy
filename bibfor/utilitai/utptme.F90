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

subroutine utptme(nomarg, valarg, iret)
    implicit none
    character(len=8), intent(in) :: nomarg
    real(kind=8), intent(in) :: valarg
    integer(kind=8), intent(out) :: iret
! person_in_charge: mathieu.courtois at edf.fr
! ----------------------------------------------------------------------
!     Affecte la valeur associée au nom nomarg du paramètre mémoire en Mo
! in  nomarg  : nom du paramètre
! in  valarg  : valeur du paramètre
! out iret    : code retour
!                =0 la valeur a été affectée
!               !=0 la valeur est invalide
!
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2)
    real(kind=8) :: vmumps, vpetsc, rlqmem, vminit
    common/msolve/vmumps, vpetsc, rlqmem, vminit
! ----------------------------------------------------------------------
    iret = 0
    if (nomarg .eq. 'LIMIT_JV') then
! ----- Limite memoire jeveux (modifiee par jermxd)
        vmxdyn = valarg*1024*1024/lois

    else if (nomarg .eq. 'MEM_TOTA') then
! --------- Limite memoire allouee lors de l'execution
        vmet = valarg*1024*1024/lois

    else if (nomarg .eq. 'RLQ_MEM') then
! -------- Reliquat memoire (consommation hors jeveux et solveur)
        rlqmem = valarg*(1024*1024)

    else if (nomarg .eq. 'MEM_MUMP') then
! --------- Consommation memoire du solveur mumps
        vmumps = valarg*(1024*1024)

    else if (nomarg .eq. 'MEM_PETS') then
! --------- Consommation memoire du solveur petsc
        vpetsc = valarg*(1024*1024)

    else if (nomarg .eq. 'MEM_INIT') then
! --------- Consommation memoire du jdc
        vminit = valarg*(1024*1024)

    else
        iret = 1
    end if
!
end subroutine
