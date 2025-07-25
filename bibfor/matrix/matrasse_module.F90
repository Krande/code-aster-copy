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

!
! Sparse Real Matrix stored in Aster format
! This module provides some tools for aster "matr_asse"
! matrices
!
module matrasse_module
!
!
! because the pointer is a result
! aslint: disable=C1310
!
    implicit none
    private
#include "asterc/asmpi_comm.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    !
    integer(kind=8), private :: ierr
    integer(kind=8), parameter, public :: lagrange1_dof = 0, physical_dof = 1
    integer(kind=8), parameter, public :: lagrange2_dof = 2
    !
    public :: get_indices_of_dofs, get_num_of_dofs
    !
contains
    !
    ! This function returns the number of degrees of freedom
    ! with a selected type (physical_dof or lagrange1_dof or lagrange2_dof)
    ! in a matrasse
    function get_num_of_dofs(type_dof, matass) result(ndof)
        ! Dummy arguments
        integer(kind=8), intent(in)                     :: type_dof
        character(len=19), intent(in)           :: matass
        integer(kind=8)                                 :: ndof
        ! Local variables
        character(len=14) :: nonu
        integer(kind=8), dimension(:), pointer :: delg => null()
        integer(kind=8) :: nbeq, iret
        aster_logical :: test_dof
        !
        call jemarq()
        !
        test_dof = (type_dof == physical_dof) .or. (type_dof == lagrange1_dof) &
                   .or. (type_dof == lagrange2_dof)
        ASSERT(test_dof)
        !
        ! La matrice existe-t-elle ?
        call jeexin(matass//'.REFA', iret)
        ASSERT(iret > 0)
        ! Quelle est sa taille ?
        call dismoi('NB_EQUA', matass, 'MATR_ASSE', repi=nbeq)
        ! Le tableau delg permet de distinguer ddls physiques/lagrange
        call dismoi('NOM_NUME_DDL', matass, 'MATR_ASSE', repk=nonu)
        call jeveuo(nonu//'.NUME.DELG', 'L', vi=delg)
        !
        ! Nombre de ddls sélectionnés
        select case (type_dof)
        case (physical_dof)
            ndof = count(delg(1:nbeq) == 0)
        case (lagrange1_dof)
            ndof = count(delg(1:nbeq) == -1)
        case (lagrange2_dof)
            ndof = count(delg(1:nbeq) == -2)
        end select
        !
        call jedema()
        !
    end function get_num_of_dofs
    !
    ! This function allocates, fills and returns an integer array with
    ! the (Fortran) indices of some selected dofs : one may select
    ! physical_dof or lagrange1_dof or lagrange2_dof
    !
    function get_indices_of_dofs(type_dof, matass) result(idof)
        !
        ! Dummy arguments
        integer(kind=8), intent(in)                     :: type_dof
        character(len=19), intent(in)           :: matass
        integer(kind=8), dimension(:), pointer           :: idof
        ! Local variables
        character(len=14) :: nonu
        integer(kind=8), dimension(:), pointer :: delg => null()
        integer(kind=8) :: nbeq, nlag1, nlag2, nphys, ndof, iret
        integer(kind=8) :: i
        aster_logical :: test_dof
        mpi_int :: mpicomm
        !
        call jemarq()
        !
        test_dof = (type_dof == physical_dof) .or. (type_dof == lagrange1_dof) &
                   .or. (type_dof == lagrange2_dof)
        ASSERT(test_dof)
        !
        idof => null()
        ! Communicateur MPI
        call asmpi_comm('GET', mpicomm)
        !
        ! La matrice existe-t-elle ?
        call jeexin(matass//'.REFA', iret)
        ASSERT(iret > 0)
        ! Quelle est sa taille ?
        call dismoi('NB_EQUA', matass, 'MATR_ASSE', repi=nbeq)
        ! Le tableau delg permet de distinguer ddls physiques/lagrange
        call dismoi('NOM_NUME_DDL', matass, 'MATR_ASSE', repk=nonu)
        call jeveuo(nonu//'.NUME.DELG', 'L', vi=delg)
        !
        ! Construction du vecteur d'indices
        !
        ! Nombre de ddls physiques
        nphys = count(delg(1:nbeq) == 0)
        ! Nombre de Lagrange 1
        nlag1 = count(delg(1:nbeq) == -1)
        ! Nombre de Lagrange 2
        nlag2 = count(delg(1:nbeq) == -2)
        ! Vérification de la cohérence des dimensions trouvées
        ASSERT(nbeq .eq. nphys+nlag1+nlag2)
        ! Nombre de ddls sélectionnés (soit lagrange1, soit physiques)
        select case (type_dof)
        case (physical_dof)
            ndof = nphys
        case (lagrange1_dof)
            ndof = nlag1
        case (lagrange2_dof)
            ndof = nlag2
        end select
        !
        ! Allocation du vecteur d'indices (Fortran)
        allocate (idof(ndof), stat=ierr)
        ASSERT(ierr == 0)
        idof(:) = 0
        !
        ndof = 0
        do i = 1, nbeq
            if (delg(i) .eq. 0) then
                if (type_dof == physical_dof) then
                    ndof = ndof+1
!          idof(ndof)= int(i, 4)
                    idof(ndof) = i
                end if
            end if
            if (delg(i) .eq. -1) then
                nlag1 = nlag1+1
                if (type_dof == lagrange1_dof) then
                    ndof = ndof+1
                    !       idof(ndof) = int(i, 4)
                    idof(ndof) = i
                end if
            end if
            if (delg(i) .eq. -2) then
                nlag2 = nlag2+1
                if (type_dof == lagrange2_dof) then
                    ndof = ndof+1
                    idof(ndof) = i
                end if
            end if
        end do
        !
        call jedema()
    end function get_indices_of_dofs
    !
end module matrasse_module
