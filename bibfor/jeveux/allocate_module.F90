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

module allocate_module
    use iso_c_binding, only: c_ptr, c_loc, C_NULL_PTR, c_associated
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
!---------------------------------------------------------------------------
! Ce module contient une variable globale  (slvec)
! permettant de garder la trace de tous les objets alloues
! par le mecanisme AS_ALLOCATE / AS_DEALLOCATE
!
! Cette variable permet de controler les fuites memoire
! et de desallouer les objets non desalloues du fait d'une interruption
! brusque (par exemple par un "try / except" dans un fichier de commande)
!---------------------------------------------------------------------------
! ne pas rendre publiques les fonctions externes
! (sinon il faudrait les définir dans un module)
    private :: utmess, assert
#include "asterf_types.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
!
!------------------------------------------------------------------------
!   -- commons jeveux :
!   --------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
!------------------------------------------------------------------------
!
    type array_1
        aster_logical, allocatable :: vl(:)
        integer(kind=8), allocatable :: vi(:)
        integer(kind=4), allocatable :: vi4(:)
        real(kind=8), allocatable :: vr(:)
        complex(kind=8), allocatable :: vc(:)
        character(len=8), allocatable :: vk8(:)
        character(len=16), allocatable :: vk16(:)
        character(len=24), allocatable :: vk24(:)
        character(len=32), allocatable :: vk32(:)
        character(len=80), allocatable :: vk80(:)
!
        aster_logical :: present
        type(c_ptr) :: ptr_ident
        character(len=3) :: tsca
    end type array_1
!
    type save_lvec
        type(array_1), allocatable :: lvec(:)
        integer(kind=8) :: nmax = 0, kfree = 0
    end type save_lvec
!
    type(save_lvec), target :: slvec
!
! -- pour mesurer le taux de remplissage de slvec
    integer(kind=8), save :: kmax = 0
!
!
contains
!
    subroutine init_slvec(slvec, nmax)
        type(save_lvec) :: slvec
        integer(kind=8) :: nmax, k
!
        allocate (slvec%lvec(nmax))
        slvec%nmax = nmax
        slvec%kfree = 1
        do k = 1, nmax
            slvec%lvec(k)%present = .false.
            slvec%lvec(k)%ptr_ident = C_NULL_PTR
            slvec%lvec(k)%tsca = ' '
        end do
    end subroutine
!
    subroutine free_slvec(slvec)
        type(save_lvec) :: slvec
!
        if (slvec%nmax > 0) then
            call deallocate_all_slvec()
            deallocate (slvec%lvec)
            slvec%nmax = 0
            slvec%kfree = 0
        end if
    end subroutine
!
!---------------------------------------------------------------------
    subroutine allocate_slvec(lon1, vl, vi, vi4, vr, &
                              vc, vk8, vk16, vk24, vk32, &
                              vk80)
        integer(kind=8) :: lon1
        aster_logical, pointer, optional :: vl(:)
        integer(kind=8), pointer, optional :: vi(:)
        integer(kind=4), pointer, optional :: vi4(:)
        real(kind=8), pointer, optional :: vr(:)
        complex(kind=8), pointer, optional :: vc(:)
        character(len=8), pointer, optional :: vk8(:)
        character(len=16), pointer, optional :: vk16(:)
        character(len=24), optional, pointer :: vk24(:)
        character(len=32), optional, pointer :: vk32(:)
        character(len=80), optional, pointer :: vk80(:)
!
        integer(kind=8) :: k, ktrou
        type(c_ptr) :: pteur_c
!
        ktrou = 0
        if (slvec%kfree < 1) then
            do k = 1, slvec%nmax
                if (.not. slvec%lvec(k)%present) then
                    ktrou = k
                    goto 1
                end if
            end do
            call utmess('F', 'DVP_7')
        else
            ktrou = slvec%kfree
            ASSERT(.not. slvec%lvec(ktrou)%present)
        end if
1       continue
        if (kmax .lt. ktrou) then
            kmax = ktrou
        end if
        slvec%kfree = -1
        if (ktrou+1 <= slvec%nmax .and. .not. slvec%lvec(ktrou+1)%present) then
            slvec%kfree = ktrou+1
        end if
!
        if (present(vl)) then
            allocate (slvec%lvec(ktrou)%vl(lon1))
            slvec%lvec(ktrou)%vl = .false.
            vl => slvec%lvec(ktrou)%vl
            pteur_c = c_loc(vl(1))
            slvec%lvec(ktrou)%tsca = 'L'
        else if (present(vi)) then
            allocate (slvec%lvec(ktrou)%vi(lon1))
            slvec%lvec(ktrou)%vi = 0
            vi => slvec%lvec(ktrou)%vi
            pteur_c = c_loc(vi(1))
            slvec%lvec(ktrou)%tsca = 'I'
        else if (present(vi4)) then
            allocate (slvec%lvec(ktrou)%vi4(lon1))
            slvec%lvec(ktrou)%vi4 = 0
            vi4 => slvec%lvec(ktrou)%vi4
            pteur_c = c_loc(vi4(1))
            slvec%lvec(ktrou)%tsca = 'S'
        else if (present(vr)) then
            allocate (slvec%lvec(ktrou)%vr(lon1))
            slvec%lvec(ktrou)%vr = 0.d0
            vr => slvec%lvec(ktrou)%vr
            pteur_c = c_loc(vr(1))
            slvec%lvec(ktrou)%tsca = 'R'
        else if (present(vc)) then
            allocate (slvec%lvec(ktrou)%vc(lon1))
            slvec%lvec(ktrou)%vc = cmplx(0.d0, 0.d0)
            vc => slvec%lvec(ktrou)%vc
            pteur_c = c_loc(vc(1))
            slvec%lvec(ktrou)%tsca = 'C'
        else if (present(vk8)) then
            allocate (slvec%lvec(ktrou)%vk8(lon1))
            slvec%lvec(ktrou)%vk8 = ' '
            vk8 => slvec%lvec(ktrou)%vk8
            pteur_c = c_loc(vk8(1))
            slvec%lvec(ktrou)%tsca = 'K8'
        else if (present(vk16)) then
            allocate (slvec%lvec(ktrou)%vk16(lon1))
            slvec%lvec(ktrou)%vk16 = ' '
            vk16 => slvec%lvec(ktrou)%vk16
            pteur_c = c_loc(vk16(1))
            slvec%lvec(ktrou)%tsca = 'K16'
        else if (present(vk24)) then
            allocate (slvec%lvec(ktrou)%vk24(lon1))
            slvec%lvec(ktrou)%vk24 = ' '
            vk24 => slvec%lvec(ktrou)%vk24
            pteur_c = c_loc(vk24(1))
            slvec%lvec(ktrou)%tsca = 'K24'
        else if (present(vk32)) then
            allocate (slvec%lvec(ktrou)%vk32(lon1))
            slvec%lvec(ktrou)%vk32 = ' '
            vk32 => slvec%lvec(ktrou)%vk32
            pteur_c = c_loc(vk32(1))
            slvec%lvec(ktrou)%tsca = 'K32'
        else if (present(vk80)) then
            allocate (slvec%lvec(ktrou)%vk80(lon1))
            slvec%lvec(ktrou)%vk80 = ' '
            vk80 => slvec%lvec(ktrou)%vk80
            pteur_c = c_loc(vk80(1))
            slvec%lvec(ktrou)%tsca = 'K80'
        else
            ASSERT(.false.)
        end if
!
        slvec%lvec(ktrou)%present = .true.
        slvec%lvec(ktrou)%ptr_ident = pteur_c
    end subroutine allocate_slvec
!
!
!---------------------------------------------------------------------
    subroutine deallocate_slvec(ierr, vl, vi, vi4, vr, &
                                vc, vk8, vk16, vk24, vk32, &
                                vk80)
        integer(kind=8), intent(out) :: ierr
        aster_logical, pointer, optional :: vl(:)
        integer(kind=8), pointer, optional :: vi(:)
        integer(kind=4), pointer, optional :: vi4(:)
        real(kind=8), pointer, optional :: vr(:)
        complex(kind=8), pointer, optional :: vc(:)
        character(len=8), pointer, optional :: vk8(:)
        character(len=16), pointer, optional :: vk16(:)
        character(len=24), optional, pointer :: vk24(:)
        character(len=32), optional, pointer :: vk32(:)
        character(len=80), optional, pointer :: vk80(:)
!
        integer(kind=8) :: k, ktrou
        type(c_ptr) :: pteur_c
!
        ierr = 1
!
        if (present(vl)) then
            if (.not. associated(vl)) goto 2
            pteur_c = c_loc(vl(1))
        else if (present(vi)) then
            if (.not. associated(vi)) goto 2
            pteur_c = c_loc(vi(1))
        else if (present(vi4)) then
            if (.not. associated(vi4)) goto 2
            pteur_c = c_loc(vi4(1))
        else if (present(vr)) then
            if (.not. associated(vr)) goto 2
            pteur_c = c_loc(vr(1))
        else if (present(vc)) then
            if (.not. associated(vc)) goto 2
            pteur_c = c_loc(vc(1))
        else if (present(vk8)) then
            if (.not. associated(vk8)) goto 2
            pteur_c = c_loc(vk8(1))
        else if (present(vk16)) then
            if (.not. associated(vk16)) goto 2
            pteur_c = c_loc(vk16(1))
        else if (present(vk24)) then
            if (.not. associated(vk24)) goto 2
            pteur_c = c_loc(vk24(1))
        else if (present(vk32)) then
            if (.not. associated(vk32)) goto 2
            pteur_c = c_loc(vk32(1))
        else if (present(vk80)) then
            if (.not. associated(vk80)) goto 2
            pteur_c = c_loc(vk80(1))
        else
            ASSERT(.false.)
        end if
!
        ktrou = 0
        if (slvec%kfree > 1) then
            if (slvec%lvec(slvec%kfree-1)%present) then
                if (c_associated(slvec%lvec(slvec%kfree-1)%ptr_ident, pteur_c)) then
                    ktrou = slvec%kfree-1
                    goto 1
                end if
            end if
        end if
        do k = 1, slvec%nmax
            if (slvec%lvec(k)%present) then
                if (c_associated(slvec%lvec(k)%ptr_ident, pteur_c)) then
                    ktrou = k
                    goto 1
                end if
            end if
        end do
        ASSERT(.false.)
1       continue
!
        if (present(vl)) then
            deallocate (slvec%lvec(ktrou)%vl, stat=ierr)
            nullify (vl)
        else if (present(vi)) then
            deallocate (slvec%lvec(ktrou)%vi, stat=ierr)
            nullify (vi)
        else if (present(vi4)) then
            deallocate (slvec%lvec(ktrou)%vi4, stat=ierr)
            nullify (vi4)
        else if (present(vr)) then
            deallocate (slvec%lvec(ktrou)%vr, stat=ierr)
            nullify (vr)
        else if (present(vc)) then
            deallocate (slvec%lvec(ktrou)%vc, stat=ierr)
            nullify (vc)
        else if (present(vk8)) then
            deallocate (slvec%lvec(ktrou)%vk8, stat=ierr)
            nullify (vk8)
        else if (present(vk16)) then
            deallocate (slvec%lvec(ktrou)%vk16, stat=ierr)
            nullify (vk16)
        else if (present(vk24)) then
            deallocate (slvec%lvec(ktrou)%vk24, stat=ierr)
            nullify (vk24)
        else if (present(vk32)) then
            deallocate (slvec%lvec(ktrou)%vk32, stat=ierr)
            nullify (vk32)
        else if (present(vk80)) then
            deallocate (slvec%lvec(ktrou)%vk80, stat=ierr)
            nullify (vk80)
        else
            ASSERT(.false.)
        end if
!
        slvec%lvec(ktrou)%present = .false.
        slvec%lvec(ktrou)%ptr_ident = C_NULL_PTR
        slvec%lvec(ktrou)%tsca = ' '
        slvec%kfree = ktrou
2       continue
    end subroutine deallocate_slvec
!
!---------------------------------------------------------------------
    subroutine deallocate_all_slvec()
! but : desallouer tous les objets de slvec qui ne n'ont pas ete
!
        integer(kind=8) :: k, n1, n2, lonty, lsic
        character(len=3) :: tsca
!
        n2 = 0
        do k = 1, slvec%nmax
            if (slvec%lvec(k)%present) then
                tsca = slvec%lvec(k)%tsca
                if (tsca .eq. 'L') then
                    n1 = size(slvec%lvec(k)%vl)
                    deallocate (slvec%lvec(k)%vl)
                    lonty = lois
                else if (tsca .eq. 'I') then
                    n1 = size(slvec%lvec(k)%vi)
                    deallocate (slvec%lvec(k)%vi)
                    lonty = lois
                else if (tsca .eq. 'S') then
                    n1 = size(slvec%lvec(k)%vi4)
                    deallocate (slvec%lvec(k)%vi4)
                    lonty = 4
                else if (tsca .eq. 'R') then
                    n1 = size(slvec%lvec(k)%vr)
                    deallocate (slvec%lvec(k)%vr)
                    lonty = 8
                else if (tsca .eq. 'C') then
                    n1 = size(slvec%lvec(k)%vc)
                    deallocate (slvec%lvec(k)%vc)
                    lonty = 16
                else if (tsca .eq. 'K8') then
                    n1 = size(slvec%lvec(k)%vk8)
                    deallocate (slvec%lvec(k)%vk8)
                    lonty = 8
                else if (tsca .eq. 'K16') then
                    n1 = size(slvec%lvec(k)%vk16)
                    deallocate (slvec%lvec(k)%vk16)
                    lonty = 16
                else if (tsca .eq. 'K24') then
                    n1 = size(slvec%lvec(k)%vk24)
                    deallocate (slvec%lvec(k)%vk24)
                    lonty = 24
                else if (tsca .eq. 'K32') then
                    n1 = size(slvec%lvec(k)%vk32)
                    deallocate (slvec%lvec(k)%vk32)
                    lonty = 32
                else if (tsca .eq. 'K80') then
                    n1 = size(slvec%lvec(k)%vk80)
                    deallocate (slvec%lvec(k)%vk80)
                    lonty = 80
                else
                    ASSERT(.false.)
                end if
                n2 = n2+n1*lonty
                slvec%lvec(k)%present = .false.
                slvec%lvec(k)%ptr_ident = C_NULL_PTR
                slvec%lvec(k)%tsca = ' '
            end if
        end do
        slvec%kfree = 1
        kmax = 0
!
!   -- il faut "rendre" la memoire a JEVEUX
!   -------------------------------------------------
        lsic = n2/lois
        mcdyn = mcdyn-lsic
        cuvtrav = cuvtrav-lsic
    end subroutine deallocate_all_slvec
!
end module
