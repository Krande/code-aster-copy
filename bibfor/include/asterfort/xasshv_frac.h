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
interface
    subroutine xasshv_frac(ds_thm,&
                           nddls, nddlm, nnop, nnops,&
                           lact, elrefp, elrefc, elc, contac,&
                           dimuel, nface, npgf, nbspg, nptf,&
                           jcohes, jptint, igeom, jbasec,&
                           nlact, cface, rinstp,&
                           rinstm, carcri, fpg, ncompv, vect,&
                           compor, jmate, ndim, idepm, idepd, pla,&
                           algocr, rela, jheavn, ncompn, ifiss,&
                           nfiss, nfh, jheafa, ncomph, pos)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: nnop
        integer(kind=8) :: nnops
        integer(kind=8) :: lact(16)
        character(len=8) :: elrefp
        character(len=8) :: elrefc
        character(len=8) :: elc
        integer(kind=8) :: contac
        integer(kind=8) :: dimuel
        integer(kind=8) :: nface
        integer(kind=8) :: npgf
        integer(kind=8) :: nbspg
        integer(kind=8) :: nptf
        integer(kind=8) :: jcohes
        integer(kind=8) :: jptint
        integer(kind=8) :: igeom
        integer(kind=8) :: jbasec
        integer(kind=8) :: nlact(2)
        integer(kind=8) :: cface(30,6)
        real(kind=8) :: rinstp
        real(kind=8) :: rinstm
        real(kind=8) :: carcri(*)
        character(len=8) :: fpg
        integer(kind=8) :: ncompv
        real(kind=8) :: vect(560)
        character(len=16) :: compor(*)
        integer(kind=8) :: jmate
        integer(kind=8) :: ndim
        integer(kind=8) :: idepm
        integer(kind=8) :: idepd
        integer(kind=8) :: pla(27)
        integer(kind=8) :: algocr
        real(kind=8) :: rela
        integer(kind=8) :: jheavn
        integer(kind=8) :: ncompn
        integer(kind=8) :: ifiss
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
        integer(kind=8) :: jheafa
        integer(kind=8) :: ncomph
        integer(kind=8) :: pos(16)
    end subroutine xasshv_frac
end interface
