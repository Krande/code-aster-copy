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
    subroutine xnmel(nnop, nfh, nfe, ddlc,&
                     ddlm, igeom, typmod, option, imate,&
                     compor, lgpg, carcri, jpintt, cnset,&
                     heavt, lonch, basloc, instam, instap, idepl, lsn,&
                     lst, sig, vi, matuu, ivectu,&
                     codret, jpmilt, nfiss, jheavn, jstno,&
                     l_line, l_nonlin, lMatr, lVect, lSigm)
        integer(kind=8) :: nfiss
        integer(kind=8) :: nnop
        integer(kind=8) :: nfh
        integer(kind=8) :: nfe
        integer(kind=8) :: ddlc
        integer(kind=8) :: ddlm
        integer(kind=8) :: igeom
        character(len=8) :: typmod(*)
        character(len=16) :: option
        integer(kind=8) :: imate
        character(len=16) :: compor(*)
        integer(kind=8) :: lgpg
        real(kind=8) :: carcri(*)
        integer(kind=8) :: jpintt
        integer(kind=8) :: cnset(128)
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: lonch(10)
        real(kind=8) :: basloc(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
        integer(kind=8) :: idepl
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop)
        real(kind=8) :: sig(*)
        real(kind=8) :: vi(*)
        real(kind=8) :: matuu(*)
        integer(kind=8) :: ivectu
        integer(kind=8) :: codret
        integer(kind=8) :: jpmilt
        integer(kind=8) :: jheavn
        integer(kind=8) :: jstno
        aster_logical, intent(in) :: l_line, l_nonlin, lMatr, lVect, lSigm
    end subroutine xnmel
end interface
