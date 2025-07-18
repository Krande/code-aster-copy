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
#include "asterf_types.h"
!
interface
    subroutine xasshm(ds_thm,&
                      nno, npg, npi, ipoids, ivf,&
                      idfde, igeom, geom, carcri, deplm,&
                      deplp, contm, contp, varim,&
                      varip, defgem, defgep, drds,&
                      drdsr, dsde, b, dfdi, dfdi2,&
                      r, sigbar, c, ck, cs,&
                      matuu, vectu, rinstm, rinstp, option,&
                      j_mater, mecani, press1, press2, tempe,&
                      dimdef, dimcon, dimuel, nbvari, nddls,&
                      nddlm, nmec, np1, ndim,&
                      compor, axi, modint, codret,&
                      nnop, nnops, nnopm, enrmec,&
                      dimenr, heavt, lonch, cnset, jpintt,&
                      jpmilt, jheavn, angmas,dimmat, enrhyd,&
                      nfiss, nfh, jfisno, work1, work2,&
                      lVect, lMatr, lVari, lSigm)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: dimmat
        integer(kind=8) :: dimenr
        integer(kind=8) :: nnops
        integer(kind=8) :: nnop
        integer(kind=8) :: ndim
        integer(kind=8) :: dimuel
        integer(kind=8) :: dimcon
        integer(kind=8) :: dimdef
        integer(kind=8) :: nno
        integer(kind=8) :: npg
        integer(kind=8) :: npi
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        integer(kind=8) :: igeom
        real(kind=8) :: geom(ndim, nnop)
        real(kind=8) :: carcri(*)
        real(kind=8) :: deplm(dimuel)
        real(kind=8) :: deplp(dimuel)
        real(kind=8) :: contm(*)
        real(kind=8) :: contp(*)
        real(kind=8) :: varim(*)
        real(kind=8) :: varip(*)
        real(kind=8) :: defgem(dimdef)
        real(kind=8) :: defgep(dimdef)
        real(kind=8) :: drds(dimenr, dimcon)
        real(kind=8) :: drdsr(dimenr, dimcon)
        real(kind=8) :: dsde(dimcon, dimenr)
        real(kind=8) :: b(dimenr, dimuel)
        real(kind=8) :: dfdi(nnop, ndim)
        real(kind=8) :: dfdi2(nnops, ndim)
        real(kind=8) :: r(dimenr)
        real(kind=8) :: sigbar(dimenr)
        real(kind=8) :: c(dimenr)
        real(kind=8) :: ck(dimenr)
        real(kind=8) :: cs(dimenr)
        real(kind=8) :: matuu(dimuel*dimuel)
        real(kind=8) :: vectu(dimuel)
        real(kind=8) :: rinstm
        real(kind=8) :: rinstp
        character(len=16) :: option
        integer(kind=8) :: j_mater
        integer(kind=8) :: mecani(5)
        integer(kind=8) :: press1(7)
        integer(kind=8) :: press2(7)
        integer(kind=8) :: tempe(5)
        integer(kind=8) :: nbvari
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: nmec
        integer(kind=8) :: np1
        character(len=16) :: compor(*)
        aster_logical :: axi
        character(len=3) :: modint
        integer(kind=8) :: codret
        integer(kind=8) :: nnopm
        integer(kind=8) :: enrmec(3)
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: lonch(10)
        integer(kind=8) :: cnset(*)
        integer(kind=8) :: jpintt
        integer(kind=8) :: jpmilt
        integer(kind=8) :: jheavn
        real(kind=8) :: angmas(3)
        integer(kind=8) :: enrhyd(3)
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
        integer(kind=8) :: jfisno
        real(kind=8) :: work1(dimcon, dimuel)
        real(kind=8) :: work2(dimenr, dimuel)
        aster_logical, intent(in) :: lVect, lMatr, lVari, lSigm
    end subroutine xasshm
end interface
