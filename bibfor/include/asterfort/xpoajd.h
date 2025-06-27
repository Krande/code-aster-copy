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
!
#include "asterf_types.h"
!
interface 
    subroutine xpoajd(elrefp, ino, nnop, lsn, lst,&
                      ninter, iainc, ncompa, typma, co, igeom,&
                      jdirno, nfiss, jheavn, ncompn, he, ndime,&
                      ndim, cmp, nbcmp, nfh, nfe,&
                      ddlc, ima, jconx1, jconx2, jcnsv1,&
                      jcnsv2, jcnsl2, nbnoc, inntot, inn,&
                      nnn, contac, lmeca, pre1, heavno,&
                      nlachm, lacthm, jbaslo, &
                      jlsn, jlst, jstno, ka, mu)
        integer(kind=8) :: nbcmp
        integer(kind=8) :: nfiss
        integer(kind=8) :: nnop
        character(len=8) :: elrefp
        integer(kind=8) :: ino
        real(kind=8) :: lsn(nfiss)
        real(kind=8) :: lst(nfiss)
        integer(kind=8) :: ninter(4)
        integer(kind=8) :: iainc
        character(len=8) :: typma
        real(kind=8) :: co(3)
        integer(kind=8) :: igeom
        integer(kind=8) :: jdirno
        integer(kind=8) :: jfisno
        integer(kind=8) :: jheavn
        integer(kind=8) :: ncompn
        integer(kind=8) :: he(nfiss)
        integer(kind=8) :: ndime
        integer(kind=8) :: ndim
        integer(kind=8) :: cmp(*)
        integer(kind=8) :: nfh
        integer(kind=8) :: nfe
        integer(kind=8) :: ddlc
        integer(kind=8) :: ima
        integer(kind=8) :: jconx1
        integer(kind=8) :: jconx2
        integer(kind=8) :: jcnsv1
        integer(kind=8) :: jcnsv2
        integer(kind=8) :: jcnsl2
        integer(kind=8) :: nbnoc
        integer(kind=8) :: inntot
        integer(kind=8) :: inn
        integer(kind=8) :: nnn
        integer(kind=8) :: contac
        aster_logical :: lmeca
        aster_logical :: pre1
        integer(kind=8) :: ncompa
        integer(kind=8) :: heavno(20, 3)
        integer(kind=8) :: nlachm(2)
        integer(kind=8) :: lacthm(16)
        integer(kind=8) :: jbaslo
        integer(kind=8) :: jlsn
        integer(kind=8) :: jlst
        integer(kind=8) :: jstno
        real(kind=8) :: ka
        real(kind=8) :: mu
    end subroutine xpoajd
end interface 
