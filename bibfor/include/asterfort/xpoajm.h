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
    subroutine xpoajm(maxfem, jtypm2, itypse, jcnse, im,&
                      n, nnose, prefno, jdirno, nnm,&
                      inm, inmtot, nbmac, he, jnivgr,&
                      iagma, ngrm, jdirgr, opmail, nfiss,&
                      ndim, ndime, jconx1, jconx2, jconq1,&
                      jconq2, ima, iad1, nnn, inn,&
                      inntot, nbnoc, nbnofi, inofi, iacoo1,&
                      iacoo2, iad9, ninter, iainc, ncompa, elrefp,&
                      jlsn, jlst, typma, igeom, jheavn, ncompn,&
                      contac, cmp, nbcmp, nfh, nfe,&
                      ddlc, jcnsv1, jcnsv2, jcnsl2, lmeca,&
                      pre1, heavno, fisco,&
                      nlachm, lacthm, jbaslo, jstno, ka, mu)
        integer(kind=8) :: nfiss
        character(len=8) :: maxfem
        integer(kind=8) :: jtypm2
        integer(kind=8) :: itypse
        integer(kind=8) :: jcnse
        integer(kind=8) :: im
        integer(kind=8) :: n
        integer(kind=8) :: nnose
        character(len=2) :: prefno(4)
        integer(kind=8) :: jdirno
        integer(kind=8) :: nnm
        integer(kind=8) :: inm
        integer(kind=8) :: inmtot
        integer(kind=8) :: nbmac
        integer(kind=8) :: he(nfiss)
        integer(kind=8) :: jnivgr
        integer(kind=8) :: iagma
        integer(kind=8) :: ngrm
        integer(kind=8) :: jdirgr
        aster_logical :: opmail
        integer(kind=8) :: ndim
        integer(kind=8) :: ndime
        integer(kind=8) :: jconx1
        integer(kind=8) :: jconx2
        integer(kind=8) :: jconq1
        integer(kind=8) :: jconq2
        integer(kind=8) :: ima
        integer(kind=8) :: iad1
        integer(kind=8) :: nnn
        integer(kind=8) :: inn
        integer(kind=8) :: inntot
        integer(kind=8) :: nbnoc
        integer(kind=8) :: nbnofi
        integer(kind=8) :: inofi
        integer(kind=8) :: iacoo1
        integer(kind=8) :: iacoo2
        integer(kind=8) :: iad9
        integer(kind=8) :: ninter(4)
        integer(kind=8) :: iainc
        character(len=8) :: elrefp
        integer(kind=8) :: jlsn
        integer(kind=8) :: jlst
        character(len=8) :: typma
        integer(kind=8) :: igeom
        integer(kind=8) :: jheavn
        integer(kind=8) :: ncompn
        integer(kind=8) :: contac
        integer(kind=8) :: cmp(*)
        integer(kind=8) :: nbcmp
        integer(kind=8) :: nfh
        integer(kind=8) :: nfe
        integer(kind=8) :: ddlc
        integer(kind=8) :: jcnsv1
        integer(kind=8) :: jcnsv2
        integer(kind=8) :: jcnsl2
        integer(kind=8) :: jstno
        aster_logical :: lmeca
        aster_logical :: pre1
        integer(kind=8) :: ncompa
        integer(kind=8) :: heavno(20,3)
        integer(kind=8) :: fisco(*)
        integer(kind=8) :: nlachm(2)
        integer(kind=8) :: lacthm(16)
        integer(kind=8) :: jbaslo
        real(kind=8) :: ka
        real(kind=8) :: mu
    end subroutine xpoajm
end interface 
