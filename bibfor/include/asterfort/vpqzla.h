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
          interface
            subroutine vpqzla(typeqz,qrn,iqrn,lqrn,qrar,qrai,qrba,qrvl, &
     &lvec,kqrn,lvalpr,nconv,omecor,ktyp,kqrnr,neqact,ilscal,irscal,    &
     &optiof,omemin,omemax,omeshi,ddlexc,dlagr,nfreq,lmasse,lraide,lamor,&
     &numedd,sigma,icscal,ivscal,iiscal,bwork,flage)
              character(len=16) :: typeqz
              integer(kind=8) :: qrn
              integer(kind=8) :: iqrn
              integer(kind=8) :: lqrn
              integer(kind=8) :: qrar
              integer(kind=8) :: qrai
              integer(kind=8) :: qrba
              integer(kind=8) :: qrvl
              integer(kind=8) :: lvec
              integer(kind=8) :: kqrn
              integer(kind=8) :: lvalpr
              integer(kind=8) :: nconv
              real(kind=8) :: omecor
              character(len=1) :: ktyp
              integer(kind=8) :: kqrnr
              integer(kind=8) :: neqact
              integer(kind=8) :: ilscal
              integer(kind=8) :: irscal
              character(len=16) :: optiof
              real(kind=8) :: omemin
              real(kind=8) :: omemax
              real(kind=8) :: omeshi
              integer(kind=8) :: ddlexc(*)
              integer(kind=8) :: dlagr(*)
              integer(kind=8) :: nfreq
              integer(kind=8) :: lmasse
              integer(kind=8) :: lraide
              integer(kind=8) :: lamor
              character(len=19) :: numedd
              complex(kind=8) :: sigma
              integer(kind=8) :: icscal
              integer(kind=8) :: ivscal
              integer(kind=8) :: iiscal
              logical(kind=4) :: bwork(*)
              aster_logical :: flage
            end subroutine vpqzla
          end interface
