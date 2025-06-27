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
#include "asterf_types.h"
interface 
subroutine renfort3d(istep,nbrenf,numr,epstf33,vecr,&
                     epsr0,epsrf,eprp0,eprpf,your,&
                     syr,sigr0,sigrf,hplr,tomr,&
                     ekr,skr,ATRref,khir,gamr,&
                     sprec,teta1,teta2,dt,ppas,&
                     theta,eprm0,eprmf,ttaref,rhor,&
                     mu00,fl3d,errr,xnr,xmuthr,&
                     eprk0,eprkf,tokr,kr,plast_seule,&
                     ann,xn,bn,ngf,ipzero,&
                     epspmf33,epspmf,eps_nl,spre0,spref)

        integer(kind=8) :: istep
        integer(kind=8) :: nbrenf
        integer(kind=8) :: numr
        real(kind=8) :: epstf33(3,3)
        real(kind=8) :: vecr(nbrenf,3)
        real(kind=8) :: epsr0
        real(kind=8) :: epsrf
        real(kind=8) :: eprp0
        real(kind=8) :: eprpf
        real(kind=8) :: your
        real(kind=8) :: syr
        real(kind=8) :: sigr0
        real(kind=8) :: sigrf
        real(kind=8) :: hplr
        real(kind=8) :: tomr
        real(kind=8) :: ekr
        real(kind=8) :: skr
        real(kind=8) :: ATRref
        real(kind=8) :: khir
        real(kind=8) :: gamr
        real(kind=8) :: sprec
        real(kind=8) :: teta1
        real(kind=8) :: teta2
        real(kind=8) :: dt
        aster_logical :: ppas
        real(kind=8) :: theta
        real(kind=8) :: eprm0
        real(kind=8) :: eprmf
        real(kind=8) :: ttaref
        real(kind=8) :: rhor
        real(kind=8) :: mu00
        aster_logical :: fl3d
        aster_logical :: errr
        real(kind=8) :: xnr
        real(kind=8) :: xmuthr
        real(kind=8) :: eprk0
        real(kind=8) :: eprkf
        real(kind=8) :: tokr
        real(kind=8) :: kr
        aster_logical :: plast_seule
        real(kind=8) :: ann(ngf,ngf+1)
        real(kind=8) :: xn(ngf)
        real(kind=8) :: bn(ngf)
        integer(kind=8) :: ngf
        integer(kind=8) :: ipzero(ngf)
        real(kind=8) :: epspmf33(3,3) 
        real(kind=8) :: epspmf
        real(kind=8) :: eps_nl
        real(kind=8) :: spre0
        real(kind=8) :: spref

    end subroutine renfort3d
end interface
