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
     subroutine breselsqp(cequi, effmy, effmz, effn, ht, bw,&
                          enrobyi, enrobys, enrobzi, enrobzs,&
                          wmaxyi, wmaxys, wmaxzi, wmaxzs,&
                          ferrcomp, precs, ferrsyme, slsyme, uc, um,&
                          kt, eys, facier, fbeton, sigelsqp,&
                          phiyi, phiys, phizi, phizs,&
                          dnsyi, dnsys, dnszi, dnszs,&
                          sigmsyi, sigmsys, sigmcyi, sigmcys,&
                          sigmszi, sigmszs, sigmczi, sigmczs,&
                          alphay, alphaz, pivoty, pivotz, etaty, etatz,&
                          wfinyi, wfinys, wfinzi, wfinzs, kvarfy, kvarfz, ierr)
        real(kind=8) :: cequi
        real(kind=8) :: effmy
        real(kind=8) :: effmz
        real(kind=8) :: effn
        real(kind=8) :: ht
        real(kind=8) :: bw
        real(kind=8) :: enrobyi
        real(kind=8) :: enrobys
        real(kind=8) :: enrobzi
        real(kind=8) :: enrobzs
        real(kind=8) :: wmaxyi
        real(kind=8) :: wmaxys
        real(kind=8) :: wmaxzi
        real(kind=8) :: wmaxzs
        integer(kind=8) :: ferrcomp
        integer(kind=8) :: precs
        integer(kind=8) :: ferrsyme
        real(kind=8) :: slsyme
        integer(kind=8) :: uc
        integer(kind=8) :: um
        real(kind=8) :: kt
        real(kind=8) :: eys
        real(kind=8) :: facier
        real(kind=8) :: fbeton
        real(kind=8) :: sigelsqp
        real(kind=8) :: phiyi
        real(kind=8) :: phiys
        real(kind=8) :: phizi
        real(kind=8) :: phizs
        real(kind=8) :: dnsyi
        real(kind=8) :: dnsys
        real(kind=8) :: dnszi
        real(kind=8) :: dnszs
        real(kind=8) :: sigmsyi
        real(kind=8) :: sigmsys
        real(kind=8) :: sigmcyi
        real(kind=8) :: sigmcys
        real(kind=8) :: sigmszi
        real(kind=8) :: sigmszs
        real(kind=8) :: sigmczi
        real(kind=8) :: sigmczs
        real(kind=8) :: alphay
        real(kind=8) :: alphaz
        integer(kind=8) :: pivoty
        integer(kind=8) :: pivotz
        integer(kind=8) :: etaty
        integer(kind=8) :: etatz
        real(kind=8) :: wfinyi
        real(kind=8) :: wfinys
        real(kind=8) :: wfinzi
        real(kind=8) :: wfinzs
        real(kind=8) :: kvarfy
        real(kind=8) :: kvarfz
        integer(kind=8) :: ierr
    end subroutine breselsqp
end interface
