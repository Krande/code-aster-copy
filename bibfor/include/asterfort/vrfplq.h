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
   subroutine vrfplq(typcmb, typco, nb, cequi, enrobi, enrobs, sigs, sigci, &
                     sigcs, alphacc, gammas, gammac, facier, eys, typdiag, fbeton, clacier, &
                     uc, ht, effrts, effref, dnsxi, dnsxs, dnsyi, dnsys, marge, theta_crit, &
                     dnsinf_crit, dnssup_crit, effn0_crit, effm0_crit, effn_crit, effm_crit, &
                     myNrd_crit, myMrd_crit, bw, tau_crit, c0c_crit, c0crd_crit)

      integer(kind=8) ::typcmb
      integer(kind=8) :: typco
      integer(kind=8) :: nb
      real(kind=8) :: cequi
      real(kind=8) :: alphacc
      real(kind=8) :: ht
      real(kind=8) :: enrobi
      real(kind=8) :: enrobs
      real(kind=8) :: sigs
      real(kind=8) :: sigci
      real(kind=8) :: sigcs
      real(kind=8) :: facier
      real(kind=8) :: fbeton
      real(kind=8) :: gammas
      real(kind=8) :: gammac
      integer(kind=8) :: clacier
      real(kind=8) :: eys
      integer(kind=8) :: typdiag
      integer(kind=8) :: uc
      real(kind=8) :: effrts(8)
      real(kind=8) :: effref(8)
      real(kind=8) :: dnsxi
      real(kind=8) :: dnsxs
      real(kind=8) :: dnsyi
      real(kind=8) :: dnsys
      real(kind=8) :: marge
      real(kind=8) :: tau_crit
      real(kind=8) :: c0c_crit
      real(kind=8) :: c0crd_crit
      real(kind=8) :: dnsinf_crit, dnssup_crit, myNrd_crit, myMrd_crit, theta
      real(kind=8) :: theta_crit, effn_crit, effm_crit, effn0_crit, effm0_crit
   end subroutine vrfplq
end interface
