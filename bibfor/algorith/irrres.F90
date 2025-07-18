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

subroutine irrres(fami, kpg, ksp, mod, nmat, &
                  materd, materf, yd, yf, deps, &
                  dy, r)
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcnrts.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
    character(len=8) :: mod
    character(len=*) :: fami
    integer(kind=8) :: nmat, kpg, ksp
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2)
    real(kind=8) :: yd(*), yf(*), deps(6), dy(*), r(*)
!
! person_in_charge: jean-luc.flejou at edf.fr
!
    real(kind=8) :: dfds(6), id3d(6), sf
    real(kind=8) :: irrad, irraf, dphi, sigd(6), sigf(6), dkooh(6, 6)
    real(kind=8) :: fkooh(6, 6), hookf(6, 6), k, n, p0, ai0, etais, p, zetaff
    real(kind=8) :: zetag
    real(kind=8) :: dp, detai, dpi, dg, etaif, seqf, epsef(6), pe, kappa, r02
    real(kind=8) :: depsa(6), depsg(6), dev(6), epsed(6)
    real(kind=8) :: rs(6), rp, re, rpi, rg, qf, r8aux
    real(kind=8) :: pk, penpe, spe, seqd
    real(kind=8) :: etaid, agdi, agfi
    integer(kind=8) :: ndt, ndi
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
!
    data id3d/1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0/
!
    sigf(1:ndt) = yf(1:ndt)
    sigd(1:ndt) = yd(1:ndt)
    call lcopil('ISOTROPE', mod, materd(1, 1), dkooh)
    call lcopil('ISOTROPE', mod, materf(1, 1), fkooh)
    call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
!
!     RECUPERATION DES CARACTERISTIQUES MATERIAUX A t+
    ai0 = materf(4, 2)
    etais = materf(5, 2)
    k = materf(7, 2)
    n = materf(8, 2)
    p0 = materf(9, 2)
    kappa = materf(10, 2)
    r02 = materf(11, 2)
    zetaff = materf(12, 2)
    penpe = materf(13, 2)
    pk = materf(14, 2)
    pe = materf(15, 2)
    spe = materf(16, 2)
    zetag = materf(17, 2)
!     RECUPERATION DES CARACTERISTIQUES MATERIAUX A t-
!     RECUPERATION GONFLEMENT DEJA INTEGRE
    agdi = materd(19, 2)
    agfi = materf(19, 2)
!
!     RECUPERATION DES VARIABLES INTERNES A t+
    p = yf(ndt+1)
    etaif = yf(ndt+2)
!     RECUPERATION DES VARIABLES INTERNES A t-
    etaid = yd(ndt+2)
!
!     RECUPERATION DES INCREMENTS DES VARIABLES INTERNES
    dp = dy(ndt+1)
    detai = dy(ndt+2)
    dpi = dy(ndt+3)
    dg = dy(ndt+4)
!
!     RECUPERATION DE L'IRRADIATION
    irrad = materd(18, 2)
    irraf = materf(18, 2)
    dphi = irraf-irrad
!
    call lcdevi(sigf, dev)
    seqf = lcnrts(dev)
    if (seqf .eq. 0.0d0) then
        dfds(:) = 0.0d0
    else
        dfds(1:ndt) = (1.5d0/seqf)*dev(1:ndt)
    end if
    epsef(1:ndt) = matmul(fkooh(1:ndt, 1:ndt), sigf(1:ndt))
    epsed(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigd(1:ndt))
    depsa(1:ndt) = ((dp+dpi))*dfds(1:ndt)
    depsg(1:ndt) = dg*id3d(1:ndt)
!
!   RESIDU EN SIGMA, HOMOGENE A DES DEFORMATIONS
    rs(1:ndt) = epsef(1:ndt)-epsed(1:ndt)
    rs(1:ndt) = rs(1:ndt)+depsa(1:ndt)
    rs(1:ndt) = rs(1:ndt)+depsg(1:ndt)
    rs(1:ndt) = rs(1:ndt)-deps(1:ndt)
    rs(1:ndt) = (-1.d0)*rs(1:ndt)
!
!  RESIDU EN DEFORMATION PLASTIQUE
    if (p .lt. pk) then
        sf = kappa*r02
    else if (p .lt. pe) then
        sf = penpe*(p-pe)+spe
    else
        sf = k*(p+p0)**n
    end if
    if (((seqf .ge. sf) .and. (dp .ge. 0.d0)) .or. (dp .gt. r8prem())) then
        rp = -(seqf-sf)/hookf(1, 1)
    else
        rp = -dp
    end if
!
!     CONTRAINTE EQUIVALENTE A T-
    call lcdevi(sigd, dev)
    seqd = lcnrts(dev)
!
!  RESIDU EN DEFORMATION D IRRADIATION
    materf(21, 2) = 0.0d0
    r8aux = 0.0d0
    if (etaid .gt. etais) then
        rpi = -(dpi-ai0*detai)
        r8aux = zetaff*(seqd+seqf)*dphi*0.5d0
    else if (etaif .le. etais) then
        rpi = -dpi
        r8aux = zetaff*(seqd+seqf)*dphi*0.5d0
    else if (detai .gt. 0.0d0) then
        rpi = -(dpi-ai0*(etaif-etais))
        r8aux = (seqd+seqf)*dphi*0.5d0
        if (seqd .gt. 0.0d0) then
            r8aux = r8aux-(seqf-seqd)*(etais-etaid)/(2.0d0*seqd)
        end if
        r8aux = r8aux*zetaff
!        INDICATEUR DE FRANCHISSEMENT DU SEUIL ETAIS
        if (etais .gt. 0.0d0) then
            materf(21, 2) = (etaif-etais)/etais
        end if
    else
        rpi = -dpi
    end if
!
!  RESIDU PAR RAPPORT A ETA, HOMOGENE A DES DEFORMATIONS : RE
    re = -(detai-r8aux)/hookf(1, 1)
!
!  RESIDU PAR RAPPORT AU GONFLEMENT
    if (dphi .gt. 0.0d0) then
        rg = -(dg-zetag*(agfi-agdi))
    else
        rg = -dg
    end if
!
!  RESIDU PAR RAPPORT A LA DEFORMATION ( C_PLAN )
    if (mod(1:6) .eq. 'C_PLAN') then
        qf = ( &
             -hookf(3, 3)*epsef(3)-hookf(3, 1)*epsef(1)-hookf(3, 2)*epsef(2)-hookf(3, 4)*epse&
             &f(4))/hookf(1, &
             1 &
             )
    end if
!
    r(1:ndt) = rs(1:ndt)
    r(ndt+1) = rp
    r(ndt+2) = re
    r(ndt+3) = rpi
    r(ndt+4) = rg
    if (mod(1:6) .eq. 'C_PLAN') r(ndt+5) = qf
!
end subroutine
