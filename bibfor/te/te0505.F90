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
subroutine te0505(option, nomte)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ntfcma.h"
#include "asterfort/projet.h"
#include "asterfort/rcfode.h"
#include "asterfort/rcfodi.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_TNL'
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
    real(kind=8) :: dfdx(9), dfdy(9), poids, r, tpg, rbid
    real(kind=8) :: uloc(2, 9), dbeta, ul(2, 10), jacob(10)
    real(kind=8) :: betaa, tpn, betai, dupgdx(9), dupgdy(9)
    real(kind=8) :: xkpt, xkptt(9), dtpgdx(9), dtpgdy(9)
    real(kind=8) :: vect(50), res(50), dbpgdx(9), dbpgdy(9)
    real(kind=8) :: xr, xrr, xaux, rr, tpg0, xk1, xk0, pn, pnp1
    integer(kind=8) :: kp, i, k, itemps, ivectt, ifon(6), igeom, imate
    integer(kind=8) :: ivite, itemp, itempi, ilagrm, ilagrp, iveres
    integer(kind=8) :: nbvf, jvalf, idim
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfde, jgano
    aster_logical :: aniso
! DEB ------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTR', 'L', itemps)
    call jevech('PTEMPER', 'L', itemp)
    call jevech('PTEMPEI', 'L', itempi)
    call jevech('PLAGRM ', 'L', ilagrm)
    call jevech('PVITESR', 'L', ivite)
    call jevech('PLAGRP ', 'E', ilagrp)
    call jevech('PVECTTR', 'E', ivectt)
    call jevech('PRESIDU', 'E', iveres)
!
    aniso = .false.
    call ntfcma(' ', zi(imate), aniso, ifon)
    nbvf = zi(ifon(1))
    jvalf = zi(ifon(1)+2)
    xr = 0.d0
    do i = 1, nbvf
        xaux = zr(jvalf+i-1)
        call rcfodi(ifon(1), xaux, rbid, xrr)
        if (xrr .gt. xr) then
            xr = xrr
        end if
    end do
    rr = 0.6d0/xr
!
    k = 0
    do i = 1, nno
        do idim = 1, 2
            k = k+1
            uloc(idim, i) = zr(ivite+k-1)
        end do
    end do
!
    do kp = 1, npg
        ul(1, kp) = 0.d0
        ul(2, kp) = 0.d0
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
        r = 0.d0
        tpg = 0.d0
        tpg0 = 0.d0
        dtpgdx(kp) = 0.d0
        dtpgdy(kp) = 0.d0
!
        do i = 1, nno
            r = r+zr(igeom+2*(i-1))*zr(ivf+k+i-1)
            tpg = tpg+zr(itempi+i-1)*zr(ivf+k+i-1)
            tpg0 = tpg0+zr(itemp+i-1)*zr(ivf+k+i-1)
            ul(1, kp) = ul(1, kp)+uloc(1, i)*zr(ivf+k+i-1)
            ul(2, kp) = ul(2, kp)+uloc(2, i)*zr(ivf+k+i-1)
            dtpgdx(kp) = dtpgdx(kp)+zr(itempi+i-1)*dfdx(i)
            dtpgdy(kp) = dtpgdy(kp)+zr(itempi+i-1)*dfdy(i)
        end do
!
        if (lteatt('AXIS', 'OUI')) poids = poids*r
        call rcfode(ifon(2), tpg, xk1, xkpt)
        call rcfode(ifon(2), tpg0, xk0, xkpt)
        pn = zr(ilagrm+kp-1)
        call rcfodi(ifon(1), pn, betaa, dbeta)
        pnp1 = pn+((tpg-betaa)*rr)
        zr(ilagrp+kp-1) = pnp1
        vect(kp) = pnp1
        jacob(kp) = poids
        xkptt(kp) = xk1-xk0
!
    end do
    call projet(2, npg, nno, vect, res)
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
        dbpgdx(kp) = 0.d0
        dbpgdy(kp) = 0.d0
        dupgdx(kp) = 0.d0
        dupgdy(kp) = 0.d0
!
        do i = 1, nno
            dupgdx(kp) = dupgdx(kp)+res(i)*dfdx(i)
            dupgdy(kp) = dupgdy(kp)+res(i)*dfdy(i)
            tpn = res(i)
            call rcfodi(ifon(1), tpn, betai, rbid)
            dbpgdx(kp) = dbpgdx(kp)+betai*dfdx(i)
            dbpgdy(kp) = dbpgdy(kp)+betai*dfdy(i)
        end do
!
    end do
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
!
        do i = 1, nno
            zr(iveres+i-1) = zr(iveres+i-1)+jacob(kp)*zr(ivf+k+i-1)*(&
                             &rr*(ul(1, kp)*dbpgdx(kp)+ul(2, kp)*dbpgdy(kp))-&
                             &   (ul(1, kp)*dupgdx(kp)+ul(2, kp)*dupgdy(kp)))+&
                             &jacob(kp)*xkptt(kp)*(dfdx(i)*dtpgdx(kp)+dfdy(i)*dtpgdy(kp))
!
        end do
!
    end do
!
! FIN ------------------------------------------------------------------
end subroutine
