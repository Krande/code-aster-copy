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

subroutine te0221(option, nomte)
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
    implicit none
#include "jeveux.h"
!
#include "asterfort/dfdm1d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moytpg.h"
#include "asterfort/rcvala.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          COQUE 1D
!                          OPTION : 'RIGI_MECA      '
!                          ELEMENT: MECXSE3
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: jgano, nbres, ndim, nnos
    real(kind=8) :: gss
!-----------------------------------------------------------------------
    parameter(nbres=2)
    character(len=16) :: nomres(nbres)
    character(len=8) :: elrefe, nompar
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: valres(nbres), valpar
    real(kind=8) :: dfdx(3), zero, un, deux, trois, douze
    real(kind=8) :: test, test2, eps, nu, h, cosa, sina, cour, r
    real(kind=8) :: coefxx, coefyy, coefxy, coeff1, coeff2
    real(kind=8) :: css, ctt, cts, dss, dts, dtt, bss, btt, bts, vfi, vfj
    real(kind=8) :: c1, c2, c3, cons, cons2, jacp, kappa, correc
    integer(kind=8) :: nno, kp, npg, imatuu, icaco
    integer(kind=8) :: ii, jj, i, j, k, ij1, ij2, ij3, kd1, kd2, kd3
    integer(kind=8) :: ipoids, ivf, idfdk, igeom, imate
    integer(kind=8) :: nbpar, iret
!
!
    data zero, un, deux, trois, douze/0.d0, 1.d0, 2.d0, 3.d0, 12.d0/
    call elref1(elrefe)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
    eps = 1.d-3
!
!
!
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icaco)
    h = zr(icaco)
    kappa = zr(icaco+1)
    correc = zr(icaco+2)
!
    call jevech('PMATERC', 'L', imate)
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call jevech('PMATUUR', 'E', imatuu)
!
    do kp = 1, npg
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, jacp, cosa, sina)
        r = zero
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
        end do
!===============================================================
!     -- RECUPERATION DE LA TEMPERATURE POUR LE MATERIAU:
!
!     -- SI LA TEMPERATURE EST CONNUE AUX POINTS DE GAUSS :
!
        call moytpg('RIGI', kp, 3, '+', valpar, &
                    iret)
        nbpar = 1
        nompar = 'TEMP'
!===============================================================
        test = abs(h*cour/deux)
        if (test .ge. un) correc = zero
        call rcvala(zi(imate), ' ', 'ELAS', nbpar, nompar, &
                    [valpar], 2, nomres, valres, icodre, 1)
!
        nu = valres(2)
        coeff1 = valres(1)/(un-(nu*nu))
        coeff2 = valres(1)/(un+nu)
!
!     CALCUL DES COEFFICIENTS RIGIDITE COQUE AXI
!     COUR : COURBURE RS
!     H    : EPAISSEUR
!     COSA : COSINUS ANGLE ALPHA: NORMALE/ HORIZONTALE
!     R    : RAYON VECTEUR
!
        c1 = -cour*h/(un-(cour*h/deux)**2)
        if (nomte .eq. 'MECXSE3') then
            c2 = -cosa*h/(r**2-(cosa*h/deux)**2)
            if (test .le. eps .or. correc .eq. zero) then
                css = coeff1*h
                bss = zero
                dss = coeff1*h**3/douze
                gss = coeff2*kappa*h/deux
            else
                cons = log((un+(h*cour/deux))/(un-(h*cour/deux)))
                css = ((c1+cons)*cosa/(r*cour**2)+cons/cour)
                bss = -((c1+deux*cons-h*cour)*cosa/(r*cour**3)+cons/(cour**2)-h/cour)*coeff1
                dss = ((cons-h*cour)/cour**3+cosa*(c1+cons*trois-deux*h*cour)/(cour**4*r) &
                       )*coeff1
                gss = css*coeff2*kappa/deux
                css = css*coeff1
            end if
            test2 = abs(h*cosa/(deux*r))
            if (test2 .ge. un) correc = zero
            if (test2 .le. eps .or. correc .eq. zero) then
                ctt = coeff1*h
                btt = zero
                dtt = coeff1*h**3/douze
            else
                cons2 = log((r+(h*cosa/deux))/(r-(h*cosa/deux)))
                c3 = r/cosa
                ctt = (cons2*c3+cour*c3*c3*(r*c2+cons2))*coeff1
                btt = -(c3**3*c2*r*cour+cour*c3*c3*(deux*cons2*c3-h)+c3*(c3*cons2-h))*coeff1
                dtt = ( &
                      cour*c3**4*r*c2-h*c3*c3*(un+deux*cour*c3)+cons2*c3**3*(un+trois*c3*cour) &
                      )*coeff1
            end if
            if (abs(cosa) .le. eps .or. abs(cour*r) .le. eps .or. abs(cosa-cour*r) .le. &
                eps .or. correc .eq. zero) then
                cts = coeff1*h
                bts = zero
                dts = coeff1*h**3/douze
            else
                c3 = r/cosa
                cts = ((cour*r**2*cons2-cosa**2*cons/cour)/(cosa*(cour*r-cosa)))*coeff1
                bts = -( &
                      -h*(un/cour+c3)+cosa*cons/(cour**2*(cosa-cour*r))+cons2*cour*cosa*c3**3&
                      &/(r*cour-cosa) &
                      )*coeff1
                dts = ( &
                      -h*(un+c3*cour+c3*c3*cour*cour)/cour**2+cons2*cour*c3**3*r/(r*cour-cosa)&
                      &-cons/(cour**3*(cour*c3-un)) &
                      )*coeff1
            end if
            jacp = jacp*r
        end if
!
        coefxx = sina*sina
        coefyy = cosa*cosa
        coefxy = -cosa*sina
        kd1 = 5
        kd2 = 3
        kd3 = 2
        do i = 1, 3*nno, 3
            kd1 = kd1+3*i-6
            kd2 = kd2+3*i-3
            kd3 = kd3+3*i
            ii = (i+2)/3
            do j = 1, i, 3
                jj = (j+2)/3
                ij1 = imatuu+kd1+j-3
                ij2 = imatuu+kd2+j-3
                ij3 = imatuu+kd3+j-3
                vfi = zr(ivf+k+ii-1)
                vfj = zr(ivf+k+jj-1)
                zr(ij1) = zr(ij1)+dfdx(ii)*dfdx(jj)*jacp*(coefxx*css+coefyy*gss)
                zr(ij2) = zr(ij2)+dfdx(ii)*dfdx(jj)*jacp*coefxy*(css-gss)
                zr(ij2+1) = zr(ij2+1)+dfdx(ii)*dfdx(jj)*jacp*(coefyy*css+coefxx*gss)
                zr(ij3) = zr(ij3)+dfdx(jj)*jacp*(cosa*gss*vfi+sina*bss*dfdx(ii))
                zr(ij3+1) = zr(ij3+1)-dfdx(jj)*jacp*(cosa*bss*dfdx(ii)-sina*gss*vfi)
                zr(ij3+2) = zr(ij3+2)+jacp*(dss*dfdx(ii)*dfdx(jj)+gss*vfi*vfj)
            end do
            do j = 1, i-3, 3
                jj = (j+2)/3
                ij1 = imatuu+kd1+j-3
                ij2 = imatuu+kd2+j-3
                vfi = zr(ivf+k+ii-1)
                vfj = zr(ivf+k+jj-1)
                zr(ij1+1) = zr(ij1+1)+dfdx(ii)*dfdx(jj)*jacp*coefxy*(css-gss)
                zr(ij1+2) = zr(ij1+2)+dfdx(ii)*jacp*(cosa*gss*vfj+sina*bss*dfdx(jj))
                zr(ij2+2) = zr(ij2+2)+dfdx(ii)*jacp*(sina*gss*vfj-cosa*bss*dfdx(jj))
            end do
        end do
!
        if (nomte .eq. 'MECXSE3') then
            kd1 = 5
            kd2 = 3
            kd3 = 2
            do i = 1, 3*nno, 3
                kd1 = kd1+3*i-6
                kd2 = kd2+3*i-3
                kd3 = kd3+3*i
                ii = (i+2)/3
                do j = 1, i, 3
                    jj = (j+2)/3
                    ij1 = imatuu+kd1+j-3
                    ij2 = imatuu+kd2+j-3
                    ij3 = imatuu+kd3+j-3
                    vfi = zr(ivf+k+ii-1)
                    vfj = zr(ivf+k+jj-1)
                    zr(ij1) = zr(ij1)+jacp*(ctt*vfi*vfj/(r*r)-nu*cts*sina*(dfdx(ii)*vfj+dfd&
                              &x(jj)*vfi)/r)
                    zr(ij2) = zr(ij2)+jacp*nu*cts*cosa*dfdx(ii)*vfj/r
                    zr(ij3) = zr(ij3)+jacp*nu*bts*deux*(coefxx*vfi*dfdx(jj)-dfdx(ii)*vfj)/r-&
                              & jacp*btt*sina*vfi*vfj/r
                    zr(ij3+1) = zr(ij3+1)+jacp*nu*bts*deux*coefxy*dfdx(jj)*vfi/r
                    zr(ij3+2) = zr(ij3+2)+jacp*sina*(dtt*sina*vfi*vfj/(r*r)+nu*dts*(vfi*df&
                                &dx(jj)+dfdx(ii)*vfj)/r)
                end do
                do j = 1, i-3, 3
                    jj = (j+2)/3
                    ij1 = imatuu+kd1+j-3
                    ij2 = imatuu+kd2+j-3
                    vfi = zr(ivf+k+ii-1)
                    vfj = zr(ivf+k+jj-1)
                    zr(ij1+1) = zr(ij1+1)+jacp*nu*cts*cosa*dfdx(jj)*vfi/r
                    zr(ij1+2) = zr(ij1+2)+jacp*nu*bts*deux*(coefxx*vfj*dfdx(ii)-dfdx(jj)*vfi)&
                                &/r-jacp*btt*sina*vfj*vfi/r
                    zr(ij2+2) = zr(ij2+2)+jacp*nu*bts*deux*coefxy*dfdx(ii)*vfj/r
                end do
            end do
        end if
!
    end do
end subroutine
