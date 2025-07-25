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
! aslint: disable=W0413
! => real zero (affect here)
!
subroutine te0234(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/defgen.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/effi.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moytpg.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 1D
!
!     OPTION : FORC_NODA (REPRISE)
!          -----------------------------------------------------------
!
!
    integer(kind=8) :: nbres, jnbspi, nbsp, itab(7)
!
    integer(kind=8) :: nbcou, npge, jvSief, jvDisp, ivectu, icou, inte, kpki, k1
!
    real(kind=8) :: cisail, zic, coef, rhos, rhot, epsx3, gsx3, sgmsx3
!
!---- DECLARATIONS LOCALES ( RAMENEES DE TE0239.F FULL_MECA )
!
    parameter(nbres=2)
    character(len=8) :: elrefe
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: valres(nbres)
    real(kind=8) :: dfdx(3), zero, un, deux
    real(kind=8) :: test, test2, eps, nu, h, cosa, sina, cour, r, tpg
    real(kind=8) :: jacp, kappa, correc
    real(kind=8) :: eps2d(4), sigtdi(5), sigmtd(5)
    real(kind=8) :: x3
    integer(kind=8) :: nno, nnos, jgano, ndim, kp, npg, i, k, icaco, iret
    integer(kind=8) :: ipoids, ivf, idfdk, igeom, imate
    aster_logical :: testl1, testl2
    real(kind=8) :: zmin, hic
!
!
    data zero, un, deux/0.d0, 1.d0, 2.d0/
!
!-- SHIFT POUR LES COURBURES
    call elref1(elrefe)
    eps = 1.d-3
!
!DEB
!
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
!
!
!-- LECTURE DU COMPORTEMENT
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    if (nbcou .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
    if (nbcou .gt. 30) then
        call utmess('F', 'ELEMENTS3_50')
    end if
!
    npge = 3
!
!---- LECTURES STANDARDS ( RAMENEES DE TE0239.F FULL_MECA )
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCACOQU', 'L', icaco)
    h = zr(icaco)
    kappa = zr(icaco+1)
    correc = zr(icaco+2)
!---- COTE MINIMALE SUR L'EPAISSEUR
!
    zmin = -h/2.d0
!---- EPAISSEUR DE CHAQUE COUCHE
!
    hic = h/nbcou
    call jevech('PMATERC', 'L', imate)
    nomres(1) = 'E'
    nomres(2) = 'NU'
    call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, &
                itab=itab)
    jvSief = itab(1)
    nbsp = itab(7)
    if (nbsp .ne. npge*nbcou) then
        call utmess('F', 'ELEMENTS_4')
    end if
!
    call jevech('PDEPLAR', 'L', jvDisp)
!---- INITIALISATION DU VECTEUR FORCE INTERNE
!
    call jevech('PVECTUR', 'E', ivectu)
    do i = 1, 3*nno
        zr(ivectu+i-1) = 0.d0
    end do
!
    kpki = 0
    do kp = 1, npg
!-- BOUCLE SUR LES POINTS D'INTEGRATION SUR LA SURFACE
!
        k = (kp-1)*nno
        call dfdm1d(nno, zr(ipoids+kp-1), zr(idfdk+k), zr(igeom), dfdx, &
                    cour, jacp, cosa, sina)
!
        call r8inir(5, 0.d0, sigmtd, 1)
        r = zero
        call moytpg('RIGI', kp, npge, '+', tpg, &
                    iret)
!
        do i = 1, nno
            r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
        end do
!
        call rcvala(zi(imate), ' ', 'ELAS', 1, 'TEMP', &
                    [tpg], 2, nomres, valres, icodre, &
                    1)
        nu = valres(2)
        cisail = valres(1)/(un+nu)
        if (nomte .eq. 'MECXSE3') jacp = jacp*r
        test = abs(h*cour/deux)
        if (test .ge. un) correc = zero
        test2 = abs(h*cosa/(deux*r))
        if (test2 .ge. un) correc = zero
!
        testl1 = (test .le. eps .or. correc .eq. zero)
        testl2 = ( &
                 test2 .le. eps .or. correc .eq. zero .or. abs(cosa) .le. eps .or. abs(cour*r) &
                 .le. eps .or. abs(cosa-cour*r) .le. eps &
                 )
!
        do icou = 1, nbcou
            do inte = 1, npge
                if (inte .eq. 1) then
                    zic = zmin+(icou-1)*hic
                    coef = 1.d0/3.d0
                else if (inte .eq. 2) then
                    zic = zmin+hic/2.d0+(icou-1)*hic
                    coef = 4.d0/3.d0
                else
                    zic = zmin+hic+(icou-1)*hic
                    coef = 1.d0/3.d0
                end if
                x3 = zic
!
                if (testl1) then
                    rhos = 1.d0
                else
                    rhos = 1.d0+x3*cour
                end if
                if (testl2) then
                    rhot = 1.d0
                else
                    rhot = 1.d0+x3*cosa/r
                end if
!
!-- CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!-- EPSSS, EPSTT, EPSSX3
!-- (EN FONCTION DES DEFORMATIONS GENERALISEES :ESS,KSS,ETT,KTT,GS)
!-- DE L'INSTANT PRECEDANT ET DES DEFORMATIONS INCREMENTALES
!-- DE L'INSTANT PRESENT
!
                call defgen(testl1, testl2, nno, r, x3, &
                            sina, cosa, cour, zr(ivf+k), dfdx, &
                            zr(jvDisp), eps2d, epsx3)
!
!
!-- CONSTRUCTION DE LA DEFORMATION GSX3 ET DE LA CONTRAINTE SGMSX3
!
                gsx3 = 2.d0*epsx3
                sgmsx3 = cisail*kappa*gsx3/2.d0
!-- JEU D'INDICES DANS LA BOUCLE SUR LES POINTS D'INTEGRATION
!                                  DE LA SURFACE MOYENNE
!
                kpki = kpki+1
                k1 = 4*(kpki-1)
!-- CALCUL DES CONTRAINTES TILDE, ON A REMPLACE ICONTP PAR ICONTM
!

                sigtdi(1) = zr(jvSief-1+k1+1)/rhos
                sigtdi(2) = x3*zr(jvSief-1+k1+1)/rhos
                sigtdi(3) = zr(jvSief-1+k1+2)/rhot
                sigtdi(4) = x3*zr(jvSief-1+k1+2)/rhot
                sigtdi(5) = sgmsx3/rhos

!
                do i = 1, 5
                    sigmtd(i) = sigmtd(i)+sigtdi(i)*0.5d0*hic*coef
                end do
!
            end do
        end do
!
!-- CALCUL DES EFFORTS INTERIEURS
!
        call effi(nomte, sigmtd, zr(ivf+k), dfdx, jacp, &
                  sina, cosa, r, zr(ivectu))
!
    end do
!
end subroutine
