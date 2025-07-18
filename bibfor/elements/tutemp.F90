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

subroutine tutemp(option, nomte, nbrddl, f, b, &
                  vout, pass, vtemp)
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/bcoudc.h"
#include "asterfort/bcoude.h"
#include "asterfort/carcou.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/moytem.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verifg.h"
#include "asterfort/vlggl.h"
#include "asterfort/vlgglc.h"
    character(len=16) :: option
! ......................................................................
!
!    - FONCTION REALISEE:  CALCUL DU SECOND MEMBRE : TRAVAIL DE LA
!                          DILATATION THERMIQUE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
! ......................................................................
!
    integer(kind=8) :: nbrddl, nbsecm, nbcoum
    parameter(nbsecm=32, nbcoum=10)
    real(kind=8) :: h, a, l, valpar, beta, r, r1
    real(kind=8) :: poicou(2*nbcoum+1), poisec(2*nbsecm+1)
    real(kind=8) :: f(nbrddl), b(4, nbrddl), vout(nbrddl), sig(4)
    real(kind=8) :: pi, deuxpi, fi, e, nu, valres(3)
    real(kind=8) :: pgl(3, 3), pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), omega
    real(kind=8) :: c(2, 2), coe1, pgl4(3, 3)
    real(kind=8) :: poids, rayon, theta, sinfi, xpg(4)
    real(kind=8) :: vtemp(nbrddl), pass(nbrddl, nbrddl)
    integer(kind=8) :: codres(3)
    character(len=8) :: nompar
    character(len=16) :: nomte, nomres(3)
    character(len=32) :: phenom
    integer(kind=8) :: nno, npg, nbcou, nbsec, m, lorien, icoude, i
    integer(kind=8) :: ipoids, ivf, icou, nbpar
    integer(kind=8) :: igeom, jout, imate, j
    integer(kind=8) :: igau, isect, icoud2, mmt, nspg
    integer(kind=8) :: jnbspi, iret, iret2
    integer(kind=8) :: ndim, nnos, jcoopg, idfdk, jdfd2, jgano
!
    integer(kind=8), parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8), parameter :: noms_cara1(nb_cara1) = (/'R1 ', 'EP1'/)
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, &
                     jdfd2=jdfd2, jgano=jgano)
!
    pi = r8pi()
    deuxpi = 2.d0*pi
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    nbsec = zi(jnbspi-1+2)
!
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCAORIE', 'L', lorien)
    call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    r1 = vale_cara1(1)
    h = vale_cara1(2)
    a = r1-h/2.d0
!
! A= RMOY, H = EPAISSEUR
! RINT = RAYON INTERIEUR
!
!
    m = 3
    if (nomte .eq. 'MET6SEG3') m = 6
!
!
    do i = 1, npg
        xpg(i) = zr(jcoopg-1+i)
    end do
!
!
!     LES POIDS POUR L'INTEGRATION DANS L'EPAISSEUR
!
    poicou(1) = 1.d0/3.d0
    do i = 1, nbcou-1
        poicou(2*i) = 4.d0/3.d0
        poicou(2*i+1) = 2.d0/3.d0
    end do
    poicou(2*nbcou) = 4.d0/3.d0
    poicou(2*nbcou+1) = 1.d0/3.d0
!
!     LES POIDS POUR L'INTEGRATION SUR LA CIRCONFERENCE
!
    poisec(1) = 1.d0/3.d0
    do i = 1, nbsec-1
        poisec(2*i) = 4.d0/3.d0
        poisec(2*i+1) = 2.d0/3.d0
    end do
    poisec(2*nbsec) = 4.d0/3.d0
    poisec(2*nbsec+1) = 1.d0/3.d0
!
!     FIN DES POIDS D'INTEGRATION
!
! CALCUL DE L = LONGUEUR DU COUDE
!
    call carcou(zr(lorien), l, pgl, rayon, theta, &
                pgl1, pgl2, pgl3, pgl4, nno, &
                omega, icoud2)
!
    if (icoud2 .ge. 10) then
        icoude = icoud2-10
        mmt = 0
    else
        icoude = icoud2
        mmt = 1
    end if
!
!---- RECUPERATION TEMPERATURE
!===============================================================
!          -- RECUPERATION DE LA TEMPERATURE :
!     -- SI LA TEMPERATURE N'EST PAS DONNEE:
    nspg = (2*nbsec+1)*(2*nbcou+1)
    iret2 = 0
    call moytem('RIGI', npg, nspg, '+', valpar, &
                iret2)
    nbpar = 1
    nompar = 'TEMP'
!===============================================================
!
!---- RECUPERATION DU COMPORTEMENT
!
    call jevech('PMATERC', 'L', imate)
    call rccoma(zi(imate), 'ELAS', 1, phenom, codres(1))
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
    call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                ' ', phenom, nbpar, nompar, [valpar], &
                2, nomres, valres, codres, 1)
    e = valres(1)
    nu = valres(2)
!
!
! DEFINITION DE LA MATRICE DE COMPORTEMENT C
! POUR LA DILATATION
!
    beta = 1.d0/(1.d0-nu**2)
!
    c(1, 1) = e*beta
    c(1, 2) = e*nu*beta
!
    c(2, 1) = e*nu*beta
    c(2, 2) = e*beta
!
!
!     FIN DE LA CONSTRUCTION DE LA MATRICE DE COMPORTEMENT C
!
    if (option .eq. 'CHAR_MECA_TEMP_R') then
!
!----- CAS DILATATION THERMIQUE
!
        do i = 1, nbrddl
            f(i) = 0.d0
        end do
!
!     DEBUT CONSTRUCTION DE B
!
! BOUCLE SUR LES POINTS DE GAUSS
        nspg = (2*nbsec+1)*(2*nbcou+1)
        do igau = 1, npg
            call verifg('RIGI', igau, nspg, '+', zi(imate), &
                        coe1, iret)
            sig(1) = (c(1, 1)+c(1, 2))*coe1
            sig(2) = (c(2, 1)+c(2, 2))*coe1
            sig(3) = 0.d0
            sig(4) = 0.d0
!
! BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
!
            do icou = 1, 2*nbcou+1
                if (mmt .eq. 0) then
                    r = a
                else
                    r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
                end if
!
!
! BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
!
                do isect = 1, 2*nbsec+1
                    if (icoude .eq. 0) then
                        call bcoude(igau, icou, isect, l, h, &
                                    a, m, nno, nbcou, nbsec, &
                                    zr(ivf), zr(idfdk), zr(jdfd2), mmt, b)
                    else if (icoude .eq. 1) then
                        fi = (isect-1)*deuxpi/(2.d0*nbsec)
!
                        sinfi = sin(fi)
                        l = theta*(rayon+r*sinfi)
                        call bcoudc(igau, icou, isect, h, a, &
                                    m, omega, xpg, nno, nbcou, &
                                    nbsec, zr(ivf), zr(idfdk), zr(jdfd2), rayon, &
                                    theta, mmt, b)
                    end if
                    do j = 1, nbrddl
                        vout(j) = b(1, j)*sig(1)+b(2, j)*sig(2)
                    end do
!
!  STOCKAGE DU VECTEUR VOUT DANS FI
!
                    poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deuxpi/(4.&
                            &d0*nbcou*nbsec)*r
!
                    do i = 1, nbrddl
                        f(i) = f(i)+vout(i)*poids
                    end do
                end do
            end do
        end do
        if (icoude .eq. 0) then
            call vlggl(nno, nbrddl, pgl, f, 'LG', &
                       pass, vtemp)
        else
            call vlgglc(nno, nbrddl, pgl1, pgl2, pgl3, &
                        pgl4, f, 'LG', pass, vtemp)
        end if
        call jevech('PVECTUR', 'E', jout)
        do i = 1, nbrddl
            zr(jout-1+i) = f(i)
        end do
    end if
end subroutine
