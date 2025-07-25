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

subroutine tusief(option, nomte, nbrddl, b, vin, &
                  mat, pass, vtemp)
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
#include "asterfort/ppgan2.h"
#include "asterfort/prmave.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verifg.h"
#include "asterfort/vlggl.h"
#include "asterfort/vlgglc.h"
    character(len=16) :: option, nomte
! ......................................................................
!
!   FONCTION REALISEE:  CALCUL DES OPTIONS EPSI_ELGA,
!                                          DEGE_ELGA,
!                                          DEGE_ELNO,
!                                          SIEF_ELGA POUR UN TUYAU
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbres, nbrddl
    parameter(nbres=9)
    character(len=8) :: nompar
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: valres(nbres), valpar, degg(24)
    real(kind=8) :: h, a, l, e, nu, beta, cisail, g, omega, dhk, r1
    real(kind=8) :: sinfi, fi, deuxpi, r, at1, at2, vpg(4), vno(4)
    real(kind=8) :: b(4, nbrddl), c(4, 4), epsthe, hk, sigth(2), xpg(4)
    real(kind=8) :: pgl(3, 3), vin(nbrddl), vout(4), mat(4, nbrddl)
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), rayon, theta
    real(kind=8) :: vtemp(nbrddl), pass(nbrddl, nbrddl), pgl4(3, 3)
    integer(kind=8) :: nno, npg, nbcou, nbsec, icoude, ndim, jcoopg, nspg
    integer(kind=8) :: ipoids, ivf, kpgs, k, nnos
    integer(kind=8) :: imate, nbpar
    integer(kind=8) :: igau, icou, isect, i, j, jin, jout, iret, indice
    integer(kind=8) :: lorien, icoud2, mmt, jnbspi
    integer(kind=8) :: nddl, m, idfdk, jdfd2, jgano
!
!
    integer(kind=8), parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8), parameter :: noms_cara1(nb_cara1) = (/'R1 ', 'EP1'/)
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, &
                     jdfd2=jdfd2, jgano=jgano)
!
    deuxpi = 2.d0*r8pi()
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    nbsec = zi(jnbspi-1+2)
!
    call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    r1 = vale_cara1(1)
    h = vale_cara1(2)
    a = r1-h/2.d0
!
    m = 3
    if (nomte .eq. 'MET6SEG3') m = 6
!
    do i = 1, npg
        xpg(i) = zr(jcoopg-1+i)
    end do
!
!     --- RECUPERATION DES ORIENTATIONS ---
!
    call jevech('PCAORIE', 'L', lorien)
    call carcou(zr(lorien), l, pgl, rayon, theta, &
                pgl1, pgl2, pgl3, pgl4, nno, &
                omega, icoud2)
    if (icoud2 .ge. 10) then
        icoude = icoud2-10
        mmt = 0
    else
        icoude = icoud2
        mmt = 1
    end if
!
    if (option .eq. 'SIEF_ELGA') then
!
        call jevech('PMATERC', 'L', imate)
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nomres(3) = 'ALPHA'
        nspg = (2*nbsec+1)*(2*nbcou+1)
        iret = 0
        call moytem('RIGI', npg, nspg, '+', valpar, &
                    iret)
        if (iret .ne. 0) valpar = 0.d0
        nbpar = 1
        nompar = 'TEMP'
        call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                    ' ', 'ELAS', nbpar, nompar, [valpar], &
                    2, nomres, valres, icodre, 1)
        e = valres(1)
        nu = valres(2)
!
!        DEFINITION DE LA MATRICE DE COMPORTEMENT C
!
        beta = e/(1.d0-nu**2)
        g = e/(2.d0*(1.d0+nu))
        cisail = 1.d0
!
        c(1, 1) = beta
        c(1, 2) = nu*beta
        c(1, 3) = 0.d0
        c(1, 4) = 0.d0
!
        c(2, 1) = nu*beta
        c(2, 2) = beta
        c(2, 3) = 0.d0
        c(2, 4) = 0.d0
!
        c(3, 1) = 0.d0
        c(3, 2) = 0.d0
        c(3, 3) = g
        c(3, 4) = 0.d0
!
        c(4, 1) = 0.d0
        c(4, 2) = 0.d0
        c(4, 3) = 0.d0
        c(4, 4) = g*cisail
!
!        FIN DE LA CONSTRUCTION DE LA MATRICE DE COMPORTEMENT C
    end if
!
    call jevech('PDEPLAR', 'L', jin)
    do i = 1, nbrddl
        vin(i) = zr(jin-1+i)
    end do
    if (icoude .eq. 0) then
        call vlggl(nno, nbrddl, pgl, vin, 'GL', &
                   pass, vtemp)
    else
        call vlgglc(nno, nbrddl, pgl1, pgl2, pgl3, &
                    pgl4, vin, 'GL', pass, vtemp)
    end if
!
    if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jout)
!
    else if (option .eq. 'DEGE_ELNO' .or. option .eq. 'DEGE_ELGA') then
!
        if (option .eq. 'DEGE_ELGA') call jevech('PDEFOPG', 'E', jout)
        if (option .eq. 'DEGE_ELNO') call jevech('PDEFOGR', 'E', jout)
        call r8inir(24, 0.d0, degg, 1)
        nbrddl = nno*(6+3+6*(m-1))
        nddl = (6+3+6*(m-1))
!
    else if (option .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jout)
!
    else
        call utmess('F', 'ELEMENTS4_49', sk=option)
!
    end if
!
    kpgs = 0
    sigth(1) = 0.d0
    sigth(2) = 0.d0
    nspg = (2*nbsec+1)*(2*nbcou+1)
!
! --- BOUCLE SUR LES POINTS DE GAUSS
! ---- BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
    do igau = 1, npg
        if (option .eq. 'SIEF_ELGA') then
            call verifg('RIGI', igau, nspg, '+', zi(imate), &
                        epsthe)
            at1 = (c(1, 1)+c(1, 2))*epsthe
            at2 = (c(2, 1)+c(2, 2))*epsthe
        end if
!
        if ((option .eq. 'SIEF_ELGA') .or. (option .eq. 'EPSI_ELGA')) then
! --- BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
            do icou = 1, 2*nbcou+1
                if (mmt .eq. 0) then
                    r = a
                else
                    r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
                end if
                do isect = 1, 2*nbsec+1
                    kpgs = kpgs+1
!
                    if (icoude .eq. 0) then
                        call bcoude(igau, icou, isect, l, h, &
                                    a, m, nno, nbcou, nbsec, &
                                    zr(ivf), zr(idfdk), zr(jdfd2), mmt, b)
                    else if (icoude .eq. 1) then
                        fi = (isect-1)*deuxpi/(2.d0*nbsec)
                        sinfi = sin(fi)
                        l = theta*(rayon+r*sinfi)
                        call bcoudc(igau, icou, isect, h, a, &
                                    m, omega, xpg, nno, nbcou, &
                                    nbsec, zr(ivf), zr(idfdk), zr(jdfd2), rayon, &
                                    theta, mmt, b)
                    end if
!
                    if (option .eq. 'EPSI_ELGA') then
                        do i = 1, 4
                            do j = 1, nbrddl
                                mat(i, j) = b(i, j)
                            end do
                        end do
                    else if (option .eq. 'SIEF_ELGA') then
                        sigth(1) = at1
                        sigth(2) = at2
                        call promat(c, 4, 4, 4, b, &
                                    4, 4, nbrddl, mat)
                    end if
                    iret = 0
                    call prmave(0, mat, 4, 4, nbrddl, &
                                vin, nbrddl, vout, 4, iret)
!
!  STOCKAGE DU VECTEUR VOUT
!
                    indice = jout-1+6*(kpgs-1)
                    zr(indice+1) = vout(1)-sigth(1)
                    zr(indice+2) = vout(2)-sigth(2)
                    zr(indice+3) = 0.d0
                    zr(indice+4) = vout(3)
                    zr(indice+5) = vout(4)
                    zr(indice+6) = 0.d0
                end do
            end do
        else if (option .eq. 'DEGE_ELNO' .or. option .eq. 'DEGE_ELGA') then
!            DEFORMATIONS GENERALISEES DE POUTRE
!
            do k = 1, nno
                if (icoude .eq. 1) then
                    l = theta*rayon
                end if
                hk = zr(ivf-1+nno*(igau-1)+k)
                dhk = zr(ivf-1+nno*npg+nno*(igau-1)+k)*(2.d0/l)
                dhk = zr(idfdk-1+nno*(igau-1)+k)*(2.d0/l)
!
!           EPXX=DU/DX
                degg(6*(igau-1)+1) = degg(6*(igau-1)+1)+dhk*vin(nddl*(k-1)+1)
!              GAXY=DV/DX -DRZ
                degg(6*(igau-1)+2) = degg(6*(igau-1)+2)+dhk*vin(nddl*(k-1)+2)- &
                                     hk*vin(nddl*(k-1)+6)
!              GAXZ=DW/DX +DRY
                degg(6*(igau-1)+3) = degg(6*(igau-1)+3)+dhk*vin(nddl*(k-1)+3)+ &
                                     hk*vin(nddl*(k-1)+5)
!              GAT=D(DRX)/DX
                degg(6*(igau-1)+4) = degg(6*(igau-1)+4)+dhk*vin(nddl*(k-1)+4)
!              KY=D(DRY)/DX
                degg(6*(igau-1)+5) = degg(6*(igau-1)+5)+dhk*vin(nddl*(k-1)+5)
!              KZ=D(DRZ)/DX
                degg(6*(igau-1)+6) = degg(6*(igau-1)+6)+dhk*vin(nddl*(k-1)+6)
            end do
        end if
    end do
!
    if (option .eq. 'DEGE_ELNO' .or. option .eq. 'DEGE_ELGA') then
        do i = 1, 6
            do igau = 1, npg
                vpg(igau) = degg(6*(igau-1)+i)
                if (option .eq. 'DEGE_ELGA') zr(jout+6*(igau-1)+i-1) = vpg(igau)
            end do
            if (option .eq. 'DEGE_ELNO') then
                call ppgan2(jgano, 1, 1, vpg, vno)
                do j = 1, nno
                    zr(jout+6*(j-1)+i-1) = vno(j)
                end do
            end if
        end do
    end if
!
end subroutine
