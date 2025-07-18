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
subroutine ef0587(nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/carcou.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ppgan2.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: nomte
!     CALCUL DE EFGE_ELNO
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbcoum, nbsecm, jnbspi
    real(kind=8) :: h, a, r1
    parameter(nbsecm=32, nbcoum=10)
    real(kind=8) :: poicou(2*nbcoum+1), poisec(2*nbsecm+1)
    real(kind=8) :: pi, deuxpi, sig(6), fno(4, 6)
    real(kind=8) :: efg(6), alphaf, betaf, alpham, betam, xa, xb, xc, xd
    real(kind=8) :: pgl(3, 3), pgl4(3, 3), vno(4), vpg(4)
    real(kind=8) :: cosfi, sinfi, hk(4, 4)
    real(kind=8) :: fi, poids, r, omega
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), rayon, theta, l
    real(kind=8) :: cp(2, 2), cv(2, 2), co(4, 4), si(4, 4), tk(4), xpg(4)
    integer(kind=8) :: nno, nnos, jgano, ndim, npg, nbcou, nbsec, lorien
    integer(kind=8) :: ipoids, ivf, icoude, ic, kp, jin, jcoopg, jdfd2
    integer(kind=8) :: i1, i2, ih, idfdk
    integer(kind=8) :: igau, icou, isect, i, jout, ino
    integer(kind=8) :: indice, k, ip, icoud2, mmt
    integer(kind=8) :: kpgs
!
    integer(kind=8) :: vali
!
    integer(kind=8), parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8) :: noms_cara1(nb_cara1)
    data noms_cara1/'R1', 'EP1'/
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jdfd2=jdfd2, &
                     jgano=jgano)
!
!
    pi = r8pi()
    deuxpi = 2.d0*pi
!
!=====RECUPERATION NOMBRE DE COUCHES ET DE SECTEURS ANGULAIRES
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    nbsec = zi(jnbspi-1+2)
    if (nbcou*nbsec .le. 0) then
        call utmess('F', 'ELEMENTS4_46')
    end if
    if (nbcou .gt. nbcoum) then
        vali = nbcoum
        call utmess('F', 'ELEMENTS5_2', si=vali)
    end if
    if (nbsec .gt. nbsecm) then
        vali = nbsecm
        call utmess('F', 'ELEMENTS5_3', si=vali)
    end if
!
!
!
!     PREMIERE FAMILLE DE POINTS DE GAUSS POUR LES CHAMPS
!
!
    do i = 1, npg
        xpg(i) = zr(jcoopg-1+i)
    end do
!
!  LES POIDS POUR L'INTEGRATION DANS L'EPAISSEUR
!
    poicou(1) = 1.d0/3.d0
    do i = 1, nbcou-1
        poicou(2*i) = 4.d0/3.d0
        poicou(2*i+1) = 2.d0/3.d0
    end do
    poicou(2*nbcou) = 4.d0/3.d0
    poicou(2*nbcou+1) = 1.d0/3.d0
!
!  LES POIDS POUR L'INTEGRATION SUR LA CIRCONFERENCE
!
    poisec(1) = 1.d0/3.d0
    do i = 1, nbsec-1
        poisec(2*i) = 4.d0/3.d0
        poisec(2*i+1) = 2.d0/3.d0
    end do
    poisec(2*nbsec) = 4.d0/3.d0
    poisec(2*nbsec+1) = 1.d0/3.d0
!
!   FIN DES POIDS D'INTEGRATION
!
!
!  CONTRUCTION DE LA MATRICE H(I,J) = MATRICE DES VALEURS DES
!  FONCTIONS DE FORMES AUX POINT DE GAUSS
!
    do k = 1, nno
        do igau = 1, npg
            hk(k, igau) = zr(ivf-1+nno*(igau-1)+k)
        end do
    end do
!
!
!
    call jevech('PCAORIE', 'L', lorien)
    call jevech('PCONTRR', 'L', jin)
    call jevech('PEFFORR', 'E', jout)
!
!       -- A= RMOY, H = EPAISSEUR
    call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    r1 = vale_cara1(1)
    h = vale_cara1(2)
    a = r1-h/2.d0
!
!       -- ORIENTATION :
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
    if (nno .eq. 3) then
        tk(1) = 0.d0
        tk(2) = theta
        tk(3) = theta/2.d0
    else if (nno .eq. 4) then
        tk(1) = 0.d0
        tk(2) = theta
        tk(3) = theta/3.d0
        tk(4) = 2.d0*theta/3.d0
    end if
!
!
    kpgs = 0
!       -- CALCUL DES EFFORTS SUR LES POINTS DE GAUSS (VFNO)
    do igau = 1, npg
!
        do i = 1, 6
            efg(i) = 0.d0
        end do
!
!         -- BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
        do icou = 1, 2*nbcou+1
            if (mmt .eq. 0) then
                r = a
            else
                r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
            end if
!
!           -- BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
            do isect = 1, 2*nbsec+1
!
                kpgs = kpgs+1
                fi = (isect-1)*deuxpi/(2.d0*nbsec)
                cosfi = cos(fi)
                sinfi = sin(fi)
!
                indice = jin-1+6*(kpgs-1)
                sig(1) = zr(indice+1)
                sig(2) = zr(indice+2)
                sig(3) = zr(indice+4)
                sig(4) = zr(indice+5)
!
                poids = poicou(icou)*poisec(isect)*h*deuxpi/(4.d0* &
                                                             nbcou*nbsec)*r
!
!
                efg(1) = efg(1)+poids*sig(1)
                efg(2) = efg(2)-poids*(sinfi*sig(4)+cosfi*sig(3))
                efg(3) = efg(3)+poids*(sinfi*sig(3)-cosfi*sig(4))
!
                efg(4) = efg(4)-poids*sig(3)*r
                efg(5) = efg(5)-poids*sig(1)*r*cosfi
                efg(6) = efg(6)+poids*sig(1)*r*sinfi
            end do
        end do
!
        do i = 1, 6
            fno(igau, i) = efg(i)
        end do
    end do
!
!
    if ((nno .eq. 3) .and. (npg .eq. 3)) then
!         -- LA BELLE PROGRAMMATION DE PATRICK
!            EST-ELLE MIEUX QUE PPGAN2 ?
        do igau = 1, npg
            do ino = 1, nno
                if (icoude .eq. 0) then
                    co(igau, ino) = 1.d0
                    si(igau, ino) = 0.d0
                else
                    co(igau, ino) = cos((1.d0+xpg(igau))*theta/2.d0-tk( &
                                        ino))
                    si(igau, ino) = sin((1.d0+xpg(igau))*theta/2.d0-tk( &
                                        ino))
                end if
            end do
        end do
        do ino = 1, nno
            if (ino .eq. 1) then
                ih = 2
                ip = 1
                i1 = 1
                i2 = 3
            else if (ino .eq. 2) then
                ih = 1
                ip = 2
                i1 = 3
                i2 = 1
            else
                do i = 1, 6
                    efg(i) = fno(2, i)
                end do
                goto 140
!
            end if
!
            cp(1, 1) = co(1, ih)*co(1, 3)+si(1, ih)*si(1, 3)
            cp(1, 2) = -co(1, ih)*si(1, 3)+si(1, ih)*co(1, 3)
            cp(2, 1) = -cp(1, 2)
            cp(2, 2) = cp(1, 1)
            cv(1, 1) = co(3, ih)*co(3, 3)+si(3, ih)*si(3, 3)
            cv(1, 2) = -co(3, ih)*si(3, 3)+si(3, ih)*co(3, 3)
            cv(2, 1) = -cp(1, 2)
            cv(2, 2) = cp(1, 1)
!
            alphaf = hk(ih, 3)*(co(1, ih)*fno(1, 1)+si(1, ih)*fno(1, 2))- &
                     hk(ih, 3)*hk(3, 1)*(cp(1, 1)*fno(2, 1)+cp(1, 2)*fno(2, 2))- &
                     hk(ih, 1)*(co(3, ih)*fno(3, 1)+si(3, ih)*fno(3, 2))+hk(ih, 1)* &
                     hk(3, 3)*(cv(1, 1)*fno(2, 1)+cv(1, 2)*fno(2, 2))
!
            betaf = hk(ih, 3)*(-si(1, ih)*fno(1, 1)+co(1, ih)*fno(1, 2))- &
                    hk(ih, 3)*hk(3, 1)*(cp(2, 1)*fno(2, 1)+cp(2, 2)*fno(2, 2))- &
                    hk(ih, 1)*(-si(3, ih)*fno(3, 1)+co(3, ih)*fno(3, 2))+hk(ih, 1)* &
                    hk(3, 3)*(cv(2, 1)*fno(2, 1)+cv(2, 2)*fno(2, 2))
!
            alpham = hk(ih, 3)*(co(1, ih)*fno(1, 4)+si(1, ih)*fno(1, 5))- &
                     hk(ih, 3)*hk(3, 1)*(cp(1, 1)*fno(2, 4)+cp(1, 2)*fno(2, 5))- &
                     hk(ih, 1)*(co(3, ih)*fno(3, 4)+si(3, ih)*fno(3, 5))+hk(ih, 1)* &
                     hk(3, 3)*(cv(1, 1)*fno(2, 4)+cv(1, 2)*fno(2, 5))
!
            betam = hk(ih, 3)*(-si(1, ih)*fno(1, 4)+co(1, ih)*fno(1, 5))- &
                    hk(ih, 3)*hk(3, 1)*(cp(2, 1)*fno(2, 4)+cp(2, 2)*fno(2, 5))- &
                    hk(ih, 1)*(-si(3, ih)*fno(3, 4)+co(3, ih)*fno(3, 5))+hk(ih, 1)* &
                    hk(3, 3)*(cv(2, 1)*fno(2, 4)+cv(2, 2)*fno(2, 5))
!
            cp(1, 1) = co(1, ih)*co(1, ip)+si(1, ih)*si(1, ip)
            cp(1, 2) = -co(1, ih)*si(1, ip)+si(1, ih)*co(1, ip)
            cp(2, 1) = -cp(1, 2)
            cp(2, 2) = cp(1, 1)
            cv(1, 1) = co(3, ih)*co(3, ip)+si(3, ih)*si(3, ip)
            cv(1, 2) = -co(3, ih)*si(3, ip)+si(3, ih)*co(3, ip)
            cv(2, 1) = -cp(1, 2)
            cv(2, 2) = cp(1, 1)
!
            xa = hk(ip, 1)*hk(ih, 3)*cp(1, 1)-hk(ip, 3)*hk(ih, 1)*cv(1, 1)
            xb = hk(ip, 1)*hk(ih, 3)*cp(1, 2)-hk(ip, 3)*hk(ih, 1)*cv(1, 2)
            xc = hk(ip, 1)*hk(ih, 3)*cp(2, 1)-hk(ip, 3)*hk(ih, 1)*cv(2, 1)
            xd = hk(ip, 1)*hk(ih, 3)*cp(2, 2)-hk(ip, 3)*hk(ih, 1)*cv(2, 2)
!
            efg(1) = (xd*alphaf-xb*betaf)/(xa*xd-xb*xc)
            efg(2) = (-xc*alphaf+xa*betaf)/(xa*xd-xb*xc)
            efg(3) = (hk(ih, i2)*fno(i1, 3)-hk(ih, i1)*fno(i2, 3)-fno(2, 3)* &
                      (hk(3, i1)*hk(ih, i2)-hk(3, i2)*hk(ih, i1)))/(hk(1, 1)*hk(2, 3) &
                                                                    -hk(1, 3)*hk(2, 1))
            efg(4) = (xd*alpham-xb*betam)/(xa*xd-xb*xc)
            efg(5) = (-xc*alpham+xa*betam)/(xa*xd-xb*xc)
            efg(6) = (hk(ih, i2)*fno(i1, 6)-hk(ih, i1)*fno(i2, 6)-fno(2, 6)* &
                      (hk(3, i1)*hk(ih, i2)-hk(3, i2)*hk(ih, i1)))/(hk(1, 1)*hk(2, 3) &
                                                                    -hk(1, 3)*hk(2, 1))
!
140         continue
!
            do i = 1, 6
                zr(jout-1+6*(ino-1)+i) = efg(i)
            end do
        end do
!
    else
!         -- UNE PROGRAMMATION STANDARD POUR MET3SEG4 :
        do ic = 1, 6
            do kp = 1, npg
                vpg(kp) = fno(kp, ic)
            end do
            call ppgan2(jgano, 1, 1, vpg, vno)
            do ino = 1, nno
                zr(jout+6*(ino-1)+ic-1) = vno(ino)
            end do
        end do
    end if
!
!
end subroutine
