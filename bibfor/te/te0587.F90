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
subroutine te0587(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/carcou.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!
!    - FONCTION REALISEE:  CALC_CHAMP POUR LES TUYAUX :
!        - EFGE_ELGA
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbcoum, nbsecm, jnbspi
    real(kind=8) :: h, a, r1
    parameter(nbsecm=32, nbcoum=10)
    real(kind=8) :: poicou(2*nbcoum+1), poisec(2*nbsecm+1)
    real(kind=8) :: pi, deuxpi, sig(6), fno(4, 6)
    real(kind=8) :: efg(6)
    real(kind=8) :: pgl(3, 3), pgl4(3, 3)
    real(kind=8) :: cosfi, sinfi
    real(kind=8) :: fi, poids, r, omega
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), rayon, theta, l
    integer(kind=8) :: nno, nnos, jgano, ndim, npg, nbcou, nbsec, lorien
    integer(kind=8) :: ipoids, ivf, ic, kp, jin, jcoopg, jdfd2
    integer(kind=8) :: idfdk
    integer(kind=8) :: igau, icou, isect, i, jout
    integer(kind=8) :: indice, icoud2, mmt
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
!  LES POIDS POUR L'INTEGRATION DANS L'EPAISSEUR
    poicou(1) = 1.d0/3.d0
    do i = 1, nbcou-1
        poicou(2*i) = 4.d0/3.d0
        poicou(2*i+1) = 2.d0/3.d0
    end do
    poicou(2*nbcou) = 4.d0/3.d0
    poicou(2*nbcou+1) = 1.d0/3.d0
!
!  LES POIDS POUR L'INTEGRATION SUR LA CIRCONFERENCE
    poisec(1) = 1.d0/3.d0
    do i = 1, nbsec-1
        poisec(2*i) = 4.d0/3.d0
        poisec(2*i+1) = 2.d0/3.d0
    end do
    poisec(2*nbsec) = 4.d0/3.d0
    poisec(2*nbsec+1) = 1.d0/3.d0
!
!
    if (option .eq. 'EFGE_ELGA') then
!     ---------------------------------
        call jevech('PCAORIE', 'L', lorien)
        call jevech('PSIEFR', 'L', jin)
        call jevech('PEFGER', 'E', jout)
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
            mmt = 0
        else
            mmt = 1
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
                    poids = poicou(icou)*poisec(isect)*h*deuxpi/ &
                            (4.d0*nbcou*nbsec)*r
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
        do ic = 1, 6
            do kp = 1, npg
                zr(jout+6*(kp-1)+ic-1) = fno(kp, ic)
            end do
        end do
    else
        call utmess('F', 'ELEMENTS4_49', sk=option)
    end if
!
end subroutine
