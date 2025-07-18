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
subroutine te0470(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/chmalg.h"
#include "asterfort/dpfch3.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matrot.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
!
    character(len=16) :: option, nomte, phenom
! .....................................................................C
! .....................................................................C
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES            C
!                          ELEMENTS 3D COEUR HOMOGENEISE               C
!                          OPTION : 'MASS_MECA      '                  C
!                                                                      C
!    - ARGUMENTS:                                                      C
!        DONNEES:      OPTION       -->  OPTION DE CALCUL              C
!                      NOMTE        -->  NOM DU TYPE ELEMENT           C
! .....................................................................C
! .....................................................................C
    integer(kind=8) :: nbres, nddl
    parameter(nbres=3, nddl=7)
    character(len=24) :: carac, ff
    character(len=8) :: elrefe, fami, poum
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: icodre(nbres), kpg, spt
    real(kind=8) :: valres(nbres), tpg, pgl(3, 3)
    integer(kind=8) :: nno1, nno2, npg1(2, 2), npg2(2, 2), npg, n, nbv
    integer(kind=8) :: imatuu, icarac, iff, imate, igeom, lorien, lsect, itype
    integer(kind=8) :: kp, k1, k2, k3, i, j, ik, ijkl, k, l, ij, ijl, lcorr
    integer(kind=8) :: ipoi2, ivf2p, ivf2g, iddx2g, iddy2g, iddz2g, ipoi3, ivf3f, iddx3f
    integer(kind=8) :: iddy3f, iddz3f, ivf3g, iddx3g, iddy3g, iddz3g, ipoi4, ivf4p
    integer(kind=8) :: ivf4f, iddx4f, iddy4f, iddz4f, ivf4g, iddx4g, iddy4g, iddz4g
    real(kind=8) :: mass(3, 3), rdp(3, 3), d(2, 2), mtor
    real(kind=8) :: dpdx(8), dpdy(8), dpdz(8), poids2, poids3, poids4, ffdpl2(8)
    real(kind=8) :: ffrot2(8), ffdpl4(8), ffrot4(8), xk(2), dfpdx3(20)
    real(kind=8) :: dfpdy3(20), dfpdz3(20), dfpdx4(20), dfpdy4(20), dfpdz4(20)
    real(kind=8) :: a(7, 7, 8, 8), b(1, 7, 9:20, 8), c(1, 1, 9:20, 9:20)
    real(kind=8) :: rhopou, rhoflu, ys, yf, xiy, xiz, rapp, ayz, b11, b12, b22
    real(kind=8) :: coord(60), ycell, xlong
!
!     ----------------------------------------
!     --- RECUPERATION FONCTIONS DE FORMES ---
!     ----------------------------------------
    if (nomte .eq. 'MECA_POHO_HEXA8') then
        elrefe = 'POHOH8'
    else
        elrefe = 'POHOH20'
    end if
    carac = '&INEL.'//elrefe//'.CARAC'
    ff = '&INEL.'//elrefe//'.FF'
! --- FAMILLES DES POINTS DE GAUSS (VOIR INI100)
    call jevete(carac, 'L', icarac)
    nno1 = zi(icarac)
    nno2 = zi(icarac+1)
    n = 1
    do i = 1, 2
        do j = 1, 2
            n = n+1
            npg1(i, j) = zi(icarac+n)
            npg2(i, j) = zi(icarac+n+4)
        end do
    end do
! --- ADRESSES DES FONCTIONS DE FORMES
    call jevete(ff, 'L', iff)
    npg = npg1(1, 1)*npg1(1, 1)*npg1(1, 2)
    iff = iff+npg+10*(npg*2*nno1)+10*npg*nno1
!
    npg = npg1(2, 1)*npg1(2, 1)*npg1(2, 2)
    ipoi2 = iff
    ivf2p = ipoi2+npg
    ivf2g = ivf2p+npg*2*nno1
    iddx2g = ivf2g+npg*nno1
    iddy2g = iddx2g+npg*nno1
    iddz2g = iddy2g+npg*nno1
!
    ipoi3 = iddz2g+npg*nno1
    npg = npg2(1, 1)*npg2(1, 1)*npg2(1, 2)
    ivf3f = ipoi3+npg
    iddx3f = ivf3f+npg*nno2
    iddy3f = iddx3f+npg*nno2
    iddz3f = iddy3f+npg*nno2
    ivf3g = iddz3f+npg*nno2
    iddx3g = ivf3g+npg*nno1
    iddy3g = iddx3g+npg*nno1
    iddz3g = iddy3g+npg*nno1
!
    ipoi4 = iddz3g+npg*nno1
    npg = npg2(2, 1)*npg2(2, 1)*npg2(2, 2)
    ivf4p = ipoi4+npg
    ivf4f = ivf4p+npg*2*nno1
    iddx4f = ivf4f+npg*nno2
    iddy4f = iddx4f+npg*nno2
    iddz4f = iddy4f+npg*nno2
    ivf4g = iddz4f+npg*nno2
    iddx4g = ivf4g+npg*nno1
    iddy4g = iddx4g+npg*nno1
    iddz4g = iddy4g+npg*nno1
!     ---------------------------------------------------
!     ---- RECUPERATION LOI DE COMPORTEMENT MATERIAU ----
!     ---------------------------------------------------
    call jevech('PMATERC', 'L', imate)
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
    if (phenom .eq. 'ELAS') then
        nomres(1) = 'RHO'
        nbv = 1
    else
        call utmess('F', 'ELEMENTS3_98')
    end if
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    tpg = 0.d0
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', phenom, 0, '   ', [tpg], &
                nbv, nomres, valres, icodre, 1)
    rhopou = valres(1)
    call rccoma(zi(imate), 'FLUIDE', 1, phenom, icodre(1))
    if (phenom .eq. 'FLUIDE') then
        nomres(1) = 'RHO'
        nbv = 1
    else
        call utmess('F', 'ELEMENTS3_98')
    end if
    tpg = 0.d0
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', phenom, 0, '   ', [tpg], &
                nbv, nomres, valres, icodre, 1)
    rhoflu = valres(1)
!     ----------------------------------------------------------------
!     --- RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS ---
!     ----------------------------------------------------------------
    call jevech('PCAGNPO', 'L', lsect)
    lsect = lsect-1
    itype = nint(zr(lsect+23))
!     --- SECTION INITIALE ---
    ayz = zr(lsect+1)
    xiy = zr(lsect+2)
    xiz = zr(lsect+3)
!
    if (itype .gt. 0) then
        call utmess('F', 'ELEMENTS3_99')
    end if
!     -------------------------------------------
!     --- RECUPERATION DES TERMES CORRECTEURS ---
!     -------------------------------------------
    call jevech('PCAPOUF', 'L', lcorr)
    lcorr = lcorr-1
    ycell = zr(lcorr+5)
    b11 = zr(lcorr+1)/ycell
    b22 = zr(lcorr+2)/ycell
    b12 = zr(lcorr+3)/ycell
    yf = zr(lcorr+4)/ycell
    ys = 1-yf
    rapp = zr(lcorr+6)
    rapp = rapp*rapp/ycell
!     --------------------------------------------------
!     --- RECUPERATION CARACTERISTIQUES GEOMETRIQUES ---
!     --------------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
! --- RECUPERATION DES ORIENTATIONS
    call jevech('PCAORIE', 'L', lorien)
! --- CALCUL DE LA MATRICE DE PASSAGE
    call matrot(zr(lorien), pgl)
! --- CHANGEMENT DE REPERE : GLOBAL --> LOCAL
    call utpvgl(nno2, 3, pgl, zr(igeom), coord)
!     -------------------------------------
!     --- RECUPERATION MATRICE DE MASSE ---
!     -------------------------------------
    call jevech('PMATUUR', 'E', imatuu)
!     ---------------------------
!     --- CALCUL DES TENSEURS ---
!     ---------------------------
    mass(1, 1) = rhoflu*b11+rhopou*ayz*rapp
    mass(1, 2) = rhoflu*b12
    mass(2, 1) = rhoflu*b12
    mass(2, 2) = rhoflu*b22+rhopou*ayz*rapp
    mass(3, 3) = rhopou*ayz*rapp
    mtor = rhopou*(xiy+xiz)*rapp
!
    rdp(1, 1) = (yf-b11)*rhoflu
    rdp(1, 2) = -b12*rhoflu
    rdp(2, 1) = -b12*rhoflu
    rdp(2, 2) = (yf-b22)*rhoflu
    rdp(3, 3) = yf*rhoflu
!
    d(1, 1) = (b11+ys)*rhoflu
    d(1, 2) = b12*rhoflu
    d(2, 1) = b12*rhoflu
    d(2, 2) = (b22+ys)*rhoflu
!     ------------------------------------------
!     --- INITIALISATION A ZERO DE A, B ET C ---
!     ------------------------------------------
    do j = 1, nno1
        do i = 1, nno1
            do l = 1, nddl
                do k = 1, nddl
                    a(k, l, i, j) = 0.d0
                end do
            end do
        end do
    end do
! --- LES MATRICES B ET C NE SONT UTILISEES QUE POUR DES MAILLES HEXA20
    do j = 1, nno1
        do i = nno1+1, nno2
            do l = 1, nddl
                b(1, l, i, j) = 0.d0
            end do
        end do
    end do
    do j = nno1+1, nno2
        do i = nno1+1, 20
            c(1, 1, i, j) = 0.d0
        end do
    end do
!     ---------------------------------------------
!     --- CHANGEMENT DE SIGNES POUR LE PASSAGE  ---
!     --- DDL ROTATION A DERIVEE DEPLACEMENT    ---
!     ---------------------------------------------
    xk(1) = 1.d0
    xk(2) = -1.d0
!     ---------------------------------------------
!     --- TERME CORRECTEUR : FONCTION D'HERMITE ---
!     --- POUR LES DDLS DE ROTATION             ---
!     ---------------------------------------------
    xlong = (coord(3*(5-1)+1)-coord(1))*0.5d0
!
!     --- CALCUL DE LA MATRICE ELEMENTAIRE DE MASSE ---
!     -------------------------------------------------
!     --- MATRICE DE MASSE : POUTRE-FLEXION ---
!     -----------------------------------------
    npg = npg1(2, 1)*npg1(2, 1)*npg1(2, 2)
    do kp = 1, npg
!
        k1 = (kp-1)*2*nno1
        k2 = (kp-1)*nno1
! --- CALCUL DES FONCTIONS DE FORME ET DE LEURS DERIVEES
        call dpfch3(nno1, nno1, zr(ipoi2+kp-1), zr(iddx2g+k2), zr(iddy2g+k2), &
                    zr(iddz2g+k2), coord(1), zr(iddx2g+k2), zr(iddy2g+k2), zr(iddz2g+k2), &
                    dpdx, dpdy, dpdz, poids2)
!
        do i = 1, nno1
            ffdpl2(i) = zr(ivf2p+k1+i-1)
            ffrot2(i) = zr(ivf2p+k1+i+nno1-1)*xlong
        end do
! --- REMPLISSAGE MATRICE ELEMENTAIRE
        do i = 1, nno1
            do j = 1, i
                do k = 1, 2
                    do l = 1, 2
                        a(k+1, l+1, i, j) = a(k+1, l+1, i, j)+mass(k, l)*poids2*ffdpl2(i)*ffdpl2&
                                         &(j)
                        a(7-k, l+1, i, j) = a(7-k, l+1, i, j)+xk(k)*mass(k, l)*poids2*ffrot2(i) &
                                           &*ffdpl2(j)
                        a(k+1, 7-l, i, j) = a(k+1, 7-l, i, j)+xk(l)*mass(k, l)*poids2*ffdpl2(i) &
                                           &*ffrot2(j)
                        a(7-k, 7-l, i, j) = a(7-k, 7-l, i, j)+xk(k)*xk(l)*mass(k, l)*poids2*ffrot&
                                         &2(i)*ffrot2(j)
                    end do
                end do
            end do
        end do
    end do
!      ---------------------------------------------------
!      --- MATRICE DE MASSE : COUPLAGE POUTRE / FLUIDE ---
!      ---------------------------------------------------
    npg = npg2(2, 2)*npg2(2, 1)*npg2(2, 1)
    do kp = 1, npg
!
        k1 = (kp-1)*2*nno1
        k2 = (kp-1)*nno2
        k3 = (kp-1)*nno1
! --- CALCUL DES FONCTIONS DE FORME ET DE LEURS DERIVEES
        do i = 1, nno1
            ffdpl4(i) = zr(ivf4p+k1+i-1)
            ffrot4(i) = zr(ivf4p+k1+i+nno1-1)*xlong
        end do
!
        call dpfch3(nno1, nno2, zr(ipoi4+kp-1), zr(iddx4f+k2), zr(iddy4f+k2), &
                    zr(iddz4f+k2), coord(1), zr(iddx4g+k3), zr(iddy4g+k3), zr(iddz4g+k3), &
                    dfpdx4, dfpdy4, dfpdz4, poids4)
! --- REMPLISSAGE DES MATRICES ELEMENTAIRES
        do i = 1, nno1
            do j = 1, i
                do k = 1, 2
                    a(7, k+1, i, j) = a(7, k+1, i, j)-poids4*ffdpl4(j)*(d(k, 1)*dfpdy4(i)+d(k,&
                                   &2)*dfpdz4(i))
                    a(7, 7-k, i, j) = a(7, 7-k, i, j)-poids4*xk(k)*ffrot4(j)*(d(k, 1)*dfpdy4(i&
                                   &)+d(k, 2)*dfpdz4(i))
                    a(k+1, 7, i, j) = a(k+1, 7, i, j)-poids4*ffdpl4(i)*(d(k, 1)*dfpdy4(j)+d(k,&
                                   &2)*dfpdz4(j))
                    a(7-k, 7, i, j) = a(7-k, 7, i, j)-poids4*xk(k)*ffrot4(i)*(d(k, 1)*dfpdy4(j&
                                   &)+d(k, 2)*dfpdz4(j))
                end do
            end do
        end do
!
        do i = nno1+1, nno2
            do j = 1, nno1
                do k = 1, 2
                    b(1, k+1, i, j) = b(1, k+1, i, j)-poids4*ffdpl4(j)*(d(k, 1)*dfpdy4(i)+d(k,&
                                   &2)*dfpdz4(i))
                    b(1, 7-k, i, j) = b(1, 7-k, i, j)-poids4*xk(k)*ffrot4(j)*(d(k, 1)*dfpdy4(i&
                                   &)+d(k, 2)*dfpdz4(i))
                end do
            end do
        end do
!
    end do
!      --------------------------------------------------------------
!      --- MATRICE DE MASSE : FLUIDE ET POUTRE-TRACTION ET TORSON ---
!      --------------------------------------------------------------
    npg = npg2(1, 1)*npg2(1, 1)*npg2(1, 2)
    do kp = 1, npg
!
        k1 = (kp-1)*nno1
        k2 = (kp-1)*nno2
! --- CALCUL DES FONCTIONS DE FORME ET DE LEURS DERIVEES
        call dpfch3(nno1, nno2, zr(ipoi3+kp-1), zr(iddx3f+k2), zr(iddy3f+k2), &
                    zr(iddz3f+k2), coord(1), zr(iddx3g+k1), zr(iddy3g+k1), zr(iddz3g+k1), &
                    dfpdx3, dfpdy3, dfpdz3, poids3)
!
        do i = 1, nno1
            ffdpl2(i) = zr(ivf3g+k1+i-1)
        end do
!
        do i = 1, nno1
            do j = 1, i
                a(1, 1, i, j) = a(1, 1, i, j)+poids3*(mass(3, 3)*ffdpl2(i)*ffdpl2(j))
!
                a(4, 4, i, j) = a(4, 4, i, j)+poids3*(mtor*ffdpl2(i)*ffdpl2(j))
!
                a(7, 7, i, j) = a(7, 7, i, j)-poids3*(rdp(1, 1)*dfpdy3(i)*dfpdy3(j)+rdp(1, 2) &
                             &*dfpdz3(i)*dfpdy3(j)+rdp(2, 1)*dfpdy3(i)*dfpdz3(j)+rdp(2&
                             &, 2)*dfpdz3(i)*dfpdz3(j)+rdp(3, 3)*dfpdx3(i)*dfpdx3(j))
            end do
        end do
        do i = nno1+1, nno2
            do j = 1, nno1
                b(1, 7, i, j) = b(1, 7, i, j)-poids3*(rdp(1, 1)*dfpdy3(i)*dfpdy3(j)+rdp(1, 2) &
                             &*dfpdz3(i)*dfpdy3(j)+rdp(2, 1)*dfpdy3(i)*dfpdz3(j)+rdp(2&
                             &, 2)*dfpdz3(i)*dfpdz3(j)+rdp(3, 3)*dfpdx3(i)*dfpdx3(j))
            end do
        end do
        do i = nno1+1, nno2
            do j = nno1+1, i
                c(1, 1, i, j) = c(1, 1, i, j)-poids3*(rdp(1, 1)*dfpdy3(i)*dfpdy3(j)+rdp(1, 2) &
                             &*dfpdz3(i)*dfpdy3(j)+rdp(2, 1)*dfpdy3(i)*dfpdz3(j)+rdp(2&
                             &, 2)*dfpdz3(i)*dfpdz3(j)+rdp(3, 3)*dfpdx3(i)*dfpdx3(j))
            end do
        end do
!
    end do
!     -------------------------------------------------
!     --- PASSAGE DU REPERE LOCAL AU REPERE GLOBAL  ---
!     -------------------------------------------------
    do i = 1, nno1
        do j = 1, i
            call chmalg(a(1, 1, i, j), pgl, nddl, nddl)
        end do
    end do
    do i = nno1+1, nno2
        do j = 1, nno1
            call chmalg(b(1, 1, i, j), pgl, 1, nddl)
        end do
    end do
! ---------------------------------------------------------------------
! --- PASSAGE DE LA MATRICE RECTANGULAIRE A LA MATRICE TRIANGULAIRE ---
! ---------------------------------------------------------------------
    do k = 1, nddl
        do l = 1, nddl
!   IL Y A ECRASEMENT SI ON INTERVERTIE L'ORDRE DES BOUCLES 400 ET 410
            do i = 1, nno1
                ik = ((i-1)*nddl+k-1)*((i-1)*nddl+k)/2
                do j = 1, i
                    ijkl = ik+(j-1)*nddl+l
                    zr(imatuu+ijkl-1) = a(k, l, i, j)
                end do
            end do
        end do
    end do
!   BOUCLE EXECUTEE QUE POUR DES MAILLES HEXA20
    imatuu = imatuu+(nddl*nno1)*(nddl*nno1+1)/2
    do i = nno1+1, nno2
        ij = (i-nno1-1)*nddl*nno1+(i-nno1-1)*(i-nno1)/2
        do j = 1, nno1
            ijl = ij+(j-1)*nddl
            do l = 1, nddl
                zr(imatuu+ijl+(l-1)) = b(1, l, i, j)
            end do
        end do
        ijl = ij+nno1*nddl
        do j = nno1+1, i
            zr(imatuu+ijl+(j-nno1-1)) = c(1, 1, i, j)
        end do
    end do
!
end subroutine
