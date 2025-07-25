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
subroutine te0471(option, nomte)
!
    implicit none
#include "jeveux.h"
#include "asterfort/chmalg.h"
#include "asterfort/dpfch3.h"
#include "asterfort/dsfch3.h"
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
!                          OPTION : 'RIGI_MECA      '                  C
!                                                                      C
!    - ARGUMENTS:                                                      C
!        DONNEES:      OPTION       -->  OPTION DE CALCUL              C
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! .....................................................................C
!......................................................................C
    integer(kind=8) :: nbres, nddl
    parameter(nbres=4, nddl=7)
    character(len=24) :: carac, ff
    character(len=8) :: elrefe, fami, poum
    character(len=16) :: nomres(nbres)
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: valres(nbres), tpg, pgl(3, 3)
    integer(kind=8) :: nno1, nno2, npg1(2, 2), npg2(2, 2), npg, n, nbv
    integer(kind=8) :: kp, k1, k2, i, j, k, l, ik, ijkl, ij, ijl, lcorr
    integer(kind=8) :: imatuu, icarac, iff, imate, igeom, lorien, lsect, itype
    integer(kind=8) :: ipoids, ivf1, idpdx1, idpdy1, idpdz1, idsdx1, idsdy1, idsdz1
    integer(kind=8) :: idsxy1, idsxz1, idsyz1, idpdx2, idpdy2, idpdz2, idsdx2, idsdy2
    integer(kind=8) :: idsdz2, idsxy2, idsxz2, idsyz2, ivf2, ivf3, ipoi3, idpdx3, idpdy3
    integer(kind=8) :: idpdz3, ivf4, idpdx4, idpdy4, idpdz4, kpg, spt
    real(kind=8) :: raid(3), a(7, 7, 8, 8), poids2, poids
    real(kind=8) :: e, nu, xiy, xiz, rtor, rapp, xjx, ayz, xk(2), coord(60)
    real(kind=8) :: ycell, xlong
    real(kind=8) :: dsdyz(16), dfpdx2(8), dfpdy2(8), dfpdz2(8), dsdxx(16)
    real(kind=8) :: dsdyy(16), dsdzz(16), dsdxy(16), dsdxz(16), d2fdpl(8)
    real(kind=8) :: d2frot(8), b(1, 7, 9:20, 8), c(1, 1, 9:20, 9:20)
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
!
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
    call jevete(ff, 'L', iff)
    npg = npg1(1, 1)*npg1(1, 1)*npg1(1, 2)
    ipoids = iff
    ivf1 = ipoids+npg
    idpdx1 = ivf1+npg*2*nno1
    idpdy1 = idpdx1+npg*2*nno1
    idpdz1 = idpdy1+npg*2*nno1
    idsdx1 = idpdz1+npg*2*nno1
    idsdy1 = idsdx1+npg*2*nno1
    idsdz1 = idsdy1+npg*2*nno1
    idsxy1 = idsdz1+npg*2*nno1
    idsxz1 = idsxy1+npg*2*nno1
    idsyz1 = idsxz1+npg*2*nno1
    ivf2 = idsyz1+npg*2*nno1
    idpdx2 = ivf2+npg*nno1
    idpdy2 = idpdx2+npg*nno1
    idpdz2 = idpdy2+npg*nno1
    idsdx2 = idpdz2+npg*nno1
    idsdy2 = idsdx2+npg*nno1
    idsdz2 = idsdy2+npg*nno1
    idsxy2 = idsdz2+npg*nno1
    idsxz2 = idsxy2+npg*nno1
    idsyz2 = idsxz2+npg*nno1
!
    npg = npg1(1, 1)*npg1(1, 1)*npg1(1, 2)
    iff = iff+npg+10*(npg*2*nno1)+10*npg*nno1
    npg = npg1(2, 1)*npg1(2, 1)*npg1(2, 2)
    iff = iff+npg+npg*2*nno1+4*npg*nno1
!
    npg = npg2(1, 1)*npg2(1, 1)*npg2(1, 2)
    ipoi3 = iff
    ivf3 = ipoi3+npg
    idpdx3 = ivf3+npg*nno2
    idpdy3 = idpdx3+npg*nno2
    idpdz3 = idpdy3+npg*nno2
    ivf4 = idpdz3+npg*nno2
    idpdx4 = ivf4+npg*nno1
    idpdy4 = idpdx4+npg*nno1
    idpdz4 = idpdy4+npg*nno1
!     -------------------------------------------------
!     --- RECUPERATION LOI DE COMPORTEMENT MATERIAU ---
!     -------------------------------------------------
    call jevech('PMATERC', 'L', imate)
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
    if (phenom .eq. 'ELAS') then
        nomres(1) = 'E'
        nomres(2) = 'NU'
        nomres(3) = 'RHO'
        nbv = 3
    else
        call utmess('F', 'ELEMENTS3_98')
    end if
    tpg = 0.d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', phenom, 0, '   ', [tpg], &
                nbv, nomres, valres, icodre, 1)
    e = valres(1)
    nu = valres(2)
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
    xjx = zr(lsect+8)
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
!     ---------------------------------------
!     --- RECUPERATION MATRICE DE RAIDEUR ---
!     ---------------------------------------
    call jevech('PMATUUR', 'E', imatuu)
!     ---------------------------
!     --- CALCUL DES TENSEURS ---
!     ---------------------------
    raid(1) = e*xiy*rapp
    raid(2) = e*xiz*rapp
    raid(3) = e*ayz*rapp
    rtor = (e/(2.d0*(1.d0+nu)))*xjx*rapp
!     ---------------------------------------------
!     --- TERME CORRECTEUR : FONCTION D'HERMITE ---
!     --- POUR LES DDLS DE ROTATION             ---
!     ---------------------------------------------
    xlong = (coord(3*(5-1)+1)-coord(1))*0.5d0
!
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
!     ---------------------------------------------------
!     --- CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR ---
!     ---------------------------------------------------
!     --- FLEXION ---
!     ---------------
    npg = npg1(1, 2)*npg1(1, 1)*npg1(1, 1)
    do kp = 1, npg
!
        k1 = (kp-1)*2*nno1
        k2 = (kp-1)*nno1
!
!        --- CALCUL DES DERIVEES SECONDES
!
        call dsfch3(nno1, 2*nno1, zr(ipoids+kp-1), zr(idpdx1+k1), zr(idpdy1+k1), &
                    zr(idpdz1+k1), zr(idsdx1+k1), zr(idsdy1+k1), zr(idsdz1+k1), zr(idsxy1+k1), &
                    zr(idsxz1+k1), zr(idsyz1+k1), coord(1), zr(idpdx2+k2), zr(idpdy2+k2), &
                    zr(idpdz2+k2), zr(idsdx2+k2), zr(idsdy2+k2), zr(idsdz2+k2), zr(idsxy2+k2), &
                    zr(idsxz2+k2), zr(idsyz2+k2), dsdxx, dsdyy, dsdzz, &
                    dsdxy, dsdyz, dsdxz, poids)
!
        do j = 1, nno1
            d2fdpl(j) = dsdxx(j)
            d2frot(j) = dsdxx(j+nno1)*xlong
        end do
!
        xk(1) = poids
        xk(2) = -poids
!
        do i = 1, nno1
            do j = 1, i
                do k = 1, 2
                    a(k+1, k+1, i, j) = a(k+1, k+1, i, j)+raid(k)*poids*d2fdpl(i)*d2fdpl(j)
                    a(7-k, k+1, i, j) = a(7-k, k+1, i, j)+raid(k)*xk(k)*d2frot(i)*d2fdpl(j)
                    a(k+1, 7-k, i, j) = a(k+1, 7-k, i, j)+raid(k)*xk(k)*d2fdpl(i)*d2frot(j)
                    a(7-k, 7-k, i, j) = a(7-k, 7-k, i, j)+raid(k)*poids*d2frot(i)*d2frot(j)
                end do
            end do
        end do
    end do
!     ---------------------------------------------
!     ---- TRACTION ET COMPRESSION ET TORSION -----
!     ---------------------------------------------
    npg = npg2(1, 1)*npg2(1, 1)*npg2(1, 2)
    do kp = 1, npg
        k1 = (kp-1)*nno1
! --- CALCUL DES FONCTIONS DE FORME ET DE LEURS DERIVEES
        call dpfch3(nno1, nno1, zr(ipoi3+kp-1), zr(idpdx4+k1), zr(idpdy4+k1), &
                    zr(idpdz4+k1), coord(1), zr(idpdx4+k1), zr(idpdy4+k1), zr(idpdz4+k1), &
                    dfpdx2, dfpdy2, dfpdz2, poids2)
!
        do i = 1, nno1
            do j = 1, i
                a(1, 1, i, j) = a(1, 1, i, j)+poids2*(raid(3)*dfpdx2(i)*dfpdx2(j))
!
                a(4, 4, i, j) = a(4, 4, i, j)+poids2*(rtor*dfpdx2(i)*dfpdx2(j))
            end do
        end do
    end do
!     -------------------------------------------------
!     --- PASSAGE DU REPERE LOCAL AU REPERE GLOBAL  ---
!     -------------------------------------------------
    do i = 1, nno1
        do j = 1, i
            call chmalg(a(1, 1, i, j), pgl, nddl, nddl)
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
!----------------------------------------------------------------------
end subroutine
