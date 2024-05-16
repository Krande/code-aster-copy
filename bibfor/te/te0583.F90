! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
subroutine te0583(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/carcou.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vlggl.h"
#include "asterfort/vlgglc.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DU SECOND MEMBRE : TRAVAIL DE LA
!                          PRESSION ET FORCES LINEIQUES TUYAUX
!     OPTIONS :'CHAR_MECA_PESA_R' 'CHAR_MECA_FR1D1D''CHAR_MECA_PRES_R'
!              'CHAR_MECA_PRES_F'
! ......................................................................
    integer :: nbrddm
    parameter(nbrddm=156)
    real(kind=8) :: h, a, l, presno(4), prespg(4), rint, r1
    real(kind=8) :: vpesan(6), fpesan(6), pesan, f(nbrddm)
    real(kind=8) :: pi, deuxpi, fi, pass(nbrddm, nbrddm)
    real(kind=8) :: fpesa1(6), fpesa2(6), fpesa3(6), vtemp(nbrddm)
    real(kind=8) :: vpesa1(6), vpesa2(6), vpesa3(6), vpesa4(6)
    real(kind=8) :: pgl(3, 3), pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), omega
    real(kind=8) :: hk, poids, rayon, theta, tk(4), ck, sk
    real(kind=8) :: cosfi, sinfi, te, pgl4(3, 3), fpesa4(6), xpg(4)
    real(kind=8) :: r8b, rext, sec, rho(1), r, time, valpar(5)
    integer :: codres(1), kpg, spt
    character(len=8) :: nompar(5), fami, poum
    character(len=16) :: phenom
    integer :: nbcou, nbsec, m, lorien, icoude
    integer :: ipoids, ivf, i, icou, ibloc, ino, nbpar, icompx, niter, iter
    integer :: igeom, lmater, jpesa, jout, lforc, iret
    integer :: igau, isect, ipres, k, ivect, nbrddl, indic0
    integer :: indic1, indic2, indic3, indic4, indic5, j
    integer :: jnbspi, nbsecm, nbcoum, itemps, ier, labsc, itab(2)
    integer :: ndim, nnos, nno, jcoopg, idfdk, jdfd2, jgano, npg
    parameter(nbsecm=32, nbcoum=10)
    real(kind=8) :: poicou(2*nbcoum+1), poisec(2*nbsecm+1), abscn(4)
    aster_logical :: normal, global
!-----------------------------------------------------------------------
    integer, parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8), parameter :: noms_cara1(nb_cara1) = ['R1 ', 'EP1']
!----------------------------------------------------------------------
    r8b = 0.d0
    call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jdfd2=jdfd2, &
                     jgano=jgano)
!
!   -- si l'abscisse curviligne est fournie, on la stocke dans abscn(:)
    labsc = 0
    if (option .eq. 'CHAR_MECA_PRES_F') then
        call tecach('ONO', 'PABSCUR', 'L', iret, iad=labsc)
        if (labsc .ne. 0) then
            ASSERT(iret .eq. 0)
            call tecach('OOO', 'PABSCUR', 'L', iret, nval=2, &
                        itab=itab)
            ASSERT(itab(1) .eq. labsc)
            ASSERT(itab(2) .eq. nno)
            do k = 1, nno
                abscn(k) = zr(labsc-1+k)
            end do
        end if
    end if
!
    if (option .eq. 'CHAR_MECA_FC1D1D') then
        icompx = 1
        niter = 2
    else
        icompx = 0
        niter = 1
    end if
    pi = r8pi()
    deuxpi = 2.d0*pi
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    nbsec = zi(jnbspi-1+2)
!
!     -- CALCUL DES POIDS DES COUCHES ET DES SECTEURS:
    poicou(1) = 1.d0/3.d0
    do i = 1, nbcou-1
        poicou(2*i) = 4.d0/3.d0
        poicou(2*i+1) = 2.d0/3.d0
    end do
    poicou(2*nbcou) = 4.d0/3.d0
    poicou(2*nbcou+1) = 1.d0/3.d0
    poisec(1) = 1.d0/3.d0
    do i = 1, nbsec-1
        poisec(2*i) = 4.d0/3.d0
        poisec(2*i+1) = 2.d0/3.d0
    end do
    poisec(2*nbsec) = 4.d0/3.d0
    poisec(2*nbsec+1) = 1.d0/3.d0
    m = 3
    if (nomte .eq. 'MET6SEG3') m = 6
!
!
    do i = 1, npg
        xpg(i) = zr(jcoopg-1+i)
    end do
    nbrddl = nno*(6+3+6*(m-1))
    if (nbrddl .gt. nbrddm) then
        call utmess('F', 'ELEMENTS4_40')
    end if
    if (nomte .eq. 'MET3SEG3') then
        if (nbrddl .ne. 63) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else if (nomte .eq. 'MET6SEG3') then
        if (nbrddl .ne. 117) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else if (nomte .eq. 'MET3SEG4') then
        if (nbrddl .ne. 84) then
            call utmess('F', 'ELEMENTS4_41')
        end if
    else
        call utmess('F', 'ELEMENTS4_42')
    end if
    call jevech('PCAORIE', 'L', lorien)
    call carcou(zr(lorien), l, pgl, rayon, theta, &
                pgl1, pgl2, pgl3, pgl4, nno, &
                omega, icoude)
    if (icoude .ge. 10) then
        icoude = icoude-10
    end if
    call jevech('PGEOMER', 'L', igeom)
    call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
    r1 = vale_cara1(1)
    h = vale_cara1(2)
    a = r1-h/2.d0
    rint = a-h/2.d0
    rext = a+h/2.d0
    sec = pi*(rext**2-rint**2)
! A= RMOY, H = EPAISSEUR RINT = RAYON INTERIEUR
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
    if (option(1:14) .eq. 'CHAR_MECA_PRES') then
!   ---------------------------------------------
!       -- LA PRESSION NE TRAVAILLE QUE SUR LE TERME WO
        if (option(15:16) .eq. '_R') then
            call jevecd('PPRESSR', ipres, 0.d0)
            do i = 1, nno
                presno(i) = zr(ipres-1+i)
            end do
        else if (option(15:16) .eq. '_F') then
            call jevech('PPRESSF', 'L', ipres)
            call jevech('PINSTR', 'L', itemps)
            valpar(4) = zr(itemps)
            nompar(4) = 'INST'
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            if (labsc .ne. 0) then
                nompar(5) = 'ABSC'
                nbpar = 5
            else
                nbpar = 4
            end if
!
            do i = 1, nno
                valpar(1) = zr(igeom+3*(i-1))
                valpar(2) = zr(igeom+3*(i-1)+1)
                valpar(3) = zr(igeom+3*(i-1)+2)
                if (labsc .ne. 0) valpar(5) = abscn(i)
                call fointe('FM', zk8(ipres), nbpar, nompar, valpar, &
                            presno(i), ier)
            end do
        end if
!
        do igau = 1, npg
            prespg(igau) = 0.d0
            do k = 1, nno
                hk = zr(ivf-1+nno*(igau-1)+k)
                prespg(igau) = hk*presno(k)+prespg(igau)
            end do
        end do
!
        call jevech('PVECTUR', 'E', ivect)
        do k = 1, nno
!           TRAVAIL SUR UX
            indic0 = ivect-1+(6+6*(m-1)+3)*(k-1)+1
!           TRAVAIL SUR UY
            indic1 = ivect-1+(6+6*(m-1)+3)*(k-1)+2
!           TRAVAIL SUR UZ
            indic2 = ivect-1+(6+6*(m-1)+3)*(k-1)+3
!           TRAVAIL SUR W0
            indic3 = ivect-1+(6+6*(m-1)+3)*(k-1)+1+6+6*(m-1)
!           TRAVAIL SUR WI1
            indic4 = ivect-1+(6+6*(m-1)+3)*(k-1)+2+6+6*(m-1)
!           TRAVAIL SUR W01
            indic5 = ivect-1+(6+6*(m-1)+3)*(k-1)+3+6+6*(m-1)
            do igau = 1, npg
!               -- boucle sur les points de simpson dans l'epaisseur
                hk = zr(ivf-1+nno*(igau-1)+k)
                if (icoude .eq. 1) then
                    ck = cos((1.d0+xpg(igau))*theta/2.d0-tk(k))
                    sk = sin((1.d0+xpg(igau))*theta/2.d0-tk(k))
                else
                    ck = 1.d0
                    sk = 0.d0
                end if
!               -- boucle sur les points de simpson sur la circonference
                do isect = 1, 2*nbsec+1
                    if (icoude .eq. 0) then
                        poids = zr(ipoids-1+igau)*poisec(isect)*(l/2.d0)*deuxpi/(2.d0*nbsec)*&
                                &rint
                        zr(indic3) = zr(indic3)+hk*poids*prespg(igau)
                    else
                        fi = (isect-1)*deuxpi/(2.d0*nbsec)
                        cosfi = cos(fi)
                        sinfi = sin(fi)
                        te = fi-omega
                        l = theta*(rayon+rint*sinfi)
                        poids = zr(ipoids-1+igau)*poisec(isect)*(l/2.d0)*deuxpi/(2.d0*nbsec)*&
                                &rint
                        zr(indic0) = zr(indic0)+hk*poids*prespg(igau)*sinfi*sk
                        zr(indic1) = zr(indic1)-hk*poids*prespg(igau)*sinfi*ck
                        zr(indic2) = zr(indic2)-hk*poids*prespg(igau)*cosfi
                        zr(indic3) = zr(indic3)+hk*poids*prespg(igau)
                        zr(indic4) = zr(indic4)+hk*poids*prespg(igau)*cos(te)
                        zr(indic5) = zr(indic5)+hk*poids*prespg(igau)*sin(te)
                    end if
                end do
            end do
        end do
        if (icoude .ne. 0) then
            call vlgglc(nno, nbrddl, pgl1, pgl2, pgl3, &
                        pgl4, zr(ivect), 'LG', pass, vtemp)
        end if
!
!
    else if ((option .eq. 'CHAR_MECA_PESA_R') .or. ( &
             option .eq. 'CHAR_MECA_FR1D1D') .or. (option .eq. 'CHAR_MECA_FC1D1D')) &
        !       -----------------------------------------------------------------------------
        !       CAS PESANTEUR ET FORCE LINEIQUE
        then
        do iter = 1, niter
            if (option .eq. 'CHAR_MECA_PESA_R') then
                call jevech('PMATERC', 'L', lmater)
                call rccoma(zi(lmater), 'ELAS', 1, phenom, codres(1))
                if (phenom .eq. 'ELAS' .or. phenom .eq. 'ELAS_ISTR' .or. phenom .eq. &
                    'ELAS_ORTH') then
                    fami = 'FPG1'
                    kpg = 1
                    spt = 1
                    poum = '+'
                    call rcvalb(fami, kpg, spt, poum, zi(lmater), &
                                ' ', phenom, 0, ' ', [r8b], &
                                1, 'RHO', rho, codres, 1)
                else
                    call utmess('F', 'ELEMENTS4_43')
                end if
                call jevech('PPESANR', 'L', jpesa)
                pesan = zr(jpesa)
                do i = 1, 3
                    vpesan(i) = rho(1)*pesan*zr(jpesa+i)
                end do
                do i = 4, 6
                    vpesan(i) = 0.d0
                end do
            else
                if (icompx .eq. 0) then
                    call jevech('PFR1D1D', 'L', lforc)
                else
                    call jevech('PFC1D1D', 'L', lforc)
                end if
                if (icompx .eq. 1) then
                    if (iter .eq. 1) then
                        do i = 1, 3
                            vpesan(i) = dble(zc(lforc-1+i))/sec
                        end do
                    else
                        do i = 1, 3
                            vpesan(i) = dimag(zc(lforc-1+i))/sec
                        end do
                    end if
                else
                    do i = 1, 3
                        vpesan(i) = zr(lforc-1+i)/sec
                    end do
                end if
                do i = 4, 6
                    vpesan(i) = 0.d0
                end do
            end if
            do i = 1, nbrddl
                f(i) = 0.d0
            end do
            if (icoude .eq. 0) then
                call utpvgl(1, 6, pgl, vpesan(1), fpesan(1))
            else
                call utpvgl(1, 6, pgl1, vpesan(1), fpesa1(1))
                call utpvgl(1, 6, pgl2, vpesan(1), fpesa2(1))
                call utpvgl(1, 6, pgl3, vpesan(1), fpesa3(1))
                if (nno .eq. 4) then
                    call utpvgl(1, 6, pgl4, vpesan(1), fpesa4(1))
                end if
            end if
!           BOUCLE SUR LES POINTS DE GAUSS DANS LA LONGUEUR
            do igau = 1, npg
!               BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
                do icou = 1, 2*nbcou+1
                    r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
!                   BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
                    do isect = 1, 2*nbsec+1
                        if (icoude .eq. 0) then
                            poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deu&
                                    &xpi/(4.d0*nbcou*nbsec)*r
                            do k = 1, nno
                                hk = zr(ivf-1+nno*(igau-1)+k)
                                ibloc = (9+6*(m-1))*(k-1)
                                f(ibloc+1) = f(ibloc+1)+poids*hk*fpesan(1)
                                f(ibloc+2) = f(ibloc+2)+poids*hk*fpesan(2)
                                f(ibloc+3) = f(ibloc+3)+poids*hk*fpesan(3)
                            end do
                        else if (icoude .eq. 1) then
                            fi = (isect-1)*deuxpi/(2.d0*nbsec)
                            cosfi = cos(fi)
                            sinfi = sin(fi)
                            l = theta*(rayon+r*sinfi)
                            poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deu&
                                    &xpi/(4.d0*nbcou*nbsec)*r
                            do k = 1, 3
                                hk = zr(ivf-1+nno*(igau-1)+1)
                                ibloc = (9+6*(m-1))*(1-1)
                                f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa1(k)
                                hk = zr(ivf-1+nno*(igau-1)+2)
                                ibloc = (9+6*(m-1))*(2-1)
                                f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa2(k)
                                hk = zr(ivf-1+nno*(igau-1)+3)
                                ibloc = (9+6*(m-1))*(3-1)
                                f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa3(k)
                                if (nno .eq. 4) then
                                    ibloc = (9+6*(m-1))*(4-1)
                                    f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa4(k)
                                end if
                            end do
                        end if
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
            if (icompx .eq. 1) then
                call jevech('PVECTUC', 'E', jout)
                if (iter .eq. 1) then
                    do j = 1, nbrddl
                        zc(jout-1+j) = f(j)
                    end do
                else
                    do j = 1, nbrddl
                        zc(jout-1+j) = dcmplx(dble(zc(jout-1+j)), dble(f(j)))
                    end do
                end if
            else
                call jevech('PVECTUR', 'E', jout)
                do i = 1, nbrddl
                    zr(jout-1+i) = f(i)
                end do
            end if
        end do
!
!
    else if ((option .eq. 'CHAR_MECA_FF1D1D')) then
!   -----------------------------------------------
!       -- CAS FORCE LINEIQUE FONCTION
!
        call jevech('PFF1D1D', 'L', lforc)
        normal = zk8(lforc+6) .eq. 'VENT'
        global = zk8(lforc+6) .eq. 'GLOBAL'
        if (normal) then
            call utmess('F', 'ELEMENTS4_44', sk=option)
        end if
        if (.not. global) then
            call utmess('F', 'ELEMENTS4_45', sk=option)
        end if
        nompar(4) = 'INST'
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
        if (itemps .ne. 0) then
            time = zr(itemps)
            valpar(4) = time
            nbpar = 4
        else
            nbpar = 3
        end if
!        NOEUDS 1 A 3
        ino = 1
        do i = 1, 3
            valpar(i) = zr(igeom-1+3*(ino-1)+i)
        end do
        do i = 1, 3
            call fointe('FM', zk8(lforc+i-1), nbpar, nompar, valpar, &
                        vpesa1(i), ier)
            vpesa1(i) = vpesa1(i)/sec
        end do
        ino = 2
        do i = 1, 3
            valpar(i) = zr(igeom-1+3*(ino-1)+i)
        end do
        do i = 1, 3
            call fointe('FM', zk8(lforc+i-1), nbpar, nompar, valpar, &
                        vpesa2(i), ier)
            vpesa2(i) = vpesa2(i)/sec
        end do
        ino = 3
        do i = 1, 3
            valpar(i) = zr(igeom-1+3*(ino-1)+i)
        end do
        do i = 1, 3
            call fointe('FM', zk8(lforc+i-1), nbpar, nompar, valpar, &
                        vpesa3(i), ier)
            vpesa3(i) = vpesa3(i)/sec
        end do
        if (nno .eq. 4) then
            ino = 4
            do i = 1, 3
                valpar(i) = zr(igeom-1+3*(ino-1)+i)
            end do
            do i = 1, 3
                call fointe('FM', zk8(lforc+i-1), nbpar, nompar, valpar, &
                            vpesa4(i), ier)
                vpesa4(i) = vpesa4(i)/sec
            end do
        end if
        do i = 4, 6
            vpesa1(i) = 0.d0
            vpesa2(i) = 0.d0
            vpesa3(i) = 0.d0
            vpesa4(i) = 0.d0
        end do
        do i = 1, nbrddl
            f(i) = 0.d0
        end do
        if (icoude .eq. 0) then
            call utpvgl(1, 6, pgl, vpesa1(1), fpesa1(1))
            call utpvgl(1, 6, pgl, vpesa2(1), fpesa2(1))
            call utpvgl(1, 6, pgl, vpesa3(1), fpesa3(1))
            if (nno .eq. 4) then
                call utpvgl(1, 6, pgl, vpesa4(1), fpesa4(1))
            end if
        else
            call utpvgl(1, 6, pgl1, vpesan(1), fpesa1(1))
            call utpvgl(1, 6, pgl2, vpesan(1), fpesa2(1))
            call utpvgl(1, 6, pgl3, vpesan(1), fpesa3(1))
            if (nno .eq. 4) then
                call utpvgl(1, 6, pgl4, vpesan(1), fpesa4(1))
            end if
        end if
!       BOUCLE SUR LES POINTS DE GAUSS DANS LA LONGUEUR
        do igau = 1, npg
!           BOUCLE SUR LES POINTS DE SIMPSON DANS L'EPAISSEUR
            do icou = 1, 2*nbcou+1
                r = a+(icou-1)*h/(2.d0*nbcou)-h/2.d0
!               BOUCLE SUR LES POINTS DE SIMPSON SUR LA CIRCONFERENCE
                do isect = 1, 2*nbsec+1
                    if (icoude .eq. 0) then
                        poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deuxpi/&
                                & (4.d0*nbcou*nbsec)*r
                        do k = 1, 3
                            hk = zr(ivf-1+nno*(igau-1)+1)
                            ibloc = (9+6*(m-1))*(1-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa1(k)
                            hk = zr(ivf-1+nno*(igau-1)+2)
                            ibloc = (9+6*(m-1))*(2-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa2(k)
                            hk = zr(ivf-1+nno*(igau-1)+3)
                            ibloc = (9+6*(m-1))*(3-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa3(k)
                            if (nno .eq. 4) then
                                ibloc = (9+6*(m-1))*(4-1)
                                f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa4(k)
                            end if
                        end do
                    else if (icoude .eq. 1) then
                        fi = (isect-1)*deuxpi/(2.d0*nbsec)
                        cosfi = cos(fi)
                        sinfi = sin(fi)
                        l = theta*(rayon+r*sinfi)
                        poids = zr(ipoids-1+igau)*poicou(icou)*poisec(isect)*(l/2.d0)*h*deuxpi/&
                                & (4.d0*nbcou*nbsec)*r
                        do k = 1, 3
                            hk = zr(ivf-1+nno*(igau-1)+1)
                            ibloc = (9+6*(m-1))*(1-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa1(k)
                            hk = zr(ivf-1+nno*(igau-1)+2)
                            ibloc = (9+6*(m-1))*(2-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa2(k)
                            hk = zr(ivf-1+nno*(igau-1)+3)
                            ibloc = (9+6*(m-1))*(3-1)
                            f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa3(k)
                            if (nno .eq. 4) then
                                ibloc = (9+6*(m-1))*(4-1)
                                f(ibloc+k) = f(ibloc+k)+poids*hk*fpesa4(k)
                            end if
                        end do
                    end if
                end do
            end do
        end do
!
!
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
