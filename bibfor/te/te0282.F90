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
subroutine te0282(option, nomte)
!
!      CALCUL DU TAUX DE RESTITUTION D'ENERGIE ELEMENTAIRE
!      BORDS ELEMENTS ISOPARAMETRIQUES 2D AVEC CHARGEMENT DE BORD
!      PRESSION-CISAILLEMENT ET FORCE REPARTIE
!
!      OPTION : 'CALC_G'    (G AVEC CHARGES REELLES)
!               'CALC_G_F'  (G AVEC CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!
! VECTEURS DIMENSIONNES POUR  NNO = 3 , NPG = 4
!
! ----------------------------------------------------------------------
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
    character(len=16) :: option, nomte
!
!
    integer(kind=8) :: nno, nnos, jgano, ndim, npg, kp, ipoids, ivf, idfdk, igeom, icode
    integer(kind=8) :: idepl, iforc, ipres, ithet, igthet, itemps, compt, i, j, k
    integer(kind=8) :: ipref, iforf
    integer(kind=8) :: jdfd2, jcoopg
!
    real(kind=8) :: xg, yg, ux, uy, fx, fy, thx, thy, the
    real(kind=8) :: tcla, tsurf, tsurp, epsi, pres, cisa, divthe, valpar(3)
    real(kind=8) :: vf, dfde, dxde, dyde, dsde, poids, dthxde, dthyde
    real(kind=8) :: dfxde, dfyde, presno, cisano, fxno, fyno
!                               2*NNO     2*NNO
    real(kind=8) :: presg(2), forcg(2), presn(6), forcn(6)
    real(kind=8) :: prod, dsde2
    real(kind=8) :: tsom
!
    character(len=8) :: nompar(3), elrefe
!
    aster_logical :: fonc, chargn, axis
!
! =====================================================================
! INITIALISATIONS
! =====================================================================
    call elref1(elrefe)
    epsi = r8prem()
    axis = lteatt('AXIS', 'OUI')
!
! RECUPERATION DES DONNEES GEOMETRIQUES LIEES AU CALCUL ELEMENTAIRE
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=jcoopg, jvf=ivf, jdfde=idfdk, jdfd2=jdfd2, &
                     jgano=jgano)
!
!
! INIT. POUR LE CALCUL DE G
    chargn = .false.
    tcla = 0.d0
    tsurf = 0.d0
    tsurp = 0.d0
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PGTHETA', 'E', igthet)
!
! TEST SUR LA NULLITE DE THETA_FISSURE
    compt = 0
    do i = 1, nno
        thx = zr(ithet+2*(i-1))
        thy = zr(ithet+2*(i-1)+1)
        if ((abs(thx) .lt. epsi) .and. (abs(thy) .lt. epsi)) then
            compt = compt+1
        end if
    end do
    if (compt .eq. nno) goto 999
!
! =====================================================================
! RECUPERATION DES CHAMPS LOCAUX
! =====================================================================
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    if (option .eq. 'CALC_G_XFEM_F') then
        fonc = .true.
        call jevech('PFF1D2D', 'L', iforf)
        call jevech('PPRESSF', 'L', ipref)
        call jevech('PINSTR', 'L', itemps)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'INST'
        valpar(3) = zr(itemps)
    else
        fonc = .false.
        call jevech('PFR1D2D', 'L', iforc)
        call jevech('PPRESSR', 'L', ipres)
    end if
!
! =====================================================================
! - SI CHARGE FONCTION RECUPERATION DES VALEURS AUX PG ET NOEUDS
! =====================================================================
!
    if (fonc) then
        do i = 1, nno
            do j = 1, 2
                valpar(j) = zr(igeom+2*(i-1)+j-1)
            end do
            do j = 1, 2
                call fointe('FM', zk8(ipref+j-1), 3, nompar, valpar, &
                            presn(2*(i-1)+j), icode)
                call fointe('FM', zk8(iforf+j-1), 3, nompar, valpar, &
                            forcn(2*(i-1)+j), icode)
            end do
        end do
    end if
!
! ======================================================================
! BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! ======================================================================
!
    do kp = 1, npg
!
! INITIALISATIONS
        k = (kp-1)*nno
        dxde = 0.d0
        dyde = 0.d0
        xg = 0.d0
        yg = 0.d0
        ux = 0.d0
        uy = 0.d0
        thx = 0.d0
        thy = 0.d0
        dfxde = 0.d0
        dfyde = 0.d0
        dthxde = 0.d0
        dthyde = 0.d0
        fx = 0.d0
        fy = 0.d0
!
! ===========================================
! CALCUL DES ELEMENTS GEOMETRIQUES
! ===========================================
!
! CALCUL DES DERIVEES PARTIELLES PREMIERES DU VECTEURS
! POSITIONS (DXDE,DYDE) AU POINT DE GAUSS,
! DU VECTEUR POSITION AU POINT DE GAUSS (XG,YG), DE SON VECTEUR
! DEPLACEMENT (UX,UY), DU CHAMP THETA FISSURE (THX,THY) ET DE SON
! GRADIENT (DTHXDE,DTHYDE).
        do i = 1, nno
            vf = zr(ivf+k+i-1)
            dfde = zr(idfdk+k+i-1)
            dxde = dxde+dfde*zr(igeom+2*(i-1))
            dyde = dyde+dfde*zr(igeom+2*(i-1)+1)
            xg = xg+vf*zr(igeom+2*(i-1))
            yg = yg+vf*zr(igeom+2*(i-1)+1)
            ux = ux+vf*zr(idepl+2*(i-1))
            uy = uy+vf*zr(idepl+2*(i-1)+1)
            thx = thx+vf*zr(ithet+2*(i-1))
            thy = thy+vf*zr(ithet+2*(i-1)+1)
            dthxde = dthxde+dfde*zr(ithet+2*(i-1))
            dthyde = dthyde+dfde*zr(ithet+2*(i-1)+1)
        end do
!
! ===========================================
! CALCUL DU CHARGEMENT ET DE SON GRADIENT
! ===========================================
!
        if (fonc) then
            valpar(1) = xg
            valpar(2) = yg
            do j = 1, 2
                call fointe('FM', zk8(ipref+j-1), 3, nompar, valpar, &
                            presg(j), icode)
                call fointe('FM', zk8(iforf+j-1), 3, nompar, valpar, &
                            forcg(j), icode)
            end do
        else
            presg(1) = 0.d0
            presg(2) = 0.d0
            forcg(1) = 0.d0
            forcg(2) = 0.d0
            do i = 1, nno
                do j = 1, 2
                    presg(j) = presg(j)+zr(ipres+2*(i-1)+j-1)*zr(ivf+k+i-1)
                    forcg(j) = forcg(j)+zr(iforc+2*(i-1)+j-1)*zr(ivf+k+i-1)
                end do
            end do
        end if
!
! VALEURS DU CHARGEMENT AUX POINTS DE GAUSS (FX,FY)
        dsde = sqrt(dxde*dxde+dyde*dyde)
        dsde2 = dsde*dsde
        pres = presg(1)
        cisa = presg(2)
        fx = forcg(1)-(dyde*pres-dxde*cisa)/dsde
        fy = forcg(2)+(dxde*pres+dyde*cisa)/dsde
!
! VALEURS DU CHARGEMENT AUX NOEUDS (FXNO,FYNO) ET DE SES DERIVEES
! AUX POINTS DE GAUSS (DFXDE,DFYDE,D2FXDE,D2FYDE)
        if (fonc) then
            do i = 1, nno
                dfde = zr(idfdk+k+i-1)
                presno = presn(2*(i-1)+1)
                cisano = presn(2*(i-1)+2)
                fxno = forcn(2*(i-1)+1)-(dyde*presno-dxde*cisano)/dsde
                fyno = forcn(2*(i-1)+2)+(dxde*presno+dyde*cisano)/dsde
                dfxde = dfxde+dfde*fxno
                dfyde = dfyde+dfde*fyno
            end do
        end if
!
! TESTS SUR LA NULLITE DES CHARGEMENTS ET DE LEURS GRADIENTS POUR EVITER
! DE FAIRE DES CALCULS INUTILES ET DETECTER LES VRAIS PROBLEMES
        if ((fx .eq. 0.d0) .and. (fy .eq. 0.d0) .and. (dfxde .eq. 0.d0) &
            .and. (dfyde .eq. 0.d0)) then
            chargn = .true.
        end if
!
! CAS PARTICULIER D'UN CALCUL SUR L'AXE
        if (xg .eq. 0.d0) then
!
! ON EST SUR L'AXE AVEC CHARGEMENTS NULS DONC G (ET DG) = 0
            if (chargn) then
                goto 799
            else if (axis) then
                call utmess('F', 'RUPTURE1_23')
            end if
        else
!
! CAS GENERAL AVEC CHARGEMENTS NULS DONC G (ET DG) = 0
            if (chargn) goto 799
        end if
!
! CALCUL DU TERME ELEMENTAIRE
        if (axis) then
            poids = zr(ipoids+kp-1)*dsde*xg
        else
            poids = zr(ipoids+kp-1)*dsde
        end if
        the = (thx*dxde+thy*dyde)/dsde2
        divthe = (dthxde*dxde+dthyde*dyde)/dsde2
!
! =======================================================
! PRISE EN COMPTE DE LA MODELISATION POUR G
! =======================================================
!
        if (axis) divthe = divthe+(thx/xg)
!
! =======================================================
! CALCUL DU TAUX DE RESTITUTION G
! =======================================================
!
        prod = (divthe*fx+dfxde*the)*ux+(divthe*fy+dfyde*the)*uy
!
        tcla = tcla+prod*poids
!
! BRANCHEMENT POUR F=0 ET DF=0
799     continue
!
! ======================================================================
! FIN DE BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! ======================================================================
    end do
!
! EXIT EN CAS DE THETA FISSURE NUL PARTOUT
999 continue
!
! ASSEMBLAGE FINAL DES TERMES DE G
    tsom = tcla+tsurf+tsurp
    zr(igthet) = tsom
!
end subroutine
