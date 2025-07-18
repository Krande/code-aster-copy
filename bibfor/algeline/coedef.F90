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
subroutine coedef(imod, fremod, nbm, young, poiss, &
                  rho, icoq, nbno, numno, nunoe0, &
                  nbnoto, coordo, iaxe, kec, geom, &
                  defm, drmax, torco, tcoef)
    implicit none
! CARACTERISATION DES DEFORMEES DES MODES PRIS EN COMPTE POUR LE
! COUPLAGE : DETERMINATION DES COEFFICIENTS DE LA DEFORMEE AXIALE
! APPELANT : MODCOQ
!-----------------------------------------------------------------------
!  IN : IMOD   : INDICE DU MODE CONSIDERE
!  IN : FREMOD : FREQUENCE DU MODE CONSIDERE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : YOUNG  : MODULE D'YOUNG DU MATERIAU STRUCTUREL
!  IN : POISS  : COEFFICIENT DE POISSON DU MATERIAU STRUCTUREL
!  IN : RHO    : MASSE VOLUMIQUE DU MATERIAU STRUCTUREL
!  IN : ICOQ   : INDICE CARACTERISANT LA COQUE SUR LAQUELLE ON TRAVAILLE
!                ICOQ=1 COQUE INTERNE   ICOQ=2 COQUE EXTERNE
!  IN : NBNO   : NOMBRE DE NOEUDS APPARTENANT AU GROUPE DE NOEUDS SUR
!                LEQUEL ON VA EXTRAIRE LES DEPLACEMENTS
!  IN : NUMNO  : LISTE DES NUMEROS DES NOEUDS APPARTENANT A CE GROUPE
!  IN : NUNOE0 : NUMERO DU NOEUD CORRESPONDANT AU DEPLACEMENT RADIAL
!                MAXIMUM SUR LE MAILLAGE
!  IN : NBNOTO : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
!  IN : COORDO : LISTE DES COORDONNEES DE TOUS LES NOEUDS DU MAILLAGE
!  IN : IAXE   : INDICE CARACTERISANT L'AXE DE REVOLUTION DES COQUES
!                IAXE = 1 : AXE X DU REPERE GLOBAL
!                IAXE = 2 : AXE Y DU REPERE GLOBAL
!                IAXE = 3 : AXE Z DU REPERE GLOBAL
!  IN : KEC    : INDICE CARACTERISTIQUE DU SENS DE L'ECOULEMENT
!                KEC =  1 : ECOULEMENT DANS LE SENS CROISSANT DU
!                PARAMETRE LE LONG DE L'AXE DE REVOLUTION DES COQUES
!                KEC = -1 : ECOULEMENT DANS LE SENS DECROISSANT
!  IN : GEOM   : VECTEUR DE GRANDEURS GEOMETRIQUES CARACTERISTIQUES
!  IN : DEFM   : TABLEAU DES DEFORMEES MODALES (DDL DE TRANSLATION DANS
!                LES DEUX DIRECTIONS ORTHOGONALES A L'AXE DE REVOLUTION)
!  IN : DRMAX  : VALEUR DU DEPLACEMENT RADIAL MAXIMUM
!  IN : TORCO  : TABLEAU CONTENANT LES ORDRES DE COQUE ET DEPHASAGES
! OUT : TCOEF  : TABLEAU CONTENANT LES COEFFICIENTS DES DEFORMEES
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mgauss.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: imod, nbm, icoq, nbno, numno(nbno), nunoe0, nbnoto, iaxe, kec
    integer(kind=8) :: iret
    real(kind=8) :: fremod, young, poiss, rho, coordo(3, nbnoto), geom(9), det
    real(kind=8) :: defm(2, nbnoto, nbm), drmax, torco(4, nbm), tcoef(10, nbm)
!
    real(kind=8) :: long, long2, long4
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ib, ic, icoea, idec, idec1, idir1, idir2
    integer(kind=8) :: idr, il, iligne, imata, ino, inogen, itab
    integer(kind=8) :: ivecb, iz, izaxe, nbnoge, numnoe
    real(kind=8) :: a, b, coshiz, cosiz, deltai, difx, dify
    real(kind=8) :: dr, drnoe0, dx, dy, epais, err, fonc
    real(kind=8) :: orcoq, orcq2, orcq4, pi, rapsi, rayo2
    real(kind=8) :: rayo4, rayon, sinhiz, siniz, somm1, somm2, thetag
    real(kind=8) :: tole, x, xgen, xno, xx, y, ygen
    real(kind=8) :: yno, z, z0, z1, zno, zz, zznoe0
!
!-----------------------------------------------------------------------
    call jemarq()
!
!
! --- 1.INITIALISATIONS
!
    pi = r8pi()
    tole = 100.d0*r8prem()
!
    if (iaxe .eq. 1) then
        idir1 = 2
        idir2 = 3
    else if (iaxe .eq. 2) then
        idir1 = 3
        idir2 = 1
    else
        idir1 = 1
        idir2 = 2
    end if
!
    long = geom(3)
    z0 = geom(4)
    z1 = geom(5)
    if (icoq .eq. 1) then
        epais = geom(6)
        rayon = geom(8)
    else
        epais = geom(7)
        rayon = geom(9)
    end if
!
    iligne = 1
    if (icoq .eq. 2) iligne = 3
    itab = 5*(icoq-1)
!
!
! --- 2.DETERMINATION DU NOMBRE D'ONDES
!
    orcoq = torco(iligne, imod)
!
!-----SI MODE DE POUTRE
    if (orcoq .eq. 1.d0) then
        rapsi = 2.d0/(rayon*rayon+epais*epais/4.d0)
        xx = dble(sqrt(rho/young*rapsi))
        tcoef(1+itab, imod) = long*dble(sqrt(2.d0*pi*fremod*xx))
!-----SINON (MODE DE COQUE)
    else
        orcq2 = orcoq*orcoq
        orcq4 = orcq2*orcq2
        rayo2 = rayon*rayon
        rayo4 = rayo2*rayo2
        long2 = long*long
        long4 = long2*long2
        x = 2.d0*pi*fremod
        x = x*x*(1.d0+1.d0/orcq2)
        a = (1.d0-poiss*poiss)*rho*x/young
        y = (orcq2-1.d0)*(orcq2-1.d0)
        b = epais*epais*y/(12.d0*rayo4)
        z = long4*orcq4*(a-b)/rayo2
        if (z .le. 0.d0) then
            call utmess('F', 'ALGELINE_24')
        else
            z = dble(sqrt(z))
            tcoef(1+itab, imod) = dble(sqrt(z))
        end if
    end if
!
!
! --- 3.EXTRACTION DE LA GENERATRICE SUR LAQUELLE ON VA CALCULER LES
! ---   DEPLACEMENTS RADIAUX
!
! --- 3.1.RECUPERATION DES DEUX COORDONNEES DEFINISSANT LA GENERATRICE
! ---     ET DE L'AZIMUT CORRESPONDANT
!
    xgen = coordo(idir1, nunoe0)
    ygen = coordo(idir2, nunoe0)
    thetag = torco(iligne+1, imod)
!
! --- 3.2.EXTRACTION DES NOEUDS APPARTENANT A LA GENERATRICE
!
    call wkvect('&&COEDEF.TEMP.NOGEN', 'V V I', nbno, inogen)
    nbnoge = 0
    do ino = 1, nbno
        numnoe = numno(ino)
        xno = coordo(idir1, numnoe)
        difx = dble(abs(xno-xgen))
        if (difx .lt. tole) then
            yno = coordo(idir2, numnoe)
            dify = dble(abs(yno-ygen))
            if (dify .lt. tole) then
                zno = coordo(iaxe, numnoe)
                if (zno .ge. z0 .and. zno .le. z1) then
                    nbnoge = nbnoge+1
                    zi(inogen+nbnoge-1) = numnoe
                end if
            end if
        end if
    end do
!
    if (nbnoge .lt. 4) then
        call utmess('F', 'ALGELINE_25')
    end if
!
!
! --- 4.CREATION D'UNE DISCRETISATION LE LONG DE LA GENERATRICE
! ---   SIMULTANEMENT ON CREE UNE DISTRIBUTION EN DR
!
    call wkvect('&&COEDEF.TEMP.ZAXE', 'V V R', nbnoge, izaxe)
    call wkvect('&&COEDEF.TEMP.DR  ', 'V V R', nbnoge, idr)
!
    do ino = 1, nbnoge
        numnoe = zi(inogen+ino-1)
        if (kec .eq. 1) then
            zr(izaxe+ino-1) = coordo(iaxe, numnoe)-z0
        else
            zr(izaxe+ino-1) = z1-coordo(iaxe, numnoe)
        end if
        dx = defm(1, numnoe, imod)
        dy = defm(2, numnoe, imod)
        zr(idr+ino-1) = dble(cos(thetag))*dx+dble(sin(thetag))*dy
    end do
!
!
! --- 5.DETERMINATION DES COEFFICIENTS DE LA DEFORMEE AXIALE
!
! --- 5.1.CALCUL DES COEFFICIENTS DU TRIANGLE SUPERIEUR DE LA MATRICE A
! ---    (PAR COLONNES DESCENDANTES) ET DES COEFFICIENTS DU SECOND
! ---     MEMBRE B
!
    call wkvect('&&COEDEF.TEMP.COEA', 'V V R', 10, icoea)
    call wkvect('&&COEDEF.TEMP.VECB', 'V V R', 4, ivecb)
!
    deltai = tcoef(1+itab, imod)
!
    do iz = 1, nbnoge
        zz = deltai*zr(izaxe+iz-1)/long
        dr = zr(idr+iz-1)
        cosiz = dble(cos(zz))
        siniz = dble(sin(zz))
        coshiz = dble(cosh(zz))
        sinhiz = dble(sinh(zz))
        zr(icoea) = zr(icoea)+cosiz*cosiz
        zr(icoea+1) = zr(icoea+1)+siniz*cosiz
        zr(icoea+2) = zr(icoea+2)+siniz*siniz
        zr(icoea+3) = zr(icoea+3)+cosiz*coshiz
        zr(icoea+4) = zr(icoea+4)+siniz*coshiz
        zr(icoea+5) = zr(icoea+5)+coshiz*coshiz
        zr(icoea+6) = zr(icoea+6)+cosiz*sinhiz
        zr(icoea+7) = zr(icoea+7)+siniz*sinhiz
        zr(icoea+8) = zr(icoea+8)+coshiz*sinhiz
        zr(icoea+9) = zr(icoea+9)+sinhiz*sinhiz
        zr(ivecb) = zr(ivecb)+dr*cosiz
        zr(ivecb+1) = zr(ivecb+1)+dr*siniz
        zr(ivecb+2) = zr(ivecb+2)+dr*coshiz
        zr(ivecb+3) = zr(ivecb+3)+dr*sinhiz
    end do
!
! --- 5.2.CREATION DE LA MATRICE A
!
    call wkvect('&&COEDEF.TEMP.MATA', 'V V R', 16, imata)
!
    do ic = 1, 4
        idec1 = ic*(ic-1)/2
        do il = 1, ic
            idec = idec1+il
            zr(imata+4*(ic-1)+il-1) = zr(icoea+idec-1)
        end do
    end do
    do ic = 1, 3
        do il = ic+1, 4
            idec = il*(il-1)/2+ic
            zr(imata+4*(ic-1)+il-1) = zr(icoea+idec-1)
        end do
    end do
!
! --- 5.3.RESOLUTION DU SYSTEME LINEAIRE
!
    call mgauss('NFVP', zr(imata), zr(ivecb), 4, 4, &
                1, det, iret)
!
! --- 5.4.CALCUL D'UNE ERREUR RELATIVE SUR LA NORME DES DR
!
    somm1 = 0.d0
    somm2 = 0.d0
!
    do iz = 1, nbnoge
        dr = zr(idr+iz-1)
        somm1 = somm1+dr*dr
        zz = deltai*zr(izaxe+iz-1)/long
        cosiz = dble(cos(zz))
        siniz = dble(sin(zz))
        coshiz = dble(cosh(zz))
        sinhiz = dble(sinh(zz))
        fonc = zr(ivecb)*cosiz+zr(ivecb+1)*siniz+zr(ivecb+2)*coshiz+zr(ivecb+3)*sinhiz
        somm2 = somm2+(dr-fonc)*(dr-fonc)
    end do
!
    err = dble(sqrt(somm2/somm1))*100.d0
!
! --- 5.5.DEDUCTION DES COEFFICIENTS DE LA DEFORMEE AXIALE
!
    if (kec .eq. 1) then
        zznoe0 = deltai*(coordo(iaxe, nunoe0)-z0)/long
    else
        zznoe0 = deltai*(z1-coordo(iaxe, nunoe0))/long
    end if
    drnoe0 = zr(ivecb)*dble(cos(zznoe0))+zr(ivecb+1)*dble(sin(zznoe0))+zr(ivecb+2)*dble(cosh&
             &(zznoe0))+zr(ivecb+3)*dble(sinh(zznoe0))
    if (dble(abs(drnoe0)) .lt. drmax*tole) then
        call utmess('F', 'ALGELINE_26')
    end if
    do ib = 1, 4
        tcoef(1+itab+ib, imod) = zr(ivecb+ib-1)*drmax/drnoe0
    end do
!
!
! --- 6.IMPRESSION DES RESULTATS DANS LE FICHIER MESSAGE
!
    call utmess('I', 'ALGELINE_28', sr=err)
    if (err .gt. 1.d0) then
        call utmess('A', 'ALGELINE_29')
    end if
!
    call jedetr('&&COEDEF.TEMP.ZAXE')
    call jedetr('&&COEDEF.TEMP.DR  ')
    call jedetr('&&COEDEF.TEMP.COEA')
    call jedetr('&&COEDEF.TEMP.VECB')
    call jedetr('&&COEDEF.TEMP.MATA')
    call jedetr('&&COEDEF.TEMP.NOGEN')
    call jedema()
!
end subroutine
