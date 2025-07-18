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
subroutine ordcoq(imod, nbm, icoq, nbno, numno, &
                  inomax, nbnoto, coordo, iaxe, defm, &
                  nunoe0, drmax, torco)
! aslint: disable=
    implicit none
! CARACTERISATION DES DEFORMEES DES MODES PRIS EN COMPTE POUR LE
! COUPLAGE : DETERMINATION DE L'ORDRE DE COQUE
! APPELANT : MODCOQ
!-----------------------------------------------------------------------
!  IN : IMOD   : INDICE DU MODE CONSIDERE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : ICOQ   : INDICE CARACTERISANT LA COQUE SUR LAQUELLE ON TRAVAILLE
!                ICOQ = 1 COQUE INTERNE  ICOQ = 2 COQUE EXTERNE
!  IN : NBNO   : NOMBRE DE NOEUDS APPARTENANT AU GROUPE DE NOEUDS SUR
!                LEQUEL ON VA EXTRAIRE LES DEPLACEMENTS
!  IN : NUMNO  : LISTE DES NUMEROS DES NOEUDS APPARTENANT A CE GROUPE
!  IN : INOMAX : INDICE DU NOEUD POUR LEQUEL ON A RELEVE LE DEPLACEMENT
!                MAXIMUM DANS LE PLAN PERPENDICULAIRE A L'AXE DE
!                REVOLUTION
!  IN : NBNOTO : NOMBRE TOTAL DE NOEUDS DU MAILLAGE
!  IN : COORDO : LISTE DES COORDONNEES DE TOUS LES NOEUDS DU MAILLAGE
!  IN : IAXE   : INDICE CARACTERISANT L'AXE DE REVOLUTION DES COQUES
!                IAXE = 1 : AXE X DU REPERE GLOBAL
!                IAXE = 2 : AXE Y DU REPERE GLOBAL
!                IAXE = 3 : AXE Z DU REPERE GLOBAL
!  IN : DEFM   : TABLEAU DES DEFORMEES MODALES (DDL DE TRANSLATION DANS
!                LES DEUX DIRECTIONS ORTHOGONALES A L'AXE DE REVOLUTION)
! OUT : NUNOE0 : NUMERO DU NOEUD CORRESPONDANT AU DEPLACEMENT RADIAL
!                MAXIMUM SUR LE MAILLAGE
! OUT : DRMAX  : VALEUR DU DEPLACEMENT RADIAL MAXIMUM
! OUT : TORCO  : TABLEAU CONTENANT LES ORDRES DE COQUE ET DEPHASAGES
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterc/r8prem.h"
#include "asterc/r8rddg.h"
#include "asterfort/dcabs2.h"
#include "asterfort/fft.h"
#include "asterfort/fointr.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: imod, nbm, icoq, nbno, numno(nbno), inomax, nbnoto, iaxe
    real(kind=8) :: coordo(3, nbnoto), defm(2, nbnoto, nbm)
    integer(kind=8) :: nunoe0
    real(kind=8) :: drmax, torco(4, nbm)
!
    integer(kind=8) :: pui
    real(kind=8) :: modmax, module
    character(len=3) :: kmod
    character(len=8) :: nompar
    aster_logical :: defini
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idir1, idir2, idisc, idr, idr2, ier
    integer(kind=8) :: ifftdr, ifm, iligne, imin, ino, inocr, ip
    integer(kind=8) :: ipmax, iprol, itheta, jno, nbnocr, nbpt
    integer(kind=8) :: nbpt2, numnoe, nunomx
    real(kind=8) :: alpha, an, beta, bn, crit, delta, dif
    real(kind=8) :: dif1, dif2, difz, dr, err, fonc, gamma
    real(kind=8) :: pas, pi, rap, rki, s11, s12
    real(kind=8) :: s22, somm, somm1, somm2, themin, theta, theta0
    real(kind=8) :: theta1, theta2, tole, x1, x2, y1, y2
    real(kind=8) :: zcr, zno
!-----------------------------------------------------------------------
    call jemarq()
!
!
! --- 1.INITIALISATIONS
!
    tole = 100.d0*r8prem()
    pi = r8pi()
!
    if (iaxe .eq. 1) then
        idir1 = 2
        idir2 = 3
        nompar = 'X'
    else if (iaxe .eq. 2) then
        idir1 = 3
        idir2 = 1
        nompar = 'Y'
    else
        idir1 = 1
        idir2 = 2
        nompar = 'Z'
    end if
!
    iligne = 1
    if (icoq .eq. 2) iligne = 3
!
!
! --- 2.DEDUCTION DE L'ALTITUDE DU CONTOUR SUR LEQUEL ON VA CALCULER
! ---   LES DEPLACEMENTS RADIAUX
!
    nunomx = numno(inomax)
    zcr = coordo(iaxe, nunomx)
!
!
! --- 3.EXTRACTION DU CONTOUR
!
    call wkvect('&&ORDCOQ.TEMP.NOCR', 'V V I', nbno, inocr)
    nbnocr = 0
    do ino = 1, nbno
        numnoe = numno(ino)
        zno = coordo(iaxe, numnoe)
        difz = dble(abs(zno-zcr))
        if (difz .lt. tole) then
            nbnocr = nbnocr+1
            zi(inocr+nbnocr-1) = numnoe
        end if
    end do
!
    if (nbnocr .le. 8) then
        call utmess('F', 'ALGELINE3_12')
    end if
!
!
! --- 4.CREATION D'UNE DISCRETISATION EN THETA SUR LE CONTOUR
! ---   SIMULTANEMENT ON CREE UNE DISTRIBUTION EN DR
!
    call wkvect('&&ORDCOQ.TEMP.THETA', 'V V R', nbnocr, itheta)
    call wkvect('&&ORDCOQ.TEMP.DR', 'V V R', nbnocr, idr)
!
    do ino = 1, nbnocr
        numnoe = zi(inocr+ino-1)
        x1 = coordo(idir1, numnoe)
        x2 = coordo(idir2, numnoe)
        defini = .true.
        if (x1 .eq. 0.d0) then
            if (x2 .gt. 0.d0) then
                theta = pi/2.d0
            else if (x2 .lt. 0.d0) then
                theta = -pi/2.d0
            else
                defini = .false.
            end if
        else if (x1 .gt. 0.d0) then
            theta = dble(atan(x2/x1))
        else if (x1 .lt. 0.d0) then
            theta = dble(atan(x2/x1))+pi
        end if
        if (defini) then
            if (theta .lt. 0.d0) theta = theta+2.d0*pi
        else
            call utmess('F', 'ALGELINE3_13')
        end if
        dr = dble(cos(theta))*defm(1, numnoe, imod)+dble(sin(theta))*defm(2, numnoe, imod)
        zr(itheta+ino-1) = theta
        zr(idr+ino-1) = dr
    end do
!
!
! --- 5.ON REORDONNE LA DISCRETISATION PAR VALEURS CROISSANTES DE THETA
! ---   SIMULTANEMENT ON REORDONNE LA DISTRIBUTION EN DR ET LA LISTE DES
! ---   NUMEROS DES NOEUDS
!
    do ino = 1, nbnocr-1
        themin = zr(itheta+ino-1)
        imin = ino
        do jno = ino+1, nbnocr
            if (zr(itheta+jno-1) .lt. themin) then
                themin = zr(itheta+jno-1)
                imin = jno
            end if
        end do
        zr(itheta+imin-1) = zr(itheta+ino-1)
        zr(itheta+ino-1) = themin
        dr = zr(idr+imin-1)
        zr(idr+imin-1) = zr(idr+ino-1)
        zr(idr+ino-1) = dr
        numnoe = zi(inocr+imin-1)
        zi(inocr+imin-1) = zi(inocr+ino-1)
        zi(inocr+ino-1) = numnoe
    end do
!
!
! --- 6.CONSTRUCTION D'UNE DISCRETISATION REGULIERE EN THETA SUR 0,2*PI
! ---   D'UN NOMBRE DE POINTS EGAL A LA PREMIERE PUISSANCE DE 2
! ---   SUPERIEURE A NBNOCR
! ---   CALCUL D'UNE NOUVELLE DISTRIBUTION EN DR PAR INTERPOLATION
!
! --- 6.1.DETERMINATION DE NBPT
!
    rap = dble(log(dble(nbnocr))/log(2.d0))
    pui = int(rap)
    dif = rap-dble(pui)
    if (dif .gt. tole) pui = pui+1
    nbpt = 1
    do i = 1, pui
        nbpt = nbpt*2
    end do
!
! --- 6.2.CONSTRUCTION DE LA DISCRETISATION REGULIERE EN THETA
!
    call wkvect('&&ORDCOQ.TEMP.DISC', 'V V R', nbpt, idisc)
    pas = 2.d0*pi/dble(nbpt)
    do ip = 1, nbpt
        zr(idisc+ip-1) = dble(ip-1)*pas
    end do
!
! --- 6.3.CALCUL D'UNE NOUVELLE DISTRIBUTION EN DR PAR INTERPOLATION
!
    call wkvect('&&ORDCOQ.TEMP.PROL', 'V V K24', 6, iprol)
    zk24(iprol) = 'FONCTION'
    zk24(iprol+1) = 'LIN LIN '
    zk24(iprol+2) = nompar
    zk24(iprol+3) = 'TOUTRESU'
    zk24(iprol+4) = 'LL      '
    zk24(iprol+5) = '&&ORDCOQ.TEMP'
!
    call wkvect('&&ORDCOQ.TEMP.DR2', 'V V R', nbpt, idr2)
    call fointr(' ', zk24(iprol), nbnocr, zr(itheta), zr(idr), &
                nbpt, zr(idisc), zr(idr2), ier)
!
!
! --- 7.FFT SUR LA DISTRIBUTION EN DR
!
    call wkvect('&&ORDCOQ.TEMP.FFTDR', 'V V C', nbpt, ifftdr)
    do ip = 1, nbpt
        zc(ifftdr+ip-1) = dcmplx(zr(idr2+ip-1), 0.d0)
    end do
    call fft(zc(ifftdr), nbpt, 1)
!
!
! --- 8.DETERMINATION DE L'ORDRE DE COQUE
! ---   CALCUL DU MODULE DES NBPT/2 PREMIERES VALEURS DE LA
! ---   TRANSFORMEE DE DR
! ---   SIMULTANEMENT ON DETECTE LE PIC ET ON CALCULE UN CRITERE
! ---   DE PRECISION (ENERGIE DU PIC/ENERGIE DU BRUIT)
!
    modmax = dcabs2(zc(ifftdr))
    somm = modmax*modmax
    ipmax = 1
    nbpt2 = nbpt/2
    do ip = 2, nbpt2
        module = dcabs2(zc(ifftdr+ip-1))
        somm = somm+module*module
        if (module .gt. modmax) then
            modmax = module
            ipmax = ip
        end if
    end do
    crit = 1.d0-modmax*modmax/somm
!
    if (ipmax .eq. 1) then
        call utmess('F', 'ALGELINE3_14')
    end if
!
    rki = dble(ipmax-1)
    torco(iligne, imod) = rki
!
!
! --- 9.DETERMINATION DU DRMAX ET DU DEPHASAGE
!
! --- 9.1.RESOLUTION DU PROBLEME DE MOINDRES CARRES
!
    s11 = 0.d0
    s12 = 0.d0
    s22 = 0.d0
    y1 = 0.d0
    y2 = 0.d0
    do ino = 1, nbnocr
        theta = zr(itheta+ino-1)
        dr = zr(idr+ino-1)
        an = dble(cos(rki*theta))
        bn = dble(sin(rki*theta))
        s11 = s11+an*an
        s12 = s12+an*bn
        s22 = s22+bn*bn
        y1 = y1+dr*an
        y2 = y2+dr*bn
    end do
!
    delta = s11*s22-s12*s12
    if (dble(abs(delta)) .lt. tole) then
        write (kmod, '(I3)') imod
        call utmess('F', 'ALGELINE3_15', sk=kmod)
    end if
    alpha = (s22*y1-s12*y2)/delta
    beta = (s11*y2-s12*y1)/delta
!
! --- 9.2.DEDUCTION DU DRMAX ET DU THETA0
!
    drmax = dble(sqrt(alpha*alpha+beta*beta))
!
    defini = .true.
    if (alpha .eq. 0.d0) then
        if (beta .gt. 0.d0) then
            theta0 = pi/2.d0
        else if (beta .lt. 0.d0) then
            theta0 = -pi/2.d0
        else
            defini = .false.
        end if
    else if (alpha .gt. 0.d0) then
        theta0 = dble(atan(beta/alpha))
    else if (alpha .lt. 0.d0) then
        theta0 = dble(atan(beta/alpha))+pi
    end if
    if (defini) then
        if (theta0 .lt. 0.d0) theta0 = theta0+2.d0*pi
    else
        write (kmod, '(I3)') imod
        call utmess('F', 'ALGELINE3_16', sk=kmod)
    end if
    theta0 = theta0/rki
    torco(iligne+1, imod) = theta0
!
! --- 9.3.CALCUL D'UNE ERREUR RELATIVE SUR LA NORME DES DR
!
    somm1 = 0.d0
    somm2 = 0.d0
!
    do ino = 1, nbnocr
        dr = zr(idr+ino-1)
        somm1 = somm1+dr*dr
        gamma = rki*(zr(itheta+ino-1)-theta0)
        fonc = drmax*dble(cos(gamma))
        somm2 = somm2+(dr-fonc)*(dr-fonc)
    end do
!
    err = dble(sqrt(somm2/somm1))*100.d0
!
!
! --- 10.RECUPERATION DU NUMERO DU NOEUD SUR LE CONTOUR DONT L'AZIMUT
! ---    EST LE PLUS PROCHE DE THETA0
!
    do ino = 1, nbnocr
        theta2 = zr(itheta+ino-1)
        if (theta2 .gt. theta0) goto 101
    end do
101 continue
    if (ino .eq. 1) then
        theta1 = zr(itheta+nbnocr-1)-2.d0*pi
        dif1 = dble(abs(theta1-theta0))
        dif2 = dble(abs(theta2-theta0))
        if (dif1 .lt. dif2) then
            nunoe0 = zi(inocr+nbnocr-1)
        else
            nunoe0 = zi(inocr)
        end if
    else if (ino .eq. nbnocr+1) then
        theta1 = zr(itheta)+2.d0*pi
        dif1 = dble(abs(theta1-theta0))
        dif2 = dble(abs(theta2-theta0))
        if (dif1 .lt. dif2) then
            nunoe0 = zi(inocr)
        else
            nunoe0 = zi(inocr+nbnocr-1)
        end if
    else
        theta1 = zr(itheta+ino-2)
        dif1 = dble(abs(theta1-theta0))
        dif2 = dble(abs(theta2-theta0))
        if (dif1 .lt. dif2) then
            nunoe0 = zi(inocr+ino-2)
        else
            nunoe0 = zi(inocr+ino-1)
        end if
    end if
!
!
! --- 11.IMPRESSION DES RESULTATS DANS LE FICHIER MESSAGE
!
500 format('DETERMINATION DE L ORDRE DE COQUE')
501 format('FFT SUR DR : LE BRUIT AUTOUR DU PIC REPRESENTE UNE ',&
     &       'ENERGIE RESIDUELLE DE ', g13.6, ' %')
502 format('KI = ', i3)
510 format('DETERMINATION DU DEPHASAGE')
511 format('ERREUR RELATIVE SUR LA NORME DES DEPLACEMENTS RADIAUX : ',&
     &        g13.6, 1x, '%')
512 format('THETA0 = ', g13.6, ' DEGRES')
!
    ifm = iunifi('MESSAGE')
    write (ifm, 500)
    write (ifm, 501) crit
    write (ifm, 502) ipmax-1
    write (ifm, 510)
    write (ifm, 511) err
    write (ifm, 512) theta0*r8rddg()
!
    call jedetr('&&ORDCOQ.TEMP.NOCR')
    call jedetr('&&ORDCOQ.TEMP.THETA')
    call jedetr('&&ORDCOQ.TEMP.DR')
    call jedetr('&&ORDCOQ.TEMP.DISC')
    call jedetr('&&ORDCOQ.TEMP.PROL')
    call jedetr('&&ORDCOQ.TEMP.DR2')
    call jedetr('&&ORDCOQ.TEMP.FFTDR')
    call jedema()
!
end subroutine
