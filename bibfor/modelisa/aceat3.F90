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

subroutine aceat3(noma, nomu, nbtuy, nbpart, nbmap, &
                  elpar, nopar, ivr, nbzk, &
                  nozk, cozk, isens, coor, epsi, &
                  crit, nno, nmmt)
    implicit none
#include "jeveux.h"
#include "asterfort/angco4.h"
#include "asterfort/angcou.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/normev.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
!
    character(len=8) :: noma, nomu, crit
    integer(kind=8) :: nbpart, nbtuy, nbmap(nbpart), elpar(nbpart, nbtuy), ivr(*)
    integer(kind=8) :: nmmt(*), nno, icmp, iavant, no4, nbcmp, icoud2
    integer(kind=8) :: nopar(nbpart, nno, nbtuy), nbzk, nozk(nbzk), isens(nbpart), ifm
    real(kind=8) :: cozk(3*nbzk), coor(*), coor3(12), zk1(3), zk2(3), zk3(3)
    real(kind=8) :: angl1(3), angl2(3), angl3(3), epsi, angl4(3)
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES TUYAUX DANS UNE CARTE
! ----------------------------------------------------------------------
! IN  NOMA   : NOM DU MAILLAGE
! IN  NOMU   : NOM DU CONCEPT RESULTAT
! IN  NBTUY  : NOMBRE D'ELEMENTS TUYAU
! IN  NBPAT  : NOMBRE DE PARTIES CONNEXES DE TUYAUX
! IN  NBMAP  : NOMBRE DE MAILLES PAR PARTIE
! IN  ELPAR  : NUMERO DES MAILLES DE CHAQUE PARTIE DANS L'ORDRE
! IN  NOPAR  : NUMERO DES NOEUDS  DE CHAQUE PARTIE DANS L'ORDRE
! IN  IVR    : (3) = INFO 2
! IN  IFM    : FICHIER MESSAGES
! IN  NBZK   : NOMBRE D'OCCURENCES DE GENE_TUYAU
! IN  NOZK   : NUMEROS DES NOEUDS OU ZK EST DONNE
! IN  COZK   : COORDONNES DES ZK
! IN  ISENS  : SENS DE PARCOURS DES MAILLES
! IN  COOR   : COORDONNES DES NOEUDS
! IN  NMMT   : INDIQUE SI MODI_METRIQUE POUR CHAQUE MAILLE
! ----------------------------------------------------------------------
    character(len=8) :: nommai, nomno1, nomno2, nomno3, nomno4
    character(len=19) :: cartor
    character(len=24) :: tmpnor, tmpvor
    integer(kind=8) :: jdcmpo, jdvlvo, izk, iok1, iok2, ipa, imfin, i, ima, nummai
    integer(kind=8) :: no1, no2, no3, icoude, im0, nbdroi, nbcoud
    real(kind=8) :: pgl1(3, 3), pgl2(3, 3), pgl3(3, 3), norme, pgl4(3, 3)
!
!-----------------------------------------------------------------------
    real(kind=8) :: dn1n2, omega, rayon, theta, vx, vy, vz
!
!-----------------------------------------------------------------------
    call jemarq()
!
!     VERIFICATION QUE LES NOEUDS DONNES SONT DES EXTREMITES
! --- AFFECTATION DES VALEURS DU TAMPON DANS LA CARTE ORIENTATION :
!     -------------------------------------------------------------
    cartor = nomu//'.CARORIEN'
    tmpnor = cartor//'.NCMP'
    tmpvor = cartor//'.VALV'
    call jeveuo(tmpnor, 'E', jdcmpo)
    call jeveuo(tmpvor, 'E', jdvlvo)
    zk8(jdcmpo) = 'ALPHA'
    zk8(jdcmpo+1) = 'BETA'
    zk8(jdcmpo+2) = 'GAMMA'
    zk8(jdcmpo+3) = 'ALPHA2'
    zk8(jdcmpo+4) = 'BETA2'
    zk8(jdcmpo+5) = 'GAMMA2'
    zk8(jdcmpo+6) = 'ALPHA3'
    zk8(jdcmpo+7) = 'BETA3'
    zk8(jdcmpo+8) = 'GAMMA3'
    icmp = 8
    if (nno .eq. 4) then
        zk8(jdcmpo+9) = 'ALPHA4'
        zk8(jdcmpo+10) = 'BETA4'
        zk8(jdcmpo+11) = 'GAMMA4'
        icmp = 11
    end if
    zk8(jdcmpo+icmp+1) = 'ICOUDE'
    zk8(jdcmpo+icmp+2) = 'DN1N2'
    zk8(jdcmpo+icmp+3) = 'RCOURB'
    zk8(jdcmpo+icmp+4) = 'ANGCOU'
    zk8(jdcmpo+icmp+5) = 'ANGZZK'
!
!     LES NOEUDS ASSOCIES A GENE_TUYAUX DOIVENT ETRE UNE
!     DES EXTREMITES
!
    do ipa = 1, nbpart
        isens(ipa) = 0
    end do
    do izk = 1, nbzk
        iok1 = 0
        iok2 = 0
        do ipa = 1, nbpart
            imfin = nbmap(ipa)
            if (nozk(izk) .eq. (nopar(ipa, 1, 1))) then
                iok1 = 1
                if (isens(ipa) .eq. 0) then
                    isens(ipa) = izk
                else
                    call utmess('F', 'MODELISA_24')
                end if
            end if
            if (nozk(izk) .eq. (nopar(ipa, 2, imfin))) then
                iok2 = 1
                if (isens(ipa) .eq. 0) then
                    isens(ipa) = -izk
                else
                    call utmess('F', 'MODELISA_24')
                end if
            end if
        end do
        if ((iok1+iok2) .ne. 1) then
            call utmess('F', 'MODELISA_25')
        end if
    end do
!
!     STOCKAGE DANS LA CARTE DES PARAMETRES GEOMETRIQUES CALCULES
!     PAR ANGCOU ET DE LA GENERATRICE CONTINUE SUR LES TUYAUX
!
    nbdroi = 0
    nbcoud = 0
!
    ifm = ivr(4)
    do ipa = 1, nbpart
        if (ivr(3) .eq. 2) write (ifm, 100) ipa, nbmap(ipa)
        izk = isens(ipa)
        if (izk .eq. 0) then
!
!        PAS DE VECTEUR FOURNI : ON EN CREE UN
!
            no1 = nopar(ipa, 1, 1)
            no2 = nopar(ipa, 2, 1)
            vx = coor(3*no2-3+1)-coor(3*no1-3+1)
            vy = coor(3*no2-3+2)-coor(3*no1-3+2)
            vz = coor(3*no2-3+3)-coor(3*no1-3+3)
            zk1(1) = -vy
            zk1(2) = vx
            zk1(3) = 0.d0
            call normev(zk1, norme)
            if (norme .lt. 1.d-4) then
                zk1(1) = 0.d0
                zk1(2) = -vz
                zk1(3) = vy
                call normev(zk1, norme)
            end if
        else
            do i = 1, 3
                zk1(i) = cozk(3*(abs(izk)-1)+i)
            end do
            call normev(zk1, norme)
        end if
        iavant = -1
        do im0 = 1, nbmap(ipa)
            if (izk .ge. 0) then
                ima = im0
            else
                ima = nbmap(ipa)-im0+1
            end if
            no1 = nopar(ipa, 1, ima)
            no2 = nopar(ipa, 2, ima)
            no3 = nopar(ipa, 3, ima)
            if (nno .eq. 4) then
                no4 = nopar(ipa, 4, ima)
            end if
            nummai = elpar(ipa, ima)
            do i = 1, 3
                coor3(i) = coor(3*no1-3+i)
                coor3(3+i) = coor(3*no2-3+i)
                coor3(6+i) = coor(3*no3-3+i)
                if (nno .eq. 4) then
                    coor3(9+i) = coor(3*no4-3+i)
                end if
            end do
            if (nno .eq. 3) then
                call angcou(coor3, zk1, izk, icoude, zk2, &
                            rayon, theta, angl1, angl2, angl3, &
                            pgl1, pgl2, pgl3, omega, dn1n2, &
                            epsi, crit, zk3)
                do i = 1, 3
                    angl4(i) = 0.d0
                end do
            else if (nno .eq. 4) then
                call angco4(coor3, zk1, izk, icoude, zk2, &
                            rayon, theta, angl1, angl2, angl3, &
                            angl4, pgl1, pgl2, pgl3, pgl4, &
                            omega, dn1n2, epsi, crit)
            end if
            do i = 1, 3
                zr(jdvlvo-1+i) = angl1(i)
                zr(jdvlvo-1+3+i) = angl2(i)
                zr(jdvlvo-1+6+i) = angl3(i)
            end do
            icmp = 9
            nbcmp = 14
            if (nno .eq. 4) then
                do i = 1, 3
                    zr(jdvlvo-1+9+i) = angl4(i)
                end do
                icmp = 12
                nbcmp = 17
            end if
!
            if (icoude .eq. 0) then
                nbdroi = nbdroi+1
            else
                nbcoud = nbcoud+1
            end if
!
!              MODI_METRIQUE
!
            icoud2 = 0
            if (nmmt(nummai) .eq. 0) then
                icoud2 = icoude+10
            else if (nmmt(nummai) .eq. 1) then
                icoud2 = icoude
            else
                nommai = int_to_char8(nummai)
                call utmess('F', 'MODELISA_26', sk=nommai)
            end if
!
            zr(jdvlvo-1+icmp+1) = icoud2
            zr(jdvlvo-1+icmp+2) = dn1n2
            zr(jdvlvo-1+icmp+3) = rayon
            zr(jdvlvo-1+icmp+4) = theta
            zr(jdvlvo-1+icmp+5) = omega
!
            call nocart(cartor, 3, nbcmp, mode='NUM', nma=1, limanu=[nummai])
!
            if (ivr(3) .eq. 2) then
                nommai = int_to_char8(nummai)
                nomno1 = int_to_char8(no1)
                nomno2 = int_to_char8(no2)
                nomno3 = int_to_char8(no3)
                if (nno .eq. 4) then
                    nomno4 = int_to_char8(no4)
                else
                    nomno4 = ' '
                end if
                if (icoude .ne. iavant) then
                    if (nno .eq. 3) then
                        if (icoude .eq. 0) then
                            write (ifm, 131)
                        else
                            write (ifm, 132)
                        end if
                    else if (nno .eq. 4) then
                        if (icoude .eq. 0) then
                            write (ifm, 131)
                        else
                            write (ifm, 132)
                        end if
                    end if
                    iavant = icoude
                end if
                if (izk .ge. 0) then
                    if (icoude .eq. 0) then
                        write (ifm, 110) nommai, nomno1, nomno2, nomno3, &
                            nomno4, (zk1(i), i=1, 3), (angl1(i), i=1, 3)
                    else
                        write (ifm, 111) nommai, nomno1, nomno2, nomno3, &
                            nomno4, (zk1(i), i=1, 3), (zk2(i), i=1, 3), rayon, &
                            theta, omega, (angl1(i), i=1, 3), (angl2(i), i=1, 3) &
                            , (angl3(i), i=1, 3)
                    end if
                else
                    if (icoude .eq. 0) then
                        write (ifm, 110) nommai, nomno1, nomno2, nomno3, &
                            nomno4, (zk1(i), i=1, 3), (angl1(i), i=1, 3)
                    else
                        write (ifm, 111) nommai, nomno1, nomno2, nomno3, &
                            nomno4, (zk2(i), i=1, 3), (zk1(i), i=1, 3), rayon, &
                            theta, omega, (angl1(i), i=1, 3), (angl2(i), i=1, 3) &
                            , (angl3(i), i=1, 3)
                    end if
                end if
            end if
            do i = 1, 3
                zk1(i) = zk2(i)
            end do
        end do
    end do
!
    write (ifm, 155) nbdroi
    write (ifm, 156) nbcoud
!
131 format(&
     &3x, 'MAILLE  NOEUD1  NOEUD2  NOEUD3  NOEUD4   TYPE   ',&
     &'Z1_X', 8x, 'Z1_Y', 8x, 'Z1_Z', 8x, 'ALPHA1', 6x, 'BETA1', 7x, 'GAMMA1')
132 format(&
     &3x, 'MAILLE  NOEUD1  NOEUD2  NOEUD3  NOEUD4   TYPE   ',&
     &'Z1_X', 8x, 'Z1_Y', 8x, 'Z1_Z', 8x, 'Z2_X', 8x, 'Z2_Y', 8x, 'Z2_Z',&
     &8x, 'RAYON', 7x, 'ANGLE', 7x, 'OMEGA',&
     &7x, 'ALPHA1', 6x, 'BETA1', 7x, 'GAMMA1',&
     &6x, 'ALPHA2', 6x, 'BETA2', 7x, 'GAMMA2',&
     &6x, 'ALPHA3', 6x, 'BETA3', 7x, 'GAMMA3')
!
100 format(3x, 'TUYAUTERIE NUMERO : ', i6, ' NOMBRE DE MAILLES : ', i6)
!
110 format(3x, 5a8, 1x, 'DROIT', 1x, 9(d11.4, 1x))
111 format(3x, 5a8, 1x, 'COUDE', 1x, 18(d11.4, 1x))
!
155 format(3x, 'NOMBRE TOTAL D ELEMENTS TUYAU DROITS ', 1x, i6)
156 format(3x, 'NOMBRE TOTAL D ELEMENTS TUYAU COUDES ', 1x, i6)
!
    call jedetr(tmpnor)
    call jedetr(tmpvor)
!
    call jedema()
end subroutine
