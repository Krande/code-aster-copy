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
subroutine raycir(jvecpg, jdtau, jvecn, nbordr, nbvec, &
                  nommet)
! person_in_charge: van-xuan.tran at edf.fr
! aslint: disable=W1501
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cer3pt.h"
#include "asterfort/dimax1.h"
#include "asterfort/dimax2.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: jvecpg, jdtau, jvecn, nbordr, nbvec
    character(len=16) :: nommet
! ---------------------------------------------------------------------
! BUT: DETERMINER LE PLUS PETIT CERCLE CIRCONSCRIT AUX POINTS
!      REPRESANTANT LE VECTEUR DE CISAILLEMENT TAU DANS LE PLAN u, v.
! ---------------------------------------------------------------------
! ARGUMENTS:
! JVECPG     IN    I : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                      LES COMPOSANTES u ET v DU VECTEUR TAU
!                      (CISAILLEMENT), POUR TOUS LES NUMEROS
!                      D'ORDRE.
! JDTAU      IN    I : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                      LES VALEURS DE DELTA_TAU_MAX POUR CHAQUE VECTEUR.
! JVECN      IN    I : ADRESSE DU VECTEUR DE TRAVAIL CONTENANT
!                      LA VALEUR DU POINTEUR PERMETTANT D'ACCEDER AU
!                      VECTEUR NORMAL ASSOCIE A DELTA_TAU_MAX.
! NBORDR     IN    I : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                      STRUCTURE DE DONNEES RESULTAT.
! NBVEC      IN    I : NOMBRE DE VECTEURS NORMAUX.
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: i, ivect, iordr, k, n1
    integer(kind=8) :: jsec1, jsec2, jsec3, jsec4, jdom1
    integer(kind=8) :: jdom2, jcoorp
    integer(kind=8) :: nbpts1, nbpts2, nbpts3, nbpts4, nbptd1, nbptd2
    integer(kind=8) :: indsec, n, ireth, iretv, nboucl, nbr, iret3p
!
    real(kind=8) :: cumin, cumax, cvmin, cvmax, cui, cvi, diamin
    real(kind=8) :: raymin, cuomi1, cvomi1, cuomi2, cvomi2, cuo1, cvo1
    real(kind=8) :: disto1, dists(4), distd1, distd2, dmaxi(6)
    real(kind=8) :: coorpt(24), dmax, cuppe1, cvppe1, cuppe2, cvppe2
    real(kind=8) :: dun, dvn, cuon, cvon, dsegn, rsegn, cupn, cvpn
    real(kind=8) :: cupn0, cvpn0, cupn1, cvpn1, cupn2, cvpn2, ray3pt
    real(kind=8) :: etir, rayon, dist, cutau, cvtau, p, epsilo, x
    real(kind=8) :: epsil1, hypot, vtest0, vtest1, cetir, cuoi, cvoi
    real(kind=8) :: rmin3p
!
    character(len=5) :: oricad
    real(kind=8), pointer :: vcer3pt(:) => null()
!
!-----------------------------------------------------------------------
!234567                                                              012
!
    call jemarq()
!
    if (nommet(1:15) .eq. 'CERCLE_APPROCHE') then
        goto 777
    end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 --------------------------
!                 |    PREMIERE METHODE    |
!                 --------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
    cetir = sqrt(3.0d0/4.0d0)
!
    epsilo = 1.0d-6
    epsil1 = 1.0d-4
!
    call wkvect('&&RAYCIR.SECT1', 'V V R', nbordr*2, jsec1)
    call wkvect('&&RAYCIR.SECT2', 'V V R', nbordr*2, jsec2)
    call wkvect('&&RAYCIR.SECT3', 'V V R', nbordr*2, jsec3)
    call wkvect('&&RAYCIR.SECT4', 'V V R', nbordr*2, jsec4)
    call wkvect('&&RAYCIR.DOM1', 'V V R', nbordr*2, jdom1)
    call wkvect('&&RAYCIR.DOM2', 'V V R', nbordr*2, jdom2)
!
    call wkvect('&&RAYCIR.COORP', 'V V R', 12, jcoorp)
    AS_ALLOCATE(vr=vcer3pt, size=6)
!
! ININTIALISATION
!
    n1 = 0
!
    do ivect = 1, nbvec
        cumin = r8maem()
        cumax = -r8maem()
        cvmin = r8maem()
        cvmax = -r8maem()
        call jerazo('&&RAYCIR.SECT1', nbordr*2, 1)
        call jerazo('&&RAYCIR.SECT2', nbordr*2, 1)
        call jerazo('&&RAYCIR.SECT3', nbordr*2, 1)
        call jerazo('&&RAYCIR.SECT4', nbordr*2, 1)
        call jerazo('&&RAYCIR.DOM1', nbordr*2, 1)
        call jerazo('&&RAYCIR.DOM2', nbordr*2, 1)
!
        do iordr = 1, nbordr
            n1 = n1+1
            cui = zr(jvecpg+(n1-1)*2)
            cvi = zr(jvecpg+(n1-1)*2+1)
!
            !
            if (cui .lt. cumin) then
                cumin = cui
            end if
            if (cui .gt. cumax) then
                cumax = cui
            end if
            if (cvi .lt. cvmin) then
                cvmin = cvi
            end if
            if (cvi .gt. cvmax) then
                cvmax = cvi
            end if
        end do
!
!-----------------------------------------------------------------------
!   ------------------------------------
!  |  TRAITEMENT DES CAS PARTICULIERS  |
!  ------------------------------------
!
        hypot = sqrt((cumax-cumin)**2+(cvmax-cvmin)**2)
        vtest0 = sqrt((cumin+((cumax-cumin)/2.0d0))**2+(cvmin+((cvmax-cvmin)/2.0d0))**2)
        if (vtest0 .gt. epsilo) then
            vtest1 = hypot/vtest0
        else
            vtest1 = 1.0d0
        end if
!
! 1.1 CAS OU TOUS LES POINTS SONT DANS UNE BOITE DONT LE CENTRE EST
!     VOISIN DE ZERO ET DONT LA NORME EST INFERIEURE A EPSILO.
!
        if (hypot .lt. epsilo) then
            zr(jdtau+(ivect-1)) = 0.0d0
            zi(jvecn+(ivect-1)) = ivect
            goto 30
!
! 1.2 CAS OU L HYPOTENUSE DE LA BOITE EST PETITE DEVANT LES
!     COORDONNEES EXTREMES DE LA BOITE => PETITE VALEUR DU CISAILLEMENT.
!
        else if (vtest1 .lt. epsil1) then
            zr(jdtau+(ivect-1)) = hypot/2.0d0
            zi(jvecn+(ivect-1)) = ivect
            goto 30
!
! 1.3 CAS OU TOUS LES POINTS SONT ALIGNES HORIZONTALEMENT
!
        else if (((cvmax-cvmin)/hypot) .lt. epsil1) then
            zr(jdtau+(ivect-1)) = hypot/2.0d0
            zi(jvecn+(ivect-1)) = ivect
            goto 30
!
! 1.4 CAS OU TOUS LES POINTS SONT ALIGNES VERTICALEMENT
!
        else if (((cumax-cumin)/hypot) .lt. epsil1) then
            zr(jdtau+(ivect-1)) = hypot/2.0d0
            zi(jvecn+(ivect-1)) = ivect
            goto 30
        end if
!
!-----------------------------------------------------------------------
!
! NOUS FAISONS UNE CORRECTION BARYCENTRIQUE
!
        cuo1 = cumin+(cumax-cumin)/2.0d0
        cvo1 = cvmin+(cvmax-cvmin)/2.0d0
!
        if ((cumax-cumin) .ge. (cvmax-cvmin)) then
            diamin = cumax-cumin
            raymin = (cumax-cumin)/2.0d0
            cuomi1 = 0.0d0
            cvomi1 = -(cvmax-cvmin)/2.0d0
            cuomi2 = 0.0d0
            cvomi2 = (cvmax-cvmin)/2.0d0
            etir = (cvmax-cvmin)/(cumax-cumin)
            oricad = 'HORIZ'
        else
            diamin = cvmax-cvmin
            raymin = (cvmax-cvmin)/2.0d0
            cuomi1 = -(cumax-cumin)/2.0d0
            cvomi1 = 0.0d0
            cuomi2 = (cumax-cumin)/2.0d0
            cvomi2 = 0.0d0
            etir = (cumax-cumin)/(cvmax-cvmin)
            oricad = 'VERTI'
        end if
!
! REDEFINITION DES POINTS EXTREMES DU CADRE APRES LA CORRECTION
! BARYCENTRIQUE (ON A TRANSLAT2 LE CADRE AUTOUR DE ZERO)
!
        cumin = -(cumax-cumin)/2.0d0
        cvmin = -(cvmax-cvmin)/2.0d0
        cumax = abs(cumin)
        cvmax = abs(cvmin)
!
! DETERMINATION DES POINTS SITUES DANS LES 4 SECTEURS ET LES 2 DOMAINES
!
        nbpts1 = 0
        nbpts2 = 0
        nbpts3 = 0
        nbpts4 = 0
        nbptd1 = 0
        nbptd2 = 0
        n1 = n1-nbordr
        do iordr = 1, nbordr
            n1 = n1+1
            cui = zr(jvecpg+(n1-1)*2)-cuo1
            cvi = zr(jvecpg+(n1-1)*2+1)-cvo1
!
            disto1 = sqrt((0.0d0-cui)**2+(0.0d0-cvi)**2)
            dists(1) = sqrt((cumin-cui)**2+(cvmax-cvi)**2)
            dists(2) = sqrt((cumax-cui)**2+(cvmax-cvi)**2)
            dists(3) = sqrt((cumax-cui)**2+(cvmin-cvi)**2)
            dists(4) = sqrt((cumin-cui)**2+(cvmin-cvi)**2)
            distd1 = sqrt((cuomi1-cui)**2+(cvomi1-cvi)**2)
            distd2 = sqrt((cuomi2-cui)**2+(cvomi2-cvi)**2)
!
            indsec = 0
            do i = 1, 4
                if ((dists(i) .gt. diamin) .and. (indsec .eq. 0)) then
                    if (cui .ge. 0.0d0) then
                        if (cvi .ge. 0.0d0) then
                            zr(jsec2+nbpts2*2) = cui
                            zr(jsec2+nbpts2*2+1) = cvi
                            nbpts2 = nbpts2+1
                            indsec = 1
                        else
                            zr(jsec3+nbpts3*2) = cui
                            zr(jsec3+nbpts3*2+1) = cvi
                            nbpts3 = nbpts3+1
                            indsec = 1
                        end if
                    else
                        if (cvi .ge. 0.0d0) then
                            zr(jsec1+nbpts1*2) = cui
                            zr(jsec1+nbpts1*2+1) = cvi
                            nbpts1 = nbpts1+1
                            indsec = 1
                        else
                            zr(jsec4+nbpts4*2) = cui
                            zr(jsec4+nbpts4*2+1) = cvi
                            nbpts4 = nbpts4+1
                            indsec = 1
                        end if
                    end if
                end if
            end do
!
            if ((distd1 .gt. raymin) .or. (disto1 .gt. raymin)) then
                zr(jdom1+nbptd1*2) = cui
                zr(jdom1+nbptd1*2+1) = cvi
                nbptd1 = nbptd1+1
            end if
!
            if ((distd2 .gt. raymin) .or. (disto1 .gt. raymin)) then
                zr(jdom2+nbptd2*2) = cui
                zr(jdom2+nbptd2*2+1) = cvi
                nbptd2 = nbptd2+1
            end if
        end do
!
! RECHERCHE DES 2 POINTS LES PLUS ELOIGNES PARMI LES POINTS DES
! SECTEURS 1, 2, 3 ET 4.
!
! EXEMPLE : CAS OU LE CADRE EST ORIENTE HORIZONTALEMENT
!
!       -----------------------------------
!       | Sect. 1                 Sect. 2 |
!       |                                 |
!       |                                 |
!       |                                 |
!       |                                 |
!       | Sect. 4                 Sect. 3 |
!       -----------------------------------
!
        do i = 1, 6
            dmaxi(i) = 0.0d0
        end do
        do i = 1, 24
            coorpt(i) = 0.0d0
        end do
        if ((nbpts1 .gt. 0) .and. (nbpts2 .gt. 0)) then
            if ((oricad .eq. 'HORIZ') .or. ((oricad .eq. 'VERTI') .and. (etir .ge. cetir))) then
                call dimax1(jsec1, jsec2, nbpts1, nbpts2, dmaxi(1), &
                            coorpt(1), coorpt(2), coorpt(3), coorpt(4))
            else
                dmaxi(1) = 0.0d0
            end if
        end if
!
        if ((nbpts1 .gt. 0) .and. (nbpts3 .gt. 0)) then
            call dimax1(jsec1, jsec3, nbpts1, nbpts3, dmaxi(2), &
                        coorpt(5), coorpt(6), coorpt(7), coorpt(8))
        else
            dmaxi(2) = 0.0d0
        end if
!
        if ((nbpts1 .gt. 0) .and. (nbpts4 .gt. 0)) then
            if ((oricad .eq. 'VERTI') .or. ((oricad .eq. 'HORIZ') .and. (etir .ge. cetir))) then
                call dimax1(jsec1, jsec4, nbpts1, nbpts4, dmaxi(3), &
                            coorpt(9), coorpt(10), coorpt(11), coorpt(12))
            else
                dmaxi(3) = 0.0d0
            end if
        end if
!
        if ((nbpts2 .gt. 0) .and. (nbpts3 .gt. 0)) then
            if ((oricad .eq. 'VERTI') .or. ((oricad .eq. 'HORIZ') .and. (etir .ge. cetir))) then
                call dimax1(jsec2, jsec3, nbpts2, nbpts3, dmaxi(4), &
                            coorpt(13), coorpt(14), coorpt(15), coorpt(16))
            else
                dmaxi(4) = 0.0d0
            end if
        end if
!
        if ((nbpts2 .gt. 0) .and. (nbpts4 .gt. 0)) then
            call dimax1(jsec2, jsec4, nbpts2, nbpts4, dmaxi(5), &
                        coorpt(17), coorpt(18), coorpt(19), coorpt(20))
        else
            dmaxi(5) = 0.0d0
        end if
!
        if ((nbpts3 .gt. 0) .and. (nbpts4 .gt. 0)) then
            if ((oricad .eq. 'HORIZ') .or. ((oricad .eq. 'VERTI') .and. (etir .ge. cetir))) then
                call dimax1(jsec3, jsec4, nbpts3, nbpts4, dmaxi(6), &
                            coorpt(21), coorpt(22), coorpt(23), coorpt(24))
            else
                dmaxi(6) = 0.0d0
            end if
        end if
!
        dmax = 0.0d0
        n = 0
        do i = 1, 6
            if (dmaxi(i) .gt. dmax) then
                dmax = dmaxi(i)
                n = i
            end if
        end do
        cuppe1 = coorpt((n-1)*4+1)
        cvppe1 = coorpt((n-1)*4+2)
        cuppe2 = coorpt((n-1)*4+3)
        cvppe2 = coorpt((n-1)*4+4)
!
! CALCUL DU CENTRE DU SEGMENT LE PLUS LONG ET DE SA LONGUEUR.
!
!
        dun = abs(cuppe1-cuppe2)/2.0d0
        dvn = abs(cvppe1-cvppe2)/2.0d0
        if (cuppe1 .lt. cuppe2) then
            cuon = cuppe1+dun
        else
            cuon = cuppe2+dun
        end if
        if (cvppe1 .lt. cvppe2) then
            cvon = cvppe1+dvn
        else
            cvon = cvppe2+dvn
        end if
!
        dsegn = sqrt((cuppe1-cuppe2)**2+(cvppe1-cvppe2)**2)
        rsegn = dsegn/2.0d0
!
! ON CHERCHE SI IL EXISTE DES POINTS Pi TELS QUE LEURS DISTANCES AU
! CENTRE 'On' SOIENT SUPERIEURES AU RAYON 'RSEGN'. SI OUI ON PREND LE
! PLUS ELOIGNE.
! SUIVANT LA POSITION DE 'On' PAR RAPPORT AU CENTRE DE LA BOITE O1,
! LA RECHERCHE EST FAITE, SOIT DANS LE DOMAINE1, SOIT DANS LE DOMAINE2.
!
        ireth = 0
        iretv = 0
!
        if (oricad(1:5) .eq. 'HORIZ') then
            if (cvon .lt. 0.0d0) then
                call dimax2(jdom1, nbptd1, cuon, cvon, rsegn, &
                            cupn, cvpn, ireth)
            else
                call dimax2(jdom2, nbptd2, cuon, cvon, rsegn, &
                            cupn, cvpn, ireth)
            end if
        else if (oricad(1:5) .eq. 'VERTI') then
            if (cuon .lt. 0.0d0) then
                call dimax2(jdom1, nbptd1, cuon, cvon, rsegn, &
                            cupn, cvpn, iretv)
            else
                call dimax2(jdom2, nbptd2, cuon, cvon, rsegn, &
                            cupn, cvpn, iretv)
            end if
        end if
!
! ON CHERCHE LE PLUS PETIT CERCLE CIRCONSCRIT PAR LA METHODE DU CERCLE
! PASSANT PAR TROIS POINTS.
!
        nboucl = 0
!
100     continue
!
        nboucl = nboucl+1
!
        if (ireth .eq. 1) then
            if (nboucl .eq. 1) then
                cupn0 = cuppe1
                cvpn0 = cvppe1
                cupn1 = cuppe2
                cvpn1 = cvppe2
                cupn2 = cupn
                cvpn2 = cvpn
                call cer3pt(cupn0, cvpn0, cupn1, cvpn1, cupn2, &
                            cvpn2, cuon, cvon, ray3pt)
!
                if (cvon .lt. 0.0d0) then
                    call dimax2(jdom1, nbptd1, cuon, cvon, ray3pt, &
                                cupn, cvpn, ireth)
                else
                    call dimax2(jdom2, nbptd2, cuon, cvon, ray3pt, &
                                cupn, cvpn, ireth)
                end if
                goto 100
!
            else if (nboucl .gt. 1) then
                zr(jcoorp) = cupn0
                zr(jcoorp+1) = cvpn0
                zr(jcoorp+2) = cupn1
                zr(jcoorp+3) = cvpn1
                zr(jcoorp+4) = cupn2
                zr(jcoorp+5) = cvpn2
                zr(jcoorp+6) = cupn
                zr(jcoorp+7) = cvpn
                zr(jcoorp+8) = cupn0
                zr(jcoorp+9) = cvpn0
                zr(jcoorp+10) = cupn1
                zr(jcoorp+11) = cvpn1
!
                ray3pt = r8maem()
                k = 0
!
                do i = 1, 3
                    cupn0 = zr(jcoorp+i*2)
                    cvpn0 = zr(jcoorp+i*2+1)
                    cupn1 = zr(jcoorp+i*2+2)
                    cvpn1 = zr(jcoorp+i*2+3)
                    cupn2 = zr(jcoorp+i*2+4)
                    cvpn2 = zr(jcoorp+i*2+5)
                    call cer3pt(cupn0, cvpn0, cupn1, cvpn1, cupn2, &
                                cvpn2, cuoi, cvoi, rmin3p)
!
                    call dimax2(jcoorp, 4, cuoi, cvoi, rmin3p, &
                                cupn, cvpn, iret3p)
!
                    if (iret3p .eq. 0) then
                        k = k+1
                    end if
!
                    if ((rmin3p .lt. ray3pt) .and. (iret3p .eq. 0)) then
                        ray3pt = rmin3p
                        cuon = cuoi
                        cvon = cvoi
                        vcer3pt(1) = cupn0
                        vcer3pt(1+1) = cvpn0
                        vcer3pt(1+2) = cupn1
                        vcer3pt(1+3) = cvpn1
                        vcer3pt(1+4) = cupn2
                        vcer3pt(1+5) = cvpn2
                    end if
!
                end do
                if (k .eq. 0) then
                    call utmess('F', 'PREPOST4_59')
                end if
!
                if (cvon .lt. 0.0d0) then
                    call dimax2(jdom1, nbptd1, cuon, cvon, ray3pt, &
                                cupn, cvpn, ireth)
                else
                    call dimax2(jdom2, nbptd2, cuon, cvon, ray3pt, &
                                cupn, cvpn, ireth)
                end if
!
                if (ireth .eq. 1) then
                    cupn0 = vcer3pt(1)
                    cvpn0 = vcer3pt(1+1)
                    cupn1 = vcer3pt(1+2)
                    cvpn1 = vcer3pt(1+3)
                    cupn2 = vcer3pt(1+4)
                    cvpn2 = vcer3pt(1+5)
                end if
!
                goto 100
            end if
!
        else if (iretv .eq. 1) then
            if (nboucl .eq. 1) then
                cupn0 = cuppe1
                cvpn0 = cvppe1
                cupn1 = cuppe2
                cvpn1 = cvppe2
                cupn2 = cupn
                cvpn2 = cvpn
                write (6, *) 'raycir 521', cupn0, cvpn0, cupn1, &
                    cvpn1, cupn2, cvpn2
!
                call cer3pt(cupn0, cvpn0, cupn1, cvpn1, cupn2, &
                            cvpn2, cuon, cvon, ray3pt)
!
                if (cuon .lt. 0.0d0) then
                    call dimax2(jdom1, nbptd1, cuon, cvon, ray3pt, &
                                cupn, cvpn, iretv)
                else
                    call dimax2(jdom2, nbptd2, cuon, cvon, ray3pt, &
                                cupn, cvpn, iretv)
                end if
!
                goto 100
!
            else if (nboucl .gt. 1) then
                zr(jcoorp) = cupn0
                zr(jcoorp+1) = cvpn0
                zr(jcoorp+2) = cupn1
                zr(jcoorp+3) = cvpn1
                zr(jcoorp+4) = cupn2
                zr(jcoorp+5) = cvpn2
                zr(jcoorp+6) = cupn
                zr(jcoorp+7) = cvpn
                zr(jcoorp+8) = cupn0
                zr(jcoorp+9) = cvpn0
                zr(jcoorp+10) = cupn1
                zr(jcoorp+11) = cvpn1
!
                ray3pt = r8maem()
                k = 0
!
                do i = 1, 3
                    cupn0 = zr(jcoorp+i*2)
                    cvpn0 = zr(jcoorp+i*2+1)
                    cupn1 = zr(jcoorp+i*2+2)
                    cvpn1 = zr(jcoorp+i*2+3)
                    cupn2 = zr(jcoorp+i*2+4)
                    cvpn2 = zr(jcoorp+i*2+5)
                    call cer3pt(cupn0, cvpn0, cupn1, cvpn1, cupn2, &
                                cvpn2, cuoi, cvoi, rmin3p)
!
                    call dimax2(jcoorp, 4, cuoi, cvoi, rmin3p, &
                                cupn, cvpn, iret3p)
!
                    if (iret3p .eq. 0) then
                        k = k+1
                    end if
!
                    if ((rmin3p .lt. ray3pt) .and. (iret3p .eq. 0)) then
                        ray3pt = rmin3p
                        cuon = cuoi
                        cvon = cvoi
                        vcer3pt(1) = cupn0
                        vcer3pt(1+1) = cvpn0
                        vcer3pt(1+2) = cupn1
                        vcer3pt(1+3) = cvpn1
                        vcer3pt(1+4) = cupn2
                        vcer3pt(1+5) = cvpn2
                    end if
!
                end do
                if (k .eq. 0) then
                    call utmess('F', 'PREPOST4_59')
                end if
!
                if (cuon .lt. 0.0d0) then
                    call dimax2(jdom1, nbptd1, cuon, cvon, ray3pt, &
                                cupn, cvpn, iretv)
                else
                    call dimax2(jdom2, nbptd2, cuon, cvon, ray3pt, &
                                cupn, cvpn, iretv)
                end if
!
                if (iretv .eq. 1) then
                    cupn0 = vcer3pt(1)
                    cvpn0 = vcer3pt(1+1)
                    cupn1 = vcer3pt(1+2)
                    cvpn1 = vcer3pt(1+3)
                    cupn2 = vcer3pt(1+4)
                    cvpn2 = vcer3pt(1+5)
                end if
!
                goto 100
            end if
!
        else
            if (nboucl .eq. 1) then
                zr(jdtau+(ivect-1)) = rsegn
                zi(jvecn+(ivect-1)) = ivect
            else if (nboucl .gt. 1) then
                zr(jdtau+(ivect-1)) = ray3pt
                zi(jvecn+(ivect-1)) = ivect
            end if
        end if
!
30      continue
    end do
!
    goto 999
!
777 continue
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 --------------------------
!                 |    DEUXIEME METHODE    |
!                 --------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!
    epsilo = 1.0d-3
    x = 5.0d-2
!
! INITIALISATION
!
    do ivect = 1, nbvec
        cuo1 = 0.0d0
        cvo1 = 0.0d0
        do iordr = 1, nbordr
            cuo1 = cuo1+zr(jvecpg+(iordr-1)*2+(ivect-1)*nbordr*2)
            cvo1 = cvo1+zr(jvecpg+(iordr-1)*2+(ivect-1)*nbordr*2+1)
        end do
!
        cuo1 = cuo1/nbordr
        cvo1 = cvo1/nbordr
        rayon = 0.0d0
!
! CALCUL RECURRENT
!
        cuon = cuo1
        cvon = cuo1
        n = 0
        nbr = 0
!
300     continue
!
        n = n+1
        if (n .gt. nbordr) then
            n = n-nbordr
        end if
!
        cutau = zr(jvecpg+(n-1)*2+(ivect-1)*nbordr*2)
        cvtau = zr(jvecpg+(n-1)*2+(ivect-1)*nbordr*2+1)
!
!
        dist = sqrt((cutau-cuon)**2+(cvtau-cvon)**2)
        p = dist-rayon
!
        if (p .gt. epsilo) then
            nbr = 0
            rayon = rayon+x*p
            cuon = cutau+rayon*((cutau-cuon)/sqrt((cutau-cuon)**2))
            cvon = cvtau+rayon*((cvtau-cvon)/sqrt((cvtau-cvon)**2))
        else
            nbr = nbr+1
        end if
!
        if (nbr .lt. nbordr) then
            goto 300
        else
            zr(jdtau+(ivect-1)) = rayon
            zi(jvecn+(ivect-1)) = ivect
        end if
!
    end do
!
!
999 continue
!
! MENAGE
!
    call jedetr('&&RAYCIR.SECT1')
    call jedetr('&&RAYCIR.SECT2')
    call jedetr('&&RAYCIR.SECT3')
    call jedetr('&&RAYCIR.SECT4')
    call jedetr('&&RAYCIR.DOM1')
    call jedetr('&&RAYCIR.DOM2')
    call jedetr('&&RAYCIR.COORP')
    AS_DEALLOCATE(vr=vcer3pt)
!
    call jedema()
end subroutine
