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

subroutine mdconf(typflu, base, noma, nbm, lnoe, &
                  nuor, iimpr, indic, veci1, vecr1, &
                  vecr2, vecr3, vecr4, vecr5)
    implicit none
! DESCRIPTION : CALCUL DES PARAMETRES DE COUPLAGE FLUIDE-STRUCTURE POUR
! -----------   DIFFERENTES CONFIGURATION
!               REMPLISSAGE DES OBJETS DE TRAVAIL DEPENDANT DU TYPE
!               DE CONFIGURATION
!
!               APPELANTS : FLUST1, FLUST2, MDITMI
!
!  IN : TYPFLU : NOM DU CONCEPT DE TYPE TYPE_FLUI_STRU DEFINISSANT LA
!                CONFIGURATION ETUDIEE
!  IN : BASE   : NOM DU CONCEPT DE TYPE MODE_MECA DEFINISSANT LA BASE
!                MODALE DU SYSTEME AVANT PRISE EN COMPTE DU COUPLAGE
!  IN : NOMA   : NOM DU CONCEPT DE TYPE MAILLAGE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : LNOE   : NOMBRE DE NOEUDS DU TUBE OU DU MAILLAGE
!  IN : IIMPR  : PARAMETRE D'IMPRESSION (IMPREESION FICHIER
!  OUT: VECI1  : VECTEUR ENTIER VARIABLE SELON LE TYPE DE CONFIGURATION
!  OUT: VECR1  : VECTEUR REEL VARIABLE SELON LE TYPE DE CONFIGURATION
!  OUT: VECR2  : VECTEUR REEL VARIABLE SELON LE TYPE DE CONFIGURATION
!  OUT: VECR3  : VECTEUR REEL VARIABLE SELON LE TYPE DE CONFIGURATION
!  OUT: VECR4  : VECTEUR REEL VARIABLE SELON LE TYPE DE CONFIGURATION
!  OUT: VECR5  : VECTEUR REEL VARIABLE SELON LE TYPE DE CONFIGURATION
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! -------------------------
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/axdipo.h"
#include "asterfort/deelpo.h"
#include "asterfort/dismoi.h"
#include "asterfort/exmano.h"
#include "asterfort/extmod.h"
#include "asterfort/extmod_sorted.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mtdscr.h"
#include "asterfort/recude.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/permnoe.h"
#include "asterfort/char8_to_int.h"
!
    character(len=24) :: valk(4)
!
! ARGUMENTS
! ---------
    character(len=8) :: typflu, base, noma
    integer(kind=8) :: nbm, lnoe, nuor(*), iimpr, indic, veci1(*)
    integer(kind=8) :: vali(2)
    real(kind=8) :: vecr1(*), vecr2(*), vecr3(*), vecr4(*), vecr5(*)
!
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: lp
    character(len=8) :: depla(3), config(4), depl, k8b, mailla, modele, nomno0
    character(len=14) :: numddl
    character(len=19) :: caelem, masse
    character(len=24) :: deeq, frhoe, frhoi, fsic, fsvi, fsvk, fsvr, mlgnma
    character(len=24) :: nomcha, pvite
!
! FONCTIONS EXTERNES
! ------------------
    real(kind=8) :: valr(3)
!     EXTERNAL      R8PI, R8PREM
!     EXTERNAL      JEXNOM
!
! DATA
! ----
!-----------------------------------------------------------------------
    integer(kind=8) :: iaxe, icoupl, ide, ideeq, idep, ifsic
    integer(kind=8) :: ifsvi, ifsvk, igrap, ik, imail, imod, ipas
    integer(kind=8) :: ipm, ipv, ireszo, iret, irhoe, irhoi
    integer(kind=8) :: irota1, irota2, itran1, itran2, itypfl, ivale, izone
    integer(kind=8) :: j, lfsvk, lfsvr, lmasg, lmasse, n1
    integer(kind=8) :: n2, nbma, nbmano, neq, numno0, numod, nuzo
    integer(kind=8) :: nzex
    real(kind=8) :: aire, alonto, cm1, cm2, cocaj, cokaj, comaj
    real(kind=8) :: comaj1, comaj2, difphi, dryi, drzi, dyi, dzi
    real(kind=8) :: ep, phi1, phi2, phie, prota, ptran, rhof
    real(kind=8) :: tolr, vmoy, vmoyto, x1, x2
!-----------------------------------------------------------------------
    data depla/'DX      ', 'DY      ', 'DZ      '/
    data config/'ASC_CEN ', 'ASC_EXC ', 'DES_CEN ', 'DES_EXC '/
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
    tolr = r8prem()
!
!
! --- 1.TYPE DE CONFIGURATION  ---
! ---   TYPE_FLUI_STRU         ---
!
    fsic = typflu//'           .FSIC'
    call jeveuo(fsic, 'L', ifsic)
    itypfl = zi(ifsic)
    icoupl = zi(ifsic+1)
!
!
! --- 2.CONFIGURATION DE TYPE FAISCEAU DE TUBE SOUS ECOULEMENT  ---
! ---   TRANSVERSE                                              ---
!
    if (itypfl .eq. 1) then
!
        fsvk = typflu//'           .FSVK'
        call jeveuo(fsvk, 'L', ifsvk)
!
        fsvi = typflu//'           .FSVI'
        call jeveuo(fsvi, 'L', ifsvi)
!
! ---    2.1.NOMBRE DE POINTS DE DISCRETISATION --> VARIABLE INDIC ---
!
        indic = lnoe
!
! ---    2.2.DIAMETRE EXTERIEUR DU TUBE  ET EPAISSEUR --> VECTEUR VECR4
! ---        (TEST SI SECTION CONSTANTE)                           ---
!
        caelem = zk8(ifsvk)
        call recude(caelem, phie, ep)
        vecr4(1) = phie
        vecr4(2) = ep
!
! ---    2.3.ABSCISSE CURVILIGNES DE VITESSE    --> VECTEUR VECR5  ---
! ---        PROFIL DE MASSE VOLUMIQUE          --> VECTEUR VECR2  ---
!
! ---    CONCEPT DE TYPE FONCTION DEFINISSANT LE PROFIL DE MASSE
! ---    VOLUMIQUE DU FLUIDE EXTERNE ET INTERNE
!
        frhoe = zk8(ifsvk+3)
        frhoi = zk8(ifsvk+2)
        frhoe = frhoe(1:19)//'.VALE'
        frhoi = frhoi(1:19)//'.VALE'
        call jeveuo(frhoe, 'L', irhoe)
        call jeveuo(frhoi, 'L', irhoi)
!
! ---    ABSCISSE CURVILIGNES
!
        do ik = 1, lnoe
            vecr5(ik) = zr(irhoe+ik-1)
            vecr2(ik) = zr(irhoe+ik+lnoe-1)
            vecr2(ik+lnoe) = zr(irhoi+ik+lnoe-1)
        end do
!
! ---    2.4.CALCUL DES NUMERO DE ZONE     --> VECTEUR VECI1 ---
! ---        CALCUL DU  PROFIL DE VITESSE  --> VECTEUR VECR1 ---
!                                              DE 1 A LNOE
! ---        CALCUL DES VITESSES MOYENNES
!            PAR ZONE                      --> VECTEUR VECR1
!                                              DE (LNOE+1) A (2*LNOE+1)
! ---        CALCUL DE LA VITESSE MOYENNE  --> VECTEUR VECR1(2*LNOE+1)
!
! ---    A PARTIR DES PROFIL DE VITESSE DE CHAQUE ZONE ON CONSTRUIT UN
! ---    PROFIL UNIQUE, AINSI QUE LES VECTEUR IRES (TYPE DE RESEAU) ET
! ---    VMOY (VITESSE MOYENNE PAR ZONE) DE CHAQUE POINT DU TUBE.
! ---    VMOYTO EST LA VITESSE MOYENNE SUR L ENSEMBLE DES ZONES.
! ---    LES PROFILS DE VITESSE NE SONT PAS NORMES.
!
        vmoyto = 0.d0
        alonto = 0.d0
        nzex = zi(ifsvi+1)
        call wkvect('&&MDCONF.TEMPO', 'V V I', 2*nzex+1, izone)
        zi(izone-1+1) = nzex
!
! ---    BOUCLE SUR LES ZONES D EXCITATION DU FLUIDE
!
!
        do nuzo = 1, nzex
            pvite = zk8(ifsvk+3+nuzo)
            ireszo = zi(ifsvi+1+nuzo)
            pvite = pvite(1:19)//'.VALE'
            call jeveuo(pvite, 'L', ipv)
! ---       RECHERCHE DES EXTREMITES DE LA ZONE 'NUZO'
            n1 = -1
            n2 = -1
            do ik = 1, lnoe
                if (abs(zr(ipv+lnoe+ik-1)) .gt. tolr) then
                    n1 = ik
                    goto 21
                end if
            end do
21          continue
!
            do ik = lnoe, 1, -1
                if (abs(zr(ipv+lnoe+ik-1)) .gt. tolr) then
                    n2 = ik
                    goto 31
                end if
            end do
31          continue
!
            zi(izone+2*(nuzo-1)-1+2) = n1
            zi(izone+2*(nuzo-1)-1+3) = n2
            if (n1 .eq. n2) then
                call utmess('F', 'ALGELINE_68', sk=zk8(ifsvk+3+nuzo))
            end if
!
            aire = 0.d0
            x1 = zr(ipv+n1-1)
            x2 = zr(ipv+n2-1)
            do ik = n1+1, n2
                aire = aire+( &
                       zr(ipv+lnoe+ik-1)+zr(ipv+lnoe+ik-2))*(zr(ipv+ik-1)-zr(ipv+ik-2) &
                                                             )/2.d0
            end do

!
            vmoy = aire/(x2-x1)
            vmoyto = vmoyto+aire
            alonto = alonto+(x2-x1)
            do ik = n1, n2
                if (veci1(ik) .ne. 0) then
                    call utmess('F', 'ALGELINE_69', sk=zk8(ifsvk+3+nuzo))
                end if
                vecr1(ik+lnoe) = vmoy
                veci1(ik) = ireszo
                vecr1(ik) = zr(ipv+lnoe+ik-1)
            end do

            !
            call jelibe(pvite)
!
! ---    FIN DE BOUCLE SUR LES ZONES D EXCITATION DU FLUIDE
        end do
!
! ---    VITESSE MOYENNE SUR L'ENSEMBLE DU TUBE
        vmoyto = vmoyto/alonto

        vecr1(1+2*lnoe) = vmoyto
!
!
! ---   2.5.DEFORMEES MODALES    --> VECTEUR VECR3  ---
!
        if (icoupl .eq. 1) then
!
! ---      DIRECTION DANS LAQUELLE AGISSENT LES FORCES FLUIDELASTIQUES
            idep = 0
            depl = zk8(ifsvk+1)
            do ide = 1, 3
                if (depla(ide) .eq. depl) idep = ide
            end do
!
! ---       DEFORMEES MODALES
!
            call dismoi('REF_MASS_PREM', base, 'RESU_DYNA', repk=masse)
            call mtdscr(masse)
            call jeveuo(masse//'.&INT', 'L', lmasse)
            call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=numddl)
            call dismoi('NOM_MAILLA', masse, 'MATR_ASSE', repk=mailla)
            call dismoi('NB_EQUA', masse, 'MATR_ASSE', repi=neq)
!
            call extmod_sorted(base, numddl, nuor, nbm, vecr3, &
                               neq, lnoe, [idep], 1)
            call permnoe(mailla, vecr3, nbm, lnoe, 1)
!
        end if
!
! ---   2.6.IMPRESSION  ---
!
        if (iimpr .eq. 1) then
            ipas = zi(ifsvi)
            vali(1) = lnoe
            vali(2) = ipas
            valk(1) = base
            valk(2) = caelem(1:8)
            valr(1) = phie
            call utmess('I', 'ALGELINE5_10', nk=2, valk=valk, ni=2, &
                        vali=vali, sr=valr(1))
            do nuzo = 1, nzex
                valk(1) = zk8(ifsvk+nuzo+3)
                vali(1) = zi(ifsvi+nuzo+1)
                call utmess('I', 'ALGELINE5_11', sk=valk(1), si=vali(1))
            end do
            call utmess('I', 'VIDE_1')
        end if
!
! --- 3.CONFIGURATION DE TYPE "GRAPPE DE COMMANDE"  ---
!
    else if (itypfl .eq. 2) then
!
        if (icoupl .eq. 1) then
!
            fsvk = typflu//'           .FSVK'
            call jeveuo(fsvk, 'L', lfsvk)
!
            fsvr = typflu//'           .FSVR'
            call jeveuo(fsvr, 'L', lfsvr)
!
!------- 3.1.CREATION D'OBJETS DE TRAVAIL
!
            mlgnma = noma//'.TYPMAIL'
            call jelira(mlgnma, 'LONMAX', nbma)
            call wkvect('&&MDCONF.TEMP.MAIL', 'V V I', nbma, imail)
!
!------- RECUPERATION DES DONNEES DANS LE CONCEPT TYPE_FLUI_STRU
!        DEDUCTION DE DONNEES COMPLEMENTAIRES
! ---    3.2.TYPE DE CONFIGURATION GRAPPE --> VARIABLE INDIC ---
!
            do igrap = 1, 4
                if (zk8(lfsvk) .eq. config(igrap)) goto 121
            end do
121         continue
            indic = igrap
            nomno0 = zk8(lfsvk+1)
            numno0 = char8_to_int(nomno0)
            caelem = zk8(lfsvk+2)
            modele = zk8(lfsvk+3)
            cm1 = zr(lfsvr)
            cm2 = cm1/3.d0
            rhof = zr(lfsvr+1)
!
! ---    3.3.DIAMETRE EXTERIEUR DU TUBE         --> VECTEUR VECR4  ---
!
            call exmano(noma, numno0, zi(imail), nbmano)
            if (nbmano .ne. 2) then
                call utmess('F', 'ALGELINE_70')
            end if
!
            call deelpo(caelem(1:8), noma, zi(imail), phi1)
            call deelpo(caelem(1:8), noma, zi(imail+1), phi2)
            difphi = dble(abs(phi1-phi2))
            if (difphi .gt. phi1*tolr) then
                call utmess('F', 'ALGELINE_71')
            else
                phie = phi1
            end if
!
            vecr4(1) = phie
!
!------- 3.4.RECUPERATION DES GRANDEURS GEOMETRIQUES CARACTERISTIQUES --
!        DEDUCTION DE COEFFICIENTS DE DIMENSIONNEMENT                ---
! ---                                            --> VECTEUR VECR2   ---
            lp = phie*986.d0/890.d0
!
            comaj = 0.5d0*rhof*phie*phie*lp
            comaj1 = comaj*cm1
            comaj2 = comaj*cm2*lp*lp
!
            cocaj = -0.5d0*rhof*phie*lp
            vecr2(1) = cocaj
            vecr2(2) = cocaj*lp*lp*0.12849663d0
!
            cokaj = -0.5d0*rhof*lp
            vecr2(3) = cokaj
            vecr2(4) = cokaj*lp*lp*0.12849663d0
!
!------- DETERMINATION DE L'AXE DIRECTEUR DE LA POUTRE
!        DEDUCTION DES DDLS A EXTRAIRE
!
            call axdipo(noma, caelem(1:8), modele, iaxe)
!
            if (iaxe .eq. 1) then
                itran1 = 2
                itran2 = 3
                irota1 = 6
                irota2 = 5
            else if (iaxe .eq. 2) then
                itran1 = 3
                itran2 = 1
                irota1 = 4
                irota2 = 6
            else
                itran1 = 1
                itran2 = 2
                irota1 = 5
                irota2 = 4
            end if
!
!------- 3.5.PONDERATIONS DUES AUX DEFORMEES MODALES
!                                                --> VECTEUR VECR3  ---
!            MASSES MODALES EN EAU               --> VECTEUR VECR1  ---
            call dismoi('REF_MASS_PREM', base, 'RESU_DYNA', repk=masse)
            call mtdscr(masse)
            call jeveuo(masse//'.&INT', 'L', lmasse)
            call dismoi('NOM_NUME_DDL', masse, 'MATR_ASSE', repk=numddl)
            call dismoi('NOM_MAILLA', masse, 'MATR_ASSE', repk=mailla)
            call dismoi('NB_EQUA', masse, 'MATR_ASSE', repi=neq)
            deeq = numddl//'.NUME.DEEQ'
            call jeveuo(deeq, 'L', ideeq)
!
            do imod = 1, nbm
!
                numod = nuor(imod)
!
                call rsadpa(base, 'L', 1, 'MASS_GENE', numod, &
                            0, sjv=lmasg, styp=k8b)
!
                call rsexch('F', base, 'DEPL', numod, nomcha, &
                            iret)
                nomcha = nomcha(1:19)//'.VALE'
                call jeveuo(nomcha, 'L', ivale)
                ipm = 0
                dyi = 0.d0
                dzi = 0.d0
                dryi = 0.d0
                drzi = 0.d0
                do j = 1, neq
                    if (zi(ideeq+(2*j)-1) .eq. 1) then
                        ipm = ipm+1
                    end if
                    if (ipm .eq. numno0) then
                        if (zi(ideeq+(2*j)-1) .eq. itran1) then
                            dyi = zr(ivale+j-1)
                        else if (zi(ideeq+(2*j)-1) .eq. itran2) then
                            dzi = zr(ivale+j-1)
                        else if (zi(ideeq+(2*j)-1) .eq. irota1) then
                            drzi = zr(ivale+j-1)
                        else if (zi(ideeq+(2*j)-1) .eq. irota2) then
                            dryi = zr(ivale+j-1)
                        end if
                    else if (ipm .eq. (numno0+1)) then
                        goto 131
                    end if
                end do
131             continue
                ptran = dyi*dyi+dzi*dzi
                prota = drzi*drzi+dryi*dryi
                vecr3(2*imod-1) = ptran
                vecr3(2*imod) = prota
                call jelibe(nomcha)
!
                vecr1(imod) = zr(lmasg)+comaj1*ptran+comaj2*prota
!
            end do
!
!
! ---   2.6.IMPRESSION  ---
!
            if (iimpr .eq. 1) then
                valk(1) = nomno0
                valk(2) = base
                valk(3) = caelem(1:8)
                valk(4) = config(indic)
                valr(1) = phie
                valr(2) = cm1
                valr(3) = rhof
                call utmess('I', 'ALGELINE5_13', nk=4, valk=valk, nr=3, &
                            valr=valr)
            end if
!
            call jedetr('&&MDCONF.TEMP.MAIL')
!
        else
! ---  PAS DE COUPLAGE
!
            if (iimpr .eq. 1) then
                call utmess('I', 'ALGELINE5_14')
            end if
!
        end if
!
! --- 4.AUTRES CONFIGURATIONS NON TRAITEES  ---
!
    else
!
        call utmess('F', 'ALGELINE_72')
!
    end if
!
    call jedema()
!
! --- FIN DE MDCONF.
end subroutine
