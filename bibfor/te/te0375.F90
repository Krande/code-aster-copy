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
subroutine te0375(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calnor.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/elref7.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ermeb3.h"
#include "asterfort/ermes3.h"
#include "asterfort/ermev3.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/intega.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jevech.h"
#include "asterfort/jexnum.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/resr3d.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/uthk.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!       CALCUL DE L'INDICATEUR D'ERREUR EN MECANIQUE 3D AVEC LA
!       METHODE DES RESIDUS EXPLICITES.
!       OPTION : 'ERME_ELEM'
!
! REMARQUE : LES PROGRAMMES SUIVANTS DOIVENT RESTER TRES SIMILAIRES
!            TE0368, TE0375, TE0377, TE0378, TE0382, TE0497
!
! ......................................................................
!
!
!
!
! DECLARATION VARIABLES LOCALES
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iadzi, iazk24
    integer(kind=8) :: ibid, iaux, iret, itab(7)
    integer(kind=8) :: noe(9, 6, 4)
    integer(kind=8) :: igeom, jtime
    integer(kind=8) :: ierr, ivois
    integer(kind=8) :: imate
    integer(kind=8) :: iad
    integer(kind=8) :: ifovr, ifovf
    integer(kind=8) :: ipes, irot
    integer(kind=8) :: iref1, iref2
    integer(kind=8) :: ndim
    integer(kind=8) :: nno, nnos, npg, ipoids, ivf, idfde, jgano
    integer(kind=8) :: ndimf
    integer(kind=8) :: nnof, nnosf, npgf, ipoidf, ivff, idfdxf, jganof
    integer(kind=8) :: nno2, nnos2, npg2, ipoid2, ivf2, idfdx2, jgano2
    integer(kind=8) :: nbcmp, kpg, spt
    integer(kind=8) :: ipg
    integer(kind=8) :: ipgf
    integer(kind=8) :: nbf
    integer(kind=8) :: tymvol, ndegre, ifa, tyv
!
    real(kind=8) :: r8bid, r8bid2, r8bid3(3), r8bid4(3)
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), hk, poids
    real(kind=8) :: fpx, fpy, fpz
    real(kind=8) :: frx(27), fry(27), frz(27)
    real(kind=8) :: fovo(3)
    real(kind=8) :: dsx, dsy, dsz
    real(kind=8) :: errest, nor, norsig, sigcal, nuest, coeff
    real(kind=8) :: ter1, ter2, ter3, hf, inte, inst
    real(kind=8) :: nx(9), ny(9), nz(9), jaco(9)
    real(kind=8) :: chx(9), chy(9), chz(9)
    real(kind=8) :: dsg11(9), dsg22(9), dsg33(9), dsg12(9), dsg13(9), dsg23(9)
    real(kind=8) :: sig11(9), sig22(9), sig12(9), sig33(9), sig13(9), sig23(9)
    real(kind=8) :: rho, valres(1)
!
    integer(kind=8) :: icodre(2)
    character(len=8) :: typmav, elrefe, fami, poum
    character(len=8) :: elreff, elrefb
    character(len=8) :: nompar(1)
    character(len=16) :: phenom
    character(len=24) :: valk(2)
!
    aster_logical :: yapr, yaro
!
! --- INITIALISATION DU TABLEAU DES NUMEROS DE NOEUDS FACE PAR FACE ----
!
!     NOE (IN,IFA,TYMVOL) : IN     : NUMERO DU NOEUD DANS LA FACE
!                           IFA    : NUMERO DE LA FACE
!                           TYMVOL : TYPE DE LA MAILLE VOLUMIQUE
!                                    1 : HEXAEDRE A 8,20 ET 27 NOEUDS
!                                    2 : PENTAEDRE A 6,15 ET 18 NOEUDS
!                                    3 : TETRAEDRE A 4 ET 10 NOEUDS
!                                    4 : PYRAMIDE A 5 ET 13 NOEUDS
!     VOIR TE003 POUR LES EXPLICATIONS DETAILLEES
!
    data noe/1, 4, 3, 2, 12, 11, 10, 9, 21, 1, 2, 6, 5, 9, 14, 17, 13, 22,&
     &         2, 3, 7, 6, 10, 15, 18, 14, 23, 3, 4, 8, 7, 11, 16, 19, 15, 24,&
     &         4, 1, 5, 8, 12, 13, 20, 16, 25, 5, 6, 7, 8, 17, 18, 19, 20, 26,&
     &         1, 3, 2, 9, 8, 7, 3*0, 4, 5, 6, 13, 14, 15, 3*0,&
     &         1, 2, 5, 4, 7, 11, 13, 10, 16, 2, 3, 6, 5, 8, 12, 14, 11, 17,&
     &         1, 4, 6, 3, 10, 15, 12, 9, 18, 9*0,&
     &         1, 3, 2, 7, 6, 5, 3*0, 2, 3, 4, 6, 10, 9, 3*0,&
     &         3, 1, 4, 7, 8, 10, 3*0, 1, 2, 4, 5, 9, 8, 3*0,&
     &         9*0, 9*0,&
     &         1, 2, 5, 6, 11, 10, 3*0, 2, 3, 5, 7, 12, 11, 3*0,&
     &         3, 4, 5, 8, 13, 12, 3*0, 4, 1, 5, 9, 10, 13, 3*0,&
     &         1, 4, 3, 2, 9, 8, 7, 6, 0, 9*0/
!
! ----------------------------------------------------------------------
1000 format(a, ' :', (6(1x, 1pe17.10)))
! ----------------------------------------------------------------------
! 1 -------------- GESTION DES DONNEES ---------------------------------
! ----------------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
! 1.1. --- LES INCONTOURNABLES
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVOISIN', 'L', ivois)
!
    call jevech('PINSTR', 'L', jtime)
    inst = zr(jtime-1+1)
!
    call jevech('PERREUR', 'E', ierr)
!
! 1.2. --- LES CARACTERISTIQUES DE LA MAILLE EN COURS
!
    call tecael(iadzi, iazk24)
    valk(1) = zk24(iazk24-1+3)
    valk(2) = option
!
    call elref1(elrefe)
!
    if (niv .ge. 2) then
        write (ifm, *) ' '
        write (ifm, *) '================================================='
        write (ifm, *) ' '
        write (ifm, *) 'MAILLE NUMERO', zi(iadzi), ', DE TYPE ', elrefe
    end if
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! 1.3. --- CHAMP DE CONTRAINTES
!
    call tecach('OOO', 'PCONTNO', 'L', iret, nval=3, &
                itab=itab)
    iad = itab(1)
    nbcmp = itab(2)/nno
!
! 1.4. --- CARTES DE PESANTEUR ET ROTATION
!
    call tecach('ONO', 'PPESANR', 'L', iret, iad=itab(1))
    if (itab(1) .ne. 0) then
        call jevech('PPESANR', 'L', ipes)
        yapr = .true.
    else
        yapr = .false.
    end if
    call tecach('ONO', 'PROTATR', 'L', iret, iad=itab(1))
    if (itab(1) .ne. 0) then
        call jevech('PROTATR', 'L', irot)
        yaro = .true.
    else
        yaro = .false.
    end if
!
! 1.5. --- FORCES VOLUMIQUES EVENTUELLES
!          VALEURS REELLES ?
    call tecach('ONO', 'PFRVOLU', 'L', iret, iad=ifovr)
!          OU FONCTIONS ?
    if (ifovr .eq. 0) then
        call tecach('ONO', 'PFFVOLU', 'L', iret, iad=ifovf)
    else
        ifovf = 0
    end if
!GN      WRITE(IFM,2000) 'IFOVR', IFOVR
!GN      WRITE(IFM,2000) 'IFOVF', IFOVF
!
! 1.6. --- FORCES ET PRESSIONS AUX BORDS
!
    call jevech('PFORCE', 'L', iref1)
!
    call jevech('PPRESS', 'L', iref2)
!
! 1.7. --- MATERIAU SI BESOIN
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    if (yapr .or. yaro) then
!
        call jevech('PMATERC', 'L', imate)
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
        nompar(1) = 'RHO'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, 1, ' ', [r8bid], &
                    1, nompar, valres, icodre, 1)
        rho = valres(1)
!GN        WRITE(IFM,1000) 'RHO', RHO
!
    end if
!
! ----------------------------------------------------------------------
! 2 -------------- CALCUL DU PREMIER TERME DE L'ERREUR -----------------
! ----------------------------------------------------------------------
!
! 2.1. --- CALCUL DU DIAMETRE HK DE LA MAILLE ----
!
    call uthk(nomte, zr(igeom), hk, ndim, niv)
!
! 2.2. --- CALCUL DE LA FORCE DE PESANTEUR ---
!
    if (yapr) then
        fpx = rho*zr(ipes)*zr(ipes+1)
        fpy = rho*zr(ipes)*zr(ipes+2)
        fpz = rho*zr(ipes)*zr(ipes+3)
    else
        fpx = 0.d0
        fpy = 0.d0
        fpz = 0.d0
    end if
!GN      WRITE(IFM,1000) 'P',FPX,FPY,FPZ
!
! 2.3. --- CALCUL DE LA FORCE DE ROTATION ---
!
    if (yaro) then
        call resr3d(zr(irot), zr(igeom), zr(ivf), rho, nno, &
                    npg, frx, fry, frz)
    else
        do 23, ipg = 1, npg
            frx(ipg) = 0.d0
            fry(ipg) = 0.d0
            frz(ipg) = 0.d0
23          continue
            end if
!GN      WRITE(IFM,1000) 'R X',(FRX(IPG),IPG = 1 , NPG)
!GN      WRITE(IFM,1000) 'R Y',(FRY(IPG),IPG = 1 , NPG)
!GN      WRITE(IFM,1000) 'R Z',(FRZ(IPG),IPG = 1 , NPG)
!
! 2.4. --- CALCUL DE LA FORCE VOLUMIQUE EVENTUELLE ---
!
            if (ifovr .ne. 0) then
                fovo(1) = zr(ifovr)
                fovo(2) = zr(ifovr+1)
                fovo(3) = zr(ifovr+2)
!
            else if (ifovf .ne. 0) then
                nompar(1) = 'INST'
                r8bid3(1) = inst
!       SI UNE COMPOSANTE N'A PAS ETE DECRITE, ASTER AURA MIS PAR
!       DEFAUT LA FONCTION NULLE &FOZERO. ON LE REPERE POUR
!       IMPOSER LA VALEUR 0 SANS FAIRE DE CALCULS INUTILES
                do 24, ibid = 1, ndim
                if (zk8(ifovf+ibid-1) (1:7) .eq. '&FOZERO') then
                    fovo(ibid) = 0.d0
                else
                    call fointe('FM', zk8(ifovf+ibid-1), 1, nompar, r8bid3, &
                                fovo(ibid), iret)
                end if
24              continue
!GN        WRITE(IFM,*) 'F X : ',ZK8(IFOVF),FOVO(1)
!GN        WRITE(IFM,*) 'F Y : ',ZK8(IFOVF+1),FOVO(2)
!GN        WRITE(IFM,*) 'F Z : ',ZK8(IFOVF+2),FOVO(3)
                end if
!
! 2.5. --- CALCUL DU TERME D'ERREUR AVEC INTEGRATION DE GAUSS ---
!
                ter1 = 0.d0
                norsig = 0.d0
!
                do 25, ipg = 1, npg
!
! ------- CALCUL DES DERIVEES DES FONCTIONS DE FORMES /X, /Y ET /Z -----
!
                    call dfdm3d(nno, ipg, ipoids, idfde, zr(igeom), &
                                poids, dfdx, dfdy, dfdz)
!
! ------- CALCUL DE LA DIVERGENCE ET DE LA NORME DE SIGMA --------------
!
                    call ermev3(nno, ipg, ivf, iad, nbcmp, &
                                dfdx, dfdy, dfdz, dsx, dsy, &
                                dsz, nor)
!
! ------- CUMUL
!
                    r8bid3(1) = fpx+frx(ipg)+dsx
                    r8bid3(2) = fpy+fry(ipg)+dsy
                    r8bid3(3) = fpz+frz(ipg)+dsz
!
! ------- PRISE EN COMPTE DE L'EFFORT VOLUMIQUE EVENTUEL ---------------
!
                    if (ifovr .ne. 0 .or. ifovf .ne. 0) then
!
!GN          WRITE(IFM,1000) 'F X', FOVO(1)
!GN          WRITE(IFM,1000) 'F Y', FOVO(2)
!GN          WRITE(IFM,1000) 'F Z', FOVO(3)
                        r8bid3(1) = r8bid3(1)+fovo(1)
                        r8bid3(2) = r8bid3(2)+fovo(2)
                        r8bid3(3) = r8bid3(3)+fovo(3)
!
                    end if
!
! ------- CUMUL DU TERME D'ERREUR
!
                    ter1 = ter1+(r8bid3(1)**2+r8bid3(2)**2+r8bid3(3)**2)*poids
                    if (niv .ge. 2) then
                        write (ifm, 1000) 'POIDS', poids
                        write (ifm, 1000) 'A2 + B2 + C2', r8bid3(1)**2+r8bid3(2)** &
                            2+r8bid3(3)**2
                        write (ifm, 1000) '==> TER1    ', ter1
                    end if
!
! ------- CALCUL DE LA NORME DE SIGMA SUR L'ELEMENT --------------------
!
                    norsig = norsig+nor*poids
!
25              end do
!
                if (ter1 .lt. 0.d0) then
                    call utmess('A', 'INDICATEUR_9', nk=2, valk=valk)
                    goto 999
                end if
!
                ter1 = hk*sqrt(ter1)
!
! ----------------------------------------------------------------------
! ------------ FIN DU CALCUL DU PREMIER TERME DE L'ERREUR --------------
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! 3. ------- CALCUL DES DEUXIEME ET TROISIEME TERMES DE L'ERREUR -------
! ----------------------------------------------------------------------
!
! 3.1. ---- INFORMATIONS SUR LA MAILLE COURANTE : ----------------------
!       TYMVOL : TYPE DE LA MAILLE VOLUMIQUE
!       NDEGRE : DEGRE DE L'ELEMENT
!       NBF    : NOMBRE DE FACES DE LA MAILLE VOLUMIQUE
!       ELREFF : DENOMINATION DE LA MAILLE FACE DE ELREFE - FAMILLE 1
!       ELREFB : DENOMINATION DE LA MAILLE FACE DE ELREFE - FAMILLE 2
!      --- REMARQUE : ON IMPOSE UNE FAMILLE DE POINTS DE GAUSS
!
                call elref7(elrefe, tymvol, ndegre, nbf, elreff, &
                            elrefb)
!GN      WRITE(6,*) 'TYPE MAILLE VOLUMIQUE COURANTE :',TYMVOL
! --- CARACTERISTIQUES DES FACES DE BORD DE LA FAMILLE 1 ---------------
                call elrefe_info(elrefe=elreff, fami='NOEU', ndim=ndimf, &
                                 nno=nnof, nnos=nnosf, npg=npgf, jpoids=ipoidf, &
                                 jvf=ivff, jdfde=idfdxf, jgano=jganof)
!GN      WRITE(IFM,2000) 'NDIMF',NDIMF
!GN      WRITE(IFM,2000) 'NNOSF,NNOF,NPGF',NNOSF,NNOF,NPGF
!GN      WRITE(IFM,1000) 'IPOIDF', (ZR(IPOIDF+IFA),IFA=0,NPGF-1)
!
! --- COMPLEMENT EVENTUEL POUR LES MAILLES QUI ONT 2 TYPES DE ---
! --- MAILLES DE BORD (PENTAEDRE, PYRAMIDE) ---
!
                if (elrefb(1:1) .ne. ' ') then
                    call elrefe_info(elrefe=elrefb, fami='NOEU', ndim=ndimf, &
                                     nno=nno2, nnos=nnos2, &
                                     npg=npg2, jpoids=ipoid2, &
                                     jvf=ivf2, jdfde=idfdx2, jgano=jgano2)
!GN       WRITE(IFM,2000) 'NDIMF,NNO2',NDIMF,NNO2
!GN       WRITE(IFM,2000) 'NNOS2,NPG2',NNOS2,NPG2
!GN       WRITE(IFM,1000) 'IPOID2', (ZR(IPOID2+IFA),IFA=0,NPG2-1)
                end if
!
! 3.2. --- BOUCLE SUR LES FACES DE LA MAILLE VOLUMIQUE --------------
!
                ter2 = 0.d0
                ter3 = 0.d0
                do 320, ifa = 1, nbf
!
! ------TEST DU TYPE DE VOISIN -----------------------------------------
!
                    tyv = zi(ivois+7+ifa)
!
                    if (tyv .ne. 0) then
!
! ------- RECUPERATION DU TYPE DE LA MAILLE VOISINE
!
                        call jenuno(jexnum('&CATA.TM.NOMTM', tyv), typmav)
                        if (niv .ge. 2) then
                            write (ifm, 1003) ifa, zi(ivois+ifa), typmav
1003                        format(i2, '-EME FACE DE NUMERO', i10, ' ==> TYPMAV = ', a)
                            write (ifm, 1000) 'TER2', ter2
                            write (ifm, 1000) 'TER3', ter3
                        end if
!
! --- QUAND ON ARRIVE AUX FACES QUAD DES PENTAEDRES OU DES PYRAMIDES ---
! --- IL FAUT REMPLACER LES CARACTERISTIQUES DE LA FAMILLE 1         ---
! --- PAR CELLES DE LA FAMILLE 2                                     ---
!
                        if ((tymvol .eq. 2 .and. ifa .ge. 3) .or. &
                            (tymvol .eq. 4 .and. ifa .ge. 5)) then
!
                            nnof = nno2
                            npgf = npg2
                            nnosf = nnos2
                            ipoidf = ipoid2
                            idfdxf = idfdx2
!
                        end if
!GN      WRITE(IFM,*) '. NPGF =', NPGF
!
! ----- CALCUL DU DIAMETRE HF DE LA FACE ----------
!
                        ibid = 0
                        call uthk(nomte, zr(igeom), hf, ibid, niv, &
                                  noe, nnosf, tymvol, ifa)
!
! ------- CALCUL DE NORMALES ET JACOBIENS AUX POINTS DE GAUSS ----------
!
                        iaux = ifa
                        call calnor('3D', zr(igeom), ibid, ibid, ibid, &
                                    r8bid, nnof, npgf, noe, iaux, &
                                    tymvol, idfdxf, jaco, nx, ny, &
                                    nz, r8bid3, r8bid4, r8bid2)
!
! ----------------------------------------------------------------------
! --------------- CALCUL DU DEUXIEME TERME DE L'ERREUR -----------------
! --------------- LE BORD VOISIN EST UN VOLUME -------------------------
! ----------------------------------------------------------------------
!
                        if (typmav(1:4) .eq. 'HEXA' .or. typmav(1:4) .eq. 'PENT' &
                            .or. typmav(1:4) .eq. &
                            'TETR' .or. typmav(1:4) .eq. 'PYRA') then
!
! ------- CALCUL DU SAUT DE CONTRAINTE ENTRE ELEMENTS ------------------
! ------- CE CHAMP DSGXX EST EXPRIME SUR LES NOEUDS DE LA FACE ---------
!
                            call ermes3(noe, ifa, tymvol, nnof, typmav, &
                                        iref1, ivois, iad, nbcmp, dsg11, &
                                        dsg22, dsg33, dsg12, dsg13, dsg23)
!
! ------- CALCUL DE L'INTEGRALE SUR LA FACE ----------------------------
! ------- ATTENTION : CELA MARCHE CAR ON A CHOISI LA FAMILLE -----------
! ------- AVEC LES POINTS DE GAUSS SUR LES NOEUDS ----------------------
!
                            do 321, ipgf = 1, npgf
                                chx(ipgf) = 0.d0
                                chy(ipgf) = 0.d0
                                chz(ipgf) = 0.d0
321                             continue
!
                                call intega(npgf, jaco, zr(ipoidf), chx, chy, &
                                            chz, dsg11, dsg22, dsg33, dsg12, &
                                            dsg13, dsg23, nx, ny, nz, &
                                            inte)
!
! ------- CALCUL DU TERME D'ERREUR -------------------------------------
!
                                if (inte .lt. 0.d0) then
                                    call utmess('A', 'INDICATEUR_9', nk=2, valk=valk)
                                    goto 999
                                end if
!
                                ter2 = ter2+0.5d0*sqrt(hf)*sqrt(inte)
                                if (niv .ge. 2) then
                                    write (ifm, 1000) 'VOLU INTE', inte
                                    write (ifm, 1000) '==> TER2 ', ter2
                                end if
!
! ----------------------------------------------------------------------
! --------------- CALCUL DU TROISIEME TERME DE L'ERREUR ----------------
! --------------- LE BORD VOISIN EST UNE FACE --------------------------
! ----------------------------------------------------------------------
!
                                elseif (typmav(1:4) .eq. 'QUAD' .or. typmav(1:4) .eq. 'TRIA' &
                                        ) then
!
! ------- CALCUL EFFORTS SURFACIQUES ET DES CONTRAINTES ----------------
!
                                call ermeb3(noe, ifa, tymvol, nnof, iref1, &
                                            iref2, ivois, igeom, iad, nbcmp, &
                                            inst, nx, ny, nz, sig11, &
                                            sig22, sig33, sig12, sig13, sig23, &
                                            chx, chy, chz)
!
! ------- CALCUL DE L'INTEGRALE SUR LA FACE ----------------------------
! ------- ATTENTION : CELA MARCHE CAR ON A CHOISI LA FAMILLE -----------
! ------- AVEC LES POINTS DE GAUSS SUR LES NOEUDS ----------------------
!
                                call intega(npgf, jaco, zr(ipoidf), chx, chy, &
                                            chz, sig11, sig22, sig33, sig12, &
                                            sig13, sig23, nx, ny, nz, &
                                            inte)
!
! ------- CALCUL DU TERME D'ERREUR -------------------------------------
!
                                if (inte .lt. 0.d0) then
                                    call utmess('A', 'INDICATEUR_9', nk=2, valk=valk)
                                    goto 999
                                end if
!
!GN       WRITE(IFM,*) '==> INTE', INTE
                                ter3 = ter3+sqrt(hf)*sqrt(inte)
                                if (niv .ge. 2) then
                                    write (ifm, 1000) 'SURF INTE', inte
                                    write (ifm, 1000) '==> TER3 ', ter3
                                end if
!
! ----------------------------------------------------------------------
! --------------- CURIEUX ----------------------------------------------
! ----------------------------------------------------------------------
!
                                else
!
                                valk(1) = typmav(1:4)
                                call utmess('F', 'INDICATEUR_10', sk=valk(1))
!
                                end if
!
                            end if
!
320                         end do
!
! ----------------------------------------------------------------------
! ------- FIN DU CALCUL DU DEUXIEME ET TROISIEME TERME DE L'ERREUR -----
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! 4. ------- MISE EN MEMOIRE DES DIFFERENTS TERMES DE L'ERREUR ---------
! ----------------------------------------------------------------------
!
                            if (ndegre .eq. 2) then
                                coeff = sqrt(96.d0)
                            else if (ndegre .eq. 1) then
                                coeff = sqrt(24.d0)
                            end if
!
                            errest = (ter1+ter2+ter3)/coeff
                            if ((errest**2+norsig) .ne. 0.d0) then
                                nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
                            else
                                nuest = 0.d0
                            end if
                            sigcal = sqrt(norsig)
!
                            zr(ierr) = errest
                            zr(ierr+1) = nuest
                            zr(ierr+2) = sigcal
!
                            errest = ter1/coeff
                            if ((errest**2+norsig) .ne. 0.d0) then
                                nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
                            else
                                nuest = 0.d0
                            end if
!       TERMRE       TERMR2
                            zr(ierr+3) = errest
                            zr(ierr+4) = nuest
!
                            errest = ter3/coeff
                            if ((errest**2+norsig) .ne. 0.d0) then
                                nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
                            else
                                nuest = 0.d0
                            end if
!       TERMNO     TERMN2
                            zr(ierr+5) = errest
                            zr(ierr+6) = nuest
!
                            errest = ter2/coeff
                            if ((errest**2+norsig) .ne. 0.d0) then
                                nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
                            else
                                nuest = 0.d0
                            end if
!       TERMSA       TERMS2
                            zr(ierr+7) = errest
                            zr(ierr+8) = nuest
!       DIAMETRE
                            zr(ierr+9) = hk
!
999                         continue
!
                            call jedema()
!
                            end subroutine
