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

subroutine te0377(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
!
!     BUT:
!       CALCUL DE L'INDICATEUR D'ERREUR EN MECANIQUE 2D AVEC LA
!       METHODE DES RESIDUS EXPLICITES.
!       OPTION : 'ERME_ELEM'
!
! REMARQUE : LES PROGRAMMES SUIVANTS DOIVENT RESTER TRES SIMILAIRES
!            TE0368, TE0375, TE0377, TE0378, TE0382, TE0497
!
! ----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/calnor.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elref1.h"
#include "asterfort/elref7.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/ermeb2.h"
#include "asterfort/ermes2.h"
#include "asterfort/ermev2.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/intenc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jevech.h"
#include "asterfort/jexnum.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/resrot.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/uthk.h"
#include "asterfort/utjac.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
!
!
! DECLARATION VARIABLES LOCALES
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iadzi, iazk24
    integer(kind=8) :: ibid, iaux, iret, itab(7), noe(9, 6, 4)
    integer(kind=8) :: igeom, jtime
    integer(kind=8) :: ierr, ivois
    integer(kind=8) :: imate
    integer(kind=8) :: iad
    integer(kind=8) :: ifovr, ifovf
    integer(kind=8) :: ipes, irot
    integer(kind=8) :: iref1, iref2, ivf2
    integer(kind=8) :: ndim
    integer(kind=8) :: nno, nnos, npg, ipoids, ivf, idfde, jgano
    integer(kind=8) :: nnof, npgf
    integer(kind=8) :: nbcmp, nbpar
    integer(kind=8) :: ipg
    integer(kind=8) :: ipgf
    integer(kind=8) :: nbf
    integer(kind=8) :: tymvol, ndegre, ifa, tyv
!
    real(kind=8) :: r8bid3(3)
    real(kind=8) :: dfdx(9), dfdy(9), hk, poids
    real(kind=8) :: fpx, fpy
    real(kind=8) :: frx(9), fry(9)
    real(kind=8) :: fovo(2)
    real(kind=8) :: dsx, dsy
    real(kind=8) :: errest, nor, norsig, sigcal, nuest, coeff
    real(kind=8) :: ter1, ter2, ter3, hf, inte, inst
    real(kind=8) :: nx(9), ny(9), nz(9), jaco(9), orien
    real(kind=8) :: chx(3), chy(3)
    real(kind=8) :: sg11(3), sg22(3), sg12(3)
    real(kind=8) :: tx(3), ty(3)
    real(kind=8) :: sig11(3), sig22(3), sig12(3)
    real(kind=8) :: e, nu, rho, valres(3)
    integer(kind=8), parameter :: nb_para = 3
    real(kind=8) :: para_vale(nb_para)
    character(len=16), parameter :: para_name(nb_para) = (/'X', 'Y', 'Z'/)
!
    integer(kind=8) :: icodre(3)
    character(len=3) :: typnor
    character(len=8) :: typmav, elrefe
    character(len=8) :: elreff, elrefb
    character(len=8) :: nompar(3)
    character(len=16) :: phenom
    character(len=24) :: valk(2)
    type(Behaviour_Integ) :: BEHinteg
!
    aster_logical :: yapr, yaro
!
! ----------------------------------------------------------------------
! ----- NORME CALCULEE : SEMI-H1 (H1) ou ENERGIE (NRJ) -----------------
! ----------------------------------------------------------------------
!
    data typnor/'NRJ'/
!
! ----------------------------------------------------------------------
111 format(a, ' :', (6(1x, 1pe17.10)))
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
    call jevech('PINSTR', 'L', jtime)
    inst = zr(jtime-1+1)
!
! - Initialisation of behaviour datastructure
!
    call behaviourInit(BEHinteg)
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
!
! 1.6. --- FORCES ET PRESSIONS AUX BORDS
!
    call jevech('PFORCE', 'L', iref1)
!
    call jevech('PPRESS', 'L', iref2)
    call elrefe_info(fami='FPG1', jvf=ivf2)
    call behaviourCoorGauss(nno, 1, ndim, &
                            ivf2, zr(igeom), BEHinteg%behavESVA%behavESVAGeom)
!
! 1.7. --- MATERIAU SI BESOIN
!
    if (yapr .or. yaro .or. typnor .eq. 'NRJ') then
!
        call jevech('PMATERC', 'L', imate)
        call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
        nbpar = 0
        if (typnor .eq. 'NRJ') then
            nbpar = nbpar+1
            nompar(nbpar) = 'E'
            nbpar = nbpar+1
            nompar(nbpar) = 'NU'
        end if
        if (yapr .or. yaro) then
            nbpar = nbpar+1
            nompar(nbpar) = 'RHO'
        end if
!
        para_vale(:) = BEHinteg%behavESVA%behavESVAGeom%coorElga(1, :)
        call rcvalb('FPG1', 1, 1, '+', zi(imate), &
                    ' ', phenom, nb_para, para_name, para_vale, &
                    nbpar, nompar, valres, icodre, 1)
!
        if (typnor .eq. 'NRJ') then
            e = valres(1)
            nu = valres(2)
        end if
        if (yapr .or. yaro) then
            rho = valres(nbpar)
        end if
!GN        WRITE(IFM,111) 'RHO, E, NU', RHO, E, NU
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
    else
        fpx = 0.d0
        fpy = 0.d0
    end if
!GN      WRITE(IFM,111) 'P',FPX,FPY,FPZ
!
! 2.3. --- CALCUL DE LA FORCE DE ROTATION ---
!
    if (yaro) then
        call resrot(zr(irot), zr(igeom), zr(ivf), rho, nno, &
                    npg, frx, fry)
    else
        do ipg = 1, npg
            frx(ipg) = 0.d0
            fry(ipg) = 0.d0
        end do
    end if
!GN      WRITE(IFM,111) 'R X',(FRX(IPG),IPG = 1 , NPG)
!GN      WRITE(IFM,111) 'R Y',(FRY(IPG),IPG = 1 , NPG)
!
! 2.4. --- CALCUL DE LA FORCE VOLUMIQUE EVENTUELLE ---
!
    if (ifovr .ne. 0) then
        fovo(1) = zr(ifovr)
        fovo(2) = zr(ifovr+1)
!
    else if (ifovf .ne. 0) then
        nompar(1) = 'INST'
        r8bid3(1) = inst
!       SI UNE COMPOSANTE N'A PAS ETE DECRITE, ASTER AURA MIS PAR
!       DEFAUT LA FONCTION NULLE &FOZERO. ON LE REPERE POUR
!       IMPOSER LA VALEUR 0 SANS FAIRE DE CALCULS INUTILES
        do ibid = 1, ndim
        if (zk8(ifovf+ibid-1) (1:7) .eq. '&FOZERO') then
            fovo(ibid) = 0.d0
        else
            call fointe('FM', zk8(ifovf+ibid-1), 1, nompar, r8bid3, &
                        fovo(ibid), iret)
        end if
        end do
!GN        WRITE(IFM,*) 'F X : ',ZK8(IFOVF),FOVO(1)
!GN        WRITE(IFM,*) 'F Y : ',ZK8(IFOVF+1),FOVO(2)
    end if
!
! 2.5. --- CALCUL DU TERME D'ERREUR AVEC INTEGRATION DE GAUSS ---
!
    ter1 = 0.d0
    norsig = 0.d0
!
    do ipg = 1, npg
!
! ------- CALCUL DES DERIVEES DES FONCTIONS DE FORMES /X ET /Y ---------
!
        call dfdm2d(nno, ipg, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy)
!
! ------- CALCUL DE L'ORIENTATION DE LA MAILLE -------------------------
!
        call utjac(.true._1, zr(igeom), ipg, idfde, 0, &
                   ibid, nno, orien)
!
! ------- CALCUL DE LA DIVERGENCE ET DE LA NORME DE SIGMA --------------
!
        iaux = ivf+(ipg-1)*nno
        ibid = 1
        call ermev2(nno, igeom, zr(iaux), zr(iad), nbcmp, &
                    dfdx, dfdy, poids, ibid, dsx, &
                    dsy, nor)
!
! ------- CUMUL
!
        r8bid3(1) = fpx+frx(ipg)+dsx
        r8bid3(2) = fpy+fry(ipg)+dsy
!
! ------- PRISE EN COMPTE DE L'EFFORT VOLUMIQUE EVENTUEL ---------------
!
        if (ifovr .ne. 0 .or. ifovf .ne. 0) then
!
!GN          WRITE(IFM,111) 'F X', FOVO(1)
!GN          WRITE(IFM,111) 'F Y', FOVO(2)
            r8bid3(1) = r8bid3(1)+fovo(1)
            r8bid3(2) = r8bid3(2)+fovo(2)
!
        end if
!
! ------- CUMUL DU TERME D'ERREUR
!
        ter1 = ter1+(r8bid3(1)**2+r8bid3(2)**2)*poids
        if (niv .ge. 2) then
            write (ifm, 111) 'POIDS', poids
            write (ifm, 111) 'A2 + B2 ', r8bid3(1)**2+r8bid3(2)**2
            write (ifm, 111) '==> TER1', ter1
        end if
!
! ------- CALCUL DE LA NORME DE SIGMA SUR L'ELEMENT --------------------
!
        norsig = norsig+nor*poids
!
    end do
!
    if (typnor(1:2) .eq. 'H1') then
!       NORME H1
        ter1 = hk*sqrt(ter1)
    else if (typnor .eq. 'NRJ') then
!       NORME EN ENERGIE
        ter1 = (hk**2)*abs(ter1)
    end if
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
! --- CARACTERISTIQUES DES FACES DE BORD -------------------------------
!     on est tente de faire l'appel a elrefe_info comme en 3d mais c'est en
!     FAIT INUTILE CAR ON N'A BESOIN QUE DE NNOF ET NPGF.
!     CELA TOMBE BIEN CAR L'APPEL MARCHE RAREMENT ...
!
    if (ndegre .eq. 1) then
        nnof = 2
    else
        nnof = 3
    end if
    npgf = nnof
!GN      CALL ELREF4 ( ELREFF,fami='RIGI',
!GN     >              NDIMF, NNOF, NNOSF, NPGF, IPOIDF, IVFF,
!GN     >              IDFDXF, JGANOF )
!GN      WRITE(IFM,222) 'NDIMF',NDIMF
!GN      WRITE(IFM,222) 'NNOSF,NNOF,NPGF',NNOSF,NNOF,NPGF
!GN      WRITE(IFM,111) 'IPOIDF', (ZR(IPOIDF+IFA),IFA=0,NPGF-1)
!
! 3.2. --- BOUCLE SUR LES FACES DE LA MAILLE VOLUMIQUE --------------
!
    ter2 = 0.d0
    ter3 = 0.d0
    do ifa = 1, nbf
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
                write (ifm, 103) ifa, zi(ivois+ifa), typmav
103             format(i2, '-EME FACE DE NUMERO', i10, ' ==> TYPMAV = ', a)
            end if
!
! ----- CALCUL DE NORMALES, TANGENTES ET JACOBIENS AUX POINTS DE GAUSS
!
            iaux = ifa
            call calnor('2D', zr(igeom), iaux, nnos, nnof, &
                        orien, ibid, ibid, noe, ibid, &
                        ibid, ibid, jaco, nx, ny, &
                        nz, tx, ty, hf)
!
! ----------------------------------------------------------------------
! --------------- CALCUL DU DEUXIEME TERME DE L'ERREUR -----------------
! --------------- LE BORD VOISIN EST UN VOLUME -------------------------
! ----------------------------------------------------------------------
!
            if (typmav(1:4) .eq. 'TRIA' .or. typmav(1:4) .eq. 'QUAD') then
!
! ------- CALCUL DU SAUT DE CONTRAINTE ENTRE ELEMENTS ------------------
!
                iaux = ifa
                call ermes2(iaux, elrefe, typmav, iref1, ivois, &
                            iad, nbcmp, sg11, sg22, sg12)
!
! ------- CALCUL DE L'INTEGRALE SUR LA FACE ----------------------------
! ------- CALCUL DU TERME D'ERREUR AVEC INTEGRATION DE NEWTON-COTES ----
! ------- ATTENTION : CELA MARCHE CAR ON A CHOISI LA FAMILLE -----------
! ------- AVEC LES POINTS DE GAUSS SUR LES NOEUDS ----------------------
!
                do ipgf = 1, npgf
                    chx(ipgf) = 0.d0
                    chy(ipgf) = 0.d0
                end do
!
                call intenc(nnof, jaco, chx, chy, sg11, &
                            sg22, sg12, nx, ny, inte)
!
! ------- CALCUL DU TERME D'ERREUR -------------------------------------
!
                if (inte .lt. 0.d0) then
                    call utmess('A', 'INDICATEUR_9', nk=2, valk=valk)
                    goto 999
                end if
!
                if (typnor(1:2) .eq. 'H1') then
!             NORME H1
                    ter2 = ter2+0.5d0*sqrt(hf)*sqrt(inte)
                else if (typnor .eq. 'NRJ') then
!             NORME EN ENERGIE
                    ter2 = ter2+0.5d0*hf*inte
                end if
                if (niv .ge. 2) then
                    write (ifm, 111) 'VOLU INTE', inte
                    write (ifm, 111) '==> TER2 ', ter2
                end if
!
! ----------------------------------------------------------------------
! --------------- CALCUL DU TROISIEME TERME DE L'ERREUR ----------------
! --------------- LE BORD VOISIN EST UNE FACE --------------------------
! ----------------------------------------------------------------------
!
            else if (typmav(1:3) .eq. 'SEG') then
!
! ------- CALCUL EFFORTS SURFACIQUES ET DES CONTRAINTES ----------------
!
                iaux = ifa
                call ermeb2(iaux, iref1, iref2, ivois, igeom, &
                            iad, elrefe, nbcmp, inst, nx, &
                            ny, tx, ty, sig11, sig22, &
                            sig12, chx, chy)
!
! ------- CALCUL DE L'INTEGRALE SUR LE BORD ----------------------------
!
                call intenc(nnof, jaco, chx, chy, sig11, &
                            sig22, sig12, nx, ny, inte)
!
! ------- CALCUL DU TERME D'ERREUR -------------------------------------
!
                if (inte .lt. 0.d0) then
                    call utmess('A', 'INDICATEUR_9', nk=2, valk=valk)
                    goto 999
                end if
!
                if (typnor(1:2) .eq. 'H1') then
!             NORME H1
                    ter3 = ter3+sqrt(hf)*sqrt(inte)
                else if (typnor .eq. 'NRJ') then
!             NORME EN ENERGIE
                    ter3 = ter3+hf*inte
                end if
                if (niv .ge. 2) then
                    write (ifm, 111) 'SURF INTE', inte
                    write (ifm, 111) '==> TER3 ', ter3
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
    end do
!
! ----------------------------------------------------------------------
! ------- FIN DU CALCUL DU DEUXIEME ET TROISIEME TERME DE L'ERREUR -----
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! 4. ------- MISE EN MEMOIRE DES DIFFERENTS TERMES DE L'ERREUR ---------
! ----------------------------------------------------------------------
!
    if (typnor(1:2) .eq. 'H1') then
!
        if (ndegre .eq. 2) then
            coeff = sqrt(96.d0)
        else
            coeff = sqrt(24.d0)
        end if
!
!      NORME H1
        errest = (ter1+ter2+ter3)/coeff
        sigcal = sqrt(norsig)
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr) = errest
        zr(ierr+1) = nuest
        zr(ierr+2) = sigcal
!
        errest = ter1/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+3) = errest
        zr(ierr+4) = nuest
!
        errest = ter3/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+5) = errest
        zr(ierr+6) = nuest
!
        errest = ter2/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+7) = errest
        zr(ierr+8) = nuest
!
    else if (typnor .eq. 'NRJ') then
!
        if (ndegre .eq. 2) then
            coeff = sqrt(96.d0*e/(1-nu))
        else
            coeff = sqrt(24.d0*e/(1-nu))
        end if
!
!      NORME EN ENERGIE
        errest = sqrt(ter1+ter2+ter3)/coeff
        sigcal = sqrt(norsig)
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr) = errest
        zr(ierr+1) = nuest
        zr(ierr+2) = sigcal
!
        errest = sqrt(ter1)/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+3) = errest
        zr(ierr+4) = nuest
!
        errest = sqrt(ter3)/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+5) = errest
        zr(ierr+6) = nuest
!
        errest = sqrt(ter2)/coeff
        if (abs(errest**2+norsig) .gt. r8prem()) then
            nuest = 100.d0*sqrt(errest**2/(errest**2+norsig))
        else
            nuest = 0.d0
        end if
!
        zr(ierr+7) = errest
        zr(ierr+8) = nuest
!
    end if
!       DIAMETRE
    zr(ierr+9) = hk
!
999 continue
!
    call jedema()
!
end subroutine
