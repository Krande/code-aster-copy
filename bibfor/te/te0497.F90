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
subroutine te0497(option, nomte)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/thmGetElemPara.h"
#include "asterfort/calnor.h"
#include "asterfort/erhmb2.h"
#include "asterfort/erhms2.h"
#include "asterfort/erhmv2.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jevech.h"
#include "asterfort/jexnum.h"
#include "asterfort/rcvalb.h"
#include "asterfort/resrot.h"
#include "asterfort/tecach.h"
#include "asterfort/uthk.h"
#include "asterfort/utjac.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: option, nomte
!
!-----------------------------------------------------------------------
!    - FONCTION REALISEE:  CALCUL DE L'ESTIMATEUR D'ERREUR EN RESIDU
!      SUR UN ELEMENT ISOPARAMETRIQUE 2D, VIA L'OPTION 'ERME_ELEM'
!      POUR LES MODELISATIONS HM SATUREES
!   -------------------------------------------------------------------
!
! REMARQUE : LES PROGRAMMES SUIVANTS DOIVENT RESTER TRES SIMILAIRES
!            TE0368, TE0375, TE0377, TE0378, TE0382, TE0497
!
!   -------------------------------------------------------------------
!     ASTER INFORMATIONS :
!       03/07/06 (SM): CREATION EN S'INSPIRANT DE TE0003.F ET DE
!                      TE0377.F .
!                      CALCUL INDICATEURS EN STATIONNAIRE
!       01/05/07 (SM): ADIMENSIONNEMENT DES INDICATEURS EN
!                      STATIONNAIRE .
!----------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: ibid, iaux, iret, itab(7)
    integer(kind=8) :: igeom
    integer(kind=8) :: ierr, ivois
    integer(kind=8) :: ierrm, imate, ifovr, ifovf
    integer(kind=8) :: ipes, irot, iref1, iref2, ndim
    integer(kind=8) :: nno, nnos, npg, jv_poids, jv_func, jv_dfunc, jv_gano
    integer(kind=8) :: jv_poids2, jv_func2, jv_dfunc2
    integer(kind=8) :: nbcmp, ipg, ifa, tyv, nbs, kpg, spt
    integer(kind=8) :: isienp, isienm, ideplp, ideplm, jkp, nbna
    integer(kind=8) :: iagd, iatyma, typ, iacmp
    integer(kind=8) :: iade2, iava2, iaptm2, igd2, ncmpm2
    integer(kind=8) :: iade3, iava3, iaptm3, igd3, ncmpm3
    integer(kind=8) :: igrdca, dimdep, dimdef, dimcon
    integer(kind=8) :: nddl_meca, npi, nddl_p1, nddl_p2, nddl_2nd, nddls, nddlm
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5), dimuel
    integer(kind=8) :: adsip, addeme, adcome, addete
    integer(kind=8) :: addep1, adcp11
    integer(kind=8) :: addep2, adde2nd, ii, noe(9, 6, 4)
!
    real(kind=8) :: ovfl
    real(kind=8) :: r8bid3(2)
    real(kind=8) :: valres(1)
    real(kind=8) :: orien, nx(9), ny(9), nz(9), tx(3), ty(3), hf
    real(kind=8) :: fpx, fpy
    real(kind=8) :: frx(9), fry(9)
    real(kind=8) :: fovo(2)
    real(kind=8) :: instpm(2)
    real(kind=8) :: biot, rholiq, unsurk, unsurm, jaco(9)
    real(kind=8) :: hk, deltat, theta
    real(kind=8) :: cyoung, rhohom, permin, viscli, porosi, poisso
    real(kind=8) :: tm2h1v(3), tm2h1b(3), tm2h1s(3)
    real(kind=8) :: tsivom, tdevom, tsivoh, tsibom, tdebom, tsibsh, tsisam, tdesam, tsissh
    real(kind=8) :: denomi
    real(kind=8) :: longc, presc, admec, adhy0, adhy1, adv1h, adhymd
!
    aster_logical :: l_axi
!
    character(len=2) :: form, noeu
    character(len=3) :: inte_type
    character(len=4) :: nompar(1)
    character(len=8) :: typema, typmav
    character(len=8) :: type_elem(2), fami, poum
!
    integer(kind=8) :: nbre1, nbre2, nbre3, nbre4
    parameter(nbre1=2, nbre2=2, nbre3=1, nbre4=2)
!
    integer(kind=8) :: nbr11, nbr12, nbr13, nbre5, nbre6
    parameter(nbr11=1, nbr12=3, nbr13=4, nbre5=4, nbre6=5)
!
    real(kind=8) :: valre1(nbre1), valre2(nbre2), valre3(nbre3), valre4(nbre4), valr11(nbr11)
    real(kind=8) :: valr12(nbr12), valre5(nbre5), valr13(nbr13), valre6(nbre6)
!
    integer(kind=8) :: codme1(nbre1), codme2(nbre2), codme3(nbre3), codme4(nbre4), codm11(nbr11)
    integer(kind=8) :: codm12(nbr12), codm13(nbr13), codme5(nbre5), codme6(nbre6)
!
    character(len=16) :: nomre1(nbre1), nomre2(nbre2), nomre3(nbre3), nomre4(nbre4), nomr11(nbr11)
    character(len=16) :: nomr12(nbr12), nomr13(nbr13), nomre5(nbre5), nomre6(nbre6)
    character(len=8) :: valk(2)
!
    aster_logical :: yapr, yaro
    type(THM_DS) :: ds_thm
!
    data nomre1/'RHO', 'BIOT_COEF'/
    data nomr13/'RHO', 'BIOT_L', 'BIOT_N', 'BIOT_T'/
    data nomr11/'PERM_IN'/
    data nomr12/'PERMIN_L', 'PERMIN_N', 'PERMIN_T'/
    data nomre2/'RHO', 'VISC'/
    data nomre3/'PORO'/
    data nomre4/'E', 'NU'/
    data nomre5/'E_L', 'E_N', 'NU_LT', 'NU_LN'/
    data nomre6/'E_L', 'E_T', 'NU_LT', 'NU_LN', 'NU_TN'/
!
! ----------------------------------------------------------------------
100 format(a, ' :', (6(1x, 1pe17.10)))
! ----------------------------------------------------------------------
! 1 -------------- GESTION DES DONNEES ---------------------------------
! ----------------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
! ------------------------------------------------------------------
!
    ovfl = r8miem()
!
! - Get all parameters for current element
!
    call thmGetElemPara(ds_thm, l_axi, &
                        type_elem, inte_type, ndim, &
                        mecani, press1, press2, tempe, second, &
                        dimdep, dimdef, dimcon, dimuel, &
                        nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                        nno, nnos, &
                        npi, npg, &
                        jv_poids, jv_func, jv_dfunc, &
                        jv_poids2, jv_func2, jv_dfunc2, &
                        jv_gano)
!
! =====================================================================
! B. --- DETERMINATION DES VARIABLES CARACTERISANT LE MILIEU ----------
! =====================================================================
    addeme = mecani(2)
    adcome = mecani(3)
    addep1 = press1(3)
    adcp11 = press1(4)
    addep2 = press2(3)
    addete = tempe(2)
    adsip = adcp11-adcome-5
    adde2nd = second(2)
!
! =====================================================================
! C. --- RECUPERATION DES DONNEES NECESSAIRES AU CALCUL ---------------
! =====================================================================
!--------------------------------------------------------------------
! 1. EVENTUELS PARAMETRES TEMPORELS
!--------------------------------------------------------------------
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itab(1))
    if (iret .eq. 0) then
        instpm(1) = zr(itab(1))
        deltat = zr(itab(1)+1)
        theta = zr(itab(1)+2)
        instpm(2) = instpm(1)-deltat
    else
        call utmess('F', 'INDICATEUR_11')
    end if
!--------------------------------------------------------------------
! 2. RECUPERATION DE LA GEOMETRIE, DU MATERIAU ET DES CHAMPS LOCAUX
!--------------------------------------------------------------------
!
! 2.1. GEOMETRIE (IGEOM), MATERIAU (IMATE) :
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
!
! 2.2. LES DEPLACEMENTS A L'INSTANT COURANT
!      1. TOUJOURS A L'INSTANT COURANT --> IDEPLP
!      2. SI TRANSITOIRE, A L'INSTANT PRECEDENT --> IDEPLM
!
    call jevech('PDEPLAR', 'L', ideplp)
!
    call jevech('PDEPLMR', 'L', ideplm)

!
! 2.3. CONTRAINTES AUX NOEUDS PAR ELEMENTS A L'INSTANT ACTUEL
!      1. TOUJOURS A L'INSTANT ACTUEL --> ISIENP
!      2. SI TRANSITOIRE, A L'INSTANT PRECEDENT --> ISIENM
!
    call tecach('ONO', 'PCONTNO', 'L', iret, nval=3, &
                itab=itab)
    isienp = itab(1)
    nbcmp = itab(2)/nno
    call tecach('ONO', 'PCONTNM', 'L', iret, nval=3, &
                itab=itab)
    isienm = itab(1)
!
! 2.4. CARTES DE PESANTEUR ET ROTATION
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
! 2.5. LES FORCES VOLUMIQUES EVENTUELLES :
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
!--------------------------------------------------------------------
! 3. RECHERCHE DES VALEURS NECESSAIRES AU CALCUL DE L'INDICATEUR
!     . COEFFICIENT DE BIOT
!     . MASSE VOLUMIQUE HOMOGENEISEE RHOHOM
!     . MASSE VOLUMIQUE DU LIQUIDE RHOLIQ
!     . CONDUCTIVITE HYDRAULIQUE
!--------------------------------------------------------------------
    nompar(1) = 'INST'
    valres(1) = instpm(1)
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_DIFFU', 1, nompar, [valres], &
                nbre1, nomre1, valre1, codme1, 0)
!
    if (codme1(1) .eq. 0 .and. codme1(2) .eq. 0) then
        rhohom = valre1(1)
        biot = valre1(2)
    else if (codme1(2) .eq. 1) then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THM_DIFFU', 1, nompar, [valres], &
                    nbr13, nomr13, valr13, codm13, 0)
        if ((codm13(1) .eq. 0) .and. (codm13(2) .eq. 0) .and. (codm13(3) .eq. 0)) then
            rhohom = valr13(1)
            biot = sqrt(valr13(2)**2+valr13(3)**2)
        elseif ((codm13(1) .eq. 0) .and. (codm13(2) .eq. 0) .and. (codm13(4) .eq. 0)) then
            rhohom = valr13(1)
            biot = sqrt(valr13(2)**2+valr13(4)**2)
        else
            ASSERT(.false.)
        end if
    else
        call utmess('F', 'ELEMENTS4_78', sk=nomre1(1)//nomre1(2))
    end if
!
! ON RECUPERE LA PERMEABILITE INTRINSEQUE
!
! => PERMIN SI ISOTROPE
! => PERMIN_X,PERMIN_Y ET PERMIN_Z SINON
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_DIFFU', 1, nompar, [valres], &
                nbr11, nomr11, valr11, codm11, 0)
!
    if (codm11(1) .eq. 0) then
        permin = valr11(1)
    else if (codm11(1) .eq. 1) then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THM_DIFFU', 1, nompar, [valres], &
                    nbr12, nomr12, valr12, codm12, 0)
        if ((codm12(1) .eq. 0) .and. (codm12(2) .eq. 0)) then
            permin = sqrt(valr12(1)**2+valr12(2)**2+valr12(1)**2)
        else if ((codm12(1) .eq. 0) .and. (codm12(3) .eq. 0)) then
            permin = sqrt(valr12(1)**2+valr12(3)**2)
        end if
    else
        call utmess('F', 'ELEMENTS4_78', sk=nomr11(1))
    end if
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_LIQU', 1, nompar, [valres], &
                nbre2, nomre2, valre2, codme2, 1)
!
    if ((codme2(1) .eq. 0) .and. (codme2(2) .eq. 0)) then
        rholiq = valre2(1)
        viscli = valre2(2)
    else
        call utmess('F', 'ELEMENTS4_69', sk=nomre2(1)//nomre2(2))
    end if
!
    if (permin .gt. ovfl) then
        unsurk = viscli/permin
    else
        call utmess('F', 'INDICATEUR_20')
    end if
!
!--------------------------------------------------------------------
! 4. SI INSTATIONNAIRE, ON RECUPERE DES COEFFICIENTS SUPPLEMENTAIRES
!     . MODULE DE BIOT
!     . MODULE DE YOUNG
!--------------------------------------------------------------------
!

!
! 4.1. RECHERCHE DE LA POROSITE INITIALE
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'THM_INIT', 1, nompar, [valres], &
                nbre3, nomre3, valre3, codme3, 1)
!
    if (codme3(1) .eq. 0) then
        porosi = valre3(1)
    else
        call utmess('F', 'ELEMENTS4_70', sk=nomre3(1))
    end if
!
! 4.2. RECHERCHE DU COEFFICIENT DE POISSON ET DU MODULE DE YOUNG
!
    call rcvalb(fami, kpg, spt, poum, zi(imate), &
                ' ', 'ELAS', 1, nompar, [valres], &
                nbre4, nomre4, valre4, codme4, 0)
!
    if ((codme4(1) .eq. 0) .and. (codme4(2) .eq. 0)) then
        cyoung = valre4(1)
        poisso = valre4(2)
    else if ((codme4(1) .eq. 1) .and. (codme4(2) .eq. 1)) then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'ELAS_ISTR', 1, nompar, [valres], &
                    nbre5, nomre5, valre5, codme5, 0)
        if ((codme5(1) .eq. 0) .and. (codme5(2) .eq. 0) .and. (codme5(3) .eq. 0) .and. &
            (codme5(4) .eq. 0)) then
            cyoung = sqrt(valre5(1)**2+valre5(2)**2)
            poisso = sqrt(valre5(3)**2+valre5(4)**2)
        else
            call rcvalb(fami, kpg, spt, poum, zi(imate), &
                        ' ', 'ELAS_ORTH', 1, nompar, [valres], &
                        nbre6, nomre6, valre6, codme6, 0)
            if ((codme6(1) .eq. 0) .and. (codme6(2) .eq. 0) .and. (codme6(3) .eq. 0) .and. &
                (codme6(4) .eq. 0)) then
                cyoung = sqrt(valre6(1)**2+valre6(2)**2)
                poisso = sqrt(valre6(3)**2+valre6(4)**2)
            end if
!
        end if
    else
        call utmess('F', 'ELEMENTS4_71', sk=nomre4(1)//nomre4(2))
    end if
!
! 4.4. ON CALCULE L'INVERSE DU MODULE DE BIOT
!
    if (cyoung .gt. ovfl) then
        unsurm = 3.d0*(biot-porosi)*(1.d0-biot)*(1.d0-2*poisso)/cyoung
!
    else
        call utmess('F', 'ELEMENTS4_67')
    end if

!
!--------------------------------------------------------------------
! 5. RECUPERATION DES GRANDEURS CARACTERISTIQUES
!--------------------------------------------------------------------
    call jevech('PGRDCA', 'L', igrdca)
    longc = zr(igrdca)
    presc = zr(igrdca+1)
!
    iaux = 0
    if (presc .le. ovfl) then
        iaux = 1
        valk(1) = 'pression'
    else if (longc .le. ovfl) then
        iaux = 1
        valk(1) = 'longueur'
    end if
!
    if (iaux .ne. 0) then
        call utmess('F', 'INDICATEUR_21', sk=valk(1))
    end if
!
! =====================================================================
! D. --- CALCUL DES INDICATEURS ---------------------------------------
! =====================================================================
!
!------------------------------------------------------------------
! 1. -------------- CALCUL DU PREMIER TERME DE L'ERREUR -----------
!------------------------------------------------------------------
!
! 2.1. --- CALCUL DU DIAMETRE HK DE LA MAILLE ----
!
    call uthk(nomte, zr(igeom), hk, ndim, niv)
!
! 2.2. --- CALCUL DE LA FORCE DE PESANTEUR ---
!
    if (yapr) then
        fpx = rhohom*zr(ipes)*zr(ipes+1)
        fpy = rhohom*zr(ipes)*zr(ipes+2)
    else
        fpx = 0.d0
        fpy = 0.d0
    end if
!GN      WRITE(IFM,100) 'P',FPX,FPY,FPZ
!
! 2.3. --- CALCUL DE LA FORCE DE ROTATION ---
!
    if (yaro) then
        call resrot(zr(irot), zr(igeom), zr(jv_func), rhohom, nno, &
                    npg, frx, fry)
    else
!
        do ipg = 1, npg
            frx(ipg) = 0.d0
            fry(ipg) = 0.d0
        end do
!
    end if
!
! 2.4. --- CALCUL DE LA FORCE VOLUMIQUE EVENTUELLE ---
!
    if (ifovr .ne. 0) then
        fovo(1) = zr(ifovr)
        fovo(2) = zr(ifovr+1)
!
    else if (ifovf .ne. 0) then
        nompar(1) = 'INST'
        r8bid3(1) = instpm(1)
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
!
    else
        fovo(1) = 0.d0
        fovo(2) = 0.d0
    end if
!
! 2.3. --- TERME VOLUMIQUE ---
!
    call erhmv2(ds_thm, l_axi, deltat, dimdep, dimdef, &
                nddl_meca, nddl_p1, nddl_p2, nddl_2nd, ndim, nno, &
                nnos, npg, nddls, nddlm, &
                dimuel, jv_poids, jv_func, jv_dfunc, jv_poids2, &
                jv_func2, jv_dfunc2, zr(igeom), fovo, zr(ideplp), &
                zr(ideplm), zr(isienp), zr(isienm), nbcmp, biot, &
                unsurm, fpx, fpy, frx, fry, &
                addeme, addep1, &
                addep2, addete, adde2nd, tm2h1v)
!
! ON ADIMENSIONNE LES INDICATEURS VOLUMIQUES
!
    admec = 1.d0/(presc**2*longc**ndim)
    tsivom = hk**2*admec*tm2h1v(1)
    tdevom = hk**2*admec*tm2h1v(2)
    adv1h = cyoung*unsurk*admec
    tsivoh = deltat*hk**2*adv1h*tm2h1v(3)

!
!------------------------------------------------------------------
! 2. CALCUL DES TERMES SURFACIQUES
!------------------------------------------------------------------
! 2.1. PHASE DE PREPARATION : ON RECUPERE LES ADRESSES NECESSAIRES
!                             AUX CALCULS
! -----------------------------------------------------------------
!
! ON RECUPERE LES ADRESSES
!    1. VOISINS
!    2. CHARGEMENTS DE TYPE FORCE_FACE
!    3. CHARGEMENTS DE TYPE PRES_REP
    call jevech('PVOISIN', 'L', ivois)
    call jevech('PFORCE', 'L', iref1)
    call jevech('PPRESS', 'L', iref2)
!
! RECHERCHE DES ADRESSES POUR LES CHARGES SUR LES SEGMENTS
!
    iagd = zi(iref1+4)
    iacmp = zi(iref1+5)
!
    iade2 = zi(iref2+4)
    iava2 = zi(iref2+5)
    iaptm2 = zi(iref2+6)
    if (iade2 .ne. 0) then
        igd2 = zi(iade2)
        ncmpm2 = zi(iacmp-1+igd2)
    end if
!
    iade3 = zi(iref2+8)
    iava3 = zi(iref2+9)
    iaptm3 = zi(iref2+10)
    if (iade3 .ne. 0) then
        igd3 = zi(iade3)
        ncmpm3 = zi(iacmp-1+igd3)
    end if
!
!------------------------------------------------------------------
! 2.2. CARACTERISATIONS DE LA MAILLE COURANTE
! -----------------------------------------------------------------
!
! TYPE DE LA MAILLE COURANTE
!
    typ = zi(ivois+7)
!
! ADRESSE DU VECTEUR TYPE MAILLE
!
    iatyma = zi(iref1+3)
    typema = zk8(iatyma-1+typ)
    form = typema(1:2)
!
! NOMBRE DE NOEUDS SOMMETS ET NOMBRE DE NOEUDS DES ARETES
!
    if (form .eq. 'TR') then
        nbs = 3
    else
        nbs = 4
    end if
!
    noeu = typema(5:5)
!
! EN THM, ON EST TOUJOURS SUR DES MAILLAGES QUADRATIQUES.
! (TRIA6 OU QUAD8 EN 2D)
! ON A DONC FORCEMENT NBNA = 3
!
    if (noeu .eq. '6' .or. noeu .eq. '8') then
        nbna = 3
    else
        ASSERT(.false.)
    end if
!
! CALCUL DE L'ORIENTATION DE LA MAILLE 2D
!     REMARQUE : ON APPELLE LE PROGRAMME GENERIQUE POUR LE PREMIER POINT
!                DE GAUSS, SACHANT QUE L'ORIENTATION NE DOIT PAS CHANGER
!
    jkp = 1
    call utjac(.true._1, zr(igeom), jkp, jv_dfunc, 0, &
               ibid, nno, orien)
!
!------------------------------------------------------------------
! 2.3. CALCUL DES TERMES LIES AUX ARETES
!------------------------------------------------------------------
! ON INITIALISE LES TERMES DE SAUT ET LES TERMES PROVENANT DES
! CONDITIONS AUX LIMITES (HYDRAULIQUE + MECANIQUE)
!
    do ii = 1, 3
        tm2h1b(ii) = 0.d0
        tm2h1s(ii) = 0.d0
    end do
!
! BOUCLE SUR LES ARETES : IMPLICITEMENT, ON NUMEROTE LOCALEMENT LES
!                         ARETES COMME LES NOEUDS SOMMETS
! . DONC LE PREMIER NOEUD DE L'ARETE A LE MEME NUMERO QUE L'ARETE : IFA
! . LE NOEUD SUIVANT EST IFA+1, SAUF SI ON EST SUR LA DERNIERE ARETE ;
!   LE NOEUD SUIVANT EST ALORS LE PREMIER, 1.
! . L'EVENTUEL NOEUD MILIEU EST LE 1ER NOEUD, DECALE DU NOMBRE DE NOEUDS
!   SOMMETS : IFA + NBS
!
    do ifa = 1, nbs
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
!
! --- CALCUL DES NORMALES, TANGENTES ET JACOBIENS AUX POINTS DE L'ARETE
!
            iaux = ifa
            call calnor('2D', zr(igeom), iaux, nbs, nbna, &
                        orien, ibid, ibid, noe, ibid, &
                        ibid, ibid, jaco, nx, ny, &
                        nz, tx, ty, hf)
!
! ------- SI L'ARRETE N'EST PAS SUR LA FRONTIERE DE LA STRUCTURE...
!         CALCUL DES TERMES DE SAUT A TRAVERS LES FACES INTERIEURES
!         DE LA MAILLE
!
            if (typmav(1:4) .eq. 'TRIA' .or. typmav(1:4) .eq. 'QUAD') then
!
                call erhms2(ifa, nbs, theta, jaco, &
                            nx, ny, zr(isienp), adsip, zr(isienm), &
                            nbcmp, typmav, zi(iref1), zi(iref2), ivois, &
                            tm2h1s)
!
! ------- SI L'ARRETE EST SUR LA FRONTIERE DE LA STRUCTURE...
!         CALCUL DES TERMES DE VERIFICATION DES CONDITIONS DE BORD
!
            else if (typmav(1:2) .eq. 'SE') then
!
                call erhmb2(ifa, nbs, ndim, theta, &
                            instpm, jaco, nx, ny, tx, &
                            ty, nbcmp, zr(igeom), ivois, zr(isienp), &
                            zr(isienm), adsip, iagd, zi(iref2), iade2, &
                            iava2, ncmpm2, iaptm2, iade3, iava3, &
                            ncmpm3, iaptm3, tm2h1b)
! ----------------------------------------------------------------
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
            if (niv .ge. 2) then
                write (ifm, 103) ifa, zi(ivois+ifa), typmav
103             format(i2, '-EME FACE DE NUMERO', i10, ' ==> TYPMAV = ', a)
                write (ifm, 100) 'TM2H1B', tm2h1b
                write (ifm, 100) 'TM2H1S', tm2h1s
            end if
!
        end if
!
    end do
!
! =====================================================================
! E. --- MISE EN FORME DES INDICATEURS --------------------------------
! =====================================================================
!
! ON ADIMENSIONNE LES INDICATEURS
!
    denomi = presc**2*longc**(ndim-2)*rholiq**2
    adhy0 = unsurk**2/denomi
    adhy1 = adhy0/(longc**2)
!
    tsisam = hk*admec*tm2h1s(1)
    tsibom = hk*admec*tm2h1b(1)
!
!
    call jevech('PERREM', 'L', ierrm)
!
    tdebom = hk*admec*tm2h1b(2)
    tdesam = hk*admec*tm2h1s(2)
!
    adhymd = deltat*hk*cyoung*permin/viscli*adhy1
    tsibsh = adhymd*tm2h1b(3)
    tsissh = adhymd*tm2h1s(3)
!
!
! ON STOCKE LES INDICATEURS
!
    call jevech('PERREUR', 'E', ierr)
!
    zr(ierr+1) = tsivom+tsisam+tsibom
    zr(ierr+2) = tdevom+tdebom+tdesam
    zr(ierr+3) = tsivoh+tsibsh+tsissh
!
    zr(ierr+4) = max(zr(ierrm+4), sqrt(zr(ierr+1)))
    zr(ierr+5) = zr(ierrm+5)+sqrt(zr(ierr+2))
    zr(ierr+6) = sqrt(zr(ierrm+6)**2+zr(ierr+3))
    zr(ierr) = zr(ierr+4)+zr(ierr+5)+zr(ierr+6)
!
    call jedema()
!
end subroutine
