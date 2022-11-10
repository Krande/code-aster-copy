! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine fluendo3d(xmat, sig0, sigf, deps, nstrs,&
                     var0, varf, nvari, nbelas3d, teta1,&
                     teta2, dt, epstf, ierr1,&
                     iso, mfr, end3d, fl3d, local,&
                     ndim, nmatbe2, iteflumax, sech,&
                     nvarbe)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!
    implicit none
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/hydramat3d.h"
#include "asterfort/thermat3d.h"
#include "asterfort/hydravar3d.h"
#include "asterfort/bwpw3d.h"
#include "asterfort/hydracomp3d.h"
#include "asterfort/dflufin3d.h"
#include "asterfort/conso3d.h"
#include "asterfort/bgpg3d.h"
#include "asterfort/dflueff3d.h"
#include "asterfort/endo3d.h"
#include "asterfort/utmess.h"
#include "asterfort/iniInt0.h"
#include "asterfort/iniReal0.h"
#include "asterfort/iniVect0.h"
#include "asterfort/iniMat0.h"
#include "asterfort/rgiRenfoStress.h"
#include "asterfort/getValVect.h"
#include "asterfort/setValVect.h"
#include "asterfort/getR6Mat6.h"
#include "asterfort/setR6Mat6.h"
#include "asterfort/getMat33Tab.h"
#include "asterfort/setMat33Tab.h"
#include "asterfort/getIntVect.h"
#include "asterfort/setIntVect.h"
#include "asterfort/setLogVect.h"
#include "asterfort/plasti3d.h"
#include "asterfort/tirViscoElas.h"
!   declaration des variables externes
    integer, intent(in) :: nstrs, nvari, nbelas3d, nmatbe2, ndim, mfr
    integer, intent(out) :: ierr1
    real(kind=8), intent(in) :: xmat(:)
    real(kind=8) :: var0(:), varf(:), epstf(6), sig0(6), sigf(:), deps(:)
    real(kind=8) :: dt, teta1, teta2, sech
    aster_logical :: iso
!   variable logique pour activer le fluage, l endo, le traitement local
    aster_logical, intent(in) :: end3d, fl3d, local
!    nombre maximal de sous-itération de fluage
    integer, intent(in) :: iteflumax, nvarbe
! ----------------------------------------------------------------------
    aster_logical :: is_ba, inputL(5)
    real(kind=8) :: beta, beta00, bg0, biotw, biotw00, ccmin0, dim3
    real(kind=8) :: cthp, cthp1, cthv, cthv1, cthvm, cwtauk0, cwtauk1
    real(kind=8) :: cwtaukm, delta, delta00, denom, denomin
    real(kind=8) :: dfin0, dfin1, dflu0, dflu1, dfmx, dhydr, dt1, dt80, dteta, dth0
    real(kind=8) :: dth00, dth1, dther, dtmaxi, dtmaxik, dtmaxim, dvrgi, dvw, ekdc
    real(kind=8) :: epc0, epc00, epleqc, epleqc0, epleqc00, vrgi
    real(kind=8) :: epleqc01, epser, epsklim1, epsklim2, epsm00, inputVR6(6,16)
    real(kind=8) :: ept, ept00, ept1, epsklim, hyd00, valr(2), inputR(25)
    real(kind=8) :: pg0, rc00, rc1, reduc1, ref00, ref1, rt00, rt1, outputR(6)
    real(kind=8) :: tauk00, tauk1, taum00, teta, treps, trepspg, gfr, epeqpc, errgf
    real(kind=8) :: umdt, vrgi00, vrgi1, vrgi2, vw, vw1, vw2, we0s
    real(kind=8) :: xflu, young00, xnsat, xnsat00, outputVR6(6,11)
    integer :: i, j, npas1, nt, ifour, iplalim, inputI(4), outputI(3)
!
    integer :: ngf
!   tableau pour la resolution des systemes lineaires
    real(kind=8), pointer :: X(:) => null()
    real(kind=8), pointer :: B(:) => null()
    real(kind=8), allocatable :: A(:, :)
    integer, pointer :: ipzero(:) => null()
!   donnes pour test fluage3d
    real(kind=8) :: epse06(6), epsk06(6), epsm06(6), sig06(6), phi0, we0, taum1
    real(kind=8) :: epse16(6), epsk16(6), epsm16(6), sig16(6)
    real(kind=8) :: we1, psik, epsm11
    real(kind=8) :: souplesse66(6, 6), raideur66(6, 6)
    real(kind=8) :: deps6(6), theta, theta1, heta, inputMat33(3,3,3), outputMat33(3,3,4)
!   theta : theta methode pour Euler, heta: taux de variation maxi
!   des deformations de fluage  lors d une sous iteration
    parameter(theta=0.5d0,heta=0.1d0)
    real(kind=8) :: young, nu, nu00, lambda, mu
    real(kind=8) :: rt33(3, 3), rtg33(3, 3), ref33(3, 3)
!   continuite pour as3d
    real(kind=8) :: rt, pglim, bg, phivg, mg, pg, ref, rc, poro, alat, kgel
!   deformation plastiques
    real(kind=8) :: epspt6(6), epspg6(6), epspc6(6)
    real(kind=8) :: epspt600(6), epspt60(6), epspg600(6), epspg60(6), epspc600(6), epspc60(6)
    real(kind=8) :: deltam, avean
!   derivees de la pression / deformation anelastiques
    real(kind=8) :: dpg_depsa6(6), dpg_depspg6(6)
!   sigf6 : contrainte effective + gel
!   deps6r : increment de deformation pour le retour radial
    real(kind=8) :: sigf6(6), deps6r(6)
!   ipla : compteur de sous iteration plastique
!   err1 : code d erreur des sous programmes
    integer :: ipla, err1
!   hydratation actuelle et seuil
    real(kind=8) :: hydr, hyds
!   deplacements imposes reduits
    real(kind=8) :: deps6r2(6)
!   complement a un de biot gel, inverse du potentiel de fluage
!   dissipation avant hydratation ( a hyd0)
    real(kind=8) :: phi00, hyd0, nrjm, sfld, mvgn
!   contrainte elastique de l etage de Kelvin
    real(kind=8) :: sigke06(6), sigke16(6)
!   porosite
    real(kind=8) :: poro2, poro1, dporo
!   endommagement capillaire du a leau
    real(kind=8) :: pw, bw 
!   resultante des contraintes intraporeuses
    real(kind=8) :: sigp
!   endommagement micro mecanique global
!   indicateur de premier pas
    aster_logical :: ppas
!   CWtauk : coeff pour la prise en compte du fluage de dessiccation sur tauk
!   eprg00 : deformation caracteristique pour l ecrouissage et l endo de rgi
!   et pour l endo de traction
!   gft :energie de fissuration par traction directe
!   pgmax : Pression gel max atteinte
!   srw : degré de saturation
    real(kind=8) :: CWtauk, eprg00, gft00, gft, pgmax, srw
!   sigf6d : contraintes endommagees
    real(kind=8) :: sigf6d(6)
!   facteur de concentration de contrainte pour le gel et module d ecrouissage de la pression
!   avec la deformation permanente de rgi et le facteur devolution de lecrouissage
    real(kind=8) :: krgi00, krgi, hplg, hpev
!   histoire des surcharges capillaires dues au variations hydriques sous charge
    real(kind=8) :: dsw6(6), bw0, pw0
!   influence de la consolidation spherique
!   coeff de consolidation anisotrope
    real(kind=8) :: cc03(3), vcc33(3, 3), vcc33t(3, 3)
!   endommagement isotrope de fluage (asymptotique et effectif)
    real(kind=8) :: ccmax0, ccmax1, dfl00, dfl0, cmp0, CWp
!   temperatures de reference et de seuil
    real(kind=8) :: tetas, tetar
!   coeff THM  / fluage debut de pas
    real(kind=8) :: CWp0, CthP0, Cthv0, dsw06(6)
!   endommagements et ouvertures de fissures dans leur base principale
    real(kind=8) :: dt3(3), dr3(3), dgt3(3), dgc3(3), dc, wl3(3)
!   traitement endommagement isotrope prepic
    real(kind=8) :: xmt, dtr
    aster_logical :: dtiso
!   avancement de la reaction de gonflement interne (rag par exple)
    real(kind=8) :: vrgi0, taar, nrjg, srsrag, aar0, aar1, trag, vrag00
    real(kind=8) :: vdef00, def1,tdef, nrjd, def0, srsdef, CNa
    real(kind=8) :: nrjp, ttrd, tfid, ttdd, tdid, exmd, exnd, cnab, cnak, ssad
    real(kind=8) :: At, St, M1, E1, M2, E2, AtF, StF, M1F, E1F, M2F, E2F, ttkf, nrjf
!   valuer fluage debut de pas
    real(kind=8) :: epsk006(6), epsm006(6)
!   Def thermiques transitoires (si Dtheta sous charge)
    real(kind=8) :: deps6r3(6), ett600(6), ett60(6)
!   14mai2015: ouvertures de fissures du pas precedent
    real(kind=8) :: wplt6(6), wplt06(6), wplt006(6),wpltx6(6), wpltx06(6), wpltx006(6)
    real(kind=8) :: wpl3(3), vwpl33(3, 3), vwpl33t(3, 3), wplx3(3), vwplx33(3, 3), vwplx33t(3, 3)
    real(kind=8) :: epstf6(6), sigmf6(6) 
!-----------------------------------------------------------------------
    if (nmatbe2 .eq. 0) then
        ngf = 22
        is_ba = ASTER_FALSE
        iplalim = 4
    else
!       RGI_BETON_BA
        ngf = 65
        is_ba = ASTER_TRUE
        iplalim = 1000
    endif
    AS_ALLOCATE(vr=X, size=ngf)
    AS_ALLOCATE(vr=B, size=ngf)
    ALLOCATE(A(ngf,(ngf+1)))
    AS_ALLOCATE(vi=ipzero, size=ngf)
!
    call iniInt0(err1, npas1, nt, ifour)
    call iniReal0(vrgi2, epleqc, epleqc0, epleqc00, deltam, avean,&
                  epleqc01, we0s, epeqpc, we0, lambda, mu, bg, mg, pg)
    call iniReal0(pw, bw, sigp, cwtauk, srw, bw0, pw0, dfl00, dfl0,&
                  cwp, dc, xmt, dtr, vrgi0, aar0, aar1, def1, def0)
    call iniReal0(At, St, M1, E1, M2, E2, AtF, StF, M1F, E1F, M2F, E2F) 
    call iniVect0(3, cc03, dt3, dr3, dgt3, dgc3, wl3)
    call iniVect0(6, deps6r2, deps6r3, deps6r, sigf6d, dsw6,&
                  epsk06, sig16, epsm16, epsm06, epse06)
    call iniVect0(6, epspg6, epspc6, epspt600, epspt60, epspg600, epspg60, &
                  epspc600,epspc60, dpg_depsa6, dpg_depspg6)
    call iniVect0(6, sigke06, sigke16,epsk006, epsm006, ett600, ett60,&
                   wplt6, wplt06, wplt006, wpltx6, wpltx06, wpltx006)
    call iniMat0(3, ref33, rt33, rtg33, vcc33, vcc33t)
    call iniMat0(6, souplesse66, raideur66)
    call iniVect0(ngf, X, B)
!
    A(:,:)=0.d0
    ipzero(:)=0
!   chargement des tailles si l endo est active
    ierr1=0
!
!    chargement des parametres materiaux : l hydratation est considéree en fin de pas  
!    pour ne pas avoir a recuperer celle du pas precedent  les exposants de De Shutter
!    sont  E   2/3  gf 1/2 rt 2/3 rc 1  ekfl 2/3  biot 1/2
!
!   hydr = hydratation, hyds = hydratation seuil, rt00 = resistances et pression limite 2/3
!   ref00 = seuil pour la refermeture des fissures
!   rc00  = resistance en compression-par cisaillement
!   delta00  = coeff drucker prager, beta00 = coeff dilatance de cisaillement
!   ept00 = deformation au pic de traction
!   hplg = taux d ecrouissage des phases effectives / RGI
!   phivg = vide accesible au gel, kgel = Rigidite gel matrice
!   gft00 = energie de fissuration en traction directe
!   epsm00 =  deformation caracteristique du potentiel de fluage
!   psik = raideur relative Kelvin / Young
!   xflu = endommagement maximum par fluage
!   tauk00 = temps caracteristique pour Kelvin
!   taum00 = temps caracteristique pour Maxwell
    call getValVect(xmat, hydr, hyds, rt00, ref00, rc00, delta00, beta00,&
                 ept00, hplg, phivg, kgel, gft00, epsm00, psik,&
                 xflu, tauk00, taum00, ind1=nbelas3d+1)
!     stockage dans une variable interne pour avoir la vitesse
    if (abs(var0(64)-1.d0) .ge. r8prem()) then
!       au 1er passage on suppose hyd0=hydr
        hyd0=hydr
!       on initialise var03d en cas de sous incrementation par fluage
        var0(48)=hyd0
    else
        hyd0=var0(48)
    end if
    dhydr=hydr-hyd0
    varf(48)=hydr
!     stockage hydratation initiale pour calcul final de pression
    hyd00=hyd0
!
!     Module d'Young et coefficient de Poisson
    if (.not. is_ba) then
        call getValVect(xmat, young00, nu00, ind1=1)
    else
        call getValVect(xmat, young00, nu00, ind1=nbelas3d+58)
    endif
!     deformation de reference  pour le fluage prise aqu 1/3 de rc
    epser=(rc00/young00)/3.d0
!     verif validite de la dilatance
    if (beta00 .gt. dsqrt(3.d0)) then
        call utmess('E', 'COMPOR3_23', sr=dsqrt(3.d0))
        ierr1=1
        go to 999
    end if
!
!   nrjm = activation du potentiel de fluage (dmax1 pour le cas endo seul, utile dans thermat3d)
!   dt80 = endommagement thermique a 80°C
!   donnees pour le calcul hydrique fin de pas
!       biotw00 =  coeff de biot pour l eau
!       xnsat00 = module de biot pour le non sature
!   poro2 = porosite ou grandeur permettant de passer du champ sech au degrès de saturation
!   vrag00 = volume maximal de rgi comprenant le volume non effectif
!   nrjf = energie de fixation des alus en HG (RSI)
!   sfld = contrainte caracteristique pour l endo capillaire
!   mvgn = exposant de Van Genuchten pour la pression capillaire
!   epc0 = deformation au pic de compression
!   ekdc = deformation cracateristique endo de compression
!   eprg00 = deformation caracteristique pour l endo de rgi
!   gfr = deformation caracteristique pour l endo de traction
!   alat = parametre induisant periode de latence initiale 
!   krgi00 = coeff de concentration de contrainte des RGI
!   tetar = temperature de reference des parametres de fluage (celsius)
!   tetas = temperature seuil pour l endo thermique (celsisu)
    call getValVect(xmat, nrjm, dt80, biotw00, xnsat00, poro2, vrag00,&
                 nrjf,  sfld, mvgn, epc0,ekdc, eprg00, gfr, alat, krgi00,&
                 tetar, tetas, ind1= nbelas3d+18)
!   volume d eau pour le non sature
    vw2=sech
!   initialisation des variables internes associee a la saturation si premier pas
    if (abs(var0(64)-1.d0) .ge. r8prem()) then
        var0(58)=vw2
        var0(59)=poro2
    end if
!   eps pic comp ne peut pas etre inferieur à rc/E
    epc00=dmax1(epc0,3.d0*epser)
!
!   dfmx = endommagement maximum par fluage
!   taar = temps caracterisqtique de la rgi à tref
!   nrjg = nrj d'activation de la rgi
!   srsrag = seuil de saturation minimal pour avoir la rgi
!   trag = temperature de reference des parametres de RAG (celsius)
!   dim3 = dimension 3 en 2D
!   tdef = temps cracteristique pour la def
!   nrjp = energie d activation de precipitation de la def
!   srsdef = seuil de saturation pour declancher la def
!   vdef00 = quantite maximale de def pouvant etre realise
!   cna = teneur en alcalin pour la def
!   ssad = rapport molaire S03 / Al2O3 du ciment
!   cnak = concentration caracteristique pour les lois de couplage de la def
!   cnab = concetration en alcalin de blocage de la def
!   exnd = exposant de la loi de couplage temperature de dissolution alcalins
!   exmd = exposant de la loi de couplage precipitation def alcalins
!   ttdd = temperature de reference pour la dissolution de la def
    call getValVect(xmat, dfmx, taar, nrjg, srsrag,trag, dim3,&
                 tdef,  nrjp, srsdef, vdef00,cna, ssad, cnak,&
                 cnab, exnd, exmd, ttdd, ind1=nbelas3d+35)
!   tdid = temps caracteristique pour la dissolution de l ettringite primaire
!   tfid = temps caracteristique pour la fixation des aluminiums en temperature
!   nrjd = energie d activation des processus de dissolution des phases primaires
!   ttrd = temperature de reference pour la precipitation de la def
!   ttkf = température seuil de fixation des alus en HG (RSI)
!   hpev = ratio pour levolution du module decrouissage
    call getValVect(xmat, tdid, tfid, nrjd, ttrd,&
                 ttkf, hpev, ind1=nbelas3d+52)
!   on peut traiter l'endommagement pre pic iso de traction si ept > Rt/E
    ept00=dmax1(rt00/young00,ept00)
!
!   indicateur de premier passage pour hydracomp3d
    if (abs(var0(64)-1.d0) .ge. r8prem()) then
        ppas=.true.
    else
        ppas=.false.
    end if
!
!   chargement de l increment de deformation imposee
    deps6(:)=0.d0
    epstf6(:)=0.d0
    deps6(1:nstrs)=deps(1:nstrs)
    epstf6(1:nstrs)=epstf(1:nstrs)
!   passage en epsilon
    deps6(4:6)=0.5d0*deps6(4:6)
    epstf6(4:6)=0.5d0*epstf6(4:6) 
!     remarque si rt rtg ref et rc dependent de  l ecrouissage
!     il faut les actualiser en fonction de l ecrouissage avant de
!     de passer dans hydramat, il faut egalement que dra_dl soit
!     parfaitement compatible avec le recalcul en fonction des deformations
!     plastiques (ou bien rajouter des vari pour stocker les resistances)
    rc1=rc00
!     l hydratation n a pas d influence sur la deofrmation plastique
!     caracteristique de cisaillement
    rt1=rt00
    ept1=ept00
    ref1=ref00
!     parametres materiau dependant eventuellement de l hydratation et de l endo de fluage
!     pour le 1er passage la reduction des resistances par fluage est
!     negligee car non utilise pour le tir visco elastique
    call hydramat3d(hyd0, hydr, hyds, young00, young,&
                    nu00, nu, rt1, rt, ref1,&
                    ref, rc1, rc, delta00, delta,&
                    beta00, beta, gft00, gft, ept1,&
                    ept, pglim, epsm00, epsm11, xnsat00,&
                    xnsat, biotw00, biotw, krgi00, krgi,&
                    iso, lambda, mu, rt33, rtg33,&
                    ref33, raideur66, souplesse66, xmt, dtiso,&
                    err1)
! - influence de la temperature sur les parametres  materiau et
!   actualisation de l endo thermique initial
    dth00=var0(49)
    call thermat3d(teta1, nrjm, tetas, tetar, DT80,&
                   dth00, dth0, CTHp0, CTHv0)
! - endommagement thermique en fin de pas
    call thermat3d(teta2, nrjm, tetas, tetar, DT80,&
                   dth0, dth1, CTHp1, CTHv1)
    varf(49)=dth1
! - chargement des variables internes du fluage (etat du squelette solide)
    do i = 1, 6
        epsk006(i)=var0(i+6)
        epsm006(i)=var0(i+12)
        sig06(i)=var0(i+18)
        sigke06(i)=var0(i+49)
        dsw06(i)=var0(73+i)
    end do
!   phi00 = dissipation visqueuse
!   dfl00 = endommagement effectif par fluage
!   dth00 = endommagement thermique
!   epleqc00 = deformation plastique equivallente de cisaillement
!   bw0, pw0 : pression capillaire
!   bg0, pg0 : pression RGI
    call getValVect(var0, phi00, dfl00, dth00, epleqc00, bw0, pw0,&
                    bg0, pg0, vectInd = [25,27,49,67,66,56,65,61])
!   tenseurs de deformations plastique et surpression capillaire
    do j = 1, 6
        epspt600(j)=var0(29+j)
        epspg600(j)=var0(35+j)
        epspc600(j)=var0(41+j)
        ett600(j)=var0(96+j)
        wplt006(j)=var0(102+j)
        wpltx006(j)=var0(67+j)
    end do
!
! - influence du degre d hydratation sur les variables internes
    call hydravar3d(hyd0, hydr, hyds, phi00, phi0,&
                    dth00, dth0, epleqc00, epleqc0, epspt600,&
                    epspt60, epspg600, epspg60, epspc600, epspc60,&
                    epsk006, epsk06, epsm006, epsm06, dfl00,&
                    dfl0, ett600, ett60, wplt006, wplt06,&
                    wpltx006, wpltx06)
! - calcul de la pression capillaire due a leau en fin de pas
    if (fl3d) then
        call bwpw3d(mfr, biotw, poro2, vw2, xnsat, mvgn, pw, bw, srw)
!       modif eventuelle des viscosites en fonction de srw
        CWtauk1=1.d0/srw
        call setValVect(varf, pw, bw, CWtauk1, vectInd=[56,66,57])
    end if
! - reevaluation de la deformation compatible avec l etat
!   de contrainte debut de pas(pour evaluer la sous incrementation)
    call hydracomp3d(we0, we0s, epse06, souplesse66, sig06,&
                     deps6, deps6r, sigke06, epsk06, psik,&
                     fl3d)
! - reevaluation des coeffs de consolidation en debut de pas
    if (fl3d) then
!       effet de l eau
        CWp0=var0(58)/var0(59)
        CWtauk0=var0(59)/var0(58)
!       effet l endo de fluage sur le potentiel de fluage
        call dflufin3d(sig06, bw0, pw0, bg0, pg0,&
                       dsw06, delta, rc, xflu, dfin0,&
                       CMp0, dfmx)
!       effet de la temperature sur le potentiel
        call conso3d(epsm11, epser, ccmin0, ccmax0, epsm06,&
                     epse06, cc03, vcc33, vcc33t, CWp0,&
                     CMp0, CthP0, Cthv0)
! ----  determination de la subdivision eventuelle du pas de temps
        dtmaxi=dt
        do i = 1, 6
!           extremes de epsk
            epsklim1=epse06(i)/psik
            epsklim2=(epse06(i)+deps6r(i))/psik
!            valeur maximale possible de l increment de deformation de kelvin
            if (dabs(epsklim1) .gt. dabs(epsklim2)) then
                epsklim=epsklim1
            else
                epsklim=epsklim2
            end if
!           cas ou epsklim est faible
            if (dabs(epsklim) .gt. 1.d-6) then
                denom=dabs(1.d0-epsk06(i)/epsklim)
            else
                denom=dabs(1.d0-epsk06(i)/1.d-6)
            end if
!           comparaison avec la valeur permettant de franchir dt en 1/heta pas
            cwtaukm=(CWtauk0+cwtauk1)/2.d0
            cthvm=(cthv0+cthv1)/2.d0
            if (abs(dt) .ge. r8prem()) then
                denomin=heta*tauk00*CWtaukm/CTHVm/dt
            else
                denomin=1.d0
            end if
!            comparaison avec les autres composantes de deformation
            if (denom .le. denomin) then
                denom=denomin
            end if
            dtmaxik=tauk00*CWtaukm/CTHVm/denom
!            cas de la deformation de maxwell
            dtmaxim=heta*(taum00*ccmin0)
!            choix entre condition maxwell et kelvin
            dtmaxi=dmin1(dtmaxi,dtmaxim,dtmaxik)
        end do
!       *** subdivision du pas si necessaire  **************************
        if (dtmaxi .lt. dt) then
            npas1=int(dt/dtmaxi)+1
            if (iteflumax .gt. 0.d0) then
                npas1=min(npas1,int(iteflumax))
            end if
        else
            npas1=1
        end if
    else
!        pas de fluage
        npas1=1
    end if
!      coeff de reduction de l increment
    reduc1=1.d0/dble(npas1)
!      reduction des pas d hydratation
    dhydr=reduc1*(hydr-hyd0)
!      pas de temps reduit
    dt1=dt*reduc1
!      increment de temperature reduit
    dteta=(teta2-teta1)*reduc1
!      initialisation de la temperature debut de pas
    teta=teta1
!      increment de volume d eau reduit
    vw1=var0(58)
    dvw=(vw2-vw1)*reduc1
!      increment de poro capillaire reduite
    poro1=var0(59)
    dporo=(poro2-poro1)*reduc1
!      increment des potentiels de rgi
    vrgi1=var0(60)
    dvrgi=(vrgi2-vrgi1)*reduc1
!      increments de deformations reduits, on reduit deps6 et non deps6r
!      car hydracomp est reapplique plus bas
    if (npas1 .ne. 1) then
        deps6r(:)=reduc1*deps6(:)
    else
        deps6r(:)=deps6(:)
    end if
!***********************************************************************
!       debut du chargement sous-discretise sur le pas de temps
!***********************************************************************
    do nt = 1, npas1
!        chargement des variables internes
        do i = 1, 6
            epsk006(i)=var0(i+6)
            epsm006(i)=var0(i+12)
            sig06(i)=var0(i+18)
            sigke06(i)=var0(i+49)
            dsw06(i)=var0(73+i)
        end do
!
!       recuperation de l endommagement par fluage dfl00 et thermique dth00
        call getValVect(var0, bw0, pw0, bg0, pg0, phi00, dth00, hyd0, vw,&
                        poro, vrgi00, epleqc00, dfl00,&
                        vectInd = [66,56,65,61,25,49,48,58,59,60,67,27])
!        actualisation du degre d hydratation
        hydr=hyd0+dhydr
!        actualisation de la temperature
        teta=teta+dteta
!        actualisation volume deau capillaire
        vw=vw+dvw
!        actualisation de la porosite capillaire
        poro=poro+dporo
!        actualisation du volume pour les rgi
        vrgi00=vrgi00+dvrgi
        call setValVect(varf, hydr, vw, poro, vrgi00, vectInd=[48,58, 59,60])
!
!        recuperation de la deformation maximale de traction
        do j = 1, 6
            epspt600(j)=var0(29+j)
            epspg600(j)=var0(35+j)
            epspc600(j)=var0(41+j)
            ett600(j)=var0(96+j)
            wplt006(j)=var0(102+j)
            wpltx006(j)=var0(67+j)
        end do
!
!       prise en compte  hydratation intermediaire si sous-increment
        call hydravar3d(hyd0, hydr, hyds, phi00, phi0,&
                        dth00, dth0, epleqc00, epleqc0, epspt600,&
                        epspt60, epspg600, epspg60, epspc600, epspc60,&
                        epsk006, epsk06, epsm006, epsm06, dfl00,&
                        dfl0, ett600, ett60, wplt006, wplt06,&
                        wpltx006, wpltx06)
!
!        stockage de la valeur de la deformation plastique cumulee
!        modifiee par l hydratation pour mise a jour incrementale
        epleqc01=epleqc0
!        effet de l endommagement de fluage et de l'hydratation sur les
!        parametres materiau
        call hydramat3d(hyd0, hydr, hyds, young00, young,&
                        nu00, nu, rt1, rt, ref1,&
                        ref, rc1, rc, delta00, delta,&
                        beta00, beta, gft00, gft, ept1,&
                        ept, pglim, epsm00, epsm11, xnsat00,&
                        xnsat, biotw00, biotw, krgi00, krgi,&
                        iso, lambda, mu, rt33, rtg33,&
                        ref33, raideur66, souplesse66, xmt, dtiso,&
                        err1)
!
!        influence de la temperature sur les parametres du materiau
!        et calcul de l endommagement thermique
        call thermat3d(teta, nrjm, tetas, tetar, DT80,&
                       dth1, DTHER, CTHP, CTHV)
        varf(49)=dth1
!        calcul de la pression capillaire due a leau
        if (fl3d) then
            call bwpw3d(mfr, biotw, poro, vw, xnsat,&
                        mvgn, pw, bw, srw)
            varf(56)=pw
            varf(66)=bw
!        reevaluation des coeffs de consolidation apres increment hydra
!        effet de l eau sur le potentiel de fluage
            Cwp=Srw
            CWtauk=1.d0/srw
            varf(57)=CWtauk
!
!        Modification eventuelle de la viscosite en fonction de Srw
            tauk1=tauk00*CWtauk/CTHV
!
!        la modif de la viscosite de Maxwell est comprise dans le coeff
!        de consolidation (cf. conso3d)
            taum1=taum00
        end if
!
!        compatibilite des anciennes contraintes avec nouveau materiau
!        deduction de la deformation de solidification de l increment
        call hydracomp3d(we0, we0s, epse06, souplesse66, sig06,&
                         deps6r, deps6r2, sigke06, epsk06, psik,&
                         fl3d)
!       coeff theta methode pour tir visco elastique
        theta1=theta
!
!***********************************************************************
!        tir visco elastique
!***********************************************************************
!
!       E.Cheignon : on devrait supprimer deps6r3 qui n'est pas différent
!       de deps6r2 (voir tirViscoElas)
        deps6r3(:) = deps6r2(:)
!
        call setValVect(inputR, delta, rc, epsm11, epser,  CWp, CthP, Cthv,&
                    dt1, theta1, tauk1, taum1, ind1=1)
!
        call setR6Mat6(inputVR6, dsw06, epsm06, sigke06, deps6r2,&
                   epsk06, epse06, sig06)
!
        call tirViscoElas(fl3d, var0, xmat, inputR, inputVR6, ngf,  &
                        deltam, avean, A, B, X, ipzero, &
                        epsk16, epsm16, epse16, &
                        sig16, sigke16, raideur66, we1)
!
!       E.Cheignon : sorti de tirViscoElas pour ne pas passer varf
        do i = 1, 6
!            varf(96+i)=ett61(i)
            varf(96+i)=0.d0
        end do
!
!       actualisation de la variation de volume total
        treps=var0(28)
        do i = 1, 3
            treps=treps+deps6r(i)
        end do
        varf(28)=treps
!
!       chargement des deformations plastiques du pas precedent
!       actualisee par l hydratation
        do j = 1, 6
!         traction
            epspt6(j)=epspt60(j)
!         rgi
            epspg6(j)=epspg60(j)
!         compression
            epspc6(j)=epspc60(j)
!         chargement des ouvertures de fissures du pas precedent
            wplt6(j)=wplt06(j)
!         chargement des ouvertures maxi de fissure          
            wpltx6(j)=wpltx06(j)
        end do
!       chargement de la def equivalente actualisee par l hydratation
        epleqc=epleqc0
!
!       recuperation de la variation volumique due a la rgi non
!       actualisee par l hydratation
        trepspg=var0(29)
!
!***********************************************************************
!       verification des criteres de plasticite et ecoulements
!***********************************************************************
!
        call setValVect(inputR, biotw, poro, vw, xnsat, pglim, teta, dt1,&
                    young00, nu00, rc, epleqc0, delta, beta, krgi,&
                    theta, epsm11, ind1=1)
        call setValVect(inputR, CWp, CthP, Cthv, tauk1, taum1, deltam, avean,&
                    phi0, epleqc01, srw, ind1=17)
!
        call setR6Mat6(inputVR6, sig0, epspt60, wplt06, wpltx06,&
                   sig06, epsm06, dpg_depsa6, dpg_depspg6, epsk16,&
                   epsm16, epse16, epspt6, sigke16, sig16, epspg6, epspc6)
        call setMat33Tab(inputMat33, ref33, rtg33, rt33)
        call setIntVect(inputI, iplalim, mfr, nstrs, ndim)
        call setLogVect(inputL, fl3d, ppas, iso, local, end3d)
        
        call plasti3d(xmat, inputR, inputVR6, inputMat33, inputI,&
                    inputL, var0, raideur66, souplesse66,&
                    A, B, X, ngf, varf, ipzero,&
                    outputR, outputVR6, outputMat33,&
                    outputI)
!
        call getValVect(var0, aar0, def0, E1, M1, E2, M2, At, St,&
                        vectInd=[62,63,97,98,99,100,101,102])
!
        call getValVect(varf, aar1, def1, E1f, M1f, E2f, M2f,&
                        Atf, Stf, vrgi, pgmax, pg, bg, trepspg,&
                        vectInd=[62,63,97,98,99,100,101,102,60,114,61,65,29])
        call getValVect(outputR, srw, epeqpc, dfin1, ccmax1, pw, bw,&
                        wpl3(1), wpl3(2), wpl3(3), wplx3(1), wplx3(2),&
                        wplx3(3), ind1=1)
        call getR6Mat6(outputVR6, dpg_depsa6, dpg_depspg6, epsk16,&
                   epsm16, epse16, epspt6, sigke16, sig16, epspg6,&
                   epspc6, dsw6)
        call getMat33Tab(outputMat33, vwpl33, vwpl33t, vwplx33, vwplx33t)
        call getIntVect(outputI, ifour, ierr1, ipla)
!
!       fin de la boucle de consistance visco-elasto-plastique
!       ****************************************************************
!       transfert des variables internes pour la sous iteration temporelle suivante
!       si le nbre de sous iteration locale le justifie 
        if (npas1 .gt. 1) var0(1:nvari)=varf(1:nvari)
!
    end do
!      fin de la boucle de discretisation du pas de temps
!***********************************************************************
!
!***********************************************************************
!       reevaluation de la pression de gel en fin de pas pour le calcul
!       des contraintes totales (effective+rgi)
    if (ipla .gt. 1) then
!          la calcul de l avancement de la rgi est inclus dans
!          le sous programme de calcul de pression
!          vrgi le volume effectif de gel pour le calcul de la pression
        vrgi0=vrgi
        pgmax=var0(114)
!          dt mis a zero pour forcer la reprise de vrgi0
        call bgpg3d(ppas, bg, pg, mg, vrgi,&
                    treps, trepspg, epspt6, epspc6, phivg,&
                    pglim, dpg_depsa6, dpg_depspg6, taar, nrjg,&
                    trag, aar0, srw, srsrag, teta,&
                    dt1, vrag00, aar1, tdef, nrjd,&
                    def0, srsdef, vdef00, def1, cna,&
                    nrjp, ttrd, tfid, ttdd, tdid,&
                    exmd, exnd, cnab, cnak, ssad,&
                    At, St, M1, E1, M2,&
                    E2, AtF, StF, M1F, E1F,&
                    M2F, E2F, vrgi0, ttkf, nrjf,&
                    alat, young00, nu00, kgel, pgmax)
    else
        call getValVect(varf, pg, bg, pgmax, vectInd=[61,65,114])
    end if
!       stockage de la pression RGI
    call setValVect(varf, pg, bg, pgmax, vectInd=[61,65,114])
!
! - endommagement de fluage
    dflu1=0.d0
    if (fl3d) then
        dflu0=var0(27)
        call dflueff3d(ccmax1, dflu0, dflu1, dfin1)
    end if
    varf(27)=dflu1
!
!
!***********************************************************************
!       contraintes dans solide et rgi en fin de pas avec
!       prise en compte de l endo thermique et de fluage
    umdt=(1.d0-dth1)*(1.d0-dflu1)
!       resultante des pressions intraporeuses RGI et Capillaire (depression)
    sigp=-bg*pg-bw*pw
!       effet sur la contrainte apparente en non sature
    do i = 1, 6
        if (i .le. 3) then
!               prise en compte de la pression rgi
            sigf6(i)=(sig16(i)+sigp+dsw6(i))*umdt
        else
            sigf6(i)=(sig16(i)+dsw6(i))*umdt
        end if
    end do
!
!***********************************************************************
!       prise en compte de l'endommagement mécanique
    if (end3d) then
!       chargement endo traction pre-pic
        dtr=var0(96)
!       chargement endo localisee pour condition de croissance
        do i = 1, 3
            dt3(i)=var0(79+i)
        end do
!       calcul des endommagements et ouvertures de fissures
        call endo3d(wpl3, vwpl33, vwpl33t, wplx3, vwplx33,&
                    vwplx33t, gft, gfr, iso, sigf6,&
                    sigf6d, rt33, ref33, souplesse66, epspg6,&
                    eprg00, a, b, x, ipzero,&
                    ngf, ekdc, epspc6, dt3, dr3,&
                    dgt3, dgc3, dc, wl3, xmt,&
                    dtiso, rt, dtr, dim3, ndim,&
                    ifour, epeqpc, young, ept, errgf)
!            stockage des endommagements de fissuration et ouverture
        varf(96)=dtr
        do i = 1, 3
!               endo de traction
            varf(79+i)=dt3(i)
!               endo refermeture
            varf(82+i)=dr3(i)
!               endo RGI en traction
            varf(85+i)=dgt3(i)
!               endo RGI en compression
            varf(88+i)=dgc3(i)
!               ouverture non visco elastique
            varf(91+i)=wl3(i)
        end do
!            endo de traction global
        varf(110)=1.d0-((1.d0-(varf(80)))*(1.d0-(varf(81)))*(1.d0-(varf(82))))
!            endo de traction de RGI global
        varf(112)=1.d0-((1.d0-(varf(86)))*(1.d0-(varf(87)))*(1.d0-(varf(88))))
!            endo de compression de RGI global
        varf(113)=1.d0-((1.d0-(varf(89)))*(1.d0-(varf(90)))*(1.d0-(varf(91))))
!            endommagement de compression
        varf(95)=dc
!            erreur de dissipation d'énergie en traction          
        varf(109)=errgf
!            traitement erreur endo
        if (err1 .eq. 1) then
            call utmess('E', 'COMPOR3_34')
            ierr1=1
            go to 999
        end if
!       contraintes totale dans la matrice apres endommagement 
        sigmf6(1:6)=sigf6d(1:6) 
    else
!       pas d endommagement
        sigf6d(1:6)=sigf6(1:6)
        sigmf6(1:6)=sigf6(1:6)
    end if
!
    if (is_ba) then
        call rgiRenfoStress(xmat, nbelas3d+nmatbe2+1, sigmf6, &
                          epstf6, epspt6, teta1, teta2, dt, ppas, theta,&
                          fl3d, end3d, wpl3, vwpl33, vwpl33t, dt3, dr3,&
                          ipzero, nvarbe, ngf, rc00, var0, varf, sigf6d, ierr1)
    endif
!
!   affectation dans le tableau de sortie des contraintes
    sigf(1:nstrs)=sigf6d(1:nstrs)
999 continue
!
    AS_DEALLOCATE(vr=X)
    AS_DEALLOCATE(vr=B)
    deallocate(A)
    AS_DEALLOCATE(vi=ipzero)
!
end subroutine
