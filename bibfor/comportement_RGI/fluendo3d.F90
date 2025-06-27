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
subroutine fluendo3d(xmat, sig0, sigf, deps, nstrs, &
                     var0, varf, nvari, nbelas3d, teta1, &
                     teta2, dt, epstf, ierr1, &
                     iso, mfr, end3d, fl3d, local, &
                     ndim, nmatbe2, iteflumax, sech, matrEndo, &
                     matdech)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!
    implicit none
#include "rgi_module.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/hydramat3d.h"
#include "asterfort/thermat3d.h"
#include "asterfort/hydravar3d.h"
#include "asterfort/bwpw3d.h"
#include "asterfort/hydracomp3d.h"
#include "asterfort/dflufin3d.h"
#include "asterfort/conso3d.h"
#include "asterfort/bgpg3d.h"
#include "asterfort/dflueff3d.h"
#include "asterfort/endotot.h"
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
    integer(kind=8), intent(in) :: nstrs, nvari, nbelas3d, nmatbe2, ndim, mfr
    integer(kind=8), intent(out) :: ierr1
    real(kind=8), intent(in) :: xmat(:)
    real(kind=8) :: var0(:), varf(:), epstf(6), sig0(6), sigf(:), deps(:)
    real(kind=8) :: dt, teta1, teta2, sech, matdech(6, 6)
    aster_logical :: iso, matrEndo
!   variable logique pour activer le fluage, l endo, le traitement local
    aster_logical, intent(in) :: end3d, fl3d, local
!    nombre maximal de sous-itération de fluage
    integer(kind=8), intent(in) :: iteflumax
! ----------------------------------------------------------------------
    aster_logical :: is_ba, inputL(5)
    real(kind=8) :: beta, beta00, bg0, biotw, biotw00, ccmin0, dim3
    real(kind=8) :: cthp, cthp1, cthv, cthv1, cthvm, cwtauk0, cwtauk1
    real(kind=8) :: cwtaukm, delta, delta00, denom, denomin
    real(kind=8) :: dfin0, dfin1, dflu0, dflu1, dfmx, dhydr, dt1, dth80, dteta, dth0
    real(kind=8) :: dth00, dth1, dther, dtmaxi, dtmaxik, dtmaxim, dvrgi, dvw, ekdc
    real(kind=8) :: epc0, epc00, epleqc, epleqc0, epleqc00, vrgi
    real(kind=8) :: epleqc01, epser, epsklim1, epsklim2, epsm00, inputVR6(6, 16)
    real(kind=8) :: ept, ept00, ept1, epsklim, hyd00, inputR(25)
    real(kind=8) :: pg0, rc00, rc1, reduc1, ref00, ref1, rt00, rt1, outputR(12)
    real(kind=8) :: tauk00, tauk1, taum00, teta, treps, trepspg, gfr, epeqpc, errgf
    real(kind=8) :: vrgi00, vrgi1, vrgi2, vw, vw1, vw2, we0s
    real(kind=8) :: xflu, young00, xnsat, xnsat00, outputVR6(6, 11)
    integer(kind=8) :: i, j, npas1, nt, ifour, iplalim, inputI(4), outputI(3), k
!
    integer(kind=8) :: ngf
    parameter(ngf=65)
!   tableau pour la resolution des systemes lineaires
    real(kind=8) :: X(ngf)
    real(kind=8) :: B(ngf)
    real(kind=8) :: A(ngf, ngf+1)
    integer(kind=8) :: ipzero(ngf)
!   donnes pour test fluage3d
    real(kind=8) :: epse06(6), epsk06(6), epsm06(6), sig06(6), phi0, we0, taum1
    real(kind=8) :: epse16(6), epsk16(6), epsm16(6), sig16(6)
    real(kind=8) :: we1, psik, epsm11
    real(kind=8) :: souplesse66(6, 6), raideur66(6, 6)
    real(kind=8) :: deps6(6), theta, theta1, heta, inputMat33(3, 3, 3), outputMat33(3, 3, 4)
!   theta : theta methode pour Euler, heta: taux de variation maxi
!   des deformations de fluage  lors d une sous iteration
    parameter(theta=0.5d0, heta=0.1d0)
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
    integer(kind=8) :: ipla, err1
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
    aster_logical :: ppas0
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
    real(kind=8) :: dt3(3), dr3(3), dgt3(3), dgc3(3), dcc, wl3(3)
!   traitement endommagement isotrope prepic
    real(kind=8) :: xmt, dtr
    aster_logical :: dtiso
!   avancement de la reaction de gonflement interne (rag par exple)
    real(kind=8) :: vrgi0, taar, nrjg, srsrag, aar0, aar1, trag, vrag00
    real(kind=8) :: vdef00, def1, tdef, nrjd, def0, srsdef, CNa
    real(kind=8) :: nrjp, ttrd, tfid, ttdd, tdid, exmd, exnd, cnab, cnak, ssad
    real(kind=8) :: At, St, M1, E1, M2, E2, AtF, StF, M1F, E1F, M2F, E2F, ttkf, nrjf
!   valuer fluage debut de pas
    real(kind=8) :: epsk006(6), epsm006(6)
!   Def thermiques transitoires (si Dtheta sous charge)
    real(kind=8) :: deps6r3(6), ett600(6), ett60(6)
!   14mai2015: ouvertures de fissures du pas precedent
    real(kind=8) :: wplt6(6), wplt06(6), wplt006(6), wpltx6(6), wpltx06(6), wpltx006(6)
    real(kind=8) :: wpl3(3), vwpl33(3, 3), vwpl33t(3, 3), wplx3(3), vwplx33(3, 3), vwplx33t(3, 3)
    real(kind=8) :: epstf6(6), sigmf6(6)
!   Objets pour matrice de decharge
    real(kind=8) ::  sigv(6), sigvd(6), matendo(6, 6), matdechac(6, 6)
    real(kind=8) ::  coef, coef1, coef2, coef3
    real(kind=8) ::  mathook(6, 6), vremp, rhov
!-----------------------------------------------------------------------
    if (nmatbe2 .eq. 0) then
        is_ba = ASTER_FALSE
        iplalim = 4
    else
!       RGI_BETON_BA
        is_ba = ASTER_TRUE
        iplalim = 1000
    end if
!
    call iniInt0(err1, npas1, nt, ifour)
    call iniReal0(vrgi2, epleqc, epleqc0, epleqc00, deltam, avean, &
                  epleqc01, we0s, epeqpc, we0, lambda, mu, bg, mg, pg)
    call iniReal0(pw, bw, sigp, cwtauk, srw, bw0, pw0, dfl00, dfl0, &
                  cwp, dcc, xmt, dtr, vrgi0, aar0, aar1, def1, def0)
    call iniReal0(At, St, M1, E1, M2, E2, AtF, StF, M1F, E1F, M2F, E2F, &
                  coef1, coef2, coef3, vremp)
    call iniVect0(3, cc03, dt3, dr3, dgt3, dgc3, wl3)
    call iniVect0(6, deps6r2, deps6r3, deps6r, sigf6d, dsw6, &
                  epsk06, sig16, epsm16, epsm06, epse06, sigv, sigvd)
    call iniVect0(6, epspg6, epspc6, epspt600, epspt60, epspg600, epspg60, &
                  epspc600, epspc60, dpg_depsa6, dpg_depspg6)
    call iniVect0(6, sigke06, sigke16, epsk006, epsm006, ett600, ett60, &
                  wplt6, wplt06, wplt006, wpltx6, wpltx06, wpltx006)
    call iniMat0(3, ref33, rt33, rtg33, vcc33, vcc33t)
    call iniMat0(6, souplesse66, raideur66, matendo, matdechac, mathook)
    call iniVect0(ngf, X, B)
!
    A(:, :) = 0.d0
    ipzero(:) = 0
!   chargement des tailles si l endo est active
    ierr1 = 0
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
    call getValVect(xmat, hydr, hyds, rt00, ref00, rc00, delta00, beta00, &
                    ept00, hplg, phivg, kgel, gft00, epsm00, psik, &
                    xflu, tauk00, taum00, ind1=HYDR)
!     stockage dans une variable interne pour avoir la vitesse
    if (abs(var0(PPAS)-1.d0) .ge. r8prem()) then
!       au 1er passage on suppose hyd0=hydr
        hyd0 = hydr
!       on initialise var03d en cas de sous incrementation par fluage
        var0(HYDF) = hyd0
    else
        hyd0 = var0(HYDF)
    end if
    dhydr = hydr-hyd0
    varf(HYDF) = hydr
!     stockage hydratation initiale pour calcul final de pression
    hyd00 = hyd0
!
!     Module d'Young et coefficient de Poisson
    if (.not. is_ba) then
        call getValVect(xmat, young00, nu00, ind1=YOUN)
    else
        call getValVect(xmat, young00, nu00, ind1=YOUM)
    end if
!     deformation de reference  pour le fluage prise aqu 1/3 de rc
    epser = (rc00/young00)/3.d0
!     verif validite de la dilatance
    if (beta00 .gt. dsqrt(3.d0)) then
        call utmess('E', 'COMPOR3_23', sr=dsqrt(3.d0))
        ierr1 = 1
        go to 999
    end if
!
!   nrjm = activation du potentiel de fluage (dmax1 pour le cas endo seul, utile dans thermat3d)
!   dth80 = endommagement thermique a 80°C
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
    call getValVect(xmat, nrjm, dth80, biotw00, xnsat00, poro2, vrag00, &
                    nrjf, sfld, mvgn, epc0, ekdc, eprg00, gfr, alat, krgi00, &
                    tetar, tetas, ind1=NRJM)
    nrjm = dmax1(nrjm, 1.d0)
!   volume d eau pour le non sature
    vw2 = sech
!   initialisation des variables internes associee a la saturation si premier pas
    if (abs(var0(PPAS)-1.d0) .ge. r8prem()) then
        var0(WSHR) = vw2
        var0(VSHR) = poro2
    end if
!   eps pic comp ne peut pas etre inferieur à rc/E
    epc00 = dmax1(epc0, 3.d0*epser)
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
    call getValVect(xmat, dfmx, taar, nrjg, srsrag, trag, dim3, &
                    tdef, nrjp, srsdef, vdef00, cna, ssad, cnak, &
                    cnab, exnd, exmd, ttdd, ind1=DFMX)
!   tdid = temps caracteristique pour la dissolution de l ettringite primaire
!   tfid = temps caracteristique pour la fixation des aluminiums en temperature
!   nrjd = energie d activation des processus de dissolution des phases primaires
!   ttrd = temperature de reference pour la precipitation de la def
!   ttkf = température seuil de fixation des alus en HG (RSI)
!   hpev = ratio pour levolution du module decrouissage
    call getValVect(xmat, tdid, tfid, nrjd, ttrd, &
                    ttkf, hpev, ind1=TDID)
!   on peut traiter l'endommagement pre pic iso de traction si ept > Rt/E
    ept00 = dmax1(rt00/young00, ept00)
!
!   indicateur de premier passage pour hydracomp3d
    if (abs(var0(PPAS)-1.d0) .ge. r8prem()) then
        ppas0 = .true.
    else
        ppas0 = .false.
    end if
!
!   chargement de l increment de deformation imposee
    deps6(:) = 0.d0
    epstf6(:) = 0.d0
    deps6(1:nstrs) = deps(1:nstrs)
    epstf6(1:nstrs) = epstf(1:nstrs)
!   passage en epsilon
    deps6(4:6) = 0.5d0*deps6(4:6)
    epstf6(4:6) = 0.5d0*epstf6(4:6)
!     remarque si rt rtg ref et rc dependent de  l ecrouissage
!     il faut les actualiser en fonction de l ecrouissage avant de
!     de passer dans hydramat, il faut egalement que dra_dl soit
!     parfaitement compatible avec le recalcul en fonction des deformations
!     plastiques (ou bien rajouter des vari pour stocker les resistances)
    rc1 = rc00
!     l hydratation n a pas d influence sur la deofrmation plastique
!     caracteristique de cisaillement
    rt1 = rt00
    ept1 = ept00
    ref1 = ref00
!     parametres materiau dependant eventuellement de l hydratation et de l endo de fluage
!     pour le 1er passage la reduction des resistances par fluage est
!     negligee car non utilise pour le tir visco elastique
    call hydramat3d(hyd0, hydr, hyds, young00, young, &
                    nu00, nu, rt1, rt, ref1, &
                    ref, rc1, rc, delta00, delta, &
                    beta00, beta, gft00, gft, ept1, &
                    ept, pglim, epsm00, epsm11, xnsat00, &
                    xnsat, biotw00, biotw, krgi00, krgi, &
                    iso, lambda, mu, rt33, rtg33, &
                    ref33, raideur66, souplesse66, xmt, dtiso, &
                    err1)
! - influence de la temperature sur les parametres  materiau et
!   actualisation de l endo thermique initial
    dth00 = var0(DTHE)
    call thermat3d(teta1, nrjm, tetas, tetar, dth80, &
                   dth00, dth0, CTHp0, CTHv0)
! - endommagement thermique en fin de pas
    call thermat3d(teta2, nrjm, tetas, tetar, dth80, &
                   dth0, dth1, CTHp1, CTHv1)
    varf(DTHE) = dth1
! - chargement des variables internes du fluage (etat du squelette solide)
    do i = 1, 6
        epsk006(i) = var0(EPK(i))
        epsm006(i) = var0(EPM(i))
        sig06(i) = var0(SIG(i))
        sigke06(i) = var0(SKE(i))
        dsw06(i) = var0(DSW(i))
    end do
!   phi00 = dissipation visqueuse
!   dfl00 = endommagement effectif par fluage
!   dth00 = endommagement thermique
!   epleqc00 = deformation plastique equivallente de cisaillement
!   bw0, pw0 : pression capillaire
!   bg0, pg0 : pression RGI
    call getValVect(var0, phi00, dfl00, dth00, epleqc00, bw0, pw0, &
                    bg0, pg0, vectInd=[PHIM, DFLU, DTHE, EPLC, BIOW, PSHR, BIOG, PRGI])
!   tenseurs de deformations plastique et surpression capillaire
    do i = 1, 6
        epspt600(i) = var0(EPTi(i))
        epspg600(i) = var0(EPGi(i))
        epspc600(i) = var0(EPCi(i))
        ett600(i) = var0(AFT1-1+i)
        wplt006(i) = var0(WID(i))
        wpltx006(i) = var0(EMT(i))
    end do
!
! - influence du degre d hydratation sur les variables internes
    call hydravar3d(hyd0, hydr, hyds, phi00, phi0, &
                    dth00, dth0, epleqc00, epleqc0, epspt600, &
                    epspt60, epspg600, epspg60, epspc600, epspc60, &
                    epsk006, epsk06, epsm006, epsm06, dfl00, &
                    dfl0, ett600, ett60, wplt006, wplt06, &
                    wpltx006, wpltx06)
! - calcul de la pression capillaire due a leau en fin de pas
    if (fl3d) then
        call bwpw3d(mfr, biotw, poro2, vw2, xnsat, mvgn, pw, bw, srw)
!       modif eventuelle des viscosites en fonction de srw
        CWtauk1 = 1.d0/srw
        call setValVect(varf, pw, bw, CWtauk1, vectInd=[PSHR, BIOW, CSHR])
    end if
! - reevaluation de la deformation compatible avec l etat
!   de contrainte debut de pas(pour evaluer la sous incrementation)
    call hydracomp3d(we0, we0s, epse06, souplesse66, sig06, &
                     deps6, deps6r, sigke06, epsk06, psik, &
                     fl3d)
! - reevaluation des coeffs de consolidation en debut de pas
    if (fl3d) then
!       effet de l eau
        CWp0 = var0(WSHR)/var0(VSHR)
        CWtauk0 = var0(VSHR)/var0(WSHR)
!       effet l endo de fluage sur le potentiel de fluage
        call dflufin3d(sig06, bw0, pw0, bg0, pg0, &
                       dsw06, delta, rc, xflu, dfin0, &
                       CMp0, dfmx)
!       effet de la temperature sur le potentiel
        call conso3d(epsm11, epser, ccmin0, ccmax0, epsm06, &
                     epse06, cc03, vcc33, vcc33t, CWp0, &
                     CMp0, CthP0, Cthv0)
! ----  determination de la subdivision eventuelle du pas de temps
        dtmaxi = dt
        do i = 1, 6
!           extremes de epsk
            epsklim1 = epse06(i)/psik
            epsklim2 = (epse06(i)+deps6r(i))/psik
!            valeur maximale possible de l increment de deformation de kelvin
            if (dabs(epsklim1) .gt. dabs(epsklim2)) then
                epsklim = epsklim1
            else
                epsklim = epsklim2
            end if
!           cas ou epsklim est faible
            if (dabs(epsklim) .gt. 1.d-6) then
                denom = dabs(1.d0-epsk06(i)/epsklim)
            else
                denom = dabs(1.d0-epsk06(i)/1.d-6)
            end if
!           comparaison avec la valeur permettant de franchir dt en 1/heta pas
            cwtaukm = (CWtauk0+cwtauk1)/2.d0
            cthvm = (cthv0+cthv1)/2.d0
            if (abs(dt) .ge. r8prem()) then
                denomin = heta*tauk00*CWtaukm/CTHVm/dt
            else
                denomin = 1.d0
            end if
!            comparaison avec les autres composantes de deformation
            if (denom .le. denomin) then
                denom = denomin
            end if
            dtmaxik = tauk00*CWtaukm/CTHVm/denom
!            cas de la deformation de maxwell
            dtmaxim = heta*(taum00*ccmin0)
!            choix entre condition maxwell et kelvin
            dtmaxi = dmin1(dtmaxi, dtmaxim, dtmaxik)
        end do
!       *** subdivision du pas si necessaire  **************************
        if (dtmaxi .lt. dt) then
            npas1 = int(dt/dtmaxi)+1
            if (iteflumax .gt. 0.d0) then
                npas1 = min(npas1, int(iteflumax))
            end if
        else
            npas1 = 1
        end if
    else
!        pas de fluage
        npas1 = 1
    end if
!      coeff de reduction de l increment
    reduc1 = 1.d0/dble(npas1)
!      reduction des pas d hydratation
    dhydr = reduc1*(hydr-hyd0)
!      pas de temps reduit
    dt1 = dt*reduc1
!      increment de temperature reduit
    dteta = (teta2-teta1)*reduc1
!      initialisation de la temperature debut de pas
    teta = teta1
!      increment de volume d eau reduit
    vw1 = var0(WSHR)
    dvw = (vw2-vw1)*reduc1
!      increment de poro capillaire reduite
    poro1 = var0(VSHR)
    dporo = (poro2-poro1)*reduc1
!      increment des potentiels de rgi
    vrgi1 = var0(PHIG)
    dvrgi = (vrgi2-vrgi1)*reduc1
!      increments de deformations reduits, on reduit deps6 et non deps6r
!      car hydracomp est reapplique plus bas
    if (npas1 .ne. 1) then
        deps6r(:) = reduc1*deps6(:)
    else
        deps6r(:) = deps6(:)
    end if
!***********************************************************************
!       debut du chargement sous-discretise sur le pas de temps
!***********************************************************************
    do nt = 1, npas1
!        chargement des variables internes
        do i = 1, 6
            epsk006(i) = var0(EPK(i))
            epsm006(i) = var0(EPM(i))
            sig06(i) = var0(SIG(i))
            sigke06(i) = var0(SKE(i))
            dsw06(i) = var0(DSW(i))
        end do
!
!       recuperation de l endommagement par fluage dfl00 et thermique dth00
        call getValVect(var0, bw0, pw0, bg0, pg0, phi00, dth00, hyd0, vw, &
                        poro, vrgi00, epleqc00, dfl00, &
                        vectInd=[BIOW, PSHR, BIOG, PRGI, PHIM, DTHE, HYDF, WSHR, &
                                 VSHR, PHIG, EPLC, DFLU])
!        actualisation du degre d hydratation
        hydr = hyd0+dhydr
!        actualisation de la temperature
        teta = teta+dteta
!        actualisation volume deau capillaire
        vw = vw+dvw
!        actualisation de la porosite capillaire
        poro = poro+dporo
!        actualisation du volume pour les rgi
        vrgi00 = vrgi00+dvrgi
        call setValVect(varf, hydr, vw, poro, vrgi00, vectInd=[HYDF, WSHR, VSHR, PHIG])
!
!        recuperation de la deformation maximale de traction
        do i = 1, 6
            epspt600(i) = var0(EPTi(i))
            epspg600(i) = var0(EPGi(i))
            epspc600(i) = var0(EPCi(i))
            ett600(i) = var0(AFT1-1+i)
            wplt006(i) = var0(WID(i))
            wpltx006(i) = var0(EMT(i))
        end do
!
!       prise en compte  hydratation intermediaire si sous-increment
        call hydravar3d(hyd0, hydr, hyds, phi00, phi0, &
                        dth00, dth0, epleqc00, epleqc0, epspt600, &
                        epspt60, epspg600, epspg60, epspc600, epspc60, &
                        epsk006, epsk06, epsm006, epsm06, dfl00, &
                        dfl0, ett600, ett60, wplt006, wplt06, &
                        wpltx006, wpltx06)
!
!        stockage de la valeur de la deformation plastique cumulee
!        modifiee par l hydratation pour mise a jour incrementale
        epleqc01 = epleqc0
!        effet de l endommagement de fluage et de l'hydratation sur les
!        parametres materiau
        call hydramat3d(hyd0, hydr, hyds, young00, young, &
                        nu00, nu, rt1, rt, ref1, &
                        ref, rc1, rc, delta00, delta, &
                        beta00, beta, gft00, gft, ept1, &
                        ept, pglim, epsm00, epsm11, xnsat00, &
                        xnsat, biotw00, biotw, krgi00, krgi, &
                        iso, lambda, mu, rt33, rtg33, &
                        ref33, raideur66, souplesse66, xmt, dtiso, &
                        err1)
!
!        influence de la temperature sur les parametres du materiau
!        et calcul de l endommagement thermique
        call thermat3d(teta, nrjm, tetas, tetar, dth80, &
                       dth1, DTHER, CTHP, CTHV)
        varf(DTHE) = dth1
!        calcul de la pression capillaire due a leau
        if (fl3d) then
            call bwpw3d(mfr, biotw, poro, vw, xnsat, &
                        mvgn, pw, bw, srw)
            varf(PSHR) = pw
            varf(BIOW) = bw
!        reevaluation des coeffs de consolidation apres increment hydra
!        effet de l eau sur le potentiel de fluage
            Cwp = Srw
            CWtauk = 1.d0/srw
            varf(CSHR) = CWtauk
!
!        Modification eventuelle de la viscosite en fonction de Srw
            tauk1 = tauk00*CWtauk/CTHV
!
!        la modif de la viscosite de Maxwell est comprise dans le coeff
!        de consolidation (cf. conso3d)
            taum1 = taum00
        end if
!
!        compatibilite des anciennes contraintes avec nouveau materiau
!        deduction de la deformation de solidification de l increment
        call hydracomp3d(we0, we0s, epse06, souplesse66, sig06, &
                         deps6r, deps6r2, sigke06, epsk06, psik, &
                         fl3d)
!       coeff theta methode pour tir visco elastique
        theta1 = theta
!
!***********************************************************************
!        tir visco elastique
!***********************************************************************
!
!       E.Cheignon : on devrait supprimer deps6r3 qui n'est pas différent
!       de deps6r2 (voir tirViscoElas)
        deps6r3(:) = deps6r2(:)
!
        call setValVect(inputR, delta, rc, epsm11, epser, CWp, CthP, Cthv, &
                        dt1, theta1, tauk1, taum1, ind1=1)
!
        call setR6Mat6(inputVR6, dsw06, epsm06, sigke06, deps6r2, &
                       epsk06, epse06, sig06)
!
        call tirViscoElas(fl3d, var0, xmat, inputR, inputVR6, ngf, &
                          deltam, avean, A, B, X, ipzero, &
                          epsk16, epsm16, epse16, &
                          sig16, sigke16, raideur66, we1)
!
!       E.Cheignon : sorti de tirViscoElas pour ne pas passer varf
        do i = 1, 6
!            varf(96+i)=ett61(i)
            varf(AFT1-1+i) = 0.d0
        end do
!
!       actualisation de la variation de volume total
        treps = var0(TEPS)
        do i = 1, 3
            treps = treps+deps6r(i)
        end do
        varf(TEPS) = treps
!
!       chargement des deformations plastiques du pas precedent
!       actualisee par l hydratation
        do j = 1, 6
!         traction
            epspt6(j) = epspt60(j)
!         rgi
            epspg6(j) = epspg60(j)
!         compression
            epspc6(j) = epspc60(j)
!         chargement des ouvertures de fissures du pas precedent
            wplt6(j) = wplt06(j)
!         chargement des ouvertures maxi de fissure
            wpltx6(j) = wpltx06(j)
        end do
!       chargement de la def equivalente actualisee par l hydratation
        epleqc = epleqc0
!
!       recuperation de la variation volumique due a la rgi non
!       actualisee par l hydratation
        trepspg = var0(TEPG)
!
!***********************************************************************
!       verification des criteres de plasticite et ecoulements
!***********************************************************************
!
        call setValVect(inputR, biotw, poro, vw, xnsat, pglim, teta, dt1, &
                        young00, nu00, rc, epleqc0, delta, beta, krgi, &
                        theta, epsm11, ind1=1)
        call setValVect(inputR, CWp, CthP, Cthv, tauk1, taum1, deltam, avean, &
                        phi0, epleqc01, srw, ind1=17)
!
        call setR6Mat6(inputVR6, sig0, epspt60, wplt06, wpltx06, &
                       sig06, epsm06, dpg_depsa6, dpg_depspg6, epsk16, &
                       epsm16, epse16, epspt6, sigke16, sig16, epspg6, epspc6)
        call setMat33Tab(inputMat33, ref33, rtg33, rt33)
        call setIntVect(inputI, iplalim, mfr, nstrs, ndim)
        call setLogVect(inputL, fl3d, ppas0, iso, local, end3d)

        call plasti3d(xmat, inputR, inputVR6, inputMat33, inputI, &
                      inputL, var0, raideur66, souplesse66, &
                      A, B, X, ngf, varf, ipzero, &
                      outputR, outputVR6, outputMat33, &
                      outputI)
!
        call getValVect(var0, aar0, def0, E1, M1, E2, M2, At, St, &
                        vectInd=[AAAR, ADEF, AFT1, AFM1, AFT2, AFM2, ATIL, STIL])
!
        call getValVect(varf, aar1, def1, E1f, M1f, E2f, M2f, &
                        Atf, Stf, vrgi, pgmax, pg, bg, trepspg, &
                        vectInd=[AAAR, ADEF, AFT1, AFM1, AFT2, AFM2, ATIL, STIL &
                                 , PHIG, PGMAX, PRGI, BIOG, TEPG])
        call getValVect(outputR, srw, epeqpc, dfin1, ccmax1, pw, bw, &
                        wpl3(1), wpl3(2), wpl3(3), wplx3(1), wplx3(2), &
                        wplx3(3), ind1=1)
        call getR6Mat6(outputVR6, dpg_depsa6, dpg_depspg6, epsk16, &
                       epsm16, epse16, epspt6, sigke16, sig16, epspg6, &
                       epspc6, dsw6)
        call getMat33Tab(outputMat33, vwpl33, vwpl33t, vwplx33, vwplx33t)
        call getIntVect(outputI, ifour, ierr1, ipla)
!
!       fin de la boucle de consistance visco-elasto-plastique
!       ****************************************************************
!       transfert des variables internes pour la sous iteration temporelle suivante
!       si le nbre de sous iteration locale le justifie
        if (npas1 .gt. 1) var0(1:nvari) = varf(1:nvari)
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
        vrgi0 = vrgi
        pgmax = var0(PGMAX)
!          dt mis a zero pour forcer la reprise de vrgi0
        call bgpg3d(ppas0, bg, pg, mg, vrgi, &
                    treps, trepspg, epspt6, epspc6, phivg, &
                    pglim, dpg_depsa6, dpg_depspg6, taar, nrjg, &
                    trag, aar0, srw, srsrag, teta, &
                    dt1, vrag00, aar1, tdef, nrjd, &
                    def0, srsdef, vdef00, def1, cna, &
                    nrjp, ttrd, tfid, ttdd, tdid, &
                    exmd, exnd, cnab, cnak, ssad, &
                    At, St, M1, E1, M2, &
                    E2, AtF, StF, M1F, E1F, &
                    M2F, E2F, vrgi0, ttkf, nrjf, &
                    alat, young00, nu00, kgel, pgmax)
    else
        call getValVect(varf, pg, bg, pgmax, vectInd=[PRGI, BIOG, PGMAX])
    end if
!       stockage de la pression RGI
    call setValVect(varf, pg, bg, pgmax, vectInd=[PRGI, BIOG, PGMAX])
!
! - endommagement de fluage
    dflu1 = 0.d0
    if (fl3d) then
        dflu0 = var0(DFLU)
        call dflueff3d(ccmax1, dflu0, dflu1, dfin1)
    end if
    varf(DFLU) = dflu1
!***********************************************************************
!       contraintes dans solide et rgi en fin de pas
!       resultante des pressions intraporeuses RGI et Capillaire (depression)
    sigp = -bg*pg-bw*pw
!       effet sur la contrainte apparente en non sature
    do i = 1, 6
        if (i .le. 3) then
!               prise en compte de la pression rgi
            sigf6(i) = sig16(i)+sigp+dsw6(i)
        else
            sigf6(i) = sig16(i)+dsw6(i)
        end if
    end do
!
!
!       prise en compte de l'endommagement mécanique
    if (end3d) then
!       chargement endo traction pre-pic
        dtr = var0(DTPP)
!       chargement endo localisee pour condition de croissance
        do i = 1, 3
            dt3(i) = var0(DTL(i))
        end do
!
!***********************************************************************
!       calcul des endommagements et ouvertures de fissures
        call endotot(dth1, dflu1, end3d, &
                     wpl3, vwpl33, vwpl33t, wplx3, vwplx33, &
                     vwplx33t, gft, gfr, iso, sigf6, &
                     sigf6d, rt33, ref33, souplesse66, epspg6, &
                     eprg00, a, b, x, ipzero, &
                     ngf, ekdc, epspc6, dt3, dr3, &
                     dgt3, dgc3, dcc, wl3, xmt, &
                     dtiso, rt, dtr, dim3, ndim, &
                     ifour, epeqpc, ept, errgf)
!            stockage des endommagements de fissuration et ouverture
        varf(DTPP) = dtr
        do i = 1, 3
!               endo de traction
            varf(DTL(i)) = dt3(i)
!               endo refermeture
            varf(DRL(i)) = dr3(i)
!               endo RGI en traction
            varf(DTG(i)) = dgt3(i)
!               endo RGI en compression
            varf(DCG(i)) = dgc3(i)
!               ouverture non visco elastique
            varf(WL(i)) = wl3(i)
        end do
!            endo de traction global
        varf(DT0) = 1.d0-((1.d0-(varf(DTL(1))))*(1.d0-(varf(DTL(2))))*(1.d0-(varf(DTL(3)))))
!            endo de traction de RGI global
        varf(DTG0) = 1.d0-((1.d0-(varf(DTG(1))))*(1.d0-(varf(DTG(2))))*(1.d0-(varf(DTG(3)))))
!            endo de compression de RGI global
        varf(DCG0) = 1.d0-((1.d0-(varf(DCG(1))))*(1.d0-(varf(DCG(2))))*(1.d0-(varf(DCG(3)))))
!            endommagement de compression
        varf(DC) = dcc
!            erreur de dissipation d'énergie en traction
        varf(ERGF) = errgf
!            traitement erreur endo
        if (err1 .eq. 1) then
            call utmess('E', 'COMPOR3_34')
            ierr1 = 1
            go to 999
        end if
!       contraintes totale dans la matrice apres endommagement
        sigmf6(1:6) = sigf6d(1:6)
    else
!       pas d endommagement
        sigf6d(1:6) = sigf6(1:6)
        sigmf6(1:6) = sigf6(1:6)
    end if

!
    if (is_ba) then
        call rgiRenfoStress(xmat, nbelas3d+nmatbe2+1, sigmf6, &
                            epstf6, epspt6, teta1, teta2, dt, ppas0, theta, &
                            fl3d, end3d, wpl3, vwpl33, vwpl33t, dt3, dr3, &
                            ipzero, ngf, rc00, var0, varf, sigf6d, &
                            matdechac, rhov, ierr1)
    end if
!***********************************************************************
!      Matrice de décharge
    if (matrEndo) then
!       Boucle de remplissage du pseudovecteurvirtuel
        do i = 1, 6
            sigv(i) = 0.d0
            do j = 1, 6
                if (i .eq. j) then
                    if (sigf6(j) .ge. (1.5d0*rt00/100.d0)) then
                        sigv(j) = sigf6(j)-(rt00/100.d0)
                    else
                        sigv(j) = 2.d0*(rt00/100.d0)
                    end if
                else
                    sigv(j) = sigf6(j)
                end if
            end do
!       Appel de la routine avec les contrainte virtuelles
            call endotot(dth1, dflu1, end3d, &
                         wpl3, vwpl33, vwpl33t, wplx3, vwplx33, &
                         vwplx33t, gft, gfr, iso, sigv, &
                         sigvd, rt33, ref33, souplesse66, epspg6, &
                         eprg00, a, b, x, ipzero, &
                         ngf, ekdc, epspc6, dt3, dr3, &
                         dgt3, dgc3, dcc, wl3, xmt, &
                         dtiso, rt, dtr, dim3, ndim, &
                         ifour, epeqpc, ept, errgf)
            if (is_ba) then
                do j = 1, 6
                    matendo(j, i) = ((sigvd(j)-sigmf6(j))/(sigv(i)-sigf6(i)))
                end do
            else
                do j = 1, 6
                    matendo(j, i) = ((sigvd(j)-sigf6d(j))/(sigv(i)-sigf6(i)))
                end do
            end if
        end do
!       Multiplication avec la matrice de rigidité
        coef = 1.d0/((1.d0+xmat(2))*(1.d0-2.d0*xmat(2)))
        coef1 = xmat(1)*(1.d0-xmat(2))*coef
        coef2 = xmat(1)*xmat(2)*coef
        coef3 = xmat(1)*(1.d0-2*xmat(2))*coef
        do i = 1, 6
            do j = 1, 6
                if (i .eq. j) then
                    if (i .le. 3) then
                        mathook(i, j) = coef1
                    else
                        mathook(i, j) = coef3
                    end if
                else if ((i .le. 3) .and. (j .le. 3)) then
                    mathook(i, j) = coef2
                else
                    mathook(i, j) = 0.d0
                end if
            end do
        end do
        do i = 1, 6
            do j = 1, 6
                vremp = 0.d0
                do k = 1, 6
                    vremp = vremp+matendo(i, k)*mathook(k, j)
                end do
                matdech(i, j) = vremp
            end do
        end do

        if (is_ba) then
!           Matrice homogeneise
            do i = 1, 6
                do j = 1, 6
                    matdech(i, j) = (1.d0-rhov)*matdech(i, j)+matdechac(i, j)
                end do
            end do
        end if
    end if
!
!   affectation dans le tableau de sortie des contraintes
    sigf(1:nstrs) = sigf6d(1:nstrs)
999 continue
!
end subroutine
