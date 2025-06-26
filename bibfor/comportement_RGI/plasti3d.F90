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
subroutine plasti3d(xmat, inputR, inputVR6, inputMat33, inputI, &
                    inputL, var0, raideur66, souplesse66, &
                    A, B, X, ngf, varf, ipzero, &
                    outputR, outputVR6, outputMat33, &
                    outputI)
! person_in_charge: etienne.grimal@edf.fr
!-----------------------------------------------------------------------
!       verification des criteres de plasticite et ecoulements
!-----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "rgi_module.h"
#include "asterfort/b3d_valp33.h"
#include "asterfort/bgpg3d.h"
#include "asterfort/bwpw3d.h"
#include "asterfort/chrep6.h"
#include "asterfort/couplagf3d.h"
#include "asterfort/couplagpf3d.h"
#include "asterfort/couplagfp3d.h"
#include "asterfort/couplagpp3d.h"
#include "asterfort/conso3d.h"
#include "asterfort/criter3d.h"
#include "asterfort/dflufin3d.h"
#include "asterfort/fludes3d.h"
#include "asterfort/getValVect.h"
#include "asterfort/gauss3d.h"
#include "asterfort/iniVect0.h"
#include "asterfort/matfluag3d.h"
#include "asterfort/majw3d.h"
#include "asterfort/transpos1.h"
#include "asterfort/utmess.h"
#include "asterfort/x6x33.h"
#include "asterfort/plastMult.h"
#include "asterfort/iniInt0.h"
#include "asterfort/setValVect.h"
#include "asterfort/getR6Mat6.h"
#include "asterfort/setR6Mat6.h"
#include "asterfort/getMat33Tab.h"
#include "asterfort/setMat33Tab.h"
#include "asterfort/getIntVect.h"
#include "asterfort/setIntVect.h"
#include "asterfort/getLogVect.h"
!   declaration des arguments
    integer(kind=8), intent(in) :: inputI(*), ngf
    aster_logical, intent(in) :: inputL(*)
    real(kind=8), intent(in) :: xmat(*), var0(*), inputR(*), inputVR6(6, *)
    real(kind=8), intent(in) :: inputMat33(3, 3, *), raideur66(6, 6), souplesse66(6, 6)
    real(kind=8), intent(inout) :: A(ngf, ngf+1), B(ngf), X(ngf), varf(*)
    real(kind=8), intent(out) :: outputR(*), outputVR6(6, *), outputMat33(3, 3, *)
    integer(kind=8), intent(out) :: outputI(*)
    integer(kind=8), intent(inout) :: ipzero(ngf)
! ----------------------------------------------------------------------
!   variables internes
!   irr : longeur de lelement dans la base des direction principales
!   controle des actions de criter3d et stockage des directions d ecoule
    integer(kind=8) :: irr, ipla2, err1, nf0, vali(2), ifour, ierr1, ipla
    integer(kind=8) :: i, j, iplalim, mfr, nstrs, ndim
!   logic1 : logique pour le controle de coherence de dissipation par fluage
    aster_logical :: fl3d, ppas, iso, local, end3d
    aster_logical :: indic2, referm3(3), limit1, logic1, goto20
    real(kind=8) :: aar0, def0, vrgi0, pgmax, hpev, rc00, mvgn
    real(kind=8) :: aar1, def1, vrgi, pg, bg, mg, treps, trepspg
    real(kind=8) :: phivg, taar, nrjg, trag, srsrag, vrag00, tdef, nrjd
    real(kind=8) :: srsdef, vdef00, cna, nrjp, ttrd, tfid, ttdd, tdid
    real(kind=8) :: exmd, exnd, cnab, cnak, ssad, ttkf, nrjf, alat, kgel
    real(kind=8) :: sfld, sig133(3, 3), sig13(3), vsig133(3, 3), vsig133t(3, 3)
    real(kind=8) :: bw0, pw0, sig16p(6), raideur66p(6, 6)
    real(kind=8) :: souplesse66p(6, 6), rt33p(3, 3), rtg33p(3, 3)
    real(kind=8) :: ref33p(3, 3), epspt6p(6), epspg6p(6), epspc6p(6)
    real(kind=8) :: depspt6p(6), depspg6p(6), depspc6p(6), epc00, epc0
!   tauc : taux de cisaillement
    real(kind=8) :: epser, epleqc, depleqc_dl, ekdc, hplg, tauc, deps6r3(6)
    real(kind=8) :: theta1, xflu, dfmx, CMp1, ccmin1, dim3, phi1, we1
    real(kind=8) :: cc13(3), vcc33(3, 3), vcc33t(3, 3), psik, srw, epeqpc
    real(kind=8) :: kveve66(6, 6), kvem66(6, 6), kmve66(6, 6), kmm66(6, 6)
    real(kind=8) :: bve6(6), bm6(6), depleqc, epsmin, wplt6(6), wpltx6(6)
    real(kind=8) :: depsk6p(6), depsm6p(6), depspt6(6), depse6(6), depsk6(6)
!   t33 : matrice pour la gestion de la taille dans castem
    real(kind=8) :: depsm6(6), depspg6(6), depspc6(6), depleqc3, t33(3, 3)
    real(kind=8) :: dphi, epse06p(6), epsk06p(6), epsm06p(6), depse6p(6)
    real(kind=8) :: biotw, poro, vw, xnsat, pglim, teta, dt1, young00, nu00
    real(kind=8) :: rc, delta, beta, krgi, theta, epsm11, CWp
    real(kind=8) :: CthP, Cthv, tauk1, taum1, deltam, avean, phi0, epleqc01
    real(kind=8) :: E1f, M1f, E2f, M2f, Atf, Stf, E1, M1, E2, M2, At, St
    real(kind=8) :: dfin1, ccmax1, pw, bw, sig0(6), epspt60(6), wplt06(6)
    real(kind=8) :: wpltx06(6), sig06(6), epsm06(6), dpg_depsa6(6), dpg_depspg6(6)
    real(kind=8) :: epsk16(6), epsm16(6), epse16(6), epspt6(6), sigke16(6), sig16(6)
    real(kind=8) :: epspg6(6), epspc6(6), dsw6(6), vwplx33(3, 3), vwplx33t(3, 3)
    real(kind=8) :: ref33(3, 3), rtg33(3, 3), rt33(3, 3), vwpl33(3, 3), vwpl33t(3, 3)
    real(kind=8) :: wpl3(3), wplx3(3)
!   nombres de criteres totaux et actifs
    integer(kind=8) :: nc
    parameter(nc=10)
!   ig(nc) : correspondance entre n° de critere actif et n° global
!   suivant ordre defini dans criter3d.f
    integer(kind=8) :: na, ig(nc), supr(nc)
    real(kind=8) :: fa(nc), dpfa_ds(nc, 6), dgfa_ds(nc, 6), dpfa_dpg(nc)
!   derivee de la resistance / multiplicateur plastique
    real(kind=8) :: dra_dl(nc), fg(nc), fglim(nc)
!   derivee fonction seuil / resistance
    real(kind=8) :: dpfa_dr(nc)
!   nombre maxi de sous iterations plastique : imax
    integer(kind=8) :: imax
    parameter(imax=1000)
    real(kind=8) :: hpla
    parameter(hpla=1.0d-3)
    integer(kind=8) :: vectind0(17), vectind1(17)
    data vectind0/VVRG, TAUG, NRJG, TRAG, SRSG, VRAG, TDEF, NRJD, SRSD, VDEF, &
        CNAD, NRJP, TTRD, TFID, TTDD, TDID, EXMD/
    data vectind1/EXND, CNAB, CNAK, SSAD, TTKF, NRJF, ALAT, KGEL, SFLD, EPC, &
        RC, EKDC, HRGI, HPEV, XFLU, DFMX, YKSY/
! ----------------------------------------------------------------------
    pw = 0.d0
    bw = 0.d0
! ----------------------------------------------------------------------

    call getValVect(inputR, biotw, poro, vw, xnsat, pglim, teta, dt1, &
                    young00, nu00, rc, epleqc, delta, beta, krgi, &
                    theta, epsm11, ind1=1)
    call getValVect(inputR, CWp, CthP, Cthv, tauk1, taum1, deltam, avean, &
                    phi0, epleqc01, srw, ind1=17)

    call getR6Mat6(inputVR6, sig0, epspt60, wplt06, wpltx06, &
                   sig06, epsm06, dpg_depsa6, dpg_depspg6, epsk16, &
                   epsm16, epse16, epspt6, sigke16, sig16, epspg6, epspc6)

    call getMat33Tab(inputMat33, ref33, rtg33, rt33)
    call getIntVect(inputI, iplalim, mfr, nstrs, ndim)
    call getLogVect(inputL, fl3d, ppas, iso, local, end3d)
    vwpl33(:, :) = 0.d0
    vwpl33t(:, :) = 0.d0
    vwplx33(:, :) = 0.d0
    vwplx33t(:, :) = 0.d0
    wpl3(:) = 0.d0
    wplx3(:) = 0.d0
!   initialisation du compteur sur la boucle de coherence plastique : ipla
!   initialisation du compteur sur la boucle de retour radial : irr
    call iniInt0(ipla, irr, err1, nf0, ifour, ierr1, na)
!   mg est à supprimer après confirmation de Pierre Morenon
    mg = 0
!   initialisation indicateur reduction nbre de criteres actifs
    indic2 = .false.
!
    dpfa_ds(:, :) = 0.d0
    dgfa_ds(:, :) = 0.d0
    ig(:) = 0
    supr(:) = 0
    call iniVect0(nc, fa, dpfa_dpg, fg, fglim, dra_dl, dpfa_dr)

    treps = varf(TEPS)
    trepspg = var0(TEPG)

    call getValVect(xmat, phivg, taar, nrjg, trag, srsrag, vrag00, tdef, &
                    nrjd, srsdef, vdef00, cna, nrjp, ttrd, tfid, ttdd, &
                    tdid, exmd, vectInd=vectind0)
    call getValVect(xmat, exnd, cnab, cnak, ssad, ttkf, nrjf, alat, kgel, &
                    sfld, epc0, rc00, ekdc, hplg, hpev, xflu, dfmx, psik, &
                    vectInd=vectind1)
    call getValVect(xmat, dim3, mvgn, vectInd=[DIM3, MVGN])

    epser = (rc00/young00)/3.d0
    epc00 = dmax1(epc0, 3.d0*epser)

!   chargement des ouvertures de fissures du pas precedent et maxi
    wplt6(1:6) = wplt06(1:6)
    wpltx6(1:6) = wpltx06(1:6)
!
10  continue
    if (ipla .lt. iplalim) then
        ipla = ipla+1
!       compteur reduction du systeme de couplage
        ipla2 = 0
        if (ipla .gt. iplalim) call utmess('A', 'COMPOR3_27', si=ipla)
!       reevaluation de la pression capillaire due a l eau si plasticite
        if (fl3d) then
            call bwpw3d(mfr, biotw, poro, vw, xnsat, &
                        mvgn, pw, bw, srw)
            varf(PSHR) = pw
            varf(BIOW) = bw
        end if
!
!   reevaluation de la pression de rgi si plasticite
!   vrgi volume effectivement formee pour ce pas
        call getValVect(var0, aar0, def0, E1, M1, E2, M2, At, St, vrgi0, &
                        pgmax, vectInd=[AAAR, ADEF, AFT1, AFM1, AFT2, AFM2, ATIL, STIL, &
                                        PHIG, PGMAX])

!   l avancement est actualisee dans bgpg
        call bgpg3d(ppas, bg, pg, mg, vrgi, &
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
!   stockage avancement et pression
        call setValVect(varf, aar1, def1, E1f, M1f, E2f, M2f, &
                        Atf, Stf, vrgi, pgmax, pg, bg, &
                        vectInd=[AAAR, ADEF, AFT1, AFM1, AFT2, AFM2, ATIL, STIL &
                                 , PHIG, PGMAX, PRGI, BIOG])
!
!   prise en compte de l'amplification des depressions capillaires
!   sous charge (l effet du chargement sur la pression capillaire
!   est traité explicitement via sig0)
        do i = 1, 6
            dsw6(i) = var0(DSW(i))
        end do
        bw0 = var0(66)
        pw0 = var0(56)
        call fludes3d(bw0, pw0, bw, pw, sfld, sig0, dsw6, nstrs)
        do i = 1, 6
            varf(DSW(i)) = dsw6(i)
        end do
!
!   possibilite ecoulement plastique couplé au fluage
!   les criteres de Rankine etant definis dans la base principale
!   des contraintes, on se place dans cette base pour la verification
!   des criteres et le retour radial le cas echeant
!
!   diagonalisation du vecteur des contraintes si 1er passage dans criter3d (irr=0)
        call x6x33(sig16, sig133)
        call b3d_valp33(sig133, sig13, vsig133)
!   construction matrice de passage inverse
        call transpos1(vsig133t, vsig133, 3)
!
!   si irr.ne.0 on est dans la boucle de sous iteration radiale
!   on garde la base du tir elastique pour le retour
!   notation pseudo vecteur pour ecoulement plastique
!   en base principale (en raison des criteres de rankine)
        call chrep6(sig16, vsig133, .false._1, sig16p)
!
!   passage de la loi de comportement dans la base principale
!   des contraintes
        if (iso) then
!      cas isotrope : rien a faire
            raideur66p(:, :) = raideur66(:, :)
            souplesse66p(:, :) = souplesse66(:, :)
        else
!      cas elasticite anisotrope : passge de la loi de comportement
!      en base principale des contraintes necessaire
            call utmess('F', 'COMPOR3_28')
        end if
!   passage des resistances et contraintes de refermeture
!   dans la base principale des contraintes
        if (iso) then
!     cas isotrope : rien a faire
            rt33p(:, :) = rt33(:, :)
            rtg33p(:, :) = rtg33(:, :)
            ref33p(:, :) = ref33(:, :)
        else
!     cas anisotrope : changement de base des resistance normales
            call utmess('F', 'COMPOR3_29')
        end if
!
!   passage des deformations plastiques de traction
!   dans la base principale du tir visco elastique
        call chrep6(epspt6, vsig133, .false._1, epspt6p)
!   passage des deformations plastiques de gel
        call chrep6(epspg6, vsig133, .false._1, epspg6p)
!   passage des deformations plastiques de cisaillement
        call chrep6(epspc6, vsig133, .false._1, epspc6p)
!   initialisation des increments de deformations
        call iniVect0(6, depspt6p, depspg6p, depspc6p)
!
!   evaluation des criteres actifs, des derivees et des directions
!   d ecoulement dans la base principale des contraintes
        call criter3d(sig16p, bg, pg, bw, pw, &
                      rt33p, rtg33p, ref33p, rc, epc00, &
                      epleqc, epspt6p, epspg6p, delta, beta, &
                      nc, ig, fg, na, fa, &
                      dpfa_ds, dgfa_ds, dpfa_dpg, dra_dl, souplesse66p, &
                      err1, depleqc_dl, irr, fglim, krgi, &
                      hpla, ekdc, hplg, dpfa_dr, tauc, &
                      epeqpc, hpev)
        indic2 = .false.
!       detection d erreur dans les criteres
        if (err1 .eq. 1) then
            ierr1 = 1
            go to 999
        end if
!   ****************************************************************
!   controle de la boucle de consitance visco elasto plastique
!   ****************************************************************
!   retour par ecoulement visco-plasticité couplée si au moins
!   un critere est actifs
        if (na .ne. 0) then
!     *** ecoulement plastique *************************************
!
!     reconstruction de la matrice de couplage fluage->fluage
!     dans la base principale actuelle
            if (fl3d) then
                call chrep6(epse16, vsig133, .false._1, epse06p)
                call chrep6(epsk16, vsig133, .false._1, epsk06p)
                call chrep6(epsm16, vsig133, .false._1, epsm06p)
            end if
!
!     pas d increments de deformation pendant la boucle de consistance
            deps6r3(:) = 0.d0
!
!     ****** couplages ds nouvelle base ****************************
            theta1 = theta
            if (fl3d) then
!        actualisation de la matrice de couplage visco plastique
!
!        effet du chargement sur le potentiel de fluage
!        la surpression hydrique est celle de debit de pas
!        si elle utilise la contrainte totale convergee
                call dflufin3d(sig16, bw, pw, bg, pg, &
                               dsw6, delta, rc, xflu, dfin1, &
                               CMp1, dfmx)
!
!               actualisation coeffs de consolidation
                call conso3d(epsm11, epser, ccmin1, ccmax1, epsm16, &
                             epse16, cc13, vcc33, vcc33t, CWp, &
                             CMp1, CthP, Cthv)
!
                call matfluag3d(epse06p, epsk06p, sig16p, psik, tauk1, &
                                taum1, deps6r3, dt1, theta1, kveve66, &
                                kvem66, kmve66, kmm66, bve6, bm6, &
                                deltam, avean, cc13, vcc33, vcc33t, &
                                vsig133, vsig133t)
!               mise a zero des seconds membres de fluage pour le retour radial
                call iniVect0(6, bve6, bm6)
!        reinitialisation et reassembalge de la matrice de couplage : a
!
!               cas du couplage fluage fluage
                call couplagf3d(A, B, ngf, kveve66, kmm66, &
                                kmve66, kvem66, bve6, bm6)
!               plasticite -> fluage
                call couplagpf3d(a, b, ngf, na, avean, &
                                 nc, dgfa_ds, deltam, kmve66)
!               fluage -> plasticite
                call couplagfp3d(a, ngf, na, nc, dpfa_ds, &
                                 dpfa_dpg, dpg_depsa6, raideur66p)
            end if
!      couplage plasticite -> plasticite
            if (fl3d) nf0 = 12
            call couplagpp3d(a, ngf, na, nc, nf0, &
                             ig, dgfa_ds, dpfa_ds, dpfa_dpg, dpg_depspg6, &
                             raideur66p, dpfa_dr, dra_dl, dpg_depsa6)
!     second membre de la plasticite : oppose des criteres actifs
            do i = 1, na
                b(nf0+i) = -fa(i)
            end do
!
!     *** resolution du retour radial ******************************
20          continue
            call gauss3d((nf0+na), A, X, B, ngf, err1, ipzero)
            if (err1 .eq. 1) then
                call utmess('E', 'COMPOR3_30')
                ierr1 = 1
                go to 999
            end if
!
!     **** verif positivite des multiplicateurs plastiques *********
            call plastMult(na, nf0, ngf, x, irr, nc, indic2, ipla2, &
                           imax, ig, a, b, dgfa_ds, goto20)
            if (goto20) goto 20
!
!     *** mise a jour des variables dependantes de la plasticite ***
!
!     initialisation de la deformation equivalente de compression
            depleqc = 0.d0
!     bouclage sur les criteres possibles
            referm3(:) = .false.
            do i = 1, na
                do j = 1, 6
                    if (ig(i) .le. 6) then
!                cas des deformations plastiques de traction
                        depspt6p(j) = depspt6p(j)+x(nf0+i)*dgfa_ds(i, j)
!                     indicateur de refermeture active
                        if (ig(i) .gt. 3) referm3(ig(i)-3) = .true.
                    else if (ig(i) .le. 9) then
!               cas des deformations plastiques de gel
                        depspg6p(j) = depspg6p(j)+x(nf0+i)*dgfa_ds(i, j)
                    else if (ig(i) .le. 10) then
!               cas des deformations plastiques de DP
                        depspc6p(j) = depspc6p(j)+x(nf0+i)*dgfa_ds(i, j)
                    else
                        call utmess('E', 'COMPOR3_32')
                        ierr1 = 1
                    end if
                end do
                if (ig(i) .eq. 10) depleqc = depleqc_dl*x(nf0+i)
            end do
!
!     *** limitation des refermetures si necessaire ****************
!     indicateur d'au moins une refermeture forcée
            limit1 = .false.
!     déformation minimale admissible
            epsmin = 0.d0
!     test des refermetures
            do j = 1, 3
                if (((epspt6p(j)+depspt6p(j)) .lt. epsmin) .and. referm3(j)) then
                    depspt6p(j) = epsmin-epspt6p(j)
                    limit1 = .true.
                    do i = 1, 6
                        depsk6p(i) = 0.d0
                        depsm6p(i) = 0.d0
                        if (i .ne. j) then
                            depspt6p(i) = 0.d0
                        end if
                        depspg6p(i) = 0.d0
                        depspc6p(i) = 0.d0
                    end do
!            on utilise pas les autres multiplicateurs
!            car la refermeture a été forcée dans une direction
                    goto 50
                end if
            end do
!
!     recuperation des increments dans la pase principale
!     cas du fluage
            depsk6p(:) = 0.d0
            depsm6p(:) = 0.d0
            if (fl3d) then
!           increment deformation kelvin
                depsk6p(1:6) = x(1:6)
!           increment deformation maxwell
                depsm6p(1:6) = x(7:12)
            end if
!
!     *** increment deformation elastique en base principale *******
50          continue
            do i = 1, 6
                depse6p(i) = -depsm6p(i)-depsk6p(i) &
                             -depspt6p(i)-depspg6p(i)-depspc6p(i)
            end do
!
!     *** retour des increments de deformation dans la base fixe ***
!     plasticite traction
            call chrep6(depspt6p, vsig133t, .false._1, depspt6)
!     elasticite
            call chrep6(depse6p, vsig133t, .false._1, depse6)
!     cas des increments remis a zero si refermeture forcee
            if (.not. limit1) then
                call chrep6(depsk6p, vsig133t, .false._1, depsk6)
                call chrep6(depsm6p, vsig133t, .false._1, depsm6)
                call chrep6(depspg6p, vsig133t, .false._1, depspg6)
                call chrep6(depspc6p, vsig133t, .false._1, depspc6)
            else
!               mise a zero des increments en cas de refermeture forcee
                call iniVect0(6, depsk6, depsm6, depspg6, depspc6)
!               actualisation variable ecrouissage isotrope cisaillement
                depleqc = 0.d0
            end if
        else
!           pas d ecoulement plastique
            call iniVect0(6, depsk6, depsm6, depse6, depspt6, &
                          depspg6, depspc6)
            depleqc = 0.d0
!
            if (fl3d) then
!     actualisation endo de fluage
                call dflufin3d(sig16, bw, pw, bg, pg, &
                               dsw6, delta, rc, xflu, dfin1, &
                               CMp1, dfmx)
!     actualisation coeffs de consolidation
                call conso3d(epsm11, epser, ccmin1, ccmax1, epsm16, &
                             epse16, cc13, vcc33, vcc33t, CWp, &
                             CMp1, CthP, Cthv)
            end if
        end if
!
!    actualisation des deformations (base fixe)
        do i = 1, 6
!       pour les visco elastique l increment du tir visco elastique
!       a deja ete comptabilise avant l ecoulement et une
!       premiere mise a jour a ete faite
!       kelvin
            epsk16(i) = epsk16(i)+depsk6(i)
            varf(EPK(i)) = epsk16(i)
!       maxwell
            epsm16(i) = epsm16(i)+depsm6(i)
            varf(i+12) = epsm16(i)
!       elastique
            epse16(i) = epse16(i)+depse6(i)
            varf(EPE(i)) = epse16(i)
!       cas des deformations plastiques
!       traction base fixe
            epspt6(i) = epspt6(i)+depspt6(i)
            varf(EPTi(i)) = epspt6(i)
!       gel dans les pores
            epspg6(i) = epspg6(i)+depspg6(i)
            varf(EPGi(i)) = epspg6(i)
!       cisaillement et dilatance
            epspc6(i) = epspc6(i)+depspc6(i)
            varf(EPCi(i)) = epspc6(i)
        end do
!    *** actualisation deformation equivalente de compression ******
!    comparaison avec calcul direct par la trace par l invariant
        depleqc3 = 0.d0
        do i = 1, 6
            depleqc3 = depleqc3+depspc6(i)**2
        end do
        depleqc3 = dsqrt(depleqc3*2.d0/3.d0)
!    actualisation et stockage
        epleqc = dmax1(epleqc+depleqc3, epleqc01)
        varf(EPLC) = epleqc
!
!    *** actualisation des ouvertures de fissures ******************
        if (end3d) then
            call majw3d(epspt60, epspt6, t33, wplt06, wplt6, &
                        wpltx06, wpltx6, wpl3, vwpl33, vwpl33t, &
                        wplx3, vwplx33, vwplx33t, local, dim3, &
                        ndim, ifour)
        else
            call iniVect0(6, wplt6, wpltx6)
        end if
!    tenseur des ouvertures
        do j = 1, 6
            varf(WID(j)) = wplt6(j)
            varf(EMT(j)) = wpltx6(j)
        end do
!    Ouverture de fissure maximale Wpl0
        varf(WPL0) = dmax1((varf(WID(1))), (varf(WID(2))), (varf(WID(3))))
!
!    *** actualisation de la dissipation et de l energie elastique *
!    pour calcul consolidation debut de pas suivant
        phi1 = phi0
        we1 = 0.d0
        do i = 1, 6
!        pas de remise a l etat initial on repart
!        du tir visco elastique
            do j = 1, 6
                sig16(i) = sig16(i)+raideur66(i, j)*depse6(j)
                sigke16(i) = sigke16(i)+raideur66(i, j)*depsk6(j)*psik
            end do
!        actualisation de la dissipation avec l increment fin de de pas
!        et non la correction radiale
            dphi = 0.5d0*(sig06(i)+sig16(i))*(epsm16(i)-epsm06(i))
            logic1 = .false.
            if (dphi .lt. 0.d0) logic1 = .true.
            phi1 = phi1+dphi
!        actualisation du potentiel elastique
            we1 = we1+0.5d0*(sig16(i)*epse16(i))
        end do
!
!   *** actualisation variation volumique plastique rgi ************
        trepspg = 0.d0
        do i = 1, 3
            trepspg = trepspg+epspg6(i)
        end do
!   actualisation variation volumique totale dans variable interne
!   dissipation et energie elastique
        call setValVect(varf, trepspg, phi1, we1, vectInd=[TEPG, PHIM, WELA])
!   actualisation des contraintes effectives (sans bgpg)
        do i = 1, 6
            varf(SIG(i)) = sig16(i)
            varf(SKE(i)) = sigke16(i)
        end do
!   *** test de consistance apres ecoulement plastique *************
        if ((na .ne. 0) .or. indic2) then
            if (ipla .le. imax) then
!            nouvelle sous iteration de consistance plastique
                goto 10
            else
!           on vient d atteindre le nbr maxi de sous iteration
                vali(1) = imax
                vali(2) = ipla
                call utmess('E', 'COMPOR3_33', ni=3, vali=vali)
                ierr1 = 1
                go to 999
            end if
        end if
    end if

    call setValVect(outputR, srw, epeqpc, dfin1, ccmax1, pw, bw, &
                    wpl3(1), wpl3(2), wpl3(3), wplx3(1), wplx3(2), &
                    wplx3(3), ind1=1)
    call setR6Mat6(outputVR6, dpg_depsa6, dpg_depspg6, epsk16, &
                   epsm16, epse16, epspt6, sigke16, sig16, epspg6, &
                   epspc6, dsw6)
    call setMat33Tab(outputMat33, vwpl33, vwpl33t, &
                     vwplx33, vwplx33t)
    call setIntVect(outputI, ifour, ierr1, ipla)

999 continue

end subroutine
