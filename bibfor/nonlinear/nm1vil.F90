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

subroutine nm1vil(fami, kpg, ksp, icdmat, materi, &
                  crit, instam, instap, tm, tp, &
                  tref, deps, sigm, vim, option, &
                  defam, defap, angmas, sigp, vip, &
                  dsidep, iret, compo, nbvalc)
!
! aslint: disable=W1504
    implicit none
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/granac.h"
#include "asterfort/nmasse.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: icdmat, kpg, ksp, iret, nbvalc
    real(kind=8) :: crit(*)
    real(kind=8) :: instam, instap
    real(kind=8) :: tm, tp, tref
    real(kind=8) :: irram, irrap
    real(kind=8) :: deps
    real(kind=8) :: sigm, vim(nbvalc)
    character(len=16) :: option, compo
    character(len=*) :: fami
    real(kind=8) :: defam, defap
    real(kind=8) :: angmas(3)
    real(kind=8) :: sigp, vip(nbvalc), dsidep, alpha
    character(len=8) :: materi
!
! ----------------------------------------------------------------------
!      VISCO_PLASTICITE FLUAGE SOUS IRRADIATION AVEC GRANDISSEMENT
!      VISC_IRRA_LOG OU GRAN_IRRA_LOG
!      LOI 1D PURE. MODIF JMP POUR ECRIRE SIMPLEMENT :
! DEPSVP=SIGMA+.EXP(-Q/T)*(A.OMEGA/(1+OMEGA*FLUENCE)+B*FLUENCE)*DFLUENCE
!
! IN  ICDMAT  : MATERIAU CODE
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  TM      : TEMPERATURE A L'INSTANT PRECEDENT
! IN  TP      : TEMPERATURE A L'INSTANT DU CALCUL
! IN  TREF    : TEMPERATURE DE REFERENCE
! IN  DEPS    : INCREMENT DE DEFORMATION-INCREMENT DEFORMATION THERMIQUE
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! IN  DEFAM   : DEFORMATIONS ANELASTIQUES A L'INSTANT PRECEDENT
! IN  DEFAP   : DEFORMATIONS ANELASTIQUES A L'INSTANT DU CALCUL
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MODULE TANGENT
! OUT IRET    : CODE RETOUR DE LA RECHERCHE DE ZERO DE F(X)=0
!                   IRET=0 => PAS DE PROBLEME
!                   IRET=1 => ECHEC
!
!
!
!     COMMON POUR LES PARAMETRES DES LOIS VISCOPLASTIQUES
    common/nmpavp/dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
    real(kind=8) :: dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
!     COMMON POUR LES PARAMETRES DES LOIS DE FLUAGE SOUS IRRADIATION
!     VISC_IRRA_LOG: A      B      CTPS    ENER
    common/nmpair/a, b, ctps, ener
    real(kind=8) :: a, b, c, ctps, ener
! PARAMETRES MATERIAUX
! ELASTIQUES
    real(kind=8) :: ep, nup, troikp, deumup
    real(kind=8) :: em, num, troikm, deumum
! AUTRES
    integer(kind=8) :: nbcgil, iret2
    parameter(nbcgil=5)
    real(kind=8) :: coegil(nbcgil)
    character(len=8) :: nomgil(nbcgil)
    integer(kind=8) :: codgil(nbcgil)
    real(kind=8) :: t1, t2
    real(kind=8) :: degran, depsan, depsim, depsgr
    real(kind=8) :: coef1, coefb, expqt
    real(kind=8) :: fluphi
    data nomgil/'A', 'B', 'CSTE_TPS', 'ENER_ACT', 'C'/
!
    iret = 0
!     PARAMETRE THETA D'INTEGRATION
!
    theta = crit(4)
    t1 = abs(theta-0.5d0)
    t2 = abs(theta-1.d0)
    prec = 0.000001d0
    if ((t1 .gt. prec) .and. (t2 .gt. prec)) then
        call utmess('F', 'ALGORITH6_55')
    end if
!
! TEMPERATURE AU MILIEU DU PAS DE TEMPS  (DANS COMMON / NMPAVP /)
    tschem = tm*(1.d0-theta)+tp*theta
! DEFORMATION PLASTIQUE CUMULEE  (DANS COMMON / NMPAVP /)
    dpc = vim(1)
! INCREMENT DE TEMPS (DANS COMMON / NMPAVP /)
    deltat = instap-instam
! CARACTERISTIQUES ELASTIQUES VARIABLES
    call nmasse(fami, kpg, ksp, '-', icdmat, &
                materi, instam, em, num, deumum, &
                troikm)
!
    call nmasse(fami, kpg, ksp, '+', icdmat, &
                materi, instap, ep, nup, deumup, &
                troikp)
!
!     IRRADIATION AU POINT CONSIDERE
!     FLUX NEUTRONIQUE
    call rcvarc('F', 'IRRA', '-', fami, kpg, &
                ksp, irram, iret2)
    if (iret2 .gt. 0) irram = 0.d0
    call rcvarc('F', 'IRRA', '+', fami, kpg, &
                ksp, irrap, iret2)
    if (iret2 .gt. 0) irrap = 0.d0
    irrap = irrap-irram+vim(2)
    irram = vim(2)
!
    fluphi = (irrap-irram)/deltat
!     RECUPERATION DES CARACTERISTIQUES DES LOIS DE FLUAGE
    call rcvalb(fami, kpg, ksp, '+', icdmat, &
                materi, compo, 0, ' ', [0.d0], &
                nbcgil, nomgil(1), coegil(1), codgil(1), 0)
!     TRAITEMENT DES PARAMETRES DE LA LOI DE FLUAGE
    if (codgil(1) .eq. 0) then
!         LOI DE TYPE VISC_IRRA_LOG
!         PARAMETRES DE LA LOI DE FLUAGE
!
        a = coegil(1)
        b = coegil(2)
        ctps = coegil(3)
        ener = coegil(4)
        if (compo(1:10) .eq. 'GRAN_IRRA_') then
            c = coegil(5)
        else
            c = 0.0d0
        end if
        if (fluphi .lt. -prec) then
            call utmess('F', 'ALGORITH6_57')
        end if
    else
        call utmess('F', 'ALGORITH6_58')
    end if
!
!     CALCUL DE LA DEFORMATION DE GRANDISSEMENT
    degran = 0.0d0
    call granac(fami, kpg, ksp, icdmat, materi, &
                compo, irrap, irram, tm, tp, &
                depsgr)
!
    if (compo(1:10) .eq. 'GRAN_IRRA_') then
        vip(3) = vim(3)+depsgr
        if (depsgr .ne. 0.0d0) then
! --- RECUPERATION DU REPERE POUR LE GRANDISSEMENT
            alpha = angmas(1)
            if (angmas(2) .ne. 0.d0) then
                call utmess('F', 'ALGORITH6_59')
            end if
!
!        INCREMENT DEFORMATION DE GRANDISSEMENT DANS LE REPERE
            degran = depsgr*cos(alpha)*cos(alpha)
        end if
    end if
!     INCREMENT DEFORMATION ANELASTIQUE
    depsan = defap-defam
!     INCREMENT DEFORMATION IMPOSEE
    depsim = depsan+degran
!
    expqt = exp(-ener/(tp+r8t0()))
!
!     coefb=expqt*((a*ctps/(1.d0+ctps*irrap))+b+c*ctps*exp(-ctps*irrap))*(irrap-irram)
    coefb = expqt*(a*(log(1.d0+ctps*irrap)-log(1.d0+ctps*irram))+ &
                   b*(irrap-irram)+ &
                   c*(exp(-ctps*irram)-exp(-ctps*irrap)))
    coef1 = ep/(1.d0+ep*coefb)
!
! CONTRAINTE ACTUALISEE
!
!
    sigp = coef1*(sigm/em+deps-depsim)
!
! DEFORMATION PLASTIQUE CUMULEE ACTUALISEE
!
    vip(1) = vim(1)+(abs(sigp)*coefb)
! MODULE TANGENT POUR MATRICE TANGENTE
!
    dsidep = coef1
    vip(2) = irrap
!
end subroutine
