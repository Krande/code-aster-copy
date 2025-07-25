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
subroutine lcmaza(fami, kpg, ksp, ndim, typmod, &
                  imate, compor, epsm, deps, vim, &
                  option, sig, vip, dsidep)
    implicit none
#include "asterf_types.h"
#include "asterc/r8nnem.h"
#include "asterfort/bptobg.h"
#include "asterfort/get_varc.h"
#include "asterfort/jacobi.h"
#include "asterfort/lcumvi.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    character(len=8) :: typmod(2)
    character(len=16) :: compor(*), option
    character(len=*) :: fami
    integer(kind=8) :: ndim, imate, kpg, ksp
    real(kind=8) :: epsm(6), deps(6), vim(4)
    real(kind=8) :: sig(6), vip(*), dsidep(6, 6)
! ----------------------------------------------------------------------
!     LOI DE COMPORTEMENT ENDOMMAGEABLE : MODELE DE MAZARS
!     POUR MAZARS  OU MAZARS_FO COMBINABLE AVEC ELAS OU ELAS_FO
!     NB. LES PARAMETRES MATERIAUX PEUVENT DEPENDRE DE LA TEMPERATURE,
!      DE L'HYDRATATION OU DU SECHAGE
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : NATURE DU MATERIAU
! IN  EPSM    : DEFORMATION EN T-
! IN  DEPS    : INCREMENT DE DEFORMATION
! IN  VIM     : VARIABLES INTERNES EN T-
! IN  TM      : TEMPERATURE EN T-
! IN  TP      : TEMPERATURE EN T+
! IN  TREF    : TEMPERATURE DE REFERENCE
! IN  OPTION  : OPTION DEMANDEE
!                 RIGI_MECA_TANG ->     DSIDEP
!                 FULL_MECA      -> SIG DSIDEP VIP
!                 RAPH_MECA      -> SIG        VIP
!
! MODIFIE 2008 M. BOTTONI
! ON AJOUTE OPTION  RAPH_COUP POUR TRAITER LE COUPLAGE AVEC UMLV
! IN OPTION  : OPTION DEMANDEE
!                 RIGI_COUP      ->     DSIDEP
!                 RAPH_COUP      -> SIG        VIP
!
! OUT SIG     : CONTRAINTE
! OUT VIP     : VARIABLES INTERNES
!                 1   -> VALEUR DE L'ENDOMMAGEMENT
!                 2   -> INDICATEUR D'ENDOMMAGEMENT
!                 3   -> TEMPERATURE MAXIMALE ATTEINTE PAR LE MATERIAU
!                 4   -> VALEUR DE EPSEQ (UTILE POUR POSTTRAITER)
! OUT DSIDEP  : MATRICE TANGENTE
! ON A BESOIN DE
!         EPSD0 = DEFORMATION SEUIL  [REEL OU FCT]
!         AT = CONSTANTE DE TRACTION     (0.7 A 1)[REEL OU FCT]
!         AC = CONSTANTE DE COMPRESSION (1 A 1.5)[REEL OU FCT]
!         BT = CONSTANTE DE TRACTION    (10 000 A 100 000)[REEL OU FCT]
!         BC = CONSTANTE DE COMPRESSION (1000 A 2000)[REEL OU FCT]
! ----------------------------------------------------------------------
    aster_logical :: rigi, resi, prog, elas, cplan, coup
    character(len=1) :: poum
    integer(kind=8) :: icodre(7)
    character(len=16) :: nomres(7)
    character(len=8) :: nompar
    integer(kind=8) :: ndimsi, nperm, nitjac, trij, ordrej
    integer(kind=8) :: i, j, l, iret
    real(kind=8) :: e, nu, epsthe, kdess, bendo, rtemp
    real(kind=8) :: ac, at, bc, bt, epsd0, tm, tp, tref
    real(kind=8) :: eps(6), epse(6), epsplu(6), epsep(3), epseq
    real(kind=8) :: sigel(6), sigelp(3)
    real(kind=8) :: temp, tmax, tmaxm, hydr, sech, sref
    real(kind=8) :: tol, toldyn, tr(6), tu(6), jacaux(3), vecpe(3, 3)
    real(kind=8) :: rac2, lambda, deuxmu, coef
    real(kind=8) :: valres(7), valpar, coplan, d, tmp1, vala, r, a, b
    real(kind=8) :: kron(6), k, y
    real(kind=8) :: epsfp(6), epscou(6), epsi(6), chi, gama, rap
    real(kind=8) :: epseqc, epsend, epsepc(3), vecpec(3, 3)
    integer(kind=8) :: idc
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
! ======================================================================
!                            INITIALISATION
! ======================================================================
!
!
!
! - Get temperatures
!
    call get_varc(fami, kpg, ksp, 'T', tm, &
                  tp, tref)
!
! -- OPTION ET MODELISATION
!
    rigi = (option(1:4) .eq. 'RIGI' .or. option(1:4) .eq. 'FULL')
    resi = (option(1:4) .eq. 'RAPH' .or. option(1:4) .eq. 'FULL')
! M.B.: NOUVELLE OPTION COUP POUR LE COUPLAGE AVEC UMLV
! MEME OPTION UTILISEE LE COUPLAGE UMLV-ENDO_ISOT_BETON
    coup = (option(6:9) .eq. 'COUP')
    cplan = (typmod(1) .eq. 'C_PLAN  ')
    prog = .false.
    elas = .true.
    ndimsi = 2*ndim
    rac2 = sqrt(2.d0)
! M.B.: INDICE POUR IDENTIFIER LES VARIABLES INTERNES DANS LES CAS:
! COUPLAGE ET ABSENCE DE COUPLAGE AVEC UMLV
    idc = 0
    if (coup) then
        idc = 21
    end if
!   DETERMINATION DE LA TEMPERATURE DE REFERENCE (TMAX) ET
!   DES CONDITIONS D HYDRATATION OU DE SECHAGE
    tmaxm = vim(3)
    call rcvarc(' ', 'SECH', 'REF', fami, kpg, &
                ksp, sref, iret)
    if (iret .ne. 0) sref = 0.d0
    if (resi) then
        temp = tp
        call rcvarc(' ', 'HYDR', '+', fami, kpg, &
                    ksp, hydr, iret)
        if (iret .ne. 0) hydr = 0.d0
        poum = '+'
        call rcvarc(' ', 'SECH', '+', fami, kpg, &
                    ksp, sech, iret)
        if (iret .ne. 0) sech = 0.d0
        if (isnan(tp)) then
            tmax = r8nnem()
            vip(idc+3) = 0.d0
        else
            tmax = max(tmaxm, tp)
            if (tmax .gt. tmaxm) vip(idc+3) = tmax
        end if
    else
        temp = tm
        call rcvarc(' ', 'HYDR', '-', fami, kpg, &
                    ksp, hydr, iret)
        if (iret .ne. 0) hydr = 0.d0
        call rcvarc(' ', 'SECH', '-', fami, kpg, &
                    ksp, sech, iret)
        if (iret .ne. 0) sech = 0.d0
        poum = '-'
        tmax = tmaxm
    end if
!  RECUPERATION DES CARACTERISTIQUES MATERIAUX QUI PEUVENT VARIER
!  AVEC LA TEMPERATURE (MAXIMALE), L'HYDRATATION OU LE SECHAGE
    nompar = 'TEMP'
    valpar = tmax
!    LECTURE DES CARACTERISTIQUES ELASTIQUES
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'ELAS', 1, nompar, [valpar], &
                2, nomres, valres, icodre, 1)
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'ELAS', 1, nompar, [valpar], &
                1, nomres(3), valres(3), icodre(3), 0)
    if ((.not. isnan(tp)) .and. (.not. isnan(tm))) then
        if ((isnan(tref)) .or. (icodre(3) .ne. 0)) then
            call utmess('F', 'CALCULEL_15')
        else
            epsthe = valres(3)*(temp-tref)
        end if
    else
        valres(3) = 0.d0
        epsthe = 0.d0
    end if
    e = valres(1)
    nu = valres(2)
    lambda = e*nu/(1.d0+nu)/(1.d0-2.d0*nu)
    deuxmu = e/(1.d0+nu)
! --- LECTURE DU RETRAIT ENDOGENE ET RETRAIT DE DESSICCATION
    nomres(1) = 'B_ENDOGE'
    nomres(2) = 'K_DESSIC'
    call rcvala(imate, ' ', 'ELAS', 0, ' ', &
                [0.d0], 1, nomres(1), valres(1), icodre(1), &
                0)
    if (icodre(1) .ne. 0) valres(1) = 0.d0
    bendo = valres(1)
    call rcvala(imate, ' ', 'ELAS', 0, ' ', &
                [0.d0], 1, nomres(2), valres(2), icodre(2), &
                0)
    if (icodre(2) .ne. 0) valres(2) = 0.d0
    kdess = valres(2)
!
!    LECTURE DES CARACTERISTIQUES D'ENDOMMAGEMENT
    nomres(1) = 'EPSD0'
    nomres(2) = 'AC'
    nomres(3) = 'BC'
    nomres(4) = 'AT'
    nomres(5) = 'BT'
    nomres(6) = 'K'
    call rcvalb(fami, kpg, ksp, poum, imate, &
                ' ', 'MAZARS', 1, nompar, [valpar], &
                6, nomres, valres, icodre, 1)
    epsd0 = valres(1)
    ac = valres(2)
    bc = valres(3)
    at = valres(4)
    bt = valres(5)
    k = valres(6)
!    M.B.: LECTURE DU PARAMETRE DE COUPLAGE AVEC UMLV
    if (coup) then
        nomres(7) = 'CHI'
        call rcvalb(fami, kpg, ksp, poum, imate, &
                    ' ', 'MAZARS', 0, ' ', [0.d0], &
                    1, nomres(7), valres(7), icodre(7), 1)
        chi = valres(7)
    end if
! ======================================================================
!       CALCUL DES GRANDEURS UTILES QUELQUE SOIT OPTION
! ======================================================================
!    1 - CALCUL DES DEFORMATIONS MECANIQUES ET THERMIQUES
!  -  MISE A JOUR DE LA DEFORMATION TOTALE
    call r8inir(6, 0.d0, eps, 1)
    if (resi) then
        do j = 1, ndimsi
            eps(j) = epsm(j)+deps(j)
        end do
    else
        do j = 1, ndimsi
            eps(j) = epsm(j)
        end do
        d = vim(1)
    end if
!    CALCUL DE LA DEFORMATION ELASTIQUE (LA SEULE QUI CONTRIBUE
!    A FAIRE EVOLUER L'ENDOMMAGEMENT)
    call r8inir(6, 0.d0, epse, 1)
    do j = 1, ndimsi
        epse(j) = eps(j)-(epsthe-kdess*(sref-sech)-bendo*hydr)*kron(j)
    end do
    if (cplan) then
        coplan = -nu/(1.d0-nu)
        epse(3) = coplan*(epse(1)+epse(2))
    end if
!    M.B.: AVEC COUPLAGE, EPSF POUR L INSTANT
!    SERT SEULEMENT AVEC  RESI A L INSTANT P, CAR LA MATRICE TANGENTE
!    N EST PAS ENCORE IMPLEMENTEE
    if (coup .and. resi) then
        call lcumvi('FT', vip, epsfp)
        call r8inir(6, 0.d0, epscou, 1)
        do j = 1, ndimsi
            epsi(j) = epse(j)
            epse(j) = epsi(j)-epsfp(j)
            epscou(j) = epsi(j)-(1-chi)*epsfp(j)
        end do
    end if
    do j = 4, ndimsi
        epse(j) = epse(j)/rac2
    end do
    if (coup .and. resi) then
        do j = 4, ndimsi
            epscou(j) = epscou(j)/rac2
        end do
    end if
!    2 - CALCUL DE EPSEQ = SQRT(TR (<EPSE>+ * <EPSE>+)  )
!        C EST EPSEQ ELASTIQUE DANS LE CAS DU COUPLAGE
!--------------------------------------------------------
!  -   ON PASSE DANS LE REPERE PROPRE DE EPS
    nperm = 12
    tol = 1.d-10
    toldyn = 1.d-2
!     MATRICE  TR = (XX XY XZ YY YZ ZZ) POUR JACOBI)
    tr(1) = epse(1)
    tr(2) = epse(4)
    tr(3) = epse(5)
    tr(4) = epse(2)
    tr(5) = epse(6)
    tr(6) = epse(3)
!     MATRICE UNITE = (1 0 0 1 0 1) (POUR JACOBI)
    tu(1) = 1.d0
    tu(2) = 0.d0
    tu(3) = 0.d0
    tu(4) = 1.d0
    tu(5) = 0.d0
    tu(6) = 1.d0
    trij = 2
    ordrej = 2
!
    call jacobi(3, nperm, tol, toldyn, tr, &
                tu, vecpe, epsep, jacaux, nitjac, &
                trij, ordrej)
    epseq = 0.d0
    do j = 1, 3
        if (epsep(j) .gt. 0.d0) then
            epseq = epseq+(epsep(j)**2)
        end if
    end do
    epseq = sqrt(epseq)
!    2BIS - CALCUL DE EPSEQC = SQRT(TR (<EPSCOU>+ * <EPSCOU>+))
!        M.B.: C EST LA DEFORMATION EQUIVALENT DANS LE CAS DU COUPLAGE
!--------------------------------------------------------
!  -   ON PASSE DANS LE REPERE PROPRE DE EPSCOU
!     MATRICE  TR = (XX XY XZ YY YZ ZZ) POUR JACOBI)
!     A CONTROLER SI LES QUANTITES SUIVANTES SERVENT AUSSI POUR LA
!      MATRICE TANGENTE  (RIGI)!
    if (coup .and. resi) then
        tr(1) = epscou(1)
        tr(2) = epscou(4)
        tr(3) = epscou(5)
        tr(4) = epscou(2)
        tr(5) = epscou(6)
        tr(6) = epscou(3)
!      MATRICE UNITE = (1 0 0 1 0 1) (POUR JACOBI)
        tu(1) = 1.d0
        tu(2) = 0.d0
        tu(3) = 0.d0
        tu(4) = 1.d0
        tu(5) = 0.d0
        tu(6) = 1.d0
        call jacobi(3, nperm, tol, toldyn, tr, &
                    tu, vecpec, epsepc, jacaux, nitjac, &
                    trij, ordrej)
        epseqc = 0.d0
        do j = 1, 3
            if (epsepc(j) .gt. 0.d0) then
                epseqc = epseqc+(epsepc(j)**2)
            end if
        end do
        epseqc = sqrt(epseqc)
    end if
! -  3     CALCUL DE <EPS>+
! ------------------------------------------------------
    call r8inir(6, 0.d0, tr, 1)
    call r8inir(6, 0.d0, epsplu, 1)
    do j = 1, 3
        if (epsep(j) .gt. 0.d0) then
            tr(j) = epsep(j)
        end if
    end do
    call bptobg(tr, epsplu, vecpe)
    do j = 4, ndimsi
        epsplu(j) = epsplu(j)*rac2
    end do
!   4 -  CALCUL DES CONTRAINTES ELASTIQUES (REPERE PRINCIPAL)
!----------------------------------------------------------------
    do j = 1, 3
        sigelp(j) = lambda*(epsep(1)+epsep(2)+epsep(3))
    end do
    do j = 1, 3
        sigelp(j) = sigelp(j)+deuxmu*epsep(j)
    end do
!
    tmp1 = 0.d0
    do j = 1, 3
        if (sigelp(j) .lt. 0.d0) then
            tmp1 = tmp1+sigelp(j)
        end if
    end do
!   5 -     CALCUL DE R
!----------------------------------------------------------------
    vala = abs(sigelp(1))+abs(sigelp(2))+abs(sigelp(3))
    r = 0.d0
    do i = 1, 3
        r = r+max(0.00000000d0, sigelp(i))
    end do
    if (vala .gt. 1.d-10) then
        r = r/(vala)
    else
        r = 1.d0
    end if
    if (r .lt. 0.00001d0) r = 0.d0
    if (r .gt. 0.99999d0) r = 1.d0
    gama = 0.d0
    rap = 0.d0
    do i = 1, 3
        rap = rap+min(0.d0, sigelp(i))
        gama = gama+(min(0.d0, sigelp(i)))**2
    end do
    if ((abs(rap) .gt. 1.d-10) .and. (r .eq. 0.d0)) then
        gama = -(sqrt(gama)/rap)
    else
        gama = 1.d0
    end if
! ======================================================================
!       CALCUL DES CONTRAINTES ET VARIABLES INTERNES
!           (OPTION FULL_MECA ET RAPH_MECA - (RESI) )
! ====================================================================
    if (resi) then
!    M.B.: EPSEND est la deformation equivalente
!    qui fait evoluer l endommagement
        if (coup) then
            epsend = epseqc
        else
            epsend = epseq
        end if
        if (gama .le. 0.d0) gama = 1.d0
        y = gama*epsend
        if (y .le. epsd0) then
!         PAS DE PROGRESSION DE L'ENDOMMAGEMENT
            d = vim(1)
        else
            a = 2.d0*r**2.d0*(at-2.d0*k*at+ac)-r*(at*(1.d0-4.d0*k)+3.d0* &
                                                  ac)+ac
            b = r**2.d0*bt+(1.d0-r**2.d0)*bc
            rtemp = b*(y-epsd0)
            if (rtemp .le. 200.0d0) then
                d = 1.d0-epsd0*(1.d0-a)/y-a*exp(-rtemp)
            else
                d = 1.d0-epsd0*(1.d0-a)/y
            end if
            d = max(vim(1), d)
            d = min(d, 0.99999d0)
            if (d .gt. vim(1)) prog = .true.
            if (d .gt. 0.d0) elas = .false.
        end if
!    2 -   MISE A JOUR DES VARIABLES INTERNES
! ------------------------------------------------------------
        vip(idc+1) = d
        if (d .eq. 0.d0) then
            vip(idc+2) = 0.d0
        else
            vip(idc+2) = 1.d0
        end if
        vip(idc+4) = epsend
!    3 - CALCUL DES CONTRAINTES
! ------------------------------------------------------------
!        ON PASSE DANS LE REPERE INITIAL LES CONTRAINTES REELLES
        call r8inir(6, 0.d0, sig, 1)
        tr(1) = sigelp(1)*(1.d0-d)
        tr(2) = sigelp(2)*(1.d0-d)
        tr(3) = sigelp(3)*(1.d0-d)
        tr(4) = 0.d0
        tr(5) = 0.d0
        tr(6) = 0.d0
        call bptobg(tr, sig, vecpe)
        do j = 4, ndimsi
            sig(j) = rac2*sig(j)
        end do
    end if
! ======================================================================
!     CALCUL  DE LA MATRICE TANGENTE DSIDEP
!         OPTION RIGI_MECA_TANG ET FULL_MECA  (RIGI)
! ======================================================================
    if (rigi) then
! - M.B.: OPTION FULL_MECA POUR LE COUPLAGE AVEC UMLV
        if (coup) d = vip(idc+1)
!   1 -  CONTRIBUTION ELASTIQUE
! -------------------------------------------------------------
        call r8inir(36, 0.d0, dsidep, 1)
        lambda = lambda*(1.d0-d)
        deuxmu = deuxmu*(1.d0-d)
        dsidep(1, 1) = lambda+deuxmu
        dsidep(2, 2) = lambda+deuxmu
        dsidep(3, 3) = lambda+deuxmu
        dsidep(1, 2) = lambda
        dsidep(2, 1) = lambda
        dsidep(1, 3) = lambda
        dsidep(3, 1) = lambda
        dsidep(2, 3) = lambda
        dsidep(3, 2) = lambda
        dsidep(4, 4) = deuxmu
        dsidep(5, 5) = deuxmu
        dsidep(6, 6) = deuxmu
!   2 -  CONTRIBUTION DUE A  L'ENDOMMAGEMENT
!             ON SYMETRISE LA MATRICE (J + Kt )/2
! ------------------------------------------------------------
        if ((.not. elas) .and. prog .and. (d .lt. 0.99999d0)) then
            if (epseq .lt. 0.0000001d0) then
                coef = 0.d0
            else
                rtemp = b*((gama*epseq)-epsd0)
                if (rtemp .le. 200.0d0) then
                    coef = epsd0*(1.d0-a)/(gama*epseq)**2+a*b/exp( &
                           rtemp)
                else
                    coef = epsd0*(1.d0-a)/(gama*epseq)**2
                end if
                coef = coef/epseq
                if (r .eq. 0.d0) coef = gama*coef
            end if
            call r8inir(6, 0.d0, sigel, 1)
            tr(1) = sigelp(1)
            tr(2) = sigelp(2)
            tr(3) = sigelp(3)
            tr(4) = 0.d0
            tr(5) = 0.d0
            tr(6) = 0.d0
            call bptobg(tr, sigel, vecpe)
            do j = 4, ndimsi
                sigel(j) = rac2*sigel(j)
            end do
            do i = 1, 6
                do j = 1, 6
                    dsidep(i, j) = dsidep(i, j)-coef*sigel(i)* &
                                   epsplu(j)
                end do
            end do
! -- CORRECTION CONTRAINTES PLANES
            if (cplan) then
                do j = 1, ndimsi
                    if (j .eq. 3) goto 300
                    do l = 1, ndimsi
                        if (l .eq. 3) goto 310
                        if (dsidep(3, 3) .ne. 0.d0) then
                            dsidep(j, l) = dsidep(j, l)-1.d0/dsidep(3, 3) &
                                           *dsidep(j, 3)*dsidep(3, l)
                        end if
310                     continue
                    end do
300                 continue
                end do
            end if
        end if
    end if
end subroutine
