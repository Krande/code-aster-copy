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
subroutine lceiab(fami, kpg, ksp, mat, option, &
                  mu, su, de, ddedt, vim, &
                  vip, r, codret)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
    integer(kind=8) :: mat, kpg, ksp, codret
    real(kind=8) :: mu(3), su(3), de(6), ddedt(6, 6)
    real(kind=8) :: vim(*), vip(*), r
    character(len=16) :: option
    character(len=*) :: fami
!
!-----------------------------------------------------------------------
!            LOI DE COMPORTEMENT COHESIVE CZM_LAB_MIX
!            POUR LES ELEMENTS D'INTERFACE 2D ET 3D.
!
! IN : FAMI,KPG,KSP,MAT,OPTION
!      MU  : LAGRANGE
!      SU  : SAUT DE U
!      VIM : VARIABLES INTERNES
!
! OUT : DE    : DELTA
!       DDEDT : DERIVEE DE DELTA
!       VIP   : VARIABLES INTERNES MISES A JOUR
!       R     : PENALISATION DU LAGRANGE
!-----------------------------------------------------------------------
!
    aster_logical :: resi, rigi, elas
    integer(kind=8) :: regime, regm, i, j, cod(6), cinema
    real(kind=8) :: sc, dc, alpha, beta, s0, d0, ka, sk, val(6)
    real(kind=8) :: t(3), pr(3, 3), tpo(3), tpon(3), tno, tnon
    real(kind=8) :: dno, ddno, dnon, dsidno
    real(kind=8) :: res, rn, deriv, bmin, bmax
    character(len=16) :: nom(6)
    character(len=1) :: poum
!
    data nom/'SIGM_C', 'GLIS_C', 'ALPHA', 'BETA', 'PENA_LAGR',&
     &          'CINEMATIQUE'/
!-----------------------------------------------------------------------
!
! ---------------------------
! -- PRINCIPALES NOTATIONS --
! ---------------------------
!
! -  CARACTERISTIQUES DE LA ZONE COHESIVE
!    SC     : CONTRAINTE CRITIQUE
!    DC     : SAUT DE DEPLACEMENT A LA CONTRAINTE CRITIQUE
!    ALPHA  : PARAMETRE DE FORME DE LA LOI (PETITS GLISSEMENTS)
!    BETA   : PARAMETRE DE FORME DE LA LOI (GRANDS GLISSEMENTS)
!    R      : PARAMETRE DE PENALISATION
!
! -  DONNEES D'ENTREE
!    MU     : LAGRANGE
!    SU     : SAUT DE U
!    VIM    : VARIABLES INTERNES
!             |1   : PLUS GRANDE NORME DU SAUT (KA)
!             |2   : REGIME DE LA LOI (REGM)
!             |      |0 : ADHERENCE INITIALE OU COURANTE
!             |      |1 : DISSIPATION
!             |      |2 : SURFACE LIBRE FINALE (RUPTURE)
!             |      |3 : SURFACE LIBRE (SOUS CONTRAINTE)
!             |3-5 : VALEURS DE DELTA
!
! -  DONNEES DE SORTIE
!    DE     : DELTA CALCULE
!    DDEDT  : DERIVEE DE DELTA
!    VIP    : VARIABLES INTERNES MISES A JOUR
!
! -  GRANDEURS LOCALES
!    REGM   : REGIME DE FONCTIONNEMENT DE LA LOI A L'INSTANT PRECEDENT
!    REGIME : NOUVEAU REGIME DE FONCTIONNEMENT
!    KA     : OUVERTURE MAXIMALE COURANTE
!    SK     : CONTRAINTE CRITIQUE COURANTE
!    T      : FORCE COHESIVE LAMBDA + R.[U]
!    PR     : MATRICE DE PROJECTION SUIVANT LA CINEMATIQUE
!    TPO    : FORCE COHESIVE PROJETEE
!    TNO    : NORME DE LA FORCE COHESIVE PROJETEE
!
! --------------------
! -- INITIALISATION --
! --------------------
!
!    OPTION CALCUL DU RESIDU OU CALCUL DE LA MATRICE TANGENTE
    resi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH'
    rigi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RIGI'
    elas = option(11:14) .eq. 'ELAS'
!
!    RECUPERATION DES PARAMETRES PHYSIQUES
    if (option .eq. 'RIGI_MECA_TANG') then
        poum = '-'
    else
        poum = '+'
    end if
!
    call rcvalb(fami, kpg, ksp, poum, mat, &
                ' ', 'CZM_LAB_MIX', 0, ' ', [0.d0], &
                6, nom, val, cod, 2)
!
    sc = val(1)
    dc = val(2)
    alpha = val(3)
    beta = val(4)
!
    d0 = dc*beta/alpha
    s0 = sc*(alpha+beta)**(alpha+beta)/(alpha**alpha*beta**beta)
    r = val(5)*sc/dc
!
!    ENTIER DECRIVANT LA CINEMATIQUE DU COMPORTEMENT DE L INTERFACE
!     (CODE DANS LA ROUTINE RCSTOC)
    cinema = nint(val(6))
!
!    LECTURE DES VARIABLES INTERNES
    ka = vim(1)
    regm = nint(vim(2))
    sk = max(0.d0, s0*(ka/d0)**alpha/(ka/d0+1.d0)**(alpha+beta))
!
! -----------------------------
! -- CALCUL DU SECOND MEMBRE --
! -----------------------------
!
!    FORCE COHESIVE AUGMENTEE : LAMBDA + R.[U]
    t(1) = mu(1)+r*su(1)
    t(2) = mu(2)+r*su(2)
    t(3) = mu(3)+r*su(3)
!
!    PROJECTEUR POUR UNE COMPOSANTE NORMALE POSITIVE
    call r8inir(9, 0.d0, pr, 1)
!
    if ((cinema .eq. 0) .and. (t(1) .ge. 0.d0)) pr(1, 1) = 1.d0
    pr(2, 2) = 1.d0
    if ((cinema .eq. 0) .or. (cinema .eq. 2)) pr(3, 3) = 1.d0
!
!    PROJECTION DE LA COMPOSANTE NORMALE POSITIVE
    tpo(1) = t(1)*pr(1, 1)
    tpo(2) = t(2)*pr(2, 2)
    tpo(3) = t(3)*pr(3, 3)
!
!    NORME DU SECOND MEMBRE PROJETE
    tno = sqrt(tpo(1)**2+tpo(2)**2+tpo(3)**2)
!
!    VECTEUR UNITE DIRECTION DE FORCE
    if (tno .gt. 0.d0) then
        tpon(1) = tpo(1)/tno
        tpon(2) = tpo(2)/tno
        tpon(3) = tpo(3)/tno
    else
!       SI LA FORCE EST NULLE, LE VECTEUR DIRECTION EST ARBITRAIRE
!                                                   (MAIS ADMISSIBLE)
        tpon(1) = 0.d0
        tpon(2) = 1.d0
        tpon(3) = 0.d0
    end if
!
! --------------------------------------------
! -- RESOLUTION DU PROBLEME 1D SUR LA NORME --
! --------------------------------------------
!
!    DETERMINATION DU REGIME DE COMPORTEMENT
    if (resi) then
!
!      SURFACE LIBRE (SOUS CONTRAINTE)
        if (tno .lt. r*ka) then
            regime = 3
!
!      ADHERENCE (INITIALE OU COURANTE)
        else if (tno .le. r*ka+sk) then
            regime = 0
!
!      ENDOMMAGEMENT
        else
            regime = 1
        end if
!
!    SINON, ON N'ACTUALISE PAS LE REGIME DE FONCTIONNEMENT DE LA LOI
    else
        regime = regm
    end if
!
    codret = 0
!
!    CALCUL DE L'ECOULEMENT 1D SELON LE REGIME DE COMPORTEMENT
    if (regime .eq. 3) then
        dno = tno/r
    else if (regime .eq. 0) then
        dno = ka
    else
!     DANS LE CAS GENERAL, UTILISATION D UN ALGORITHME DE NEWTON
!     1 - ADIMENSIONNALISATION DES VARIABLES
        tnon = tno/s0
        rn = r*d0/s0
!
!     2 - DETERMINATION DE BORNES BMIN ET BMAX POUR NEWTON, AINSI
!         QUE D UN POINT D INITIALISATION JUDICIEUX
        dnon = alpha/beta
100     continue
        res = dnon**alpha/(dnon+1.d0)**(alpha+beta)-(tnon-rn*dnon)
        if (res .lt. 0) goto 110
        dnon = dnon/100.d0
        goto 100
110     continue
        bmin = dnon
        bmax = tnon/rn
        dnon = bmax
!
!     3 - BOUCLE DE CONVERGENCE DE L ALGORITHME DE NEWTON
        i = 0
200     continue
!         TEST DU CRITERE
        res = dnon**alpha/(dnon+1.d0)**(alpha+beta)-(tnon-rn*dnon)
        if (abs(res/rn) .lt. 1.d-12) goto 210
!
!         DIAGNOSTIC DE NON-CONVERGENCE
        if (i .ge. 20) then
            codret = 1
            goto 999
        end if
        i = i+1
!
!         NOUVEL ESTIMATEUR
        deriv = (alpha-beta*dnon)*dnon**(alpha-1.d0)/(dnon+1.d0)**(alpha+beta+1.d0)+rn
        dnon = dnon-res/deriv
!
!         PROJECTION SUR LES BORNES DE L'INTERVALLE
        if (dnon .lt. bmin) dnon = bmin
        if (dnon .gt. bmax) dnon = bmax
        goto 200
!
210     continue
        dno = dnon*d0
    end if
!
! ------------------------------------
! -- CONSTRUCTION DE LA SOLUTION 3D --
! ------------------------------------
!
!    CALCUL DU SAUT DE DEPLACEMENT 3D
    if (resi) then
        call r8inir(6, 0.d0, de, 1)
        de(1) = dno*tpon(1)
        de(2) = dno*tpon(2)
        de(3) = dno*tpon(3)
    end if
!
!    MISE A JOUR DES VARIABLES INTERNES
    if (resi) then
        vip(1) = max(ka, dno)
        vip(2) = regime
        vip(3) = de(1)
        vip(4) = de(2)
        vip(5) = de(3)
    end if
!
! ----------------------
! -- MATRICE TANGENTE --
! ----------------------
!
    if (rigi) then
!
!      AJUSTEMENT POUR PRENDRE EN COMPTE *_MECA_ELAS
        if (elas) then
            if (regime .eq. 1) regime = 0
        end if
!
!      CALCUL DU COEFFICIENT 1D DE LA MATRICE TANGENTE
        if (regime .eq. 3) then
            ddno = 1.d0/r
        else if (regime .eq. 0) then
            ddno = 0.d0
        else
            if (dno .gt. 0.d0) then
                dnon = dno/d0
                dsidno = s0/d0*(alpha-beta*dnon)*dnon**(alpha-1.d0)/(dnon+1.d0)**(alpha+beta+1.d&
                         &0)
                ddno = 1.d0/(r+dsidno)
            else
                ddno = 0.d0
            end if
        end if
!
!      MATRICE TANGENTE 3D
        call r8inir(36, 0.d0, ddedt, 1)
        if (tno .gt. 0.d0) then
            do i = 1, 3
                do j = 1, 3
                    ddedt(i, j) = ddno*tpon(i)*tpon(j)+dno/tno*(pr(i, j)-tpon(i)*tpon(j))
                end do
            end do
        else
!         CAS OU TNO EST RIGOUREUSEMENT NUL
            do i = 1, 3
                do j = 1, 3
                    ddedt(i, j) = ddno*tpon(i)*tpon(j)
                end do
            end do
        end if
!
    end if
!
999 continue
!
end subroutine
