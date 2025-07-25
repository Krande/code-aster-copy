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
subroutine lcumsf(sigi, sigf, nstrs, vari, nvari, &
                  cmat, nmat, isph, tdt, hini, &
                  hfin, varf)
!
!
! ROUTINE APPELE DANS FLU
! LCUMSF     SOURCE    BENBOU   01/03/26
!
!_______________________________________________________________________
!
! ROUTINE QUI CALCULE LA DEFORMATION DE FLUAGE DE DESSICCATION
! ET DE FLUAGE PROPRE A LA FIN DU PAS DE TEMPS
! ELLES SONT SAUVEGARDES DANS VARF(NVARI)
!
! IN  SIGI     : CONTRAINTES INITIALES
! IN  SIGF     : CONTRAINTES FINALES
! IN  NSTRS    : DIMENSION DES VECTEURS CONTRAINTE ET DEFORMATION
! IN  VARI     : VARIABLES INTERNES INITIALES
! IN  NVARI    : DIMENSION DES VECTEURS VARIABLES INTERNES
! IN  CMAT     : VECTEUR DE PARAMETRES (MATERIAU ET AUTRE)
! IN  NMAT     : DIMENSION DE CMAT
! IN  ISPH     : MODE DE CALCUL DE LA PARTIE SPHERIQUE
! IN  TDT      : PAS DE TEMPS
! IN  HINI     : HUMIDITE INITIALE
! IN  HFIN     : HUMIDITE FINALE
!
! OUT VARF     : VARIABLES INTERNES FINALES
!_______________________________________________________________________
!
    implicit none
#include "asterc/r8prem.h"
#include "asterfort/lcumfb.h"
#include "asterfort/lcumfd.h"
#include "asterfort/lcumfs.h"
    integer(kind=8) :: i, icou, ifou, ifpo, ides, isph, j, nmat, nstrs, nvari
    real(kind=8) :: ersp, eisp
    real(kind=8) :: varf(20), cmat(15), sigf(nstrs), sigi(nstrs)
! MODIFI DU 6 JANVIER 2003 - YLP SUPPRESSION DES DECLARATIONS
! IMPLICITES DES TABLEAUX
!      REAL*8 AFD(NSTRS),BFD(NSTRS,NSTRS),CFD(NSTRS,NSTRS)
!      REAL*8 AFPD(NSTRS),SIGFI(NSTRS),VARI(20)
    real(kind=8) :: afd(6), bfd(6, 6), cfd(6, 6)
    real(kind=8) :: afpd(6), sigfi(6), vari(20)
    real(kind=8) :: efde(6), epsdvr(6), epsdvi(6), sigidv(6), sigfdv(6)
    real(kind=8) :: afps, bfps, bfpd, cfpd, cfps, deisp, hfin, hini
    real(kind=8) :: sigfsp, sigisp, tdt, visp
!
! RECUPERATION DES VALEURS DES PARAMETRES MATERIAU
!
!
    ifou = nint(cmat(12))
    ifpo = nint(cmat(13))
    ides = nint(cmat(14))
    icou = nint(cmat(15))
!
! --- VISCOSITE SPHERIQUE IRREVERSIBLE
    visp = cmat(6)
!
! RECUPERATION DES VARIABLES INTERNES INITIALES
!
!   FLUAGE PROPRE
!
!     - SPHERIQUE  : ERSP EISP
!     - DEVIATOIRE : EPSDVR EPSDVI
!
!   FLUAGE DE DESSICCATION : EFDE
!
    ersp = vari(1)
    eisp = vari(2)
!
    epsdvr(1) = vari(3)
    epsdvi(1) = vari(4)
    epsdvr(2) = vari(5)
    epsdvi(2) = vari(6)
    epsdvr(3) = vari(7)
    epsdvi(3) = vari(8)
    epsdvr(4) = vari(12)
    epsdvi(4) = vari(13)
    epsdvr(5) = vari(14)
    epsdvi(5) = vari(15)
    epsdvr(6) = vari(16)
    epsdvi(6) = vari(17)
!
    efde(1) = vari(9)
    efde(2) = vari(10)
    efde(3) = vari(11)
    efde(4) = vari(18)
    efde(5) = vari(19)
    efde(6) = vari(20)
!
! INITIALISATION DES VARIABLES
!
! MOFI DU 6 JANVIER 2003 - YLP NSTRS --> 6
!      DO 11 I=1,NSTRS
!        AFPD(I)  = 0.D0
!        AFD(I)   = 0.D0
!        SIGFI(I) = 0.D0
!        DO 12 J=1,NSTRS
!          BFD(I,J)  = 0D0
!          CFD(I,J)  = 0D0
!12      CONTINUE
!11    CONTINUE
    do i = 1, 6
        afpd(i) = 0.d0
        afd(i) = 0.d0
        sigfi(i) = 0.d0
        do j = 1, 6
            bfd(i, j) = 0d0
            cfd(i, j) = 0d0
        end do
    end do
    do i = 1, 6
        sigidv(i) = 0.d0
        sigfdv(i) = 0.d0
    end do
!
! TEST DE COUPLAGE AVEC ICOU
!
    if (icou .eq. 0) then
        do i = 1, nstrs
            sigfi(i) = sigi(i)
        end do
    else
        do i = 1, nstrs
            sigfi(i) = sigf(i)
        end do
    end if
!
! CALCUL DES CONTRAINTES SPHERIQUES INI ET FIN
!
!  => EQUATION (3.11-1)
!
    sigisp = (sigi(1)+sigi(2))/3.d0
    sigfsp = (sigfi(1)+sigfi(2))/3.d0
!
    if (ifou .ne. -2) then
        sigisp = sigi(3)/3.d0+sigisp
        sigfsp = sigfi(3)/3.d0+sigfsp
    end if
!
! CALCUL DES CONTRAINTES DEVIATOIRES INI ET FIN
!
!  => EQUATION (3.11-2)
!
    sigidv(1) = sigi(1)-sigisp
    sigidv(2) = sigi(2)-sigisp
!
    sigfdv(1) = sigfi(1)-sigfsp
    sigfdv(2) = sigfi(2)-sigfsp
!
    if (ifou .ne. -2) then
        sigidv(3) = sigi(3)-sigisp
        sigidv(4) = sigi(4)
!
        sigfdv(3) = sigfi(3)-sigfsp
        sigfdv(4) = sigfi(4)
!
        if (ifou .eq. 2) then
            sigidv(5) = sigi(5)
            sigidv(6) = sigi(6)
!
            sigfdv(5) = sigfi(5)
            sigfdv(6) = sigfi(6)
        end if
    else
        sigidv(3) = sigi(3)
!
        sigfdv(3) = sigfi(3)
    end if
!_______________________________________________________________________
!
!  FLUAGE DE DESSICCATION : CALCUL DES MATRICES
!
! --> DEBRANCHE le 11 octobre 2002 - ylp.
! --> REBRANCHE le 25 aout 2004 - ylp.
!_______________________________________________________________________
!
    if ((ides .eq. 1) .or. (ides .eq. 2)) then
        call lcumfb(sigi, nstrs, vari, nvari, cmat, &
                    nmat, tdt, hini, hfin, afd, &
                    bfd, cfd)
    end if
!  CALCUL DES DEFORMATIONS FINALES FLUAGE DE DESSICCATION
!
!  => EQUATION (3.11-3)
!
    do i = 1, nstrs
        efde(i) = afd(i)+efde(i)+bfd(i, i)*sigi(i)+cfd(i, i)*sigfi(i)
    end do
!
!  FLUAGE PROPRE
!
    if (ifpo .ne. 0) then
!
!   FLUAGE PROPRE SPHERIQUE REVERSIBLE ET IRREVERSIBLE
!
!  => EQUATION (3.11-4)
!
        call lcumfs(vari, nvari, cmat, nmat, 1, &
                    isph, tdt, hini, hfin, afps, &
                    bfps, cfps)
        ersp = ersp+afps+bfps*sigisp+cfps*sigfsp
!
        call lcumfs(vari, nvari, cmat, nmat, 2, &
                    isph, tdt, hini, hfin, afps, &
                    bfps, cfps)
        deisp = afps+bfps*sigisp+cfps*sigfsp
!
!     TEST SUR LA CROISSANCE DE LA DEFORMATION DE FLUAGE
!          PROPRE SPHERIQUE
!
        if (isph .ne. 0) then
            if (deisp .gt. 0.d0) then
                isph = 2
                goto 60
            end if
            if ((sigfsp/visp .gt. -r8prem()) .and. (sigisp/visp .gt. -r8prem())) then
                isph = 2
                goto 60
            end if
        end if
!
        eisp = eisp+deisp
!
!   FLUAGE PROPRE DEVIATORIQUE REVERSIBLE ET IRREVERSIBLE
!
!  => EQUATION (3.11-5)
!
        call lcumfd(vari, nvari, nstrs, cmat, nmat, &
                    1, tdt, hini, hfin, afpd, &
                    bfpd, cfpd)
        do i = 1, nstrs
            epsdvr(i) = epsdvr(i)+afpd(i)+bfpd*sigidv(i)+cfpd*sigfdv(i)
        end do
        call lcumfd(vari, nvari, nstrs, cmat, nmat, &
                    2, tdt, hini, hfin, afpd, &
                    bfpd, cfpd)
        do i = 1, nstrs
            epsdvi(i) = epsdvi(i)+afpd(i)+bfpd*sigidv(i)+cfpd*sigfdv(i)
        end do
    end if
!
!  SAUVEGARDE DES DEFORMATIONS DE FLUAGE
!
!    PROPRE
!
!      - SPHERIQUE
!
    varf(1) = ersp
    varf(2) = eisp
!
!      - DEVIATORIQUE
!
    varf(3) = epsdvr(1)
    varf(4) = epsdvi(1)
    varf(5) = epsdvr(2)
    varf(6) = epsdvi(2)
    varf(7) = epsdvr(3)
    varf(8) = epsdvi(3)
    varf(12) = epsdvr(4)
    varf(13) = epsdvi(4)
    varf(14) = epsdvr(5)
    varf(15) = epsdvi(5)
    varf(16) = epsdvr(6)
    varf(17) = epsdvi(6)
!
!    DESSICCATION
!
    varf(9) = efde(1)
    varf(10) = efde(2)
    varf(11) = efde(3)
    varf(18) = efde(4)
    varf(19) = efde(5)
    varf(20) = efde(6)
!
60  continue
!
end subroutine
