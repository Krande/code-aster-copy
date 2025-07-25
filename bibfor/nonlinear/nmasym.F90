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

subroutine nmasym(fami, kpg, ksp, icodma, option, &
                  xlong0, a, tmoins, tplus, dlong0, &
                  effnom, vim, effnop, vip, klv, &
                  fono)
    implicit none
#include "asterfort/nm1das.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
    integer(kind=8) :: kpg, ksp, neq, nbt, nvar, icodma
    parameter(neq=6, nbt=21, nvar=4)
!
    character(len=*) :: fami, option
    real(kind=8) :: xlong0, a, syc, syt, etc, ett, cr
    real(kind=8) :: e, dlong0, tmoins, tplus
    real(kind=8) :: effnom, vim(nvar)
    real(kind=8) :: effnop, vip(nvar), fono(neq), klv(nbt)
! -------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ASYMETRIQUE LINEAIRE - VON MISES-
!    POUR UN MODELE BARRE ELEMENT MECA_BARRE
!
! -------------------------------------------------------------------
! IN  :
!       XLONG0 : LONGUEUR DE L'ELEMENT DE BARRE AU REPOS
!       A      : SECTION DE LA BARRE
!       TMOINS : INSTANT PRECEDENT
!       TPLUS  : INSTANT COURANT
!       XLONGM : LONGEUR DE L'ELEMENT AU TEMPS MOINS
!       DLONG0 : INCREMENT D'ALLONGEMENT DE L'ELEMENT
!       EFFNOM : EFFORT NORMAL PRECEDENT
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
!
! OUT : EFFNOP : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       FONO   : FORCES NODALES COURANTES
!       KLV    : MATRICE TANGENTE
!
!----------VARIABLES LOCALES
!
    real(kind=8) :: sigm, deps, dsdem, dsdep, sigp, xrig
    integer(kind=8) :: nbpar, nbres, kpg1, spt
!
    real(kind=8) :: valpar, valres(4)
    integer(kind=8) :: icodre(4)
    character(len=8) :: nompar, famil, poum
    character(len=16) :: nomela, nomasl(4)
    data nomela/'E'/
    data nomasl/'SY_C', 'DC_SIGM_EPSI', 'SY_T', 'DT_SIGM_EPSI'/
!
!----------INITIALISATIONS
!
    call r8inir(nbt, 0.d0, klv, 1)
    call r8inir(neq, 0.d0, fono, 1)
!
!----------RECUPERATION DES CARACTERISTIQUES
!
    deps = dlong0/xlong0
    sigm = effnom/a
!
! --- CARACTERISTIQUES ELASTIQUES
!
    nbres = 2
    nbpar = 0
    nompar = '  '
    valpar = 0.d0
    famil = 'FPG1'
    kpg1 = 1
    spt = 1
    poum = '+'
    call rcvalb(famil, kpg1, spt, poum, icodma, &
                ' ', 'ELAS', 0, nompar, [valpar], &
                1, nomela, valres, icodre, 1)
    e = valres(1)
!
! --- CARACTERISTIQUES ECROUISSAGE LINEAIRE ASYMETRIQUE
!
!
!JMP  NBRES = 5
    nbres = 4
    nbpar = 0
    call rcvalb(fami, 1, 1, '+', icodma, &
                ' ', 'ECRO_ASYM_LINE', nbpar, nompar, [valpar], &
                nbres, nomasl, valres, icodre, 1)
    syc = valres(1)
    etc = valres(2)
    syt = valres(3)
    ett = valres(4)
!JMP    CR     = VALRES(5) MODELE DE RESTAURATION PAS AU POINT
!
    cr = 0.d0
!
!
    call nm1das(fami, kpg, ksp, e, syc, &
                syt, etc, ett, cr, tmoins, &
                tplus, icodma, sigm, deps, vim, &
                sigp, vip, dsdem, dsdep)
    effnop = sigp*a
!
! --- CALCUL DU COEFFICIENT NON NUL DE LA MATRICE TANGENTE
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (option(11:14) .eq. 'ELAS') then
            xrig = e*a/xlong0
        else
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                xrig = dsdem*a/xlong0
            else
                xrig = dsdep*a/xlong0
            end if
        end if
        klv(1) = xrig
        klv(7) = -xrig
        klv(10) = xrig
    end if
!
! --- CALCUL DES FORCES NODALES
!
    fono(1) = -effnop
    fono(4) = effnop
!
! -------------------------------------------------------------
!
end subroutine
