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

subroutine te0141(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!     CALCULE LA MATRICE DE MASSE ELEMENTAIRE DES ELEMENTS DE POUTRE
!     D'EULER ET DE TIMOSHENKO
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!       'MASS_MECA'      : CALCUL DE LA MATRICE DE MASSE COHERENTE
!       'MASS_MECA_DIAG' : CALCUL DE LA MATRICE DE MASSE CONCENTREE
!       'MASS_MECA_EXPLI': ......
!       'MASS_FLUI_STRU' : CALCUL DE LA MATRICE DE MASSE AJOUTEE
!       'M_GAMMA'        : CALCUL DU VECTEUR M_GAMMA
!
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!       'MECA_POU_D_E'  : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!       'MECA_POU_D_T'  : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!       'MECA_POU_D_EM' : POUTRE DROITE MULTIFIBRE D EULER (SECT. CONST)
!       'MECA_POU_D_TG' : POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!       'MECA_POU_D_TGM': POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!                         MULTI-FIBRES SECTION CONSTANTE
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmavec.h"
#include "asterfort/pmfmas.h"
#include "asterfort/pomass.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rhoequ.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/vecma.h"
!
    character(len=*) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: lmater, iret, nbpar, labsc
    integer(kind=8) :: lorien, iacce, ivect, lmat
    integer(kind=8) :: nno, nc, ntc, nbv, kanl, kpg, spt
    integer(kind=8) :: inbf, nbgf
    real(kind=8) :: valpar, xl
    real(kind=8) :: absmoy
    real(kind=8) :: e, g, xnu, rho, rhos, rhofi, rhofe, cm, phie, phii
    real(kind=8) :: pgl(3, 3), mlv(105)
    real(kind=8) :: matv(105), matp(14, 14)
    character(len=8) :: nompar, fami, poum
    character(len=16) :: ch16
    character(len=24) :: mator
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbres = 6
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
    data nomres/'E', 'NU', 'RHO', 'PROF_RHO_F_INT', 'PROF_RHO_F_EXT', 'COEF_MASS_AJOU'/
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara1 = 2
    real(kind=8) :: vale_cara1(nb_cara1)
    character(len=8) :: noms_cara1(nb_cara1)
    data noms_cara1/'R1', 'EP1'/
!
! --------------------------------------------------------------------------------------------------
!
!   caracteristiques des elements
    nno = 2
    nc = 6
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    mator = ' '
    if ((nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
        nno = 2
        nc = 7
    end if
    ntc = nc*nno
    nbv = ntc*(ntc+1)/2
!   recuperation des caracteristiques materiaux
    rho = 0.d0
    valres(:) = 0.0d0
    call jevech('PMATERC', 'L', lmater)
    call moytem('RIGI', 2, 1, '+', valpar, iret)
    nompar = 'TEMP'
    nbpar = 1
!
    if (option .eq. 'MASS_FLUI_STRU') then
        call poutre_modloc('CAGEP1', noms_cara1, nb_cara1, lvaleur=vale_cara1)
        call jevech('PABSCUR', 'L', labsc)
        absmoy = (zr(labsc-1+1)+zr(labsc-1+2))/2.0d0
        if (nomte .eq. 'MECA_POU_D_TGM') then
            call jevech('PNBSP_I', 'L', inbf)
            nbgf = zi(inbf+1)
            if (nbgf .ne. 1) call utmess('F', 'ELEMENTS3_3')
        end if
        call rcvalb(fami, kpg, spt, poum, zi(lmater), mator, 'ELAS_FLUI', 1, 'ABSC', [absmoy], &
                    nbres, nomres, valres, codres, 1)
        e = valres(1)
        xnu = valres(2)
        rhos = valres(3)
        rhofi = valres(4)
        rhofe = valres(5)
        cm = valres(6)
        phie = vale_cara1(1)*2.0d0
        g = e/(2.0d0*(1.0d0+xnu))
        if (phie .eq. 0.d0) then
            call utmess('F', 'ELEMENTS3_26')
        end if
        phii = (phie-2.0d0*vale_cara1(2))
        call rhoequ(rho, rhos, rhofi, rhofe, cm, phii, phie)
!
    else if ((option .eq. 'MASS_MECA') .or. (option .eq. 'MASS_MECA_DIAG') .or. &
             (option .eq. 'MASS_MECA_EXPLI') .or. (option .eq. 'M_GAMMA')) then
        if ((nomte .ne. 'MECA_POU_D_EM') .and. (nomte .ne. 'MECA_POU_D_TGM')) then
            call rcvalb(fami, kpg, spt, poum, zi(lmater), ' ', 'ELAS', nbpar, nompar, [valpar], &
                        3, nomres, valres, codres, 1)
            e = valres(1)
            xnu = valres(2)
            rho = valres(3)
            g = e/(2.0d0*(1.0d0+xnu))
        end if
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!   coordonnees des noeuds
    xl = lonele()
!   recuperation des orientations
    call jevech('PCAORIE', 'L', lorien)
!   calcul de la matrice de masse locale
    kanl = 1
    if (option .eq. 'MASS_MECA_DIAG' .or. option .eq. 'MASS_MECA_EXPLI') then
        kanl = 0
    end if
    if ((nomte .eq. 'MECA_POU_D_EM') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
        call pmfmas(nomte, option, rho, zi(lmater), kanl, mlv)
    else
        call pomass(nomte, e, xnu, rho, kanl, mlv)
    end if
!
    if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        call matrot(zr(lorien), pgl)
        call utpslg(nno, nc, pgl, mlv, matv)
        call vecma(matv, nbv, matp, ntc)
        call pmavec('ZERO', ntc, matp, zr(iacce), zr(ivect))
    else
        call jevech('PMATUUR', 'E', lmat)
        call matrot(zr(lorien), pgl)
        call utpslg(nno, nc, pgl, mlv, zr(lmat))
    end if
!
end subroutine
