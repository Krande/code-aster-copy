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

subroutine te0151(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL
!       - ENERGIE DE DEFORMATION
!       - ENERGIE CINETIQUE
!     POUR LES ELEMENTS DE POUTRE D'EULER ET DE TIMOSHENKO.
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'EPOT_ELEM' : ENERGIE DE DEFORMATION
!        'ECIN_ELEM' : ENERGIE CINETIQUE
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!     'MECA_POU_D_E'   : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!     'MECA_POU_D_T'   : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!     'MECA_POU_D_TG'  : POUTRE DROITE DE TIMOSHENKO(SECTION CONSTANTE)
!                        AVEC GAUCHISSEMENT
!     'MECA_POU_D_EM'  : POUTRE DROITE D'EULER MULTI-FIBRE
!                        (SECTION CONSTANTE)
!     'MECA_POU_D_TGM' :POUTRE DROITE DE TIMOSHENKO(SECTION CONSTANTE)
!                        MULTIFIBRE AVEC GAUCHISSEMENT
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmfrig.h"
#include "asterfort/pomass.h"
#include "asterfort/porigi.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/ptenci.h"
#include "asterfort/ptenpo.h"
#include "asterfort/ptenth.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/verifm.h"
!
    character(len=*) :: option, nomte
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: i, idis, iret, istruc, itype, jdepl, jende
    integer(kind=8) :: jfreq, jmasd, jvite, kanl, lmater, lorien
    integer(kind=8) :: nbpar, nbres, nc, nno, npg
    real(kind=8) :: tvar, e, enerth, g, rho, valpar, xl, xnu
! --------------------------------------------------------------------------------------------------
    character(len=3) :: stopz
    character(len=4) :: fami
    character(len=8) :: nompar, famil, poum
    character(len=16) :: ch16
    real(kind=8) :: ul(14), ug(14), pgl(3, 3), klc(14, 14), klv(105)
    real(kind=8) :: epsthe
    integer(kind=8) :: kpg, spt, nklv
! --------------------------------------------------------------------------------------------------
    parameter(nbres=3)
    integer(kind=8) :: codres(nbres)
    real(kind=8) :: valres(nbres)
    character(len=16) :: nomres(nbres)
    data nomres/'E', 'NU', 'RHO'/
! --------------------------------------------------------------------------------------------------
!
    call jevech('PMATERC', 'L', lmater)
    valres(:) = 0.0d0
!
    fami = 'RIGI'
    npg = 3
    istruc = 1
    nno = 2
    if (nomte .eq. 'MECA_POU_D_EM') npg = 2
!
    call moytem(fami, npg, 1, '+', valpar, iret)
    call verifm(fami, npg, 1, '+', zi(lmater), epsthe, iret)
    nbpar = 1
    nompar = 'TEMP'
    famil = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    call rcvalb(famil, kpg, spt, poum, zi(lmater), ' ', 'ELAS', nbpar, nompar, [valpar], &
                2, nomres, valres, codres, 1)
    call rcvalb(famil, kpg, spt, poum, zi(lmater), ' ', 'ELAS', nbpar, nompar, [valpar], &
                1, nomres(3), valres(3), codres(3), 1)
!
    e = valres(1)
    xnu = valres(2)
    rho = valres(3)
    g = e/(2.0d0*(1.0d0+xnu))
!
!   recuperation des caracteristiques generales des sections
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
    xl = lonele()
    call poutre_modloc('CAGNPO', ['TVAR'], 1, valeur=tvar)
    itype = nint(tvar)
!
    nc = 6
    if (nomte(1:13) .eq. 'MECA_POU_D_TG') nc = 7
!
    nklv = 2*nc*(2*nc+1)/2
    if (option .ne. 'ECIN_ELEM') then
        call jevech('PDEPLAR', 'L', jdepl)
        do i = 1, 2*nc
            ug(i) = zr(jdepl+i-1)
        end do
    else
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iret, iad=jvite)
!       iret ne peut valoir que 0 (tout est ok) ou 2 (champ non fourni)
        if (iret .eq. 0) then
            do i = 1, 2*nc
                ug(i) = zr(jvite+i-1)
            end do
        else
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=jdepl)
            if (iret .eq. 0) then
                do i = 1, 2*nc
                    ug(i) = zr(jdepl+i-1)
                end do
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
    end if
!
!   matrice de rotation PGL. vecteur deplacement ou vitesse local  ul = pgl * ug
    call utpvgl(nno, nc, pgl, ug, ul)

!
!   energie de deformation
    if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jende)
!       calcul de la matrice de rigidite locale
        if ((nomte .eq. 'MECA_POU_D_EM') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
            call pmfrig(nomte, zi(lmater), klv)
        else
            call porigi(nomte, e, xnu, xl, klv)
        end if
!       matrice rigidite ligne > matrice rigidite carre
        call vecma(klv, nklv, klc, 2*nc)
!        energie de deformation
        idis = 1
        call ptenpo(nc*2, ul, klc, zr(jende), itype, idis)
        if (epsthe .ne. 0.0d0) then
            call ptenth(ul, xl, epsthe, 2*nc, klc, enerth)
            zr(jende) = zr(jende)-enerth
        end if
!
    else if (option .eq. 'ECIN_ELEM') then
!       energie cinetique
        call jevech('PENERCR', 'E', jende)
        call jevech('PMASDIA', 'L', jmasd)
        call jevech('POMEGA2', 'L', jfreq)
        kanl = zi(jmasd)
!        calcul de la matrice de masse locale
        if (rho .ne. 0.0d0) then
!           KANL = 0 ou 1. masses concentrees ou coherentes
            call pomass(nomte, e, xnu, rho, kanl, klv)
!           matrice masse ligne > matrice masse carre
            call vecma(klv, 78, klc, 12)
!           energie cinetique
            idis = 1
            call ptenci(12, ul, klc, zr(jfreq), zr(jende), itype, kanl, idis)
        else if (codres(3) .ne. 0) then
            call utmess('F', 'ELEMENTS3_31')
        end if
!
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!
end subroutine
