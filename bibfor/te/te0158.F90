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

subroutine te0158(option, nomte)
!
!
    implicit none
    character(len=*) :: option, nomte
!
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jpd1ff.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matela.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/pmfdge.h"
#include "asterfort/pmfpti.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/verifm.h"
#include "asterfort/pmfmats.h"
!
! --------------------------------------------------------------------------------------------------
!
!   OPTION :
!       DEGE_ELNO           : DEFORMATIONS GENERALISEES DE POUTRE
!   NOMTE  :
!       'MECA_POU_D_E'      : POUTRE DROITE D'EULER       (SECTION VARIABLE)
!       'MECA_POU_D_T'      : POUTRE DROITE DE TIMOSHENKO (SECTION VARIABLE)
!       'MECA_POU_D_TG'     : POUTRE DROITE DE TIMOSHENKO (GAUCHISSEMENT)
!       'MECA_POU_D_EM'     : POUTRE DROITE MULTIFIBRE D EULER (SECT. CONST)
!       'MECA_POU_D_TGM'    : POUTRE DROITE MULTIFIBRE DE TIMOSHENKO (GAUCHISSEMENT)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jeffg, lmater, iret, lorien, jdepl, i, j, kp, nc
    integer(kind=8) :: itemp, jtab(7), istrxr
!
    character(len=4) :: fami
    character(len=16) :: ch16
    integer(kind=8) :: ipos, in
    integer(kind=8) :: npg, nno, ivf, ipoids
    real(kind=8) :: b(4), gg, xi, wi
    real(kind=8) :: ul(14), pgl(3, 3), d1b(6, 12), dege(3, 7), d1btg(7, 14)
    real(kind=8) :: degem(6), alpha
    real(kind=8) :: temp, e, xnu, epsthe, g, xl
    real(kind=8) :: a, xiy, xiz, alfay, alfaz, phiy, phiz
    real(kind=8) :: ksi1, d1b3(2, 3), ey, ez
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara = 9
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara), nomat
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1', 'EY2', 'EZ2'/
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, nno=nno, npg=npg, jpoids=ipoids, jvf=ivf)
!
    if (option .eq. 'DEGE_ELNO') then
        call jevech('PDEFOGR', 'E', jeffg)
    else
        ch16 = option
        call utmess('F', 'ELEMENTS2_47', sk=ch16)
    end if
!
!   RECUPERATION DES COORDONNEES DES NOEUDS
    xl = lonele()
!
!   RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS
    call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)
    ey = (vale_cara(6)+vale_cara(8))/2.d0
    ez = (vale_cara(7)+vale_cara(9))/2.d0
!
!   RECUPERATION DES ORIENTATIONS
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
!
!   RECUPERATION DES DEPLACEMENTS
    call jevech('PDEPLAR', 'L', jdepl)
!
    if (nomte .ne. 'MECA_POU_D_EM') then
!       CARACTERISTIQUES MATERIAUX
        call jevech('PMATERC', 'L', lmater)
        call pmfmats(nomat)
!
        call verifm(fami, npg, 1, '+', zi(lmater), epsthe, iret)
        itemp = 0
        if (iret .eq. 0) itemp = 1
!
        call moytem(fami, npg, 1, '+', temp, iret)
        call matela(zi(lmater), nomat, itemp, temp, e, xnu)
!
        g = e/(2.0d0*(1.0d0+xnu))
!
        if (nomte .eq. 'MECA_POU_D_E') then
            nc = 6
            phiy = 0.0d0
            phiz = 0.0d0
        else if (nomte .eq. 'MECA_POU_D_T') then
            nc = 6
            phiy = e*xiz*12.d0*alfay/(xl*xl*g*a)
            phiz = e*xiy*12.d0*alfaz/(xl*xl*g*a)
        else if ((nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
            nc = 7
            phiy = e*xiz*12.d0*alfay/(xl*xl*g*a)
            phiz = e*xiy*12.d0*alfaz/(xl*xl*g*a)
        else
            ch16 = nomte
            call utmess('F', 'ELEMENTS2_42', sk=ch16)
        end if
!
!       PASSAGE DES DEPLACEMENTS DANS LE REPERE LOCAL
        call utpvgl(nno, nc, pgl, zr(jdepl), ul)
!
        if ((nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
!           PASSAGE DE G (CENTRE DE GRAVITE) A C (CENTRE DE TORSION)
            do i = 1, 2
                ul(7*(i-1)+2) = ul(7*(i-1)+2)-ez*ul(7*(i-1)+4)
                ul(7*(i-1)+3) = ul(7*(i-1)+3)+ey*ul(7*(i-1)+4)
            end do
        end if
!
!       BOUCLE SUR LES POINTS DE GAUSS
        if ((nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_TGM')) then
            do kp = 1, npg
                call jsd1ff(kp, xl, phiy, phiz, d1btg)
                do i = 1, nc
                    dege(kp, i) = 0.0d0
                    do j = 1, 2*nc
                        dege(kp, i) = dege(kp, i)+d1btg(i, j)*ul(j)
                    end do
                end do
                dege(kp, 1) = dege(kp, 1)-epsthe
            end do
        else
            do kp = 1, npg
                call jpd1ff(kp, xl, phiy, phiz, d1b)
                do i = 1, nc
                    dege(kp, i) = 0.0d0
                    do j = 1, 2*nc
                        dege(kp, i) = dege(kp, i)+d1b(i, j)*ul(j)
                    end do
                end do
                dege(kp, 1) = dege(kp, 1)-epsthe
            end do
        end if
!       POUR LE POINT 1
        ksi1 = -sqrt(5.d0/3.d0)
        d1b3(1, 1) = ksi1*(ksi1-1.d0)/2.0d0
        d1b3(1, 2) = 1.d0-ksi1*ksi1
        d1b3(1, 3) = ksi1*(ksi1+1.d0)/2.0d0
!       POUR LE POINT 2
        ksi1 = sqrt(5.d0/3.d0)
        d1b3(2, 1) = ksi1*(ksi1-1.d0)/2.0d0
        d1b3(2, 2) = 1.d0-ksi1*ksi1
        d1b3(2, 3) = ksi1*(ksi1+1.d0)/2.0d0
!
        do i = 1, nc
            do kp = 1, 3
                zr(jeffg+i-1) = zr(jeffg+i-1)+dege(kp, i)*d1b3(1, kp)
                zr(jeffg+nc+i-1) = zr(jeffg+nc+i-1)+dege(kp, i)*d1b3(2, kp)
            end do
        end do
    else
!
!       poutre multifibres MECA_POU_D_EM
        nc = 6
!       passage des deplacements dans le repere local
        call utpvgl(nno, nc, pgl, zr(jdepl), ul)
!       alpha modes incompatibles
        call tecach('ONO', 'PSTRXRR', 'L', iret, nval=7, itab=jtab)
        istrxr = jtab(1)
        alpha = zr(istrxr-1+15)
!
        do in = 1, 2
            call pmfpti(-in, zr(ipoids), zr(ivf), xl, xi, wi, b, gg)
            call pmfdge(b, gg, ul, alpha, degem)
            ipos = jeffg+nc*(in-1)
            do i = 1, nc
                zr(ipos+i-1) = degem(i)
            end do
        end do
    end if
!
end subroutine
