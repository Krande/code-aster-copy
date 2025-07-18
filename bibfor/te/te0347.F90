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
subroutine te0347(option, nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!     ELEMENTS :
!        POU_D_TG
!        POU_D_T
!        POU_D_E
!
!     CALCUL DES OPTIONS :
!        SIEF_ELNO
!        FORC_NODA
!        REFE_FORC_NODA
!        VARI_ELNO
!
!     POUR LES CONTRAINTES ET LES FORC_NODA
!       RECOPIE DES POINTS 1 ET 2 SI NPG=2
!       RECOPIE DES POINTS 1 ET 3 SI NPG=3
!     QUI CONTIENNENT DEJA LES EFFORTS AUX NOEUDS
!
! IN  OPTION : OPTION DE CALCUL
! IN  NOMTE  : NOM DU TYPE ELEMENT
!
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/moytem.h"
#include "asterfort/porea2.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/utpvlg.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jtab(7), nno, nc, ichg, icompo, ichn, lgpg, nbvar, i, k, npg
    integer(kind=8) :: lorien, icgp, icontn, jvSief, ivectu, in, iret(2)
    integer(kind=8) :: igeom, kp, kk, imate
    integer(kind=8) :: istrxm, iretc
!
    real(kind=8) :: pgl(3, 3), fs(14), d1b3(2, 3), ksi1, forref, momref
    real(kind=8) :: sigp(7), d1b(7, 14), co(3), ey, ez, xl, temp
    real(kind=8) :: valres(2), e, nu, g, aa, xiy, xiz, alfay, alfaz
    real(kind=8) :: phiy, phiz, gamma
    character(len=2) :: nomres(2)
    character(len=4) :: fami
    character(len=8) :: peffor
!
    aster_logical :: lefgno, reactu, okelem
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nb_cara = 7
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1'/
!
! --------------------------------------------------------------------------------------------------
!
    okelem = (nomte .eq. 'MECA_POU_D_TG') .or. (nomte .eq. 'MECA_POU_D_T') .or. &
             (nomte .eq. 'MECA_POU_D_E')
    ASSERT(okelem)
!
    nno = 2
    fami = 'RIGI'
!   nombre de points de gauss
    call elrefe_info(fami=fami, npg=npg)
    ASSERT((npg .eq. 2) .or. (npg .eq. 3))
    if (nomte .eq. 'MECA_POU_D_TG') then
        nc = 7
    else
        nc = 6
    end if
    lefgno = (option .eq. 'SIEF_ELNO')
    if (lefgno) peffor = 'PSIEFNOR'
!
    if (option .eq. 'REFE_FORC_NODA  ') then
        call jevech('PVECTUR', 'E', ivectu)
        call terefe('EFFORT_REFE', 'MECA_POUTRE', forref)
        call terefe('MOMENT_REFE', 'MECA_POUTRE', momref)
        do in = 1, nno
            do i = 1, 3
                zr(ivectu+(in-1)*nc+i-1) = forref
            end do
            do i = 4, nc
                zr(ivectu+(in-1)*nc+i-1) = momref
            end do
        end do
!
    else if (option .eq. 'VARI_ELNO  ') then
        call jevech('PVARIGR', 'L', ichg)
        call jevech('PCOMPOR', 'L', icompo)
        call jevech('PVARINR', 'E', ichn)
!
        call tecach('OOO', 'PVARIGR', 'L', iret(1), nval=7, &
                    itab=jtab)
        lgpg = max(jtab(6), 1)*jtab(7)
        read (zk16(icompo+1), '(I16)') nbvar
!       pour les variables internes, on projette avec les fonctions de forme sur les
!       noeuds debut et fin de l'element
!       pour le point 1
        ksi1 = -sqrt(5.d0/3.d0)
        d1b3(1, 1) = ksi1*(ksi1-1.d0)/2.0d0
        d1b3(1, 2) = 1.d0-ksi1*ksi1
        d1b3(1, 3) = ksi1*(ksi1+1.d0)/2.0d0
!       pour le point 2
        ksi1 = sqrt(5.d0/3.d0)
        d1b3(2, 1) = ksi1*(ksi1-1.d0)/2.0d0
        d1b3(2, 2) = 1.d0-ksi1*ksi1
        d1b3(2, 3) = ksi1*(ksi1+1.d0)/2.0d0
        do i = 1, nbvar
            do k = 1, 3
                zr(ichn+i-1) = zr(ichn+i-1)+zr(ichg+lgpg*(k-1)+i-1)*d1b3(1, k)
                zr(ichn+lgpg+i-1) = zr(ichn+lgpg+i-1)+zr(ichg+lgpg*(k-1)+i-1)*d1b3(2, k)
            end do
        end do
!
!
! --- ------------------------------------------------------------------
    else if (lefgno .or. option .eq. 'FORC_NODA') then
!        recopie des valeurs au point gauss 1 et [2|3] qui contiennent deja les efforts aux noeuds
!           npg=2 : recopie des points 1 et 2
!           npg=3 : recopie des points 1 et 3
        if (lefgno) then
            call jevech('PCONTRR', 'L', icgp)
            call jevech(peffor, 'E', icontn)
            if (npg .eq. 2) then
                do i = 1, nc
                    zr(icontn-1+i) = zr(icgp-1+i)
                    zr(icontn-1+i+nc) = zr(icgp-1+i+nc)
                end do
            else
                do i = 1, nc
                    zr(icontn-1+i) = zr(icgp-1+i)
                    zr(icontn-1+i+nc) = zr(icgp-1+i+nc+nc)
                end do
            end if
        else if (option .eq. 'FORC_NODA') then
            call jevech('PSIEFR', 'L', jvSief)
            call jevech('PCAORIE', 'L', lorien)
            call jevech('PVECTUR', 'E', ivectu)
            call jevech('PGEOMER', 'L', igeom)
!
            call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=icompo)
            reactu = .false.
            if (iretc .eq. 0) reactu = (zk16(icompo+2) .eq. 'GROT_GDEP')
!
            if (nomte .eq. 'MECA_POU_D_TG') then
                call jevech('PMATERC', 'L', imate)
                xl = lonele()
                call r8inir(2*nc, 0.d0, fs, 1)
                co(1) = 5.d0/9.d0
                co(2) = 8.d0/9.d0
                co(3) = 5.d0/9.d0
!
!               THERMIQUE A T+
                call moytem(fami, npg, 1, '+', temp, &
                            iret(1))
                nomres(1) = 'E'
                nomres(2) = 'NU'
                call rcvalb(fami, 1, 1, '+', zi(imate), &
                            ' ', 'ELAS', 1, 'TEMP', [temp], &
                            2, nomres, valres, iret, 1)
                e = valres(1)
                nu = valres(2)
                g = e/(2.d0*(1.d0+nu))
!
                call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
                aa = vale_cara(1)
                xiy = vale_cara(2)
                xiz = vale_cara(3)
                alfay = vale_cara(4)
                alfaz = vale_cara(5)
                ey = vale_cara(6)
                ez = vale_cara(7)
!
                phiy = e*xiz*12.d0*alfay/(xl*xl*g*aa)
                phiz = e*xiy*12.d0*alfaz/(xl*xl*g*aa)
!
                do kp = 1, npg
                    call jsd1ff(kp, xl, phiy, phiz, d1b)
                    do i = 1, nc
                        sigp(i) = zr(jvSief-1+nc*(kp-1)+i)
                    end do
                    do k = 1, 2*nc
                        do kk = 1, nc
                            fs(k) = fs(k)+xl*sigp(kk)*d1b(kk, k)*co(kp)*0.50d0
                        end do
                    end do
                end do
!               prendre en compte centre de torsion
                fs(4) = fs(4)-ez*fs(2)+ey*fs(3)
                fs(11) = fs(11)-ez*fs(9)+ey*fs(10)
            else
                if (npg .eq. 2) then
                    do in = 1, nc
                        fs(in) = -zr(jvSief+in-1)
                        fs(in+nc) = zr(jvSief+in+nc-1)
                    end do
                else
                    do in = 1, nc
                        fs(in) = -zr(jvSief+in-1)
                        fs(in+nc) = zr(jvSief+in+nc+nc-1)
                    end do
                end if
            end if
!
            if (reactu) then
                call jevech('PSTRXMR', 'L', istrxm)
                gamma = zr(istrxm+3-1)
                call porea2(nno, nc, zr(igeom), gamma, pgl, &
                            xl, "PDEPLAR")
            else
                call matrot(zr(lorien), pgl)
            end if
            call utpvlg(nno, nc, pgl, fs, zr(ivectu))
        end if
!
    end if
end subroutine
