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
subroutine te0156(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/angvx.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/utpvlg.h"
    character(len=16) :: option, nomte
!-----------------------------------------------------------------------
! REALISE LES OPTIONS :
!     SIEF_ELNO
!                POUR  LES CONTRAINTES DE L'ELEMENT MECA_BARRE
!     FORC_NODA      : FORCES NODALE DE L'ELEMENT MECA_BARRE
!
! ----------------------------------------------------------------------
! IN OPTION    : K16 :  OPTION DE CALCUL
!                       'FORC_NODA' OU  'SIEF_ELNO'
!                       OU 'REFE_FORC_NODA'
! IN NOMTE     : K16 : NOM DU TYPE ELEMENT
!                      'MECA_BARRE'
!                      'MECA_2D_BARRE'
!
!
!
    integer(kind=8) :: ivectu, icontg, lorien, nno, nc, ino, i
    integer(kind=8) :: jvCompor, jvDisp, jvGeom, iretc
    real(kind=8) :: fs(6), pgl(3, 3), vect(6), forref
    real(kind=8) :: w(6), ang1(3), xd(3)
    aster_logical :: reactu
!
!     ------------------------------------------------------------------
!
    if (option .eq. 'REFE_FORC_NODA') then
        nno = 2
        if (nomte .eq. 'MECA_2D_BARRE') then
            nc = 2
        else if (nomte .eq. 'MECA_BARRE') then
            nc = 3
        end if
        call jevech('PVECTUR', 'E', ivectu)
        call terefe('EFFORT_REFE', 'MECA_BARRE', forref)
        do ino = 1, nno
            do i = 1, nc
                zr(ivectu+(ino-1)*nc+i-1) = forref
            end do
        end do
!
    else if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', icontg)
        call jevech('PCAORIE', 'L', lorien)
        call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=jvCompor)
        reactu = .false.
        if (iretc .eq. 0) reactu = (zk16(jvCompor+2) .eq. 'PETIT_REAC')
!
!        PARAMETRES EN SORTIE
        call jevech('PVECTUR', 'E', ivectu)
        nno = 2
        nc = 3
        do i = 1, nno*nc
            fs(i) = 0.d0
        end do
        fs(1) = -zr(icontg)
        fs(4) = zr(icontg)
!
!
        if (reactu) then
            call jevech('PGEOMER', 'L', jvGeom)
            call jevech('PDEPLAR', 'L', jvDisp)

            if (nomte .eq. 'MECA_BARRE') then
                do i = 1, 3
                    w(i) = zr(jvGeom-1+i)+zr(jvDisp-1+i)
                    w(i+3) = zr(jvGeom+2+i)+zr(jvDisp+2+i)
                    xd(i) = w(i+3)-w(i)
                end do
            else if (nomte .eq. 'MECA_2D_BARRE') then
                w(1) = zr(jvGeom-1+1)+zr(jvDisp-1+1)
                w(2) = zr(jvGeom-1+2)+zr(jvDisp-1+2)
                w(3) = 0.d0
                w(4) = zr(jvGeom-1+3)+zr(jvDisp-1+3)
                w(5) = zr(jvGeom-1+4)+zr(jvDisp-1+4)
                w(6) = 0.d0
                xd(1) = w(4)-w(1)
                xd(2) = w(5)-w(2)
                xd(3) = 0.d0
            end if
            call angvx(xd, ang1(1), ang1(2))
            ang1(3) = zr(lorien+2)
            call matrot(ang1, pgl)
        else
            call matrot(zr(lorien), pgl)
        end if
        call utpvlg(nno, nc, pgl, fs, vect)
!
!        ECRITURE DANS LE VECTEUR VECTU SUIVANT L'ELEMENT
!
        if (nomte .eq. 'MECA_BARRE') then
            do i = 1, 6
                zr(ivectu+i-1) = vect(i)
            end do
        else if (nomte .eq. 'MECA_2D_BARRE') then
            zr(ivectu) = vect(1)
            zr(ivectu+1) = vect(2)
            zr(ivectu+2) = vect(4)
            zr(ivectu+3) = vect(5)
        end if
    end if
end subroutine
