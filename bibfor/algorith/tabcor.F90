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

subroutine tabcor(model, mate, mateco, ma1, ma2, moint, &
                  num, ndble, icor)
    implicit none
#include "jeveux.h"
#include "asterfort/calflu.h"
#include "asterfort/vtcreb.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvr8.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    character(len=*) :: moint, mate, mateco
!     AUTEUR : G. ROUSSEAU
!     BUT : CREER UNE TABLE DE CORRESPONDANCE ENTRE NOEUDS FLUIDES
!           D'INTERFACE ET NOEUDS DE STRUCTURES POUR DES MAILLAGES
!           STRUCTURES ET FLUIDE DISTINCTS
!     IN: MATE: NOM DU MATERIAU FLUIDE
!     IN: MA1 : MAILLAGE DE LA STRUCTURE
!     IN: MA2 : MAILLAGE DU FLUIDE
!     IN: MOINT: NOM DU MODELE THERMIQUE D INTERFACE
!     IN: NUM: NUMEDDL ASSOCIE A L INTERFACE
!     IN: CMP: COMPOSANTE DE DEPL_R A PLONGER
!     IN: NDBLE: INDICATEUR DE NOEUDS DOUBLES
!     OUT:ICOR(2):ADRESSES JEVEUX DES TABLEAUX DE CORRESPONDANCES
!              NOEUDS SIMPLES
!
!
!---------------------------------------------------------------------
    integer(kind=8) :: nbvale, nbrefe, ibid, nbno1
    integer(kind=8) :: ino1, ino2, icor(2), itb1, itb2, ncmp2, nbno2, ichnul
    integer(kind=8) :: nec2, iprn2
    integer(kind=8) :: ndble, nbptr
    real(kind=8) :: x1, y1, z1, x2, y2, z2
    real(kind=8) :: tailmi, dista2, tailm2
    character(len=2) :: model
    character(len=8) :: gd2, ma1, ma2
    character(len=14) :: num
    character(len=19) :: chnul, cn2, pchno2
    real(kind=8), pointer :: geom1(:) => null()
    real(kind=8), pointer :: geom2(:) => null()
    real(kind=8), parameter :: epsi = 1.d-2
! -----------------------------------------------------------------
!
! RECUPERATION DE LA TAILLE DE REFERENCE
!
    call getvr8(' ', 'DIST_REFE', scal=tailmi)
    call dismoi('NB_NO_MAILLA', ma2, 'MAILLAGE', repi=nbno2)
    tailm2 = (epsi*tailmi)**2
! ON CREE UN CHAMNO BIDON SUR L INTERFACE THERMIQUE
!
    chnul = '&&TABCOR.CHNUL'
    call vtcreb(chnul, 'V', 'R', &
                nume_ddlz=num)
    call jeveuo(chnul//'.VALE', 'E', ichnul)
!
!
    cn2 = '&&TABCOR.BIDON'
    call calflu(chnul, moint, mate, mateco, num, cn2, &
                nbrefe, nbvale, 'X')
!
!
!
! PARCOURS DU MAILLAGE THERMIQUE D INTERFACE
!
!
    call dismoi('NOM_GD', cn2, 'CHAM_NO', repk=gd2)
!
    call dismoi('NB_NO_MAILLA', ma1, 'MAILLAGE', repi=nbno1)
!
! CREATION D'UNE TABLE DE CORRESPONDANCES NOEUDS DE STRUCTURE
! AVEC NOEUDS DU FLUIDE DE L'INTERFACE (NOEUDS SIMPLES)
!
    call jecreo('&&TABCOR.CORRE1', 'V V I')
    call jeecra('&&TABCOR.CORRE1', 'LONMAX', nbno1)
    call jeecra('&&TABCOR.CORRE1', 'LONUTI', nbno1)
    call jeveut('&&TABCOR.CORRE1', 'E', itb1)
    icor(1) = itb1
!
! SI IL Y A DES NOEUDS DOUBLES
!
    if (ndble .eq. 1) then
        call jecreo('&&TABCOR.CORRE2', 'V V I')
        call jeecra('&&TABCOR.CORRE2', 'LONMAX', nbno1)
        call jeecra('&&TABCOR.CORRE2', 'LONUTI', nbno1)
        call jeveut('&&TABCOR.CORRE2', 'E', itb2)
        icor(2) = itb2
    end if
!
!
!
!      DO 2 I=1,NBRF
!
!2     CONTINUE
!
!
! RECUPERATION DES CARACTERISTIQUES DU CHAMP AUX NOEUDS DE
! L'INTERFACE
!
    call dismoi('NUME_EQUA', cn2, 'CHAM_NO', repk=pchno2)
    call jenonu(jexnom(pchno2//'.LILI', '&MAILLA'), ibid)
    call jeveuo(jexnum(pchno2//'.PRNO', ibid), 'L', iprn2)
    call dismoi('NB_EC', gd2, 'GRANDEUR', repi=nec2)
!
! RECUPERATION DES COORDONNEES DU MAILLAGE
!
    call jeveuo(ma2//'.COORDO    .VALE', 'L', vr=geom2)
    call jeveuo(ma1//'.COORDO    .VALE', 'L', vr=geom1)
!
! PARCOURS SUR LES NOEUDS DU MAILLAGE STRUCTURE
! ON REPERE LES NOEUDS COINCIDENTS GEOMETRIQUEMENT
! AVEC LES NOEUDS DE L INTERFACE FLUIDE
!
    do ino1 = 1, nbno1
!
        x1 = geom1((ino1-1)*3+1)
        y1 = geom1((ino1-1)*3+2)
        if (model .eq. '3D') then
            z1 = geom1((ino1-1)*3+3)
        end if
!
! PARCOURS SUR LES NOEUDS DU MAILLAGE FLUIDE
! ON REPERE LES NOEUDS COINCIDENTS GEOMETRIQUEMENT
! AVEC LES NOEUDS DE LA STRUCTURE ET ON RECOPIE
! LE NOEUD COINCIDENT DS LE TABLEAU DE CORRESPONDANCE
!
        nbptr = 0
        do ino2 = 1, nbno2
!
            ncmp2 = zi(iprn2-1+(ino2-1)*(nec2+2)+2)
            if (ncmp2 .eq. 0) goto 20
!
! CRITERE GEOMETRIQUE DE PROXIMITE
!
            x2 = geom2((ino2-1)*3+1)
            y2 = geom2((ino2-1)*3+2)
            if (model .eq. '3D') then
                z2 = geom2((ino2-1)*3+3)
                dista2 = (x1-x2)**2+(y1-y2)**2+(z1-z2)**2
            else
                dista2 = (x1-x2)**2+(y1-y2)**2
            end if
!
            if (dista2 .lt. (tailm2)) then
!
!
                nbptr = nbptr+1
                if (nbptr .eq. 1) then
                    zi(itb1+ino1-1) = ino2
                    if (ndble .eq. 0) goto 10
                else
                    zi(itb2+ino1-1) = ino2
                    goto 10
                end if
!
! PAS DE RECHERCHE DE NOEUDS DOUBLES
!
!                 GOTO 10
            end if
!
20          continue
        end do
!
10      continue
    end do
!
! --- MENAGE
    call detrsd('CHAM_NO', chnul)
    call detrsd('CHAM_NO', cn2)
!
end subroutine
