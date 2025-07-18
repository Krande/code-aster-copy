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

subroutine chtpcn(chno1, tgeom, tailmi, tmin, epsi, &
                  base, chno2)
    implicit none
!
#include "jeveux.h"
#include "asterfort/antece.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/nueq_chck.h"
!
    character(len=*) :: chno1, base, chno2
    real(kind=8) :: tgeom(6), tmin, epsi, tailmi, val
!
!----------------------------------------------------------------------
! AUTEUR: G.ROUSSEAU
!
! BUT:
! ----
! TRANSPORTER UN CHAMNO DE TEMP_R DEFINI SUR UNE PARTIE
! DU MAILLAGE D INTERFACE SUR UNE AUTRE PARTIE DE L INTERFACE
! CORRESPONDANT AUX CONTOURS IMMERGES D UNE
! SOUS-STRUCTURE NON MAILLEES PAR UNE TRANSFORMATION GEOMETRIQUE
!----------------------------------------------------------------------
!
! ARGUMENTS:
! ----------
! IN/JXIN  CHNO1: K19 : CHAM_NO DONT ON VA RECUPERER LES VALEURS
! IN       BASE   : K1  : NOM DE LA BASE SUR LAQUELLE LE CHAM_NO DOIT
!                         ETRE CREE
! IN       TAILMI  : R   : TAILLE DE MAILLE MIN
! IN       TGEOM : L_R8: TABLE DES COMPOSANTES DE LA TRANSFORMATION
!                        GEOMETRIQUE
!          3 COMPOSANTES DE TRANSL PUIS 3 ANGLES NAUTIQUES
!          DE ROTATION
! IN       TMIN : R8 : TEMP MINIMALE EN DECA DE LAQUELLE ON PEUT
!          AFFECTER
!          AU NOEUD UNE VALEUR DU CHAMNO A TRANSPORTER
!
! IN/JXOUT CHNO2: K19 : NOM DU CHAM_NO A CREER
!
!-----------------------------------------------------------------------
!
!
!
!
!
!
    integer(kind=8) :: nbante, ino1
    character(len=8) :: gd1, ma
    character(len=24) :: valk(2)
    character(len=8) :: diff, chnaff
    character(len=19) :: cn1, cn2, pchno1, pchno2
!
!-----------------------------------------------------------------------
    integer(kind=8) :: ieq1, ieq2, ino2, i_ligr_mesh
    integer(kind=8) :: iprn1, iprn2, ival1, ival2, nbcn1
    integer(kind=8) :: nbnaff, nbno, nbnrcp, ncmp1, ncmp2, nec
    real(kind=8), pointer :: val1(:) => null()
    real(kind=8), pointer :: val2(:) => null()
    integer(kind=8), pointer :: nueq1(:) => null()
    integer(kind=8), pointer :: nueq2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    cn1 = chno1
    cn2 = chno2
    call copisd('CHAMP_GD', base, cn1, cn2)
!
!
! ------------------------------ VERIFICATIONS -------------------------
!
    call dismoi('NOM_GD', cn1, 'CHAM_NO', repk=gd1)
    call dismoi('NUME_EQUA', cn1, 'CHAM_NO', repk=pchno1)
    call dismoi('NUME_EQUA', cn2, 'CHAM_NO', repk=pchno2)
    call jeveuo(cn1//'.VALE', 'L', vr=val1)
    call jelira(cn1//'.VALE', 'LONMAX', nbcn1)
!
!
    call jeveuo(cn2//'.VALE', 'E', vr=val2)
!

!
! - Protection: no matrix shrinking
!
    call nueq_chck(pchno1)
    call nueq_chck(pchno2)
!
    call jenonu(jexnom(pchno1//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(pchno1//'.PRNO', i_ligr_mesh), 'L', iprn1)
    call jenonu(jexnom(pchno2//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(pchno2//'.PRNO', i_ligr_mesh), 'L', iprn2)
    call jeveuo(pchno1//'.NUEQ', 'L', vi=nueq1)
    call jeveuo(pchno2//'.NUEQ', 'L', vi=nueq2)
!
    call dismoi('NOM_MAILLA', cn1, 'CHAM_NO', repk=ma)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
!
    call dismoi('NB_EC', gd1, 'GRANDEUR', repi=nec)
!
! NOMBRE DE NOEUDS A AFFECTER
!
    nbnaff = 0
    do ino1 = 1, nbno
        ncmp1 = zi(iprn1-1+(ino1-1)*(nec+2)+2)
        if (ncmp1 .eq. 0) goto 10
        ival1 = zi(iprn1-1+(ino1-1)*(nec+2)+1)
        ieq1 = nueq1(ival1-1+1)
        val = val1(ieq1)
        if (abs(val) .lt. tmin) nbnaff = nbnaff+1
10      continue
    end do
!
!
!
    nbnrcp = 0
    do ino2 = 1, nbno
!
        ival2 = zi(iprn2-1+(ino2-1)*(nec+2)+1)
        ncmp2 = zi(iprn2-1+(ino2-1)*(nec+2)+2)
        ieq2 = nueq2(ival2-1+1)
!
        if (ncmp2 .eq. 0) goto 1
!
        call antece(ino2, ma, tgeom, tailmi, epsi, &
                    nbante, ino1)
!
!
        if (nbante .gt. 1) then
!
            call utmess('F', 'CALCULEL2_7')
!
        else
!
            if (nbante .eq. 0) then
!
!
                val2(ieq2) = 0.0d0
!
!
            else
!
                if (nbante .eq. 1) then
!
!
                    ival1 = zi(iprn1-1+(ino1-1)*(nec+2)+1)
                    ncmp1 = zi(iprn1-1+(ino1-1)*(nec+2)+2)
                    ieq1 = nueq1(ival1-1+1)
                    val = val1(ieq1)
!
                    if (abs(val) .gt. tmin) then
                        val2(ieq2) = val
                        nbnrcp = nbnrcp+1
                    else
                        val2(ieq2) = 0.0d0
                    end if
!
!
                end if
!
            end if
        end if
!
!
1       continue
    end do
!
    if ((nbnrcp .lt. nbnaff) .and. (nbnrcp .gt. (nbnaff/2))) then
!
        call codent((nbnaff-nbnrcp), 'D0', diff(1:8))
        call codent((nbnaff), 'D0', chnaff(1:8))
!
        valk(1) = diff
        valk(2) = chnaff
        call utmess('A', 'CALCULEL2_8', nk=2, valk=valk)
!
    else
        if (nbnrcp .lt. (nbnaff/2)) then
!
            call utmess('A', 'CALCULEL2_9')
        end if
!
    end if
!
!
    if (nbnrcp .gt. nbnaff) then
!
        call utmess('F', 'CALCULEL2_10')
!
    end if
!
!
!
    call jedema()
end subroutine
