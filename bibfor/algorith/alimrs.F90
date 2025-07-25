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

subroutine alimrs(mate, mateco, ma1, ma2, moint, ndble, &
                  nume_ddl, cn1, chno, cmp, icor)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calflu.h"
#include "asterfort/vtcreb.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nueq_chck.h"
!
    character(len=*) :: moint, cmp, chno, mate, mateco
!     AUTEUR : G. ROUSSEAU
!     BUT : IMMERGER UN CHAMP AUX NOEUDS DE DEPLACEMENT
!           D UNE STRUCTURE APPARTENANT A UN MAILLAGE
!           DANS UN MODELE D INTERFACE THERMIQUE APPARTENANT
!           A UN AUTRE MAILLAGE
!     IN: MATE: NOM DU MATERIAU FLUIDE
!     IN: MA1 : MAILLAGE DE LA STRUCTURE
!     IN: MA2 : MAILLAGE DU FLUIDE
!     IN: MOINT: NOM DU MODELE THERMIQUE D INTERFACE
!     IN: NDBLE: INDICATEUR DE NOEUDS DOUBLES
!     IN: NUM: NUMEDDL ASSOCIE A L INTERFACE
!     IN: CN1: CHAMP AUX NOEUDS DE DEPL_R MODAL DE LA STRUCTURE
!     IN: CHNO: NOM DU CHAMP THERMIQUE A CREER
!     IN: CMP: COMPOSANTE DE DEPL_R A PLONGER
!     IN: ICOR(1):  ADRESSE JEVEUX D UN TABLEAU DE CORRESPONDANCE
!               NOEUDS FLUIDES - NOEUDS STRUCTURES (NOEUDS SIMPLES)
!     IN: ICOR(2):  ADRESSE JEVEUX D UN TABLEAU DE CORRESPONDANCE
!               NOEUDS FLUIDES - NOEUDS STRUCTURES (NOEUDS DOUBLES)
!
!
!
!---------------------------------------------------------------------
    aster_logical :: test
    integer(kind=8) :: nbvale, nbrefe
    integer(kind=8) :: ino1, nocmp, icor(2), ichnul, ndble, i_ligr_mesh
    real(kind=8) :: tailmi
    character(len=8) :: gd1, gd2, ma1, ma2
    character(len=14) :: nume_ddl
    character(len=19) :: chnul, cn1, cn2, pchno1, pchno2
! -----------------------------------------------------------------
!---------------------------------------------------------------------
! TESTS PRELIMINAIRES : NUMERO DE COMPOSANTE A TRAITER
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iadg1, ieq1, ieq2, igeom1
    integer(kind=8) :: igeom2, iprn1, iprn2, ival1, ival2
    integer(kind=8) :: nbno1, ncmp2, nec1, nec2
    real(kind=8), pointer :: val1(:) => null()
    real(kind=8), pointer :: val2(:) => null()
    integer(kind=8), pointer :: nueq1(:) => null()
    integer(kind=8), pointer :: nueq2(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    if (cmp .eq. 'DX') nocmp = 1
    if (cmp .eq. 'DY') nocmp = 2
    if (cmp .eq. 'DZ') nocmp = 3
!
! RECUPERATION DE LA TAILLE DE REFERENCE
!
    call getvr8(' ', 'DIST_REFE', scal=tailmi)
!
!
! ON CREE UN CHAMNO BIDON SUR L INTERFACE THERMIQUE
!
    chnul = '&&ALIMRS.CHNUL'
    call vtcreb(chnul, 'V', 'R', &
                nume_ddlz=nume_ddl)
    call jeveuo(chnul//'.VALE', 'E', ichnul)
!
    cn2 = chno
    call calflu(chnul, moint, mate, mateco, nume_ddl, cn2, &
                nbrefe, nbvale, 'X')
!
!
!
! PARCOURS DU MAILLAGE STRUCTURE
!
    call dismoi('NOM_GD', cn2, 'CHAM_NO', repk=gd2)
    call dismoi('NOM_GD', cn1, 'CHAM_NO', repk=gd1)
    call dismoi('NB_NO_MAILLA', ma1, 'MAILLAGE', repi=nbno1)
    call dismoi('NUME_EQUA', cn2, 'CHAM_NO', repk=pchno2)
    call dismoi('NUME_EQUA', cn1, 'CHAM_NO', repk=pchno1)
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
!
    call dismoi('NB_EC', gd1, 'GRANDEUR', repi=nec1)
    call dismoi('NB_EC', gd2, 'GRANDEUR', repi=nec2)
!
    call jeveuo(pchno1//'.NUEQ', 'L', vi=nueq1)
    call jeveuo(pchno2//'.NUEQ', 'L', vi=nueq2)
!
    call jeveuo(cn1//'.VALE', 'L', vr=val1)
    call jeveuo(cn2//'.VALE', 'E', vr=val2)
!
!
    call jeveuo(ma2//'.COORDO    .VALE', 'L', igeom2)
    call jeveuo(ma1//'.COORDO    .VALE', 'L', igeom1)
!
! PARCOURS SUR LES NOEUDS DU MAILLAGE STRUCTURE
! ON REPERE LES NOEUDS COINCIDENTS GEOMETRIQUEMENT
! AVEC LES NOEUDS DE L INTERFACE FLUIDE
!
    do ino1 = 1, nbno1
!
        if (zi(icor(1)+ino1-1) .eq. 0) goto 10
!
        ncmp2 = zi(iprn2-1+(zi(icor(1)+ino1-1)-1)*(nec2+2)+2)
        if (ncmp2 .eq. 0) goto 20
!
        ival2 = zi(iprn2-1+(zi(icor(1)+ino1-1)-1)*(nec2+2)+1)
        ieq2 = nueq2(ival2-1+1)
!
        ival1 = zi(iprn1-1+(ino1-1)*(nec1+2)+1)
        iadg1 = iprn1-1+(ino1-1)*(nec1+2)+3
!
! EXISTENCE DE LA COMPOSANTE A PLONGER DANS LE CHAMNO
! DE DEPL_R DE LA STRUCTURE : SI ELLE EXISTE, ON LA
! RECOPIE DANS LE CHAMNO DE TEMP_R CREE SUR
! L INTERFACE
        test = exisdg(zi(iadg1), nocmp)
!
        if (test) then
!
            ieq1 = nueq1(ival1-1+nocmp)
!
! RECOPIE DE LA VALEUR DE DX, DY
! OU DZ DANS LA CMP DE TEMP DE CN2 DU NOEUD
! GEOMETRIQUEMENT COINCIDENT
!
            val2(ieq2) = val1(ieq1)
            goto 20
!
        end if
!
! PRECHERCHE DE NOEUDS DOUBLES
!
20      continue
        if (ndble .eq. 1) then
!
            if (zi(icor(2)+ino1-1) .eq. 0) goto 10
!
            ncmp2 = zi(iprn2-1+(zi(icor(2)+ino1-1)-1)*(nec2+2)+2)
            if (ncmp2 .eq. 0) goto 10
!
            ival2 = zi(iprn2-1+(zi(icor(2)+ino1-1)-1)*(nec2+2)+1)
            ieq2 = nueq2(ival2-1+1)
!
            ival1 = zi(iprn1-1+(ino1-1)*(nec1+2)+1)
            iadg1 = iprn1-1+(ino1-1)*(nec1+2)+3
!
! EXISTENCE DE LA COMPOSANTE A PLONGER DANS LE CHAMNO
! DE DEPL_R DE LA STRUCTURE : SI ELLE EXISTE, ON LA
! RECOPIE DANS LE CHAMNO DE TEMP_R CREE SUR
! L INTERFACE
            test = exisdg(zi(iadg1), nocmp)
!
            if (test) then
!
                ieq1 = nueq1(ival1-1+nocmp)
!
! RECOPIE DE LA VALEUR DE DX, DY
! OU DZ DANS LA CMP DE TEMP DE CN2 DU NOEUD
! GEOMETRIQUEMENT COINCIDENT
!
                val2(ieq2) = val1(ieq1)
                goto 10
!
            end if
        end if
!
!
10      continue
    end do
!
! --- MENAGE
!
    call detrsd('CHAM_NO', chnul)
!
    call jedema()
end subroutine
