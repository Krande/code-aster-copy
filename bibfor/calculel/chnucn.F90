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

subroutine chnucn(chno1, numdd2, ncorr, tcorr, base, &
                  chno2)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/nueq_chck.h"
!
    character(len=*) :: chno1, numdd2, base, chno2, tcorr(*)
    integer(kind=8) :: ncorr
!
!-----------------------------------------------------------------------
! BUT:
! ----
! CREER UN CHAM_NO S'APPUYANT SUR UN NUME_DDL
! ET CONTENANT LES VALEURS D'UN AUTRE CHAM_NO
!-----------------------------------------------------------------------
! ARGUMENTS:
! ----------
! IN/JXIN  CHNO1: K19 : CHAM_NO DONT ON VA RECUPERER LES VALEURS
! IN/JXIN  NUMDD2 : K14 : NUME_MECA DU CHAM_NO A CREER
! IN       BASE   : K1  : NOM DE LA BASE SUR LAQUELLE LE CHAM_NO DOIT
!                         ETRE CREE
! IN       NCORR  : I   : DIMENSION DE TCORR
! IN       TCORR  : L_K8: TABLE DE CORRESPONDANCE DES COMPOSANTES
!
! IN/JXOUT CHNO2: K19 : NOM DU CHAM_NO A CREER
!
!-----------------------------------------------------------------------
! USAGE:
! ------
! CETTE ROUTINE RECOPIE LES VALEURS DU CHAM_NO (CHNO1) DANS UN NOUVEAU
! CHAM_NO (CHNO2) QUI S'APPUIE SUR LA NUMEROTATION (NUMDD2).CE CHAM_NO
! EST CREE SUR LA BASE (BASE).
!
!   CHNO2 DOIT ETRE DIFFERENT DE CHNO1.
!   CHNO2 EST ECRASE S'IL EXISTE DEJA.
!
! ON NE TRAITE QUE LES NOEUDS DU MAILLAGE (PAS LES NOEUDS DE LAGRANGE)
! ON NE TRAITE POUR L'INSTANT QUE LES CHAM_NO DE TYPE R8
!
! SI UNE COMPOSANTE DE CHNO2 N'EST PAS AFFECTEE DANS CHNO1, ON LA
! MET A ZERO.
!
! LA GRANDEUR ASSOCIEE A CHNO2 PEUT ETRE DIFFERENTE DE CELLE DE
! CHNO1. DANS CE CAS, ON UTILISE LA TABLE DE CORRESPONDANCE DES
! COMPOSANTES (TCORR)
! (ON SE SERVIRA SYSTEMATIQUEMENT DE TCORR LORSQUE NCORR /= 0)
!
!   NCORR EST UN ENTIER PAIR.
!   SI NCORR = 0 , LA GRANDEUR ASSOCIEE A CHNO1 DOIT ETRE IDENTIQUE A
!                  CELLE DE NUMDD2.
!   LA CORRESPONDANCE : TCORR(2*(I-1)+1) --> TCORR(2*(I-1)+2) DOIT ETRE
!   INJECTIVE.
!
!
! EXEMPLE 1 :
! -----------
! ON DISPOSE D'UN CHAM_NO_DEPL_R SUR L'ENSEMBLE DU MODELE GLOBAL : CHG
! ON DISPOSE D'UN SOUS-MODELE ET D'UN NUME_DDL ASSOCIE : NUL
! ON PEUT ALORS CREER LE CHAM_NO_DEPL_R "PROJETE" SUR LE SOUS-MODELE:CHL
!
!
!
! INVERSEMENT, ON POURRAIT "PROLONGER" (PAR DES ZEROS) UN CHAMP LOCAL :
!
!
! EXEMPLE 2 :
! -----------
! ON DISPOSE D'UN CHAM_NO_DEPL_R ASSOCIE A 1 MODELE MECANIQUE (CHDEPL)
! ON DISPOSE D'UN NUME_DDL ASSOCIE A 1 MODELE THERMIQUE (NUTH)
! ON PEUT ALORS CREER LE CHAM_NO_TEMP_R (CHTEMP)
! QUI CONTIENDRA COMME COMPOSANTE: 'TEMP' LES VALEURS DE CHDEPL POUR
! LA COMPOSANTE 'DY'
!
!       TCORR(1)='DY'
!       TCORR(2)='TEMP'
!
!
!
! EXEMPLE 3 :
! -----------
! ON DISPOSE D'UN CHAM_NO_DEPL_R (CHDEPL)
! ON DISPOSE DU NUME_DDL ASSOCIE A CE CHAM_NO : NUDEPL
! ON PEUT ALORS CREER LE CHAM_NO_DEPL_R (CHDEPL2)
! AVEC LES CORRESPONDANCES SUIVANTES :
!
!  CHNO2_DX = 0.
!  CHNO2_DY = CHNO1_DX
!  CHNO2_DZ = CHNO1_DY
!
!
!       TCORR(1)='DX'
!       TCORR(2)='DY'
!       TCORR(3)='DY'
!       TCORR(4)='DZ'
!
!
!
!  LA COMPOSANTE 'DX' DE CHNO2 N'ETANT PAS DANS LA TABLE DE
!  CORRESPONDANCE, ON LUI AFFECTERA LA VALEUR : 0.
!
!
!
!
!
    character(len=1) :: base2
    character(len=8) :: gd1, gd2, repk, tysca1, tysca2, ma, cmp1, cmp2
    character(len=14) :: nu2
    character(len=19) :: cn1, cn2, pchno1, pchno2
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, iacmp1, iacmp2, iadg1, iadg2, i_ligr_mesh
    integer(kind=8) :: iaval2, ico1, ico2, ieq1
    integer(kind=8) :: ieq2, ino, iprn1, iprn2
    integer(kind=8) :: iret, ival1, ival2, j1, j2, nbno, ncmmx1
    integer(kind=8) :: ncmmx2, ncmp1, ncmp2, nec1, nec2, nval1
    integer(kind=8) :: nval2
    integer(kind=8), pointer :: corr2(:) => null()
    real(kind=8), pointer :: vale(:) => null()
    integer(kind=8), pointer :: nueq1(:) => null()
    integer(kind=8), pointer :: nueq2(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    base2 = base
    cn1 = chno1
    cn2 = chno2
    nu2 = numdd2
!
! ------------------------------ VERIFICATIONS -------------------------
!
    call dismoi('NOM_GD', cn1, 'CHAM_NO', repk=gd1)
    call dismoi('NUME_EQUA', cn1, 'CHAM_NO', repk=pchno1)
    pchno2 = nu2//'.NUME'
    call dismoi('NOM_GD', nu2, 'NUME_DDL', repk=gd2)
!
    call dismoi('TYPE_SCA', gd1, 'GRANDEUR', repk=tysca1)
    call dismoi('TYPE_SCA', gd2, 'GRANDEUR', repk=tysca2)
    if (tysca1 .ne. 'R') then
        call utmess('F', 'CALCULEL_92', sk=cn1)
    end if
    if (tysca2 .ne. 'R') then
        call utmess('F', 'CALCULEL_93', sk=nu2)
    end if
!
!
    call dismoi('NB_EQUA', cn1, 'CHAM_NO', repi=nval1)
!
! ------------------------------- REFE --------------------------------
!
!     -- SI CN2 EXISTE DEJA, ON LE DETRUIT :
    call jeexin(cn2//'.REFE', iret)
    if (iret .gt. 0) call detrsd('CHAMP_GD', cn2)
!
    call wkvect(cn2//'.REFE', base2//' V K24', 4, i1)
    call jeecra(cn2//'.REFE', 'DOCU', cval='CHNO')
    call jeveuo(nu2//'.NUME.REFN', 'L', i2)
    zk24(i1+1) = nu2//'.NUME'
!
! ------------------------------- VALE --------------------------------
!
    call dismoi('NB_EQUA', nu2, 'NUME_DDL', repi=nval2)
    call wkvect(cn2//'.VALE', base2//' V R', nval2, iaval2)
    call jeveuo(cn1//'.VALE', 'L', vr=vale)
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
    call dismoi('NOM_MAILLA', nu2, 'NUME_DDL', repk=repk)
    ASSERT(ma .eq. repk)
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nbno)
!
    call dismoi('NB_EC', gd1, 'GRANDEUR', repi=nec1)
    call dismoi('NB_EC', gd2, 'GRANDEUR', repi=nec2)
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', gd1), 'L', iacmp1)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', gd2), 'L', iacmp2)
    call jelira(jexnom('&CATA.GD.NOMCMP', gd1), 'LONMAX', ncmmx1)
    call jelira(jexnom('&CATA.GD.NOMCMP', gd2), 'LONMAX', ncmmx2)
!
!     -- REMPLISSAGE DE L'OBJET '.CORR2' :
!     ------------------------------------
    AS_ALLOCATE(vi=corr2, size=ncmmx2)
    if (ncorr .eq. 0) then
!       LES GRANDEURS G1 ET G2 DOIVENT ETRE IDENTIQUES
        ASSERT(gd1 .eq. gd2)
        do i2 = 1, ncmmx2
            corr2(i2) = i2
        end do
    else
        ASSERT(ncorr .eq. 2*(ncorr/2))
        do i = 1, ncorr/2
            cmp1 = tcorr(2*(i-1)+1)
            cmp2 = tcorr(2*(i-1)+2)
            j1 = indik8(zk8(iacmp1), cmp1, 1, ncmmx1)
            j2 = indik8(zk8(iacmp2), cmp2, 1, ncmmx2)
            if (j2 .ne. 0) corr2(j2) = j1
        end do
    end if
!
    do ino = 1, nbno
        ival1 = zi(iprn1-1+(ino-1)*(nec1+2)+1)
        ival2 = zi(iprn2-1+(ino-1)*(nec2+2)+1)
        ncmp1 = zi(iprn1-1+(ino-1)*(nec1+2)+2)
        ncmp2 = zi(iprn2-1+(ino-1)*(nec2+2)+2)
        iadg1 = iprn1-1+(ino-1)*(nec1+2)+3
        iadg2 = iprn2-1+(ino-1)*(nec2+2)+3
        if (ncmp1*ncmp2 .eq. 0) goto 1
        ico2 = 0
        do i2 = 1, ncmmx2
            if (exisdg(zi(iadg2), i2)) then
                ico2 = ico2+1
                i1 = corr2(i2)
!
                if (.not. (exisdg(zi(iadg1), i1))) then
                    ico1 = 0
                else
                    ico1 = 0
                    do j1 = 1, i1
                        if (exisdg(zi(iadg1), j1)) ico1 = ico1+1
                    end do
                end if
!
                if (ico1 .gt. 0) then
!             --RECOPIE D'UNE VALEUR :
                    ieq1 = nueq1(ival1-1+ico1)
                    ieq2 = nueq2(ival2-1+ico2)
                    zr(iaval2-1+ieq2) = vale(ieq1)
                end if
!
            end if
        end do
1       continue
    end do
!
    AS_DEALLOCATE(vi=corr2)
!
    call jedema()
end subroutine
