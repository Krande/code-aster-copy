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

subroutine op0113()
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
!
! ----------------------------------------------------------------------
!
! OPERATEUR MODI_MODELE_XFEM
!
!
! ----------------------------------------------------------------------
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/initel.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedup1.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xcodec.h"
#include "asterfort/xcpmod.h"
#include "asterfort/xmolig.h"
#include "asterfort/xtyele.h"
#include "asterfort/xverm2.h"
#include "asterfort/xvermo.h"
!
    integer(kind=8) :: iret, iel, ima, nmoth
    integer(kind=8) :: i, j2
    integer(kind=8) :: jmofis
    integer(kind=8) :: nbma, nelt
    integer(kind=8) :: nb1
    integer(kind=8) :: nfiss, jnfis
    integer(kind=8) :: ndim
    character(len=16) :: motfac, k16bid, line_quad, face
    character(len=19) :: ligr1, ligr2
    character(len=24) :: liel1, liel2
    character(len=24) :: mail2
    character(len=24) :: trav
    integer(kind=8) :: jmail2, jtab, jxc
    character(len=8) :: modelx, mod1, modthx, noma, k8cont, k8condi, decou
    aster_logical :: linter
    character(len=8), pointer :: lgrf1(:) => null()
    character(len=8), pointer :: lgrf2(:) => null()
!
    data motfac/' '/
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! - This functionnality is NOT qualified for nuclear safety studies
!
    call utmess('A', 'QUALITY1_1')
!
! --- NOM DU MODELE MODIFIE
!
    call getres(modelx, k16bid, k16bid)
    ligr2 = modelx//'.MODELE'
    liel2 = ligr2//'.LIEL'
!
! --- NOM DU MODELE INITIAL
!
    call getvid(motfac, 'MODELE_IN', iocc=1, scal=mod1, nbret=iret)
    ligr1 = mod1//'.MODELE'
    liel1 = ligr1//'.LIEL'
!
! --- ACCES AU MAILLAGE INITIAL
!
    call jeveuo(ligr1//'.LGRF', 'L', vk8=lgrf1)
    noma = lgrf1(1)
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
!
! --- MOT-CLE MODELE_THER : CAS PARTICULIER (X-FEM THERMO MECA)
!     QUE L'ON TRAITE SEPAREMENT DU CAS GENERAL
!
    call getvid(motfac, 'MODELE_THER', iocc=1, scal=modthx, nbret=nmoth)
    if (nmoth .eq. 1) then
        call xcpmod(mod1, modthx, modelx)
        goto 999
    end if
!
! --- RECUPERER LE NOMBRE DE FISSURES
!
    call getvid(motfac, 'FISSURE', iocc=1, nbval=0, nbret=nfiss)
    nfiss = -nfiss
!
! --- CREATION DES OBJETS POUR MULTIFISSURATION DANS MODELE MODIFIE
!
    call wkvect(modelx//'.NFIS', 'G V I', 1, jnfis)
    call wkvect(modelx//'.FISS', 'G V K8', nfiss, jmofis)
    zi(jnfis) = nfiss
!
! --- RECUPERER LES FISSURES ET REMPLISSAGE DE MODELX//'.FISS'
!
    call getvid(motfac, 'FISSURE', iocc=1, nbval=nfiss, vect=zk8(jmofis), &
                nbret=iret)
!
!     VERIFICATION DE LA COHERENCE DES MOT-CLES FISSURE ET MODELE_IN
!     (COHERENCE DES MAILLAGES SOUS-JACENTS AUX FISSURES ET MODELE)
    call xvermo(nfiss, zk8(jmofis), noma)
!
!     VERIFS POUR LES MODELISATIONS "EXOTIQUES" (multi-h, thermique, HM)
    call xverm2(nfiss, zk8(jmofis), mod1)
!
!
! --- CONTACT ?
!
    call getvtx(motfac, 'CONTACT', iocc=1, scal=k8cont, nbret=iret)
    call wkvect(modelx//'.XFEM_CONT', 'G V I', 1, jxc)
    if (k8cont .eq. 'SANS') then
        zi(jxc) = 0
    else if (k8cont .eq. 'MORTAR') then
        zi(jxc) = 2
    else if (k8cont .eq. 'STANDARD') then
        call dismoi('LINE_QUAD', ligr1, 'LIGREL', repk=line_quad)
        if (line_quad .eq. 'LINE') then
            zi(jxc) = 1
        else if (line_quad .eq. 'QUAD') then
            zi(jxc) = 3
        else
            call utmess('F', 'XFEM2_3')
        end if
    else
        ASSERT(.false.)
    end if
!
    call getvtx(motfac, 'DECOUPE_FACETTE', iocc=1, scal=face, nbret=iret)
    if (iret .eq. 0) then
        decou = ' '
    else
        decou = face(1:8)
    end if
!
! --- CREATION DU TABLEAU DE TRAVAIL
!
    trav = '&&OP0113.TAB'
    call wkvect(trav, 'V V I', nbma*5, jtab)
!
    do i = 1, nbma
        zi(jtab-1+5*(i-1)+4) = 1
    end do
!
! ---------------------------------------------------------------------
!     1)  REMPLISSAGE DE TAB : NBMA X 5 : GR1 | GR2 | GR3 | GR0 | ITYP
! ---------------------------------------------------------------------
!
    call xtyele(mod1, trav, nfiss, zk8(jmofis), zi(jxc), &
                ndim, linter)
!
! ---------------------------------------------------------------------
!       2)  MODIFICATION DE TAB EN FONCTION DE L'ENRICHISSEMENT
! ---------------------------------------------------------------------
!
    call xmolig(liel1, trav)
!
! --- ON COMPTE LE NB DE MAILLES DU LIGREL1 (= NB DE GREL DE LIEL2)
!
    nelt = 0
    do ima = 1, nbma
        if (zi(jtab-1+5*(ima-1)+5) .ne. 0) then
            nelt = nelt+1
        end if
    end do
    if (nelt .eq. 0) then
        call utmess('F', 'XFEM2_51')
    end if
!
!-----------------------------------------------------------------------
!     3)  CONSTRUCTION DU .LIEL2
!-----------------------------------------------------------------------
!
    call jecrec(liel2, 'G V I', 'NU', 'CONTIG', 'VARIABLE', &
                nelt)
    call jeecra(liel2, 'LONT', 2*nelt)
!
    iel = 0
    do ima = 1, nbma
        if (zi(jtab-1+5*(ima-1)+5) .eq. 0) goto 300
        iel = iel+1
        call jecroc(jexnum(liel2, iel))
        call jeecra(jexnum(liel2, iel), 'LONMAX', 2)
        call jeveuo(jexnum(liel2, iel), 'E', j2)
        zi(j2-1+1) = ima
        zi(j2-1+2) = zi(jtab-1+5*(ima-1)+5)
300     continue
    end do
!
    call jelira(liel2, 'NUTIOC', nb1)
    ASSERT(nb1 .eq. nelt)
!
!-----------------------------------------------------------------------
!     4)  CONSTRUCTION DU .MAILLE
!-----------------------------------------------------------------------
!
    mail2 = modelx//'.MAILLE'
    call wkvect(mail2, 'G V I', nbma, jmail2)
!        write(6,*)'****** KORUPTION : VERIFICATION DU MODEL.MAILLE ******'
    do ima = 1, nbma
!        write(6,*)ima,':',zi(jtab-1+5*(ima-1)+4),zi(jtab-1+5*(ima-1)+5)
        zi(jmail2-1+ima) = zi(jtab-1+5*(ima-1)+5)
    end do
!        write(6,*)'******************************************************'
!
!-----------------------------------------------------------------------
!     5) DUPLICATION DU .NOMA, .NBNO
!                ET DES .NEMA, .SSSA S'ILS EXISTENT
!        PUIS .REPE, .PRNM ET .PRNS AVEC CALL ADALIG CORMGI ET INITEL
!-----------------------------------------------------------------------
!
    call jedupo(ligr1//'.NBNO', 'G', ligr2//'.NBNO', .false._1)
    call jedupo(ligr1//'.LGRF', 'G', ligr2//'.LGRF', .false._1)
    call jeveuo(ligr2//'.LGRF', 'E', vk8=lgrf2)
    lgrf2(2) = modelx
!
    call jedup1(mod1//'.NEMA', 'G', modelx//'.NEMA')
    call jedup1(mod1//'.SSSA', 'G', modelx//'.SSSA')
!
    call adalig(ligr2)
    call cormgi('G', ligr2)
    call initel(ligr2)
!
!-----------------------------------------------------------------------
!     6)  CALCUL DU DÉCOUPAGE EN SOUS-TETRAS, DES FACETTES DE CONTACT
!         ET VERIFICATION DES CRITERES DE CONDITIONNEMENT
!-----------------------------------------------------------------------
!
    call getvtx(motfac, 'PRETRAITEMENTS', iocc=1, scal=k8condi, nbret=iret)
!
    call xcodec(noma, modelx, k8condi, linter, decou)
!
! --- MENAGE
!
    call jedetr(trav)
!
999 continue
    call jedema()
end subroutine
