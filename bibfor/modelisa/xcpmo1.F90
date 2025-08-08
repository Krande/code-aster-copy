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

subroutine xcpmo1(modmes, modthx, modmex)
!
    implicit none
#include "jeveux.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/cormgi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/indk16.h"
#include "asterfort/initel.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/utmess.h"
#include "asterfort/xtmafi.h"
!
    character(len=8) :: modmes, modthx, modmex
!
! ----------------------------------------------------------------------
!
! routine xfem : MODI_MODELE_XFEM
!
!     Traitement particulier realise dans le cas ou le mot-cle
!     MODELE_THER est present. Cette routine est appelee par xcpmod
!     une fois que la copie de modmes dans modmex a ete effectuee
!
! --> on modifie dans cette routine certains objets de modmex en fonction
!     de mothx :
!       - le '.TYFE'   modmex//.TYFE'
!       - le ligrel    modmex//.MODELE'
!
! ----------------------------------------------------------------------
!
! in     modmes : nom du modele mecanique sain (mot-cle modele_in)
! in     modthx : nom du modele thermique x-fem (mot-cle modele_ther)
! in/out modmex : nom du modele mecanique x-fem produit par l'operateur
!
! ----------------------------------------------------------------------
!
    character(len=24) :: lieltp, mesmai, lismai
    character(len=19) :: ligthx, ligmes, ligmex, ligrtp
    character(len=16) :: ktyelt, ktyelm
    character(len=8) :: noma, valk8(3), nommax
    character(len=1) :: k1bid
    integer(kind=8) :: igr, jeltp, imx
    integer(kind=8) :: ima, iexi
    integer(kind=8) :: nmamex
    integer(kind=8) :: nfiss, nbmx, nutyelt, nutyelm
    integer(kind=8) :: nbelmx, nbel2, cpt
!
    integer(kind=8), pointer :: tabmx(:) => null()
    character(len=8), pointer :: fiss(:) => null()
!
    integer(kind=8) :: nma3d, nma2d, nma1d, nenrch
    parameter(nma3d=4)
    parameter(nma2d=2)
    parameter(nma1d=1)
    parameter(nenrch=3)
!   nma3d : HEXA8, PENTA6, PYRAM5, TETRA4 -> 4
!   nma2d : TRIA3, QUAD4 -> 2
!   nma1d : SEG2 -> 1
!   nenrch : XH, XT, XHT
!
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!                        THERMIQUE X-FEM
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!
    integer(kind=8) :: nel3dthx, nelplthx, nelaxthx
    parameter(nel3dthx=(nma3d+nma2d)*nenrch)
    parameter(nelplthx=(nma2d+nma1d)*nenrch)
    parameter(nelaxthx=(nma2d+nma1d)*nenrch)
    character(len=16) :: ele3dthx(nel3dthx)
    character(len=16) :: eleplthx(nelplthx)
    character(len=16) :: eleaxthx(nelaxthx)
!
!   elements 3D lineaires thermiques X-FEM
!   -------------------------------------------------------------------
    data ele3dthx/ &
        !   principaux
        'THER_XH_HEXA8   ', 'THER_XT_HEXA8   ', 'THER_XHT_HEXA8  ', &
        'THER_XH_PENTA6  ', 'THER_XT_PENTA6  ', 'THER_XHT_PENTA6 ', &
        'THER_XH_PYRAM5  ', 'THER_XT_PYRAM5  ', 'THER_XHT_PYRAM5 ', &
        'THER_XH_TETRA4  ', 'THER_XT_TETRA4  ', 'THER_XHT_TETRA4 ', &
        !   de bord
        'THER_XH_FACE4   ', 'THER_XT_FACE4   ', 'THER_XHT_FACE4  ', &
        'THER_XH_FACE3   ', 'THER_XT_FACE3   ', 'THER_XHT_FACE3  '/
!
!   elements PLAN lineaires thermiques X-FEM
!   -------------------------------------------------------------------
    data eleplthx/ &
        !   principaux
        'THPLQU4_XH      ', 'THPLQU4_XT      ', 'THPLQU4_XHT     ', &
        'THPLTR3_XH      ', 'THPLTR3_XT      ', 'THPLTR3_XHT     ', &
        !   de bord
        'THPLSE2_XH      ', 'THPLSE2_XT      ', 'THPLSE2_XHT     '/
!
!   elements AXIS lineaires thermiques X-FEM
!   -------------------------------------------------------------------
    data eleaxthx/ &
        !   principaux
        'THAXQU4_XH      ', 'THAXQU4_XT      ', 'THAXQU4_XHT     ', &
        'THAXTR3_XH      ', 'THAXTR3_XT      ', 'THAXTR3_XHT     ', &
        !   de bord
        'THAXSE2_XH      ', 'THAXSE2_XT      ', 'THAXSE2_XHT     '/
!
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!                        MECANIQUE X-FEM
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!
    integer(kind=8) :: nel3dmex, nelplmex, nelaxmex
    parameter(nel3dmex=(nma3d+nma2d)*nenrch)
    parameter(nelplmex=(2*nma2d+nma1d)*nenrch)
    parameter(nelaxmex=(nma2d+nma1d)*nenrch)
    character(len=16) :: ele3dmex(nel3dmex)
    character(len=16) :: eleplmex(nelplmex)
    character(len=16) :: eleaxmex(nelaxmex)
!
!   elements 3D lineaires mecaniques X-FEM
!   -------------------------------------------------------------------
    data ele3dmex/ &
        !   principaux
        'MECA_XH_HEXA8   ', 'MECA_XT_HEXA8   ', 'MECA_XHT_HEXA8  ', &
        'MECA_XH_PENTA6  ', 'MECA_XT_PENTA6  ', 'MECA_XHT_PENTA6 ', &
        'MECA_XH_PYRAM5  ', 'MECA_XT_PYRAM5  ', 'MECA_XHT_PYRAM5 ', &
        'MECA_XH_TETRA4  ', 'MECA_XT_TETRA4  ', 'MECA_XHT_TETRA4 ', &
        !   de bord
        'MECA_XH_FACE4   ', 'MECA_XT_FACE4   ', 'MECA_XHT_FACE4  ', &
        'MECA_XH_FACE3   ', 'MECA_XT_FACE3   ', 'MECA_XHT_FACE3  '/
!
!   elements C_PLAN/D_PLAN lineaires mecaniques X-FEM
!   -------------------------------------------------------------------
    data eleplmex/ &
        !   principaux
        'MECPQU4_XH      ', 'MECPQU4_XT      ', 'MECPQU4_XHT     ', &
        'MECPTR3_XH      ', 'MECPTR3_XT      ', 'MECPTR3_XHT     ', &
        'MEDPQU4_XH      ', 'MEDPQU4_XT      ', 'MEDPQU4_XHT     ', &
        'MEDPTR3_XH      ', 'MEDPTR3_XT      ', 'MEDPTR3_XHT     ', &
        !   de bord
        'MEPLSE2_XH      ', 'MEPLSE2_XT      ', 'MEPLSE2_XHT     '/
!
!   elements AXIS lineaires mecaniques X-FEM
!   -------------------------------------------------------------------
    data eleaxmex/ &
        !   principaux
        'MEAXQU4_XH      ', 'MEAXQU4_XT      ', 'MEAXQU4_XHT     ', &
        'MEAXTR3_XH      ', 'MEAXTR3_XT      ', 'MEAXTR3_XHT     ', &
        !   de bord
        'MEAXSE2_XH      ', 'MEAXSE2_XT      ', 'MEAXSE2_XHT     '/
!
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!                        MECANIQUE CLASSIQUE
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!
    integer(kind=8) :: nel3dmec, nelplmec, nelaxmec
    parameter(nel3dmec=nma3d+nma2d)
    parameter(nelplmec=2*nma2d+nma1d)
    parameter(nelaxmec=nma2d+nma1d)
    character(len=16) :: ele3dmec(nel3dmec)
    character(len=16) :: eleplmec(nelplmec)
    character(len=16) :: eleaxmec(nelaxmec)
    character(len=8) :: partsd
    integer(kind=8), pointer :: mmes(:) => null()
    integer(kind=8), pointer :: mmex(:) => null()
    integer(kind=8), pointer :: mthx(:) => null()
!
!   elements 3D lineaires mecaniques classiques
! ---------------------------------------------------------------------
    data ele3dmec/ &
        !   principaux
        'MECA_HEXA8      ', &
        'MECA_PENTA6     ', &
        'MECA_PYRAM5     ', &
        'MECA_TETRA4     ', &
        !   de bord
        'MECA_FACE4      ', &
        'MECA_FACE3      '/
!
!   elements C_PLAN/D_PLAN lineaires mecaniques classiques
! ---------------------------------------------------------------------
    data eleplmec/ &
        !   principaux
        'MECPQU4         ', &
        'MECPTR3         ', &
        'MEDPQU4         ', &
        'MEDPTR3         ', &
        !   de bord
        'MEPLSE2         '/
!
!   elements AXIS lineaires mecaniques classiques
! ---------------------------------------------------------------------
    data eleaxmec/ &
        !   principaux
        'MEAXQU4         ', &
        'MEAXTR3         ', &
        !   de bord
        'MEAXSE2         '/
!
! ---------------------------------------------------------------------
! - Debut code
! ---------------------------------------------------------------------
!
    call jemarq()
!
! ---------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modmex, 'MODELE', repk=ligmex)
    call dismoi('NOM_LIGREL', modthx, 'MODELE', repk=ligthx)
    call dismoi('NOM_LIGREL', modmes, 'MODELE', repk=ligmes)
!
! - recuperation de la liste de toutes les mailles fissurees
! - (de dimension n et n-1 car appel a xtmafi avec ndim == 0)
!
    lismai = '&&XMOT2M.NUM_MAILLES'
    mesmai = '&&XMOT2M.MES_MAILLES'
!
    call dismoi('NOM_MAILLA', modthx, 'MODELE', repk=noma)
    call dismoi('NB_FISS_XFEM', modthx, 'MODELE', repi=nfiss)
    call jeveuo(modthx//'.FISS', 'L', vk8=fiss)
!
    call xtmafi(0, fiss, nfiss, lismai, mesmai, nbmx, model=modthx)
    call jeveuo(lismai, 'L', vi=tabmx)
!
! - recuperation du '.TYFE' de modthx et modmes
!
    call jeveuo(ligthx//'.TYFE', 'L', vi=mthx)
    call jeveuo(ligmes//'.TYFE', 'L', vi=mmes)
!
! - on s'assure que toute maille affectee par un element thermique
! - enrichi dans modthx est bien affectee par un element mecanique
! - sain dans modmes
!
    do ima = 1, nbmx
!
        imx = tabmx(ima)
!
        nutyelt = mthx(imx)
        nutyelm = mmes(imx)
!
        if (nutyelm .eq. 0) then
            nommax = int_to_char8(imx)
            valk8(1) = nommax
            valk8(2) = modthx
            valk8(3) = modmes
            call utmess('F', 'XFEM_85', nk=3, valk=valk8)
        end if
!
    end do
!
! - copie integrale du contenu de modmes dans modmex
!
    call copisd('MODELE', 'G', modmes, modmex)
!
! - on supprime modmex//'.PARTSD' s'il existe car MODI_MODELE_XFEM
! - avec mot-cle FISSURE ne recree pas cet objet s'il existe dans
! - le modele sain renseigne dans MODELE_IN
!
    call dismoi('PARTITION', modmex, 'MODELE', repk=partsd)
    call exisd('PARTITION', partsd, iexi)
    if (iexi .ne. 0) then
        call detrsd('PARTITION', partsd)
    end if
!
! - recuperation du '.TYFE' de modmex
!
    call jeveuo(ligmex//'.TYFE', 'E', vi=mmex)
    call jelira(ligmex//'.TYFE', 'LONMAX', nmamex, k1bid)
!
! - modification du '.TYFE' de modmex pour les mailles fissurees
!
    do ima = 1, nbmx
!
        imx = tabmx(ima)
!
        nutyelt = mthx(imx)
        call jenuno(jexnum('&CATA.TE.NOMTE', nutyelt), ktyelt)
!
        nutyelm = mmex(imx)
        call jenuno(jexnum('&CATA.TE.NOMTE', nutyelm), ktyelm)
!
!       MECANIQUE 3D
        if (indk16(ele3dmec, ktyelm, 1, nel3dmec) .gt. 0) then
            if (indk16(ele3dthx, ktyelt, 1, nel3dthx) .eq. 0) then
                ASSERT(.false.)
            end if
            ktyelm = ktyelm(1:4)//ktyelt(5:16)
            if (indk16(ele3dmex, ktyelm, 1, nel3dmex) .eq. 0) then
                ASSERT(.false.)
            end if
!
!       MECANIQUE C_PLAN/D_PLAN
        else if (indk16(eleplmec, ktyelm, 1, nelplmec) .gt. 0) then
            if (indk16(eleplthx, ktyelt, 1, nelplthx) .eq. 0) then
                ASSERT(.false.)
            end if
            ktyelm = ktyelm(1:4)//ktyelt(5:16)
            if (indk16(eleplmex, ktyelm, 1, nelplmex) .eq. 0) then
                ASSERT(.false.)
            end if
!
!       MECANIQUE AXIS
        else if (indk16(eleaxmec, ktyelm, 1, nelaxmec) .gt. 0) then
            if (indk16(eleaxthx, ktyelt, 1, nelaxthx) .eq. 0) then
                ASSERT(.false.)
            end if
            ktyelm = ktyelm(1:4)//ktyelt(5:16)
            if (indk16(eleaxmex, ktyelm, 1, nelaxmex) .eq. 0) then
                ASSERT(.false.)
            end if
!
        else
            ASSERT(.false.)
!
        end if
!
        call jenonu(jexnom('&CATA.TE.NOMTE', ktyelm), nutyelm)
        mmex(imx) = nutyelm
!
    end do
!
! - ligrel temporaire ligrtp (remplacera celui de modmex : ligmex)
!
    ligrtp = '&&XMOT2M'//'.MODELE'
    lieltp = ligrtp//'.LIEL'
!
    nbelmx = 0
    do igr = 1, nbgrel(ligmex)
        nbelmx = nbelmx+nbelem(ligmex, igr)
    end do
!
    call jecrec(lieltp, 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbelmx)
    call jeecra(lieltp, 'LONT', 2*nbelmx)
!
    cpt = 0
    do ima = 1, nmamex
        if (mmex(ima) .eq. 0) cycle
        cpt = cpt+1
        call jecroc(jexnum(lieltp, cpt))
        call jeecra(jexnum(lieltp, cpt), 'LONMAX', 2)
        call jeveuo(jexnum(lieltp, cpt), 'E', jeltp)
        zi(jeltp-1+1) = ima
        zi(jeltp-1+2) = mmex(ima)
    end do
!
    call jelira(lieltp, 'NUTIOC', nbel2)
    ASSERT(nbel2 .eq. nbelmx)
!
    call jedupo(ligmex//'.NBNO', 'G', ligrtp//'.NBNO', .false._1)
    call jedupo(ligmex//'.LGRF', 'G', ligrtp//'.LGRF', .false._1)
    call jedupo(ligmex//'.TYFE', 'G', ligrtp//'.TYFE', .false._1)
!
! - on ecrase ligmex avec ligrtp
!
    call adalig(lieltp)
    call cormgi('V', lieltp)
    call initel(lieltp)
!
    call detrsd('LIGREL', ligmex)
    call copisd('LIGREL', 'G', ligrtp, ligmex)
    call detrsd('LIGREL', ligrtp)
!
! - menage final
!
    call jedetr(mesmai)
    call jedetr(lismai)
!
    call jedema()
!
end subroutine
