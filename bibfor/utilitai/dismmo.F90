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
subroutine dismmo(questi, nomobz, repi, repkz, ierd)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismlg.h"
#include "asterfort/dismma.h"
#include "asterfort/dismqu.h"
#include "asterfort/dismte.h"
#include "asterfort/dismzc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lteatt.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi, nomobz, repkz
!
! --------------------------------------------------------------------------------------------------
!
!     --     DISMOI(MODELE)
!    IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!    OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ico, igrel
    integer(kind=8) :: iret, elemTypeNume, nbgrel, lielSize
    character(len=4) :: tytm
    character(len=8) :: mesh, model
    character(len=16) :: elemTypeName, nomodl, nomod2
    character(len=19) :: modelLigrel
    character(len=32) :: repk
    character(len=8), pointer :: lgrf(:) => null(), k8cond(:) => null()
    integer(kind=8), pointer :: nfis(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    repk = ' '
    repi = 0
    ierd = 0
!
    model = nomobz
    modelLigrel = model//'.MODELE'
!
    if (questi .eq. 'NOM_LIGREL') then
        repk = modelLigrel

    else if (questi .eq. 'Z_CST') then
        call dismzc(questi, modelLigrel, repi, repk, ierd)

    elseif (questi .eq. 'PARTITION') then
        call dismlg(questi, modelLigrel, repi, repk, ierd)

    elseif ((questi .eq. 'DIM_GEOM') .or. (questi .eq. 'NB_SM_MAILLA') .or. &
            (questi .eq. 'NB_SS_ACTI') .or. (questi .eq. 'NB_NL_MAILLA')) then
        call dismlg(questi, modelLigrel, repi, repk, ierd)

    elseif ((questi .eq. 'AXIS') .or. (questi .eq. 'CALC_RIGI') .or. &
            (questi .eq. 'PHENOMENE')) then
        call dismlg(questi, modelLigrel, repi, repk, ierd)

    elseif ((questi .eq. 'EXI_AMOR') .or. (questi .eq. 'EXI_ELEM') .or. &
            (questi .eq. 'EXI_RDM') .or. (questi .eq. 'EXI_COQUE') .or. &
            (questi .eq. 'EXI_GRILLE') .or. (questi .eq. 'EXI_COQ3D') .or. &
            (questi .eq. 'EXI_PLAQUE') .or. (questi .eq. 'EXI_TUYAU') .or. &
            (questi .eq. 'EXI_POUX') .or. (questi .eq. 'EXI_STRX') .or. &
            (questi .eq. 'EXI_STR2') .or. (questi .eq. 'EXI_THM') .or. &
            (questi .eq. 'EXI_HHO') .or. (questi .eq. 'EXI_HHO_CSTE') .or. &
            (questi .eq. 'EXI_HHO_LINE') .or. &
            (questi .eq. 'EXI_HHO_QUAD') .or. (questi .eq. 'EXI_HHO_CUBI') .or. &
            (questi .eq. 'EXI_HHO_QUAR') .or. &
            (questi .eq. 'EXI_NO_HHO') .or. &
            (questi .eq. 'EXI_AXIS') .or. (questi .eq. 'EXI_COQSOL') .or. &
            (questi .eq. 'EXI_IMPE_ABSO') .or. (questi .eq. 'EXI_CABLE') .or. &
            (questi .eq. 'EXI_POUTRE') .or. (questi .eq. 'EXI_INCO') .or. &
            (questi .eq. 'EXI_SECH') .or. (questi .eq. 'EXI_NON_SECH')) then
        call dismlg(questi, modelLigrel, repi, repk, ierd)

    else if (questi .eq. 'ELEM_VOLU_QUAD') then
        call dismqu(questi, modelLigrel, repi, repk, ierd)

    else if (questi .eq. 'NOM_MAILLA') then
        call jeveuo(modelLigrel//'.LGRF', 'L', vk8=lgrf)
        mesh = lgrf(1)
        repk = mesh

    else if (questi .eq. 'MODELISATION') then
        call jeexin(modelLigrel//'.LIEL', iret)
        if (iret .eq. 0) goto 20
        call jelira(modelLigrel//'.LIEL', 'NUTIOC', nbgrel)
        if (nbgrel .le. 0) goto 20
        ico = 0
        nomodl = ' '
!
        do igrel = 1, nbgrel
            call jeveuo(jexnum(modelLigrel//'.LIEL', igrel), 'L', vi=liel)
            call jelira(jexnum(modelLigrel//'.LIEL', igrel), 'LONMAX', lielSize)
            elemTypeNume = liel(lielSize)
            call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
            call dismte('MODELISATION', elemTypeName, repi, repk, ierd)
            nomod2 = repk(1:16)
!           -- ON ESPERE QUE LES NOMTE '#PLUSIEURS' SONT DES ELEMENTS
!              DE BORD ET QUE L'ON PEUT LES IGNORER ET QU'IL EN RESTE
!              D'AUTRES PLUS SIGNIFICATIFS :
            if (nomod2 .ne. '#PLUSIEURS') then
                if (nomodl .ne. nomod2) then
                    ico = ico+1
                    nomodl = nomod2
                end if
            end if
        end do
        ASSERT(ico .ge. 1)
!
        if (ico .eq. 1) then
            repk = nomodl
        else if (ico .gt. 1) then
            repk = '#PLUSIEURS'
        end if
        goto 30
!
20      continue
        repk = '#AUCUNE'
!
30      continue

    elseif ((questi .eq. 'NB_NO_MAILLA') .or. (questi .eq. 'NB_MA_MAILLA') .or. &
            (questi .eq. 'NB_NO_SS_MAX')) then
        call jeveuo(modelLigrel//'.LGRF', 'L', vk8=lgrf)
        mesh = lgrf(1)
        call dismma(questi, mesh, repi, repk, ierd)

    else if (questi .eq. 'NB_FISS_XFEM') then
        call jeexin(model//'.NFIS', iret)
        if (iret .gt. 0) then
            call jeveuo(model//'.NFIS', 'L', vi=nfis)
            repi = nfis(1)
        else
            repi = 0
        end if

    else if (questi .eq. 'PRE_COND_XFEM') then
        call jeexin(model//'.PRE_COND', iret)
        if (iret .gt. 0) then
            call jeveuo(model//'.PRE_COND', 'L', vk8=k8cond)
            repk = k8cond(1)
        else
            repk = 'NON'
        end if

    else if (questi .eq. 'BESOIN_MATER') then
        call jeexin(modelLigrel//'.LIEL', iret)
        if (iret .gt. 0) then
            call jelira(modelLigrel//'.LIEL', 'NUTIOC', nbgrel)
            repk = 'NON'
            do igrel = 1, nbgrel
                call jeveuo(jexnum(modelLigrel//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(modelLigrel//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                call dismte('MODELISATION', elemTypeName, repi, repk, ierd)
                if (lteatt('TYPMOD', '0D', typel=elemTypeName)) then
                    repk = 'OUI'
                    goto 70
                end if
            end do
        else
            repk = 'NON'
        end if

    else if (questi .eq. 'EXI_ELTVOL') then
!          (EXISTENCE D'ELEMENTS DONT LA MAILLE EST VOLUMIQUE)
!
        call jeexin(modelLigrel//'.LIEL', iret)
        if (iret .gt. 0) then
            call jelira(modelLigrel//'.LIEL', 'NUTIOC', nbgrel)
            repk = 'NON'
            do igrel = 1, nbgrel
                call jeveuo(jexnum(modelLigrel//'.LIEL', igrel), 'L', vi=liel)
                call jelira(jexnum(modelLigrel//'.LIEL', igrel), 'LONMAX', lielSize)
                elemTypeNume = liel(lielSize)
                call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
                call dismte('TYPE_TYPMAIL', elemTypeName, repi, tytm, ierd)
                if (tytm .eq. 'VOLU') then
                    repk = 'OUI'
                    goto 70
                end if
            end do
        else
            repk = 'NON'
        end if
    else if (questi .eq. 'EXI_XFEM') then
        call jeexin(model//'.FISS', iret)
        if (iret .gt. 0) then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
    else
        goto 60
    end if
!
    goto 70
!
!     -- SORTIE ERREUR :
!     ------------------
60  continue
    ierd = 1
    goto 80
!
!     -- SORTIE NORMALE :
!     ------------------
70  continue
    ierd = 0
    repkz = repk
!
!
80  continue
    call jedema()
end subroutine
