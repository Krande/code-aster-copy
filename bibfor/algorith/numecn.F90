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
subroutine numecn(modelZ, champZ, numeEquaZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/idenob.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/nueffe.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: modelZ, champZ
    character(len=*), intent(out):: numeEquaZ
!
! --------------------------------------------------------------------------------------------------
!
!  IN/JXIN   : MODELE : MODELE
!  IN/JXIN   : CHAMP  : CHAMP "MODELE" POUR LA NUMEROTATION
!  VAR/JXOUT : NUME   : NUME_EQUA
!
! --------------------------------------------------------------------------------------------------
!
! BUT CREER UN NUME_EQUA (SANS STOCKAGE)
!
! CETTE ROUTINE ETANT APPELEE DANS UNE BOUCLE SUR LES NUMEROS D'ORDRE
! ON CHERCHE A LIMITER LE NOMBRE DE NUME DIFFERENTS CREES
! EN COMPARANT 2 APPELS SUCCESSIFS
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: renumSans = "SANS"
    character(len=24), parameter :: listLigrJv = '&&NUMECN.LIST_LIGR'
    character(len=24), parameter :: listLigrSave = '&&NUMECN.LIST_LIGRS'
    character(len=24) :: noojb, modelLigrel
    character(len=24), pointer :: listLigr(:) => null()
    character(len=19) :: prfchn, ligrelName
    character(len=24) :: numeEquaSave
    integer(kind=8) :: nbLili, iLili, i2, iret, nb2, iexi, nbLigr
    character(len=14) :: numeDof
    aster_logical :: newnum
    save numeEquaSave
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    numeEquaZ = " "
    call dismoi('NUME_EQUA', champZ, 'CHAM_NO', repk=prfchn)
    call dismoi("NOM_LIGREL", modelZ, "MODELE", repk=modelLigrel)

! - Get number of LILI object
    call jelira(prfchn//'.LILI', 'NOMMAX', nbLili)

! - CALCUL DE LLIGR : LISTE DES LIGRELS
    if (nbLili .eq. 1) then
! ----- Only "&MAILLAGE" => only model in LIGREL
        call wkvect(listLigrJv, 'V V K24', 1, vk24=listLigr)
        listLigr(1) = modelLigrel
        nbLigr = 1
    else
!       ON N'AJOUTE QUE LES LIGRELS QUI EXISTENT ENCORE :
        nb2 = 0
        do iLili = 2, nbLili
            call jenuno(jexnum(prfchn//'.LILI', iLili), ligrelName)
            call jeexin(ligrelName//'.LIEL', iret)
            if (iret .ne. 0) then
                if (ligrelName .ne. modelLigrel) then
                    nb2 = nb2+1
                end if
            end if
        end do
        call wkvect(listLigrJv, 'V V K24', nb2+1, vk24=listLigr)
        nbLigr = nb2+1
        i2 = 1
        listLigr(i2) = modelLigrel
        do iLili = 2, nbLili
            call jenuno(jexnum(prfchn//'.LILI', iLili), ligrelName)
            call jeexin(ligrelName//'.LIEL', iret)
            if (iret .ne. 0) then
                if (ligrelName .ne. modelLigrel) then
                    i2 = i2+1
                    listLigr(i2) = ligrelName
                end if
            end if
        end do
    end if

! - Create new numbering or not ?
    newnum = ASTER_TRUE
    call jeexin(listLigrSave, iexi)
    if (iexi .gt. 0) then
        if (idenob(listLigrJv, listLigrSave)) then
            newnum = ASTER_FALSE
        end if
        call jedetr(listLigrSave)
    end if
    call jedupo(listLigrJv, 'V', listLigrSave, ASTER_FALSE)

! - Create new numbering if required or get previous one
    if (newnum) then
! ----- Generate new name of numeDof
        noojb = '12345678.00000.NUME.PRNO'
        call gnomsd(' ', noojb, 10, 14)
        numeDof = noojb(1:14)
        call nueffe(nbLigr, listLigr, 'VG', numeDof, renumSans, modelZ)
        call jedetr(numeDof//'     .ADLI')
        call jedetr(numeDof//'     .ADNE')
        numeEquaZ = numeDof//'.NUME'
        numeEquaSave = numeEquaZ
    else
        numeEquaZ = numeEquaSave
    end if
!
    call jedetr(listLigrJv)
    call jedema()
end subroutine
