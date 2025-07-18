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

subroutine mgutdm(mdgenz, nmsstz, nusst, questi, repi, &
                  repkz)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 13/11/92
!-----------------------------------------------------------------------
!  BUT:      < MODELE GENERALISE UTILITAIRE DIS MOI >
!
!  UTILITAIRE PERMETTANT D'ACCEDER AUX CONCEPT RELATIFS AUX
!  SOUS-STRUCTURES D'UN MODELE GENERALISE
!
!  LISTE DES QUESTIONS POSSIBLES:
!    NOM_MACR_ELEM
!    NOM_BASE_MODALE
!    NOM_MAILLAGE
!    NOM_MODELE
!    NOM_NUME_DDL
!    NOM_LIST_INTERF
!    NB_CMP_MAX
!
!-----------------------------------------------------------------------
!
! MDGENZ   /I/: NOM UTILISATEUR DU MODELE GENERALISE
! NMSSTZ   /I/: NOM K8 DE LA SOUS-STRUCTURE
! NUSST    /I/: NUMERO DE LA SOUS-STRUCTURE
! QUESTI   /I/: QUESTION
! REPI     /O/: REPONSE ENTIERE
! REPKZ    /O/: REPONSE CARACTERE
!
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
!
!
    integer(kind=8) :: repi, iret, llmcl, llref, nuss, nusst
    integer(kind=8) :: vali
    character(len=*) :: questi
    character(len=24) :: repk
    character(len=24) :: valk(2), nume
    character(len=8) :: modgen, nommcl, basmod, nomsst
    character(len=*) :: mdgenz, nmsstz, repkz
    integer(kind=8), pointer :: desc(:) => null()
!
!-----------------------------------------------------------------------
!
!------------RECUPERATION NUMRERO DE SOUS-STRUCTURE ET VERIFS-----------
!
    call jemarq()
    modgen = mdgenz
    nomsst = nmsstz
    repk = repkz
!
    if (nomsst(1:1) .ne. ' ') then
        call jenonu(jexnom(modgen//'      .MODG.SSNO', nomsst), nuss)
        if (nuss .eq. 0) then
            valk(1) = modgen
            valk(2) = nomsst
            call utmess('F', 'ALGORITH13_49', nk=2, valk=valk)
        end if
    else
        nuss = nusst
        call jeexin(jexnum(modgen//'      .MODG.SSME', nuss), iret)
        if (nuss .eq. 0) then
            valk(1) = modgen
            vali = nuss
            call utmess('F', 'ALGORITH13_50', sk=valk(1), si=vali)
        end if
        call jenuno(jexnum(modgen//'      .MODG.SSNO', nuss), nomsst)
    end if
!
!
    if (questi(1:13) .eq. 'NOM_MACR_ELEM') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        repk(1:8) = zk8(llmcl)
    else if (questi(1:15) .eq. 'NOM_BASE_MODALE') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        nommcl = zk8(llmcl)
        call jeveuo(nommcl//'.MAEL_REFE', 'L', llref)
        repk(1:8) = zk24(llref)
    else if (questi(1:12) .eq. 'NOM_MAILLAGE') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        nommcl = zk8(llmcl)
        call jeveuo(nommcl//'.MAEL_REFE', 'L', llref)
        repk(1:8) = zk24(llref+1)
    else if (questi(1:12) .eq. 'NOM_NUME_DDL') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        nommcl = zk8(llmcl)
        call jeveuo(nommcl//'.MAEL_REFE', 'L', llref)
        basmod(1:8) = zk24(llref)
        call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=repk)
    else if (questi(1:12) .eq. 'NOM_MODELE  ') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        nommcl = zk8(llmcl)
        call jeveuo(nommcl//'.MAEL_REFE', 'L', llref)
        basmod(1:8) = zk24(llref)
        call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=nume)
        call dismoi('NOM_MODELE', nume, 'NUME_DDL', repk=repk)
    else if (questi(1:15) .eq. 'NOM_LIST_INTERF') then
        call jeveuo(jexnum(modgen//'      .MODG.SSME', nuss), 'L', llmcl)
        nommcl = zk8(llmcl)
        call jeveuo(nommcl//'.MAEL_REFE', 'L', llref)
        basmod(1:8) = zk24(llref)
!       call utimsd(6, 2, .false._1, .true._1,basmod(1:8)//'           .REFD', 1, ' ')
        call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=repk, arret='C', &
                    ier=iret)
    else if (questi(1:10) .eq. 'NB_CMP_MAX') then
        call jeveuo(modgen//'      .MODG.DESC', 'L', vi=desc)
        repi = desc(2)
    else
        repk = questi
        call utmess('F', 'UTILITAI_49', sk=repk)
        goto 999
    end if
!
999 continue
    repkz = repk
    call jedema()
end subroutine
