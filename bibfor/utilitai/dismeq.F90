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

subroutine dismeq(questi, nomobz, repi, repkz, ierd)
    implicit none
!     --     DISMOI(NUME_EQUA)
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/dismlg.h"
#include "asterfort/dismgd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexnom.h"
!
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi
    character(len=*) :: nomobz, repkz
    character(len=32) :: repk
! ----------------------------------------------------------------------
!     IN:
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE NUME_EQUA
!     OUT:
!       REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, 1 --> PB)
!
!     LISTE DES QUESTIONS ADMISSIBLES:
!        'NB_DDLACT'
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=19) :: noligr
    character(len=19) :: nomob
    aster_logical :: isLagr, isDbLagr
    integer(kind=8), pointer :: nequ(:) => null()
    integer(kind=8), pointer :: delg(:) => null()
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, nbddlb, nbnos, neq, nlili
    character(len=24), pointer :: refn(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    nomob = nomobz
    repk = ' '
    repi = 0
    ierd = 0
!
    if (questi .eq. 'NB_DDLACT') then
!     --------------------------------
        call jelira(nomob//'.NUEQ', 'LONMAX', neq)
        call jelira(nomob//'.LILI', 'NUTIOC', nlili)
        nbddlb = 0
        do i = 2, nlili
            call jenuno(jexnum(nomob//'.LILI', i), noligr)
            call dismlg('NB_NO_SUP', noligr, nbnos, repk, ierd)
            nbddlb = nbddlb+nbnos
        end do
        repi = neq-3*(nbddlb/2)
!
    else if (questi .eq. 'JOINTS') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        repk = refn(5) (1:19)
!
    else if (questi .eq. 'NB_EQUA') then
!     --------------------------------
        call jeveuo(nomob//'.NEQU', 'L', vi=nequ)
        repi = nequ(1)
!
!
    else if (questi .eq. 'NOM_GD') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        repk = refn(2) (1:8)
!
    else if (questi .eq. 'NUM_GD') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        call jenonu(jexnom('&CATA.GD.NOMGD', refn(2) (1:8)), repi)
!
    else if (questi(1:9) .eq. 'NUM_GD_SI') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        call dismgd(questi, refn(2) (1:8), repi, repk, ierd)
!
    else if (questi .eq. 'NOM_MODELE') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        repk = refn(3)
!
    else if (questi .eq. 'NOM_MAILLA') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        repk = refn(1)
    else if (questi .eq. 'TYPE_SUPERVIS' .or. questi .eq. 'TYPE_SCA') then
        call jeveuo(nomob//'.REFN', 'L', vk24=refn)
        call dismgd(questi, refn(2) (1:8), repi, repk, ierd)
!
    else if (questi .eq. 'EXIS_LAGR') then
        call jeveuo(nomob//'.DELG', 'L', vi=delg)
        call jeveuo(nomob//'.NEQU', 'L', vi=nequ)
        neq = nequ(1)
        repk = 'NON'
        do i = 1, neq
            if (delg(i) .lt. 0) then
                REPK = 'OUI'
                goto 10
            end if
        end do
10      continue
    else if (questi .eq. 'SIMP_LAGR') then
        call jeveuo(nomob//'.DELG', 'L', vi=delg)
        call jeveuo(nomob//'.NEQU', 'L', vi=nequ)
        neq = nequ(1)
        isLagr = ASTER_FALSE
        isDbLagr = ASTER_FALSE
        do i = 1, neq
            if (delg(i) .lt. 0) then
                isLagr = ASTER_TRUE
                if (delg(i) .eq. -2) then
                    isDbLagr = ASTER_TRUE
                    goto 20
                end if
            end if
        end do
20      continue
        repk = 'NON'
        if (isLagr .and. .not. isDbLagr) repk = 'OUI'
    else
        ierd = 1
    end if
!
    repkz = repk
    call jedema()
end subroutine
