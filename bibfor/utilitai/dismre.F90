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

subroutine dismre(questi, nomobz, repi, repkz, ierd)
    implicit none
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/dismgd.h"
#include "asterfort/dismlg.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/zerosd.h"
    integer(kind=8) :: repi, ierd
    character(len=*) :: questi, nomobz, repkz
!
!     --     DISMOI(RESUELEM)
!
!       QUESTI : TEXTE PRECISANT LA QUESTION POSEE
!       NOMOBZ : NOM D'UN OBJET DE TYPE LIGREL
!
! OUT : REPI   : REPONSE ( SI ENTIERE )
!       REPKZ  : REPONSE ( SI CHAINE DE CARACTERES )
!       IERD   : CODE RETOUR (0--> OK, -1 --> CHAMP INEXISTANT)
!
! ----------------------------------------------------------------------
!
    integer(kind=8) ::  iret, gd, jnoli, l
    character(len=8) :: k8bid, nogd
    character(len=19) :: nomob
    character(len=24) :: questl, k24
    character(len=32) :: repk
    integer(kind=8), pointer :: desc(:) => null()
! DEB-------------------------------------------------------------------
!
    call jemarq()
!
    ierd = 0
    nomob = nomobz
    repk = ' '
    questl = questi
!
    call jeexin(nomob//'.DESC', iret)
    if (iret .eq. 0) then
        ierd = 1
        goto 9999
    end if
!
    call jeveuo(nomob//'.DESC', 'L', vi=desc)
    gd = desc(1)
    call jenuno(jexnum('&CATA.GD.NOMGD', gd), nogd)
!
    if (questi .eq. 'TYPE_CHAMP') then
        call jelira(nomob//'.DESC', 'DOCU', cval=k8bid)
        ASSERT(k8bid(1:4) .eq. 'RESL')
        repk = 'RESL'
!
    else if (questi .eq. 'TYPE_SUPERVIS') then
        call jelira(nomob//'.DESC', 'DOCU', cval=k8bid)
        ASSERT(k8bid(1:4) .eq. 'RESL')
        repk = '????'
!
    else if (questi .eq. 'TYPE_MATRICE') then
        call dismgd(questi, nogd, repi, repk, ierd)
!
    else if (questi .eq. 'NOM_OPTION') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        repk = zk24(jnoli-1+2) (1:16)
!
    else if (questi .eq. 'NOM_MAILLA') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        call dismlg(questi, zk24(jnoli), repi, repk, ierd)
!
    else if (questl(1:6) .eq. 'NUM_GD') then
        repi = gd
!
    else if (questl(1:6) .eq. 'NOM_GD') then
        repk = nogd
!
    else if (questi .eq. 'NOM_LIGREL') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        repk = zk24(jnoli)
!
    else if (questi .eq. 'PARTITION') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        call dismlg(questi, zk24(jnoli), repi, repk, ierd)
!
    else if (questi .eq. 'NOM_MODELE') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        call dismlg(questi, zk24(jnoli), repi, repk, ierd)
!
    else if (questi .eq. 'TYPE_SCA') then
        l = lxlgut(nogd)
        repk = nogd(l:l)
!
    else if (questi .eq. 'MPI_COMPLET') then
        call jeveuo(nomob//'.NOLI', 'L', jnoli)
        k24 = zk24(jnoli-1+3)
        ASSERT(k24 .eq. 'MPI_COMPLET' .or. k24 .eq. 'MPI_INCOMPLET')
        if (k24 .eq. 'MPI_COMPLET') then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
!
    else if (questi .eq. 'ZERO') then
        if (zerosd('RESUELEM', nomob)) then
            repk = 'OUI'
        else
            repk = 'NON'
        end if
!
    else
        ierd = 1
    end if
!
9999 continue
    repkz = repk
!
    call jedema()
end subroutine
