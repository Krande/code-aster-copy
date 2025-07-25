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

subroutine rslipa(nomsd, nopara, nomobj, jadd, nbval)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: jadd, nbval, n1, j1
    character(len=*) :: nomsd, nopara, nomobj
! person_in_charge: jacques.pellet at edf.fr
!
! ----------------------------------------------------------------------
!   extraire d'une sd_resultat, la liste des valeurs d'un parametre
!   et recopier ces valeurs dans l'objet nomobj dont on rend l'adresse.
! ----------------------------------------------------------------------
! in  : nomsd  : nom de la structure "resultat".
! in  : nopara : nom du parametre ('inst','freq', ...)
! in  : nomobj : nom de l'objet jeveux a creer (k24)
! out : jadd   : adresse de l'objet nomobj
! out : nbval  : longueur de l'objet nomobj
!-----------------------------------------------------------------------
! remarques :
!  - l'objet retourne (nomobj) contient les valeurs du parametre dans
!    l'ordre des numeros de rangement.
!    il est "parallele" a l'objet .ordr :
!    do k=1,lonuti(.ordr) :
!       iordr=.ordr(k)
!       nomobj(k) == "rsadpa(nopara,iordr)"
!  - si le parametre est absent pour certains numeros d'ordre, la valeur
!    stockee est "NaN".
!  - cette routine ne fait pas jemarq/jedema pour ne pas
!    invalider l'adresse jeveux jadd
! ----------------------------------------------------------------------
    integer(kind=8) :: kk, jpara, i1, jtava, l1
    character(len=8) :: k8b, tsca
    character(len=5) :: nom1
    character(len=24) :: nomk24
    character(len=16) :: nompar
    character(len=19) :: noms2
    integer(kind=8), pointer :: ordr(:) => null()
! ----------------------------------------------------------------------
!
    noms2 = nomsd
    nompar = nopara
    nomk24 = nomobj
!
    call jenonu(jexnom(noms2//'.NOVA', nompar), i1)
    ASSERT(i1 .gt. 0)
    call jeveuo(jexnum(noms2//'.TAVA', i1), 'L', jtava)
    nom1 = zk8(jtava-1+1) (1:5)
    call jelira(noms2//nom1, 'TYPE', cval=tsca)
    if (tsca .eq. 'K') then
        call jelira(noms2//nom1, 'LTYP', l1)
        if (l1 .eq. 8) then
            tsca = 'K8'
        else if (l1 .eq. 16) then
            tsca = 'K16'
        else if (l1 .eq. 24) then
            tsca = 'K24'
        else if (l1 .eq. 32) then
            tsca = 'K32'
        else if (l1 .eq. 80) then
            tsca = 'K80'
        else
            ASSERT(.false.)
        end if
    end if
!
    call jeveuo(noms2//'.ORDR', 'L', vi=ordr)
    call jelira(noms2//'.ORDR', 'LONUTI', n1)
!
    call wkvect(nomk24, 'V V '//tsca, n1, j1)
!
    do kk = 1, n1
        call rsadpa(noms2, 'L', 1, nompar, ordr(kk), &
                    0, sjv=jpara, styp=k8b, istop=0)
        if (tsca .eq. 'R') then
            zr(j1-1+kk) = zr(jpara)
        else if (tsca .eq. 'C') then
            zc(j1-1+kk) = zc(jpara)
        else if (tsca .eq. 'I') then
            zi(j1-1+kk) = zi(jpara)
        else if (tsca .eq. 'K8') then
            zk8(j1-1+kk) = zk8(jpara)
        else if (tsca .eq. 'K16') then
            zk16(j1-1+kk) = zk16(jpara)
        else if (tsca .eq. 'K24') then
            zk24(j1-1+kk) = zk24(jpara)
        else if (tsca .eq. 'K32') then
            zk32(j1-1+kk) = zk32(jpara)
        else if (tsca .eq. 'K80') then
            zk80(j1-1+kk) = zk80(jpara)
        else
            ASSERT(.false.)
        end if
    end do
!
!     -- pour eviter les effets de bord (,ibid,ibid):
    jadd = j1
    nbval = n1
!
end subroutine
