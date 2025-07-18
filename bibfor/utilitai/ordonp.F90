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

subroutine ordonp(nomfon)
    implicit none
#include "jeveux.h"
!
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedupo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
    character(len=19) :: nomfon
! person_in_charge: mathieu.courtois at edf.fr
! ----------------------------------------------------------------------
!     APPELE PAR ORDONN POUR REORDONNER LES FONCTIONS D UNE NAPPE
!     PAR ORDRE CROISSANT DES PARAMETRES
    character(len=19) :: fonc0
    character(len=24) :: chval, chpar, sfval, sfpar
    integer(kind=8) :: ipar, lpar, nbpara, ior, ival, lval
    integer(kind=8) :: i, j, it, nbp, nbpt
    real(kind=8) :: xt
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     OBJET INITIAL RECOPIE DANS FONC0
    fonc0 = '&&ORDONP.FONC      '
    chval = nomfon//'.VALE'
    chpar = nomfon//'.PARA'
    sfval = fonc0//'.VALE'
    sfpar = fonc0//'.PARA'
!
    call jelira(chpar, 'LONUTI', nbpara)
!     RECUPERE LES PARAMETRES
    call jedupo(chpar, 'V', sfpar, .false._1)
    call jeveuo(sfpar, 'E', ipar)
    call jelira(sfpar, 'LONUTI', nbpara)
!
    call jedupo(chval, 'V', sfval, .false._1)
!
    call jedetr(chpar)
    call jedetr(chval)
!
!     TABLEAU D'ORDRE
    call wkvect(fonc0//'.ORDR', 'V V I', nbpara, ior)
    do i = 1, nbpara
        zi(ior-1+i) = i
    end do
!
!     TRI DES PARAMETRES
    do i = 1, nbpara-1
        do j = i+1, nbpara
            if (zr(ipar-1+i) .gt. zr(ipar-1+j)) then
                xt = zr(ipar-1+i)
                it = zi(ior-1+i)
                zr(ipar-1+i) = zr(ipar-1+j)
                zi(ior-1+i) = zi(ior-1+j)
                zr(ipar-1+j) = xt
                zi(ior-1+j) = it
            end if
        end do
    end do
!
!     CALCULE LA TAILLE CUMULEE DE LA COLLECTION
    nbpt = 0
    do i = 1, nbpara
        call jelira(jexnum(sfval, i), 'LONMAX', nbp)
        nbpt = nbpt+nbp
    end do
!
!     --- CREATION DE L'OBJET NOMFON.PARA ---
    call wkvect(chpar, 'G V R', nbpara, lpar)
!     --- CREATION DE LA COLLECTION NOMFON.VALE ---
    call jecrec(chval, 'G V R', 'NU', 'CONTIG', 'VARIABLE', &
                nbpara)
    call jeecra(chval, 'LONT', nbpt)
    do i = 1, nbpara
!        REMPLISSAGE DU .PARA
        zr(lpar-1+i) = zr(ipar-1+i)
!        REMPLISSAGE DES .VALE EN FONCTION DE L'ORDRE
        call jelira(jexnum(sfval, zi(ior-1+i)), 'LONMAX', nbp)
        call jeveuo(jexnum(sfval, zi(ior-1+i)), 'E', ival)
        call jecroc(jexnum(chval, i))
        call jeecra(jexnum(chval, i), 'LONMAX', nbp)
        call jeecra(jexnum(chval, i), 'LONUTI', nbp)
        call jeveuo(jexnum(chval, i), 'E', lval)
        do j = 1, nbp
            zr(lval+j-1) = zr(ival+j-1)
        end do
    end do
!
!     DESTRUCTION DES OBJETS DE TRAVAIL
    call jedetr(sfpar)
    call jedetr(sfval)
    call jedetr(fonc0//'.ORDR')
!
    call jedema()
end subroutine
