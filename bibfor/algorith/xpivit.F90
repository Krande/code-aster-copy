! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine xpivit(jcesd, jcesv, jcesl, ifiss,&
                  ndim, nummae, iface, xpc, ypc,&
                  nvit)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/xxmmvd.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
    integer :: jcesd(10), jcesv(10), jcesl(10)
    integer :: ndim, nummae, iface, ifiss
    real(kind=8) :: xpc, ypc
    integer :: nvit
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (CONTACT - GRANDS GLISSEMENTS)
!
! SERT À DEFINIR SI LE POINT D'INTEGRATION EST VITAL OU NON
!
! TRAVAIL EFFECTUE EN COLLABORATION AVEC L'I.F.P.
!
! ----------------------------------------------------------------------
!
!  JCES*(2)  : POINTEURS DE LA SD SIMPLE DES INFOS SUR ARETES COUPEES
!  JCES*(4)  : POINTEURS DE LA SD SIMPLE DE CONNECTIVITÉ DES FACETTES
! IN NOMA   : NOM DU MAILLAGE (utile pour la numérotation globale)
! IN NDIM  : DIMENSION DU MODELE
! IN NUMMAE : POSITION DE LA MAILLE ESCLAVE
! IN IFACE  : NUMERO LOCAL DE LA FACETTE ESCLAVE
! IN XPC    : COORDONNEE X DU POINT D'INTEGRATION DE CONTACT SUR
!             LA MAILLE ESCLAVE
! IN YPC    : COORDONNEE Y DU POINT D'INTEGRATION DE CONTACT SUR
!             LA MAILLE ESCLAVE
!
! OUT NVIT  : VAUT 1 SI LE POINT D'INTEGRATION EST SUR UNE ARETE
!             VITALE (0 SINON)
! OUT GROUP : NUMERO DE GROUPE SI LE POINT EST SUR UNE ARETE
!             APPARTENANT À UN GROUPE D'ARETES CONNECTÉES (0 SINON)
! OUT NARET : NUMERO D'ARETE DU GROUPE SI LE POINT DE CONTACT EST
!             SUR UNE ARETE APPARTENANT À UN GROUPE (0 SINON)
!
!
!
!
    integer :: zxain, pint, aret, iad
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    zxain = xxmmvd('ZXAIN')
    pint = 0
    aret = 0
    nvit = 0
!
! --- ON REGARDE SI LE POINT D'INTEGRATION EST SUR UNE ARETE
!
    if (ndim .eq. 2) then
        if (xpc .eq. -1) then
            pint = 1
        else if (xpc.eq.1) then
            pint = 2
        endif
    else if (ndim.eq.3) then
        if ((xpc.eq.0) .and. (ypc.eq.0)) then
            pint = 1
        else if ((xpc.eq.1).and.(ypc.eq.0)) then
            pint = 2
        else if ((xpc.eq.0).and.(ypc.eq.1)) then
            pint = 3
        endif
    endif
!
! --- SI IL EST SUR UNE ARETE
!
    if (pint .ne. 0) then
!
! --- ON RECUPERE LE NUMERO DU POINT D'INTERSECTION COUPEE CORESPONDANT
!
        call cesexi('S', jcesd(4), jcesl(4), nummae, 1,&
                    ifiss, ndim*( iface-1)+pint, iad)
        ASSERT(iad.gt.0)
        pint = zi(jcesv(4)-1+iad)
!
! --- ON RECUPERE LE NUMERO D'ARETE COUPÉE ET L'INFO VITAL OU PAS DE
! --- CETTE ARETE COUPEE
!
        call cesexi('S', jcesd(2), jcesl(2), nummae, 1,&
                    ifiss, zxain*(pint- 1)+1, iad)
        ASSERT(iad.gt.0)
! --- Numéro local de l'arête
        aret = nint(zr(jcesv(2)-1+iad))
        call cesexi('S', jcesd(2), jcesl(2), nummae, 1,&
                    ifiss, zxain*(pint- 1)+5, iad)
        ASSERT(iad.gt.0)
! --- 5ème composante du .AI: information arête vitale ou non (0 ou 1)
        nvit = nint(zr(jcesv(2)-1+iad))
!
    endif
!
    call jedema()
end subroutine
