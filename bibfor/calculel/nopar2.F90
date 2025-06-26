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
subroutine nopar2(nomopt, nomgd, statut, nompar, istop, &
                  iret)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: nomopt, nomgd, statut
    character(len=8), intent(out) :: nompar
    integer(kind=8), intent(in), optional :: istop
    integer(kind=8), intent(out), optional :: iret
!
! ----------------------------------------------------------------------
!     ENTREES:
!     NOMOPT   : NOM D'1 OPTION
!     NOMGD    : NOM D'1 GRANDEUR
!     STATUT : 'IN'/'OUT'/'INOUT'
!     ISTOP : si absent, on émet les messages d'erreurs le cas échéant
!             si 0, on sort sans erreur (mais avec un code retour)
!
!     SORTIES:
!     NOMPAR : NOM DU PARAMETRE DE L'OPTION NOMOPT
!              QUI CORRESPOND A LA GRANDEUR NOMGD ET AU STATUT STATUT
!              RQUE : ERREUR <F> SI ON EN TROUVE : 0,2,3,...
!              RQUE : SI INOUT, ON CHERCHE D'ABORD DANS OUT ET SI ON NE
!                     TROUVE PAS, ON CHERCHE DANS IN
!
!     SI NOMGD=' ' ET SI STATUT='OUT' :
!         - SI L'OPTION N'A QU'UN PARAMETRE 'OUT', ON LE REND
!         - SI L'OPTION A PLUSIEURS PARAMETRES 'OUT' => ERREUR <F>
!
!     IRET : 0 si paramètre trouvé, 1 sinon
!
! ----------------------------------------------------------------------
!
!
!     VARIABLES LOCALES:
!     ------------------
    integer(kind=8) :: opt, nbin, nbout, nbtrou, itrou, gd, gd2, iadesc, iaoppa, kk
    integer(kind=8) :: istop2
    character(len=16) :: nomop2
    character(len=8) :: nomgd2
    character(len=8) :: statu2, outrou
    character(len=24) :: valk(3)
!
! DEB-------------------------------------------------------------------
    nomop2 = nomopt
    nomgd2 = nomgd
    statu2 = statut
    istop2 = 1
    if (present(istop)) then
        istop2 = istop
    end if
!
    call jenonu(jexnom('&CATA.OP.NOMOPT', nomop2), opt)
    if (nomgd2 .ne. ' ') call jenonu(jexnom('&CATA.GD.NOMGD', nomgd2), gd)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iadesc)
    call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iaoppa)
    nbin = zi(iadesc-1+2)
    nbout = zi(iadesc-1+3)
!
    nbtrou = 0
    outrou = ' '
!
!
    if (statu2 .eq. 'OUT') then
        if (nomgd2 .eq. ' ') then
            ASSERT(nbout .eq. 1)
            nbtrou = nbtrou+1
            itrou = 1
            outrou = 'OUT'
        else
            do kk = 1, nbout
                gd2 = zi(iadesc-1+4+nbin+kk)
                if (gd .eq. gd2) then
                    nbtrou = nbtrou+1
                    itrou = kk
                    outrou = 'OUT'
                end if
            end do
        end if
!
    else if (statu2 .eq. 'IN') then
        do kk = 1, nbin
            gd2 = zi(iadesc-1+4+kk)
            if (gd .eq. gd2) then
                nbtrou = nbtrou+1
                itrou = kk
                outrou = 'IN'
            end if
        end do
!
    else if (statu2 .eq. 'INOUT') then
        do kk = 1, nbout
            gd2 = zi(iadesc-1+4+nbin+kk)
            if (gd .eq. gd2) then
                nbtrou = nbtrou+1
                itrou = kk
                outrou = 'OUT'
            end if
        end do
!
        if (nbtrou .eq. 0) then
            do kk = 1, nbin
                gd2 = zi(iadesc-1+4+kk)
                if (gd .eq. gd2) then
                    nbtrou = nbtrou+1
                    itrou = kk
                    outrou = 'IN'
                end if
            end do
        end if
!
    else
        ASSERT(.false.)
    end if
!
    if (nbtrou .eq. 1) then
        if (outrou .eq. 'OUT') then
            nompar = zk8(iaoppa-1+nbin+itrou)
!
        else if (outrou .eq. 'IN') then
            nompar = zk8(iaoppa-1+itrou)
!
        else
            ASSERT(.false.)
        end if
!
    else if (istop2 .eq. 1) then
        if (nbtrou .eq. 0) then
            valk(1) = statu2
            valk(2) = nomgd2
            valk(3) = nomop2
            call utmess('F', 'CALCULEL3_84', nk=3, valk=valk)
        end if
        if (nbtrou .gt. 1) then
            valk(1) = statu2
            valk(2) = nomgd2
            valk(3) = nomop2
            call utmess('F', 'CALCULEL3_85', nk=3, valk=valk)
        end if
    end if
!
    if (present(iret)) then
        iret = 0
        if (nbtrou .ne. 1) iret = 1
    else
        ASSERT(istop2 .eq. 1)
    end if
!
end subroutine
