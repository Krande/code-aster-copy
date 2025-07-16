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

subroutine checkNormals(model, slave, master)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/chbord.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/orilma.h"
#include "asterfort/ornorm.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/existGrpMa.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: slave, master
!
!      OPERATEURS :     AFFE_CHAR_MECA ET AFFE_CHAR_MECA_C
!                                      ET AFFE_CHAR_MECA_F
!                       DEFI_CONTACT
!
!     VERIFICATION DES NORMALES AUX MAILLES SURFACIQUES EN 3D
!     ET LINEIQUES EN 2D
!     V1 : ON VERIFIE QUE LES NORMALES SONT HOMOGENES
!     V2 : ON VERIFIE QUE LES NORMALES SONT SORTANTES
!
!-----------------------------------------------------------------------
    integer(kind=8), parameter :: nbobj = 2
    integer(kind=8) :: ier
    integer(kind=8) :: ndim, ndim1, vali
    integer(kind=8) :: iobj, iCell, nbCell, nconex
    integer(kind=8) :: cellNume, cellTypeNume, nbmapr, nbmabo, ntrait
    integer(kind=8) ::  jmab, norien, jlima, nbmamo
    aster_logical, parameter :: reorie = ASTER_FALSE
    aster_logical :: l_exi, l_exi_p
    character(len=8) ::  cellTypeName, mesh
    character(len=24) :: grmama, nogr
    character(len=24) :: valk(2)
    character(len=24) :: objet(nbobj)
    integer(kind=8), pointer :: listCellNume(:) => null(), typmail(:) => null()
!
!     INITIALISATIONS
    ier = 0
!
    ndim = 0
    mesh = ' '
    call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
    grmama = mesh//'.GROUPEMA       '
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)
!
! ---     RECUPERATION DE LA DIMENSION DU PROBLEME
!
    objet(1:nbobj) = [slave, master]
    do iobj = 1, nbobj
        nogr = objet(iobj)
!
! ---   RECUPERATION DU NOMBRE DE MAILLES DU GROUP_MA :
! ---------------------------------------------
        ! check si GROUP_MA existe
        call existGrpMa(mesh, nogr, l_exi, l_exi_p)

        if ((.not. l_exi) .and. (l_exi_p)) then
            goto 211
        end if

        call jelira(jexnom(grmama, nogr), 'LONUTI', nbCell)
        call jeveuo(jexnom(grmama, nogr), 'L', vi=listCellNume)
!
        do iCell = 1, nbCell
            cellNume = listCellNume(iCell)
            cellTypeNume = typmail(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)

            if (cellTypeName(1:3) .eq. 'POI') then
                goto 211
            else if (cellTypeName(1:3) .eq. 'SEG') then
                ndim1 = 2
                if (ndim .ne. ndim1) then
                    goto 211
                end if
!
            end if
        end do
!
! ---   FIN DE BOUCLE SUR LES MAILLES DU GROUP_MA
!
        norien = 0
!
        if (nbCell .gt. 0) then
!
            call wkvect('&&CHCKNO.MAILLE_BORD', 'V V I', nbCell, jmab)
            call chbord(model, nbCell, listCellNume, zi(jmab), nbmapr, nbmabo)
            if (nbmapr .eq. nbCell .and. nbmabo .eq. 0) then
                call ornorm(mesh, listCellNume, nbCell, reorie, norien, nconex)
                if (nconex .gt. 1) then
                    call utmess('F', 'MESH3_99')
                end if
            else
                nbmamo = 0
                jlima = 1
                call orilma(mesh, ndim, listCellNume, nbCell, norien, &
                            ntrait, reorie, nbmamo, zi(jlima))
                if ((ntrait .ne. 0)) then
                    call utmess('A', 'CONTACT2_20')
                end if
            end if
            call jedetr('&&CHCKNO.MAILLE_BORD')

            if (norien .ne. 0) then
                ier = ier+1
                valk(1) = nogr
                vali = norien
                call utmess('E', 'MODELISA8_56', sk=valk(1), si=vali)
            end if
        end if
211     continue
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA4_24')
    end if
!
end subroutine
