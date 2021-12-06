! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterc/r8prem.h"
#include "asterfort/chbord.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/orilma.h"
#include "asterfort/ornorm.h"
#include "asterfort/utmamo.h"
#include "asterfort/utmess.h"
#include "asterfort/utmotp.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
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
    integer, parameter :: nbobj = 2
    integer :: ier
    integer :: ndim, ndim1, vali
    integer :: iobj, ima, nbmail
    integer :: numa, idtyma, nutyma, nbmapr, nbmabo, ntrait
    integer :: jgro, jmab, norien,  jlima, nbmamo
    aster_logical, parameter :: reorie = ASTER_FALSE
    character(len=8) ::  typel, mesh
    character(len=24) :: grmama, mailma, nogr
    character(len=24) :: valk(2)
    character(len=24) :: objet(nbobj)
!
!     INITIALISATIONS
    ier    = 0
!
    ndim = 0
    mesh = ' '
    call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
!
    grmama = mesh//'.GROUPEMA'
    mailma = mesh//'.NOMMAI'
    call jeveuo(mesh//'.TYPMAIL', 'L', idtyma)
!
! ---     RECUPERATION DE LA DIMENSION DU PROBLEME
!
    objet(1:nbobj) = [slave, master]
    do iobj = 1, nbobj
        nogr = objet(iobj)
!
! ---   RECUPERATION DU NOMBRE DE MAILLES DU GROUP_MA :
! ---------------------------------------------
        call jelira(jexnom(grmama, nogr), 'LONUTI', nbmail)
        call jeveuo(jexnom(grmama, nogr), 'L', jgro)
!
        do ima = 1, nbmail
            numa = zi(jgro-1+ima)
            nutyma = zi(idtyma+numa-1)
!
! ---  TYPE DE LA MAILLE :
! -----------------
            call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
!
! ---  CAS D'UNE MAILLE POINT
! ----------------------
            if (typel(1:3) .eq. 'POI') then
                goto 211
!
! ---   CAS D'UNE MAILLE SEG
! --------------------
            else if (typel(1:3) .eq. 'SEG') then
                ndim1 = 2
                if (ndim .ne. ndim1) then
                    goto 211
                endif
!
            endif
        end do
!
! ---   FIN DE BOUCLE SUR LES MAILLES DU GROUP_MA
!
        norien = 0
!
        if (nbmail.gt.0) then
!
            call wkvect('&&CHCKNO.MAILLE_BORD', 'V V I', nbmail, jmab)
            call chbord(model, nbmail, zi(jgro), zi( jmab), nbmapr, nbmabo)
            if (nbmapr .eq. nbmail .and. nbmabo .eq. 0) then
                call ornorm(mesh, zi(jgro), nbmail, reorie, norien)
            else
                nbmamo = 0
                jlima = 1
                call orilma(mesh, ndim, zi(jgro), nbmail, norien,&
                            ntrait, reorie, nbmamo, zi(jlima ))
                if ((ntrait .ne. 0)) then
                    call utmess('A', 'CONTACT2_20')
                endif
            endif
            call jedetr('&&CHCKNO.MAILLE_BORD')

            if (norien .ne. 0) then
                ier = ier + 1
                valk(1) = nogr
                vali = norien
                call utmess('E', 'MODELISA8_56', sk=valk(1), si=vali)
            endif
        end if
211 continue
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'MODELISA4_24')
    endif
!
end subroutine
