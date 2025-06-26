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

subroutine arlini(mail, base, dime, nbno, &
                  nbma, nbnoma)

    implicit none

#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/assert.h"
#include "asterfort/jexnom.h"
#include "asterfort/wkvect.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"

!     ARGUMENTS:
!     ----------
    character(len=8) :: mail
    character(len=1) :: base
    integer(kind=8) :: dime
    integer(kind=8) :: nbma, nbno, nbnoma

! ----------------------------------------------------------------------

! UTILITAIRE - SD MAILLAGE
! PREPARE LES OBJETS DE BASE DE LA SD MAILLAGE

! ----------------------------------------------------------------------

! IN  MAIL   : NOM DU MAILLAGE
! IN  BASE   : BASE DE CREATION 'G' OU 'V'
! IN  DIME   : DIMENSION DE L'ESPACE (2 OU 3)
! IN  NBNO   : NOMBRE DE NOEUDS DU MAILLAGE
! IN  NBMA   : NOMBRE DE MAILLES DU MAILLAGE
! IN  NBNOMA : LONGUEUR DU VECTEUR CONNECTIVITE DES MAILLES

    character(len=24) :: maidim, cooval, coodsc
    integer(kind=8) :: jdime, jcoor, jcods
    character(len=24) :: typmai
    integer(kind=8) :: jtypm
    character(len=24) :: connex
    integer(kind=8) :: ntgeo

! ----------------------------------------------------------------------

    call jemarq()

! --- NOM DES OBJETS JEVEUX

    maidim = mail(1:8)//'.DIME'
    cooval = mail(1:8)//'.COORDO    .VALE'
    coodsc = mail(1:8)//'.COORDO    .DESC'
    typmai = mail(1:8)//'.TYPMAIL        '
    connex = mail(1:8)//'.CONNEX         '

! --- VERIFICATIONS

    if ((dime < 2) .or. (dime > 3)) then
        ASSERT(.false.)
    end if

    if ((nbma <= 0) .or. (nbno <= 0) .or. (nbnoma <= 0)) then
        ASSERT(.false.)
    end if

! --- RECUPERATION DU NUMERO IDENTIFIANT LE TYPE DE CHAM_NO GEOMETRIE

    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)

! --- DIMENSIONS

    call wkvect(maidim, base//' V I', 6, jdime)
    zi(jdime) = nbno
    zi(jdime+1) = 0
    zi(jdime+2) = nbma
    zi(jdime+3) = 0
    zi(jdime+4) = 0
    zi(jdime+5) = dime

! --- CHAM_NO DES COORDONNEES DES NOEUDS

    call wkvect(cooval, base//' V R', 3*nbno, jcoor)

    call jecreo(coodsc, base//' V I')
    call jeecra(coodsc, 'LONMAX', 3, ' ')
    call jeecra(coodsc, 'DOCU', 0, 'CHGO')
    call jeveuo(coodsc, 'E', jcods)
    zi(jcods) = ntgeo
    zi(jcods+1) = -3
    zi(jcods+2) = 14

! --- TYPE DES MAILLES

    call wkvect(typmai, base//' V I', nbma, jtypm)

! --- CONNECTIVITES DES MAILLES

    call jecrec(connex, base//' V I', 'NU', 'CONTIG', 'VARIABLE', nbma)
    call jeecra(connex, 'LONT', nbnoma, ' ')

    call jedema()

end subroutine
