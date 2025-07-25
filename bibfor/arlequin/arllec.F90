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

subroutine arllec(motcle, iocc, modele, noma, nomb, &
                  model, cine, dime)

    implicit none

#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/infniv.h"
#include "asterfort/getvtx.h"
#include "asterfort/wkvect.h"
#include "asterfort/arlver.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedema.h"

    character(len=16) :: motcle
    integer(kind=8) ::      iocc
    character(len=8) ::  modele
    character(len=10) :: noma, nomb
    character(len=8) ::  model(3), cine(3)
    integer(kind=8) ::      dime

! ----------------------------------------------------------------------

! ROUTINE ARLEQUIN

! LECTURE ET VERIFICATION DES MAILLES DES MODELES

! ----------------------------------------------------------------------

! IN  MOTCLE : MOT-CLEF FACTEUR POUR ARLEQUIN
! IN  IOCC   : OCCURRENCE DU MOT CLEF-FACTEUR ARLEQUIN
! IN  MODELE : NOM DU MODELE
! I/O NOMA   : NOM DE LA SD POUR STOCKAGE MAILLES GROUP_MA_1
! I/O NOMB   : NOM DE LA SD POUR STOCKAGE MAILLES GROUP_MA_2
! OUT MODEL  : MODELISATION ASSOCIEE AUX MAILLES '3D',
!              UN POUR CHAQUE GROUPE + ZONE DE COLLAGE
! OUT CINE   : CINEMATIQUE ASSOCIEE AUX MAILLES
!              'SOLIDE' OU 'POUTRE'
!              UN POUR CHAQUE GROUPE + ZONE DE COLLAGE
! OUT DIME   : DIMENSION DE L'ESPACE GLOBAL 2 OU 3

    integer(kind=8) ::      nbev1, nbev2
    character(len=8) ::  k8bid
    character(len=16) ::  option
    integer(kind=8) ::      jgrm1, jgrm2
    integer(kind=8) ::      ifm, niv, iop
    character(len=6) :: nompro
    parameter(nompro='ARLLEC')

! ----------------------------------------------------------------------

    call jemarq()
    call infniv(ifm, niv)

! --- LECTURE MAILLE GROUPE_MA_1

    call getvtx(motcle, 'OPTION', iocc=iocc, scal=option, nbret=iop)
    if (option .eq. '3D_POU_ARLEQUIN') then
        call getvtx(motcle, 'GROUP_MA_1', iocc=iocc, nbval=0, &
                    nbret=nbev1, scal=k8bid)
        call wkvect('&&'//nompro//'.GMA1', 'V V K8', -nbev1, jgrm1)
        call getvtx(motcle, 'GROUP_MA_1', iocc=iocc, nbval=-nbev1, &
                    vect=zk8(jgrm1), nbret=nbev1)
        call arlver(modele, zk8(jgrm1), nbev1, noma, model(1), cine(1))

! --- LECTURE MAILLE GROUPE_MA_2

        call getvtx(motcle, 'GROUP_MA_2', iocc, 0, k8bid, nbret=nbev2)
        call wkvect('&&'//nompro//'.GMA2', 'V V K8', -nbev2, jgrm2)
        call getvtx(motcle, 'GROUP_MA_2', iocc, -nbev2, zk8(jgrm2), &
                    nbret=nbev2)
        call arlver(modele, zk8(jgrm2), nbev2, nomb, model(2), cine(2))
    end if

! --- DIMENSION DE L'ESPACE GLOBAL

    dime = 3

! --- MENAGE

    call jedetr('&&'//nompro//'.GMA1')
    call jedetr('&&'//nompro//'.GMA2')

    call jedema()
end subroutine
