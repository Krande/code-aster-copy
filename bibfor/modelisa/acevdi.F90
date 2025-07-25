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

subroutine acevdi(nbocc, nomaz, nomoz, mcf, nlm, &
                  nlg, nln, nlj, ier)
!
!
    implicit none
    integer(kind=8) :: nbocc, nlm, nlg, nln, nlj, ier
    character(len=*) :: nomaz, nomoz, mcf
!
! --------------------------------------------------------------------------------------------------
!
!        AFFE_CARA_ELEM
!           VERIFICATION DES MOTS CLES POUR L'ELEMENT DISCRET
!
! --------------------------------------------------------------------------------------------------
!
! IN
!     NBOCC    :  NOMBRE D'OCCURENCE
!     NOMAZ    :  NOM DU MAILLAGE
!     NOMOZ    :  NOM DU MODELE
!     MCF      :  MOT CLEF
! OUT
!     NLM      :  NOMBRE TOTAL DE MAILLE
!     NLG      :  NOMBRE TOTAL DE GROUPE DE MAILLE
!     NLN      :  NOMBRE TOTAL DE NOEUD
!     NLJ      :  NOMBRE TOTAL DE GROUP_NO
!     IER      :  ERREUR
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/acevd2.h"
#include "asterfort/assert.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/verdis.h"
#include "asterfort/verima.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
! --------------------------------------------------------------------------------------------------
    integer(kind=8) ::  i3d, i2d, ndim1, ioc, nc, ng, nm, nsom, nbmail, un
    integer(kind=8) :: n1, ima, nbgrm, ig, jmail, numa, nutyma, lmax2
    character(len=4) :: type
    character(len=8) :: k8b, nomu, noma, nomo, nomail, typel
    character(len=16) :: concep, cmd
    character(len=24) :: grmama, cara, nogrm
    character(len=24) :: valk(4)
    aster_logical :: l_parallel_mesh
! --------------------------------------------------------------------------------------------------
    integer(kind=8), pointer :: typmail(:) => null()
    character(len=24), pointer :: group_ma(:) => null()
! --------------------------------------------------------------------------------------------------
    call getres(nomu, concep, cmd)
!
    noma = nomaz
    nomo = nomoz
    nlm = 0
    nlg = 0
    nln = 0
    nlj = 0
    un = 1
    grmama = noma//'.GROUPEMA'
    l_parallel_mesh = isParallelMesh(noma)
!
!   Vecteur du type des mailles du maillage :
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
!   Vérification des dimensions / modélisations
    call verdis(nomo, noma, 'E', i3d, i2d, ndim1, ier)
    ASSERT((mcf .eq. 'DISCRET_2D') .or. (mcf .eq. 'DISCRET'))
!
!   boucle sur les occurences
    do ioc = 1, nbocc
        call getvtx(mcf, 'CARA', iocc=ioc, scal=cara, nbret=nc)
!
        call getvem(noma, 'GROUP_MA', mcf, 'GROUP_MA', ioc, 0, k8b, ng)
        call getvem(noma, 'MAILLE', mcf, 'MAILLE', ioc, 0, k8b, nm)
!
        nsom = ng+nm
        if ((nsom .eq. ng) .or. (nsom .eq. nm)) then
            nlm = max(nlm, -nm)
            nlg = max(nlg, -ng)
        end if
!
!       Vérification du bon type de maille en fonction de cara
        if ((cara(2:7) .eq. '_T_D_N') .or. (cara(2:8) .eq. '_TR_D_N') .or. &
            (cara(2:5) .eq. '_T_N') .or. (cara(2:6) .eq. '_TR_N')) then
            type = 'POI1'
        else
            type = 'SEG2'
        end if
!
        if (nm .ne. 0) then
            if (l_parallel_mesh) then
                call utmess('F', 'MODELISA7_86')
            end if
            nbmail = -nm
            call wkvect('&&ACEVDI.MAILLE', 'V V K8', nbmail, jmail)
            call getvtx(mcf, 'MAILLE', iocc=ioc, nbval=nbmail, vect=zk8(jmail), nbret=n1)
            do ima = 1, nbmail
                nomail = zk8(jmail+ima-1)
                call verima(noma, nomail, un, 'MAILLE')
                numa = char8_to_int(nomail)
                nutyma = typmail(numa)
                call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
                if (typel(1:4) .ne. type) then
                    valk(1) = nomail
                    valk(2) = type
                    valk(3) = typel
                    valk(4) = cara
                    call utmess('F', 'MODELISA_56', nk=4, valk=valk)
                end if
            end do
            call jedetr('&&ACEVDI.MAILLE')
        end if
!
        if (ng .ne. 0) then
            nbgrm = -ng
            AS_ALLOCATE(vk24=group_ma, size=nbgrm)
            call getvtx(mcf, 'GROUP_MA', iocc=ioc, nbval=nbgrm, vect=group_ma, nbret=n1)
            do ig = 1, nbgrm
                nogrm = group_ma(ig)
                un = 1
                call verima(noma, nogrm, un, 'GROUP_MA')
                if (un .eq. 0 .and. l_parallel_mesh) cycle
                call jelira(jexnom(grmama, nogrm), 'LONUTI', nbmail)
                call jeveuo(jexnom(grmama, nogrm), 'L', jmail)
                do ima = 1, nbmail
                    numa = zi(jmail+ima-1)
                    nutyma = typmail(numa)
                    call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
                    if (typel(1:4) .ne. type) then
                        nomail = int_to_char8(numa)
                        valk(1) = nomail
                        valk(2) = type
                        valk(3) = typel
                        valk(4) = cara
                        call utmess('F', 'MODELISA_56', nk=4, valk=valk)
                    end if
                end do
            end do
            AS_DEALLOCATE(vk24=group_ma)
        end if
!
    end do
!
    lmax2 = max(1, nlm, nlg, nln, nlj)
    call acevd2(noma, nomo, mcf, lmax2, nbocc)
end subroutine
