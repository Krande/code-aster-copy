! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine mdidisvisc(sd_nl_, nbnoli, nomres, nbsauv, temps)
!
    implicit none
    character(len=*) :: sd_nl_
    integer          :: nbnoli
    character(len=8) :: nomres
    integer          :: nbsauv
    real(kind=8)     :: temps(*)
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/nlget.h"
#include "asterfort/ulexis.h"
#include "asterfort/ulopen.h"
!
! --------------------------------------------------------------------------------------------------
!
!                       IMPRESSION DES RESULTATS SUR LES DIS_VISC
!
! --------------------------------------------------------------------------------------------------
!
! IN
!   sd_nl_  : nom de structure de données pour les calculs non linéaires
!   nomres  : nom du concept résultat
!   nbnoli  : nombre de liaison non-linéaire
!   nbsauv  : nombre de pas sauvegardé
!   temps   : instant de sauvegarde
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ific, nbocc, iocc, iret, ii
    integer :: nbvint, nbvdisc, it, indx, jvint
    integer :: vv, nltype_i, start, finish, jvindx, jdesc
    character(len=8) :: noeud1, noeud2, sd_nl
!
    call jemarq()
    sd_nl = sd_nl_

    call getfac('IMPRESSION', nbocc)
    if (nbocc .eq. 0) then
        goto 999
    end if

!   Récupération des variables internes
!   Longueur maximum d'un bloc de variables internes
    call jeveuo(nomres//'           .DESC', 'L', jdesc)
    nbvint = zi(jdesc-1+4)
    if (nbvint .eq. 0) goto 999

    call jeveuo(nomres//'        .NL.VINT', 'L', jvint)
    call jeveuo(nomres//'        .NL.VIND', 'L', jvindx)
!
    do iocc = 1, nbocc
        call getvis('IMPRESSION', 'UNITE_DIS_VISC', iocc=iocc, scal=ific, nbret=iret)
        if (iret .ne. 1) cycle
!       Impression des informations sur les DIS_VISC
        if (.not. ulexis(ific)) then
            call ulopen(ific, ' ', ' ', 'NEW', 'O')
        end if
!
        do ii = 1, nbnoli
            call nlget(sd_nl, _NL_TYPE, iocc=ii, iscal=nltype_i)
            if (nltype_i .ne. NL_DIS_VISC) cycle

            start = zi(jvindx-1+ii)+7
            finish = zi(jvindx-1+ii+1)

            nbvdisc = finish-start
            ASSERT(nbvdisc .eq. 4)

!           Noeuds du discret
            call nlget(sd_nl, _NO1_NAME, iocc=ii, kscal=noeud1)
            call nlget(sd_nl, _NO2_NAME, iocc=ii, kscal=noeud2)

!           Impressions des variables internes dans l'ordre de stockage
            write (ific, 100) '#'
            write (ific, 100) '#--------------------------------------------------'
            write (ific, 100) '#RESULTAT '//nomres
            write (ific, 101) '#DIS_VISC ', ii, ' '//noeud1//' '//noeud2
            write (ific, 102) 'INST', 'FORCE', 'DEPLVISC', 'DEPL', 'PUISS'
            do it = 1, nbsauv
                indx = jvint-1+(it-1)*nbvint+start
                write (ific, 103) temps(it), (zr(indx+vv-1), vv=1, nbvdisc)
            end do
        end do
!       On ferme le fichier pour être sûr que le flush soit fait
        call ulopen(-ific, ' ', ' ', ' ', ' ')
    end do
100 format(A)
101 format(A, I5, A)
102 format(5(1X, A18))
103 format(5(1X, 1pE18.10E3))
!
999 continue
    call jedema()
end subroutine
