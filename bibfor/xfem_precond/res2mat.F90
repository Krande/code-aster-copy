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

subroutine res2mat(resu, inst, chmat, nommat, mu, ka, lvarc, varcns, cplan)
!
! person_in_charge: samuel.geniaut at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jelira.h"
#include "asterfort/jeexin.h"
#include "asterfort/codent.h"
#include "asterfort/vrcins.h"
#include "asterfort/celces.h"
#include "asterfort/detrsd.h"
!
    character(len=*) :: resu, chmat
    real(kind=8) :: inst
    real(kind=8), optional :: mu, ka
    aster_logical, optional :: lvarc
    aster_logical, optional :: cplan
    character(len=19), optional :: varcns
    character(len=8), optional :: nommat
!*****************************************************************
! BUT :: ROUTINE GENERIQUE POUR L EXTRACTION DU CHAMP MATERIAU
!           DANS UN CONCEPT RESULTAT
!   IN:
!      * RESU  : CHAMP RESULTAT ASTER
!   OUT:
!      * CHMAT : CHAMP MATERIAU ASTER
!      * MU    : PARAMETRE MATERIAU ELASTIQUE (OPT.)
!      * KA    : PARAMETRE MATERIAU ELASTIQUE (OPT.)
!*****************************************************************
    character(len=2) :: codret
    aster_logical :: cplan2
    character(len=6) :: k6
    character(len=8) :: nommatz, model, is_varc, carael
    character(len=16) :: kmodl
    character(len=19) :: carmat, varcel
    integer(kind=8) :: jvalm, ianorc, nbrc, irc, jvalk, jvale, ncmpa, j, ndim, nmat, imat
    integer(kind=8) :: iret
    character(len=32) :: nomrc
    real(kind=8) :: nu, ka2, e, mu2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    call jemarq()
!
    if (present(mu)) mu = 0.d0
    if (present(ka)) ka = 0.d0
    if (present(lvarc)) lvarc = .false.
    chmat = ' '
    cplan2 = .false.
!
    call dismoi('CHAM_MATER', resu(1:8), 'RESULTAT', repk=chmat)
    call dismoi('NOM_MODELE', resu, 'RESULTAT', repk=model)
    call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
    call dismoi('PHENOMENE', model, 'MODELE', repk=kmodl)
!   ON NE COMPREND POURQUOI ON DOIT LIRE LE CHAMP MATER :: SORTIE
    if (kmodl .eq. 'THERMIQUE') goto 99
    if (ndim .lt. 2 .or. chmat .eq. ' ') goto 99
    if (ndim .eq. 2) then
        call dismoi('MODELISATION', model, 'MODELE', repk=kmodl)
        if (kmodl .eq. 'C_PLAN') cplan2 = .true.
    else
        cplan2 = .false.
    end if
!
    carmat = chmat(1:8)//'.CHAMP_MAT'
    call jeveuo(carmat(1:19)//'.VALE', 'L', jvalm)
    call jelira(carmat(1:19)//'.VALE', 'LONMAX', nmat)
    k6 = ' '
    do imat = 1, nmat
        nommatz = zk8(jvalm-1+imat)
        if (nommatz .eq. ' ') goto 10
        call jeveuo(nommatz//'.MATERIAU.NOMRC', 'L', ianorc)
        call jelira(nommatz//'.MATERIAU.NOMRC', 'LONUTI', nbrc)
        do irc = 1, nbrc
            nomrc = zk32(ianorc-1+irc)
            if (nomrc .eq. 'ELAS') then
                call codent(irc, 'D0', k6)
                exit
            end if
        end do
        if (k6 .ne. ' ') exit
10      continue
    end do
!   SI ON NE TROUVE PAS DE MATER ELAS :: SORTIE
!   ATTENTION : L EXTRACTION DES PARAMETRES ELESTIQUES EST DOUTEUSE
!               EN PRESENCE DE VARIABLES DE COMMANDE MATERIAU
!               IL CONVIENT DE CHOISIR LA BONNE VARIABLE DE COMMANDE SINON
!               RISQUE DE RESU FAUX !!!
    ASSERT(k6 .ne. ' ')
    call jeveuo(nommatz//'.CPT.'//k6//'.VALK', 'L', jvalk)
    call jelira(nommatz//'.CPT.'//k6//'.VALK', 'LONMAX', ncmpa)
    call jeveuo(nommatz//'.CPT.'//k6//'.VALR', 'L', jvale)
    nu = -1.d0
    e = -1.d0
    ka2 = 0.
    mu2 = 0.
    do j = 1, ncmpa
        if (zk16(jvalk-1+j) .eq. 'NU') then
            nu = zr(jvale-1+j)
        end if
        if (zk16(jvalk-1+j) .eq. 'E') then
            e = zr(jvale-1+j)
        end if
    end do
    if (e .gt. 0 .and. nu .ge. 0) goto 50
!   SI LES MATERIAUX ELASTIQUES NE SONT PAS TROUVES ON PASSE PAR LES VARIABLES
!     DE COMMANDE
    call dismoi('EXI_VARC', chmat, 'CHAM_MATER', repk=is_varc)
    if (is_varc .eq. 'OUI') then
        ASSERT(present(lvarc))
        lvarc = .true.
        call jeexin(chmat(1:8)//'.CVRCVARC', iret)
        ASSERT(iret .gt. 0)
        call dismoi('CARA_ELEM', resu, 'RESULTAT', repk=carael)
        ASSERT(carael .ne. '#PLUSIEURS')
        if (carael .eq. '#AUCUN') then
            carael = ' '
        end if
        ASSERT(present(varcns))
        varcel = '&&RES2MAT.VARCEL'
        call vrcins(model, chmat, carael, inst, varcel, codret, nompaz='PVARCNO')
        ASSERT(codret .eq. 'OK')
        call exisd('CHAM_ELEM', varcel, iret)
        ASSERT(iret .eq. 1)
        call celces(varcel, 'V', varcns)
        call detrsd('CHAM_ELEM_S', varcel)
        goto 99
    end if
!
!    print*,' *** KOR : E,NU=',e,nu
50  continue
    ka2 = 3.d0-4.d0*nu
    mu2 = e/(2.d0*(1.d0+nu))
    if (cplan2) ka2 = (3.d0-nu)/(1.d0+nu)
!
99  continue
!
!   ECRITURES ET SORTIE
    if (present(mu)) mu = mu2
    if (present(ka)) ka = ka2
    if (present(cplan)) cplan = cplan2
    if (present(nommat)) nommat = nommatz
!
    call jedema()
!
end subroutine
