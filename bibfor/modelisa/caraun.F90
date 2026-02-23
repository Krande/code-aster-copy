! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine caraun(sdcont, zoneKeyword, nbUnilZone, &
                  nbgdcuJv, coefcuJv, &
                  compcuJv, multcuJv, penacuJv, ntCmp)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/cazouu.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "Contact_type.h"
!
    character(len=8), intent(in) :: sdcont
    character(len=16), intent(in) :: zoneKeyword
    integer(kind=8), intent(in) :: nbUnilZone
    character(len=24), intent(in) :: nbgdcuJv, coefcuJv, compcuJv, multcuJv, penacuJv
    integer(kind=8), intent(out) :: ntCmp
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_CONTACT
!
! Read informations for LIAISON_UNILATER in command
!
! --------------------------------------------------------------------------------------------------
!
! In  sdcont           : name of contact concept (DEFI_CONTACT)
! In  nbUnilZone       : number of zones with LIAISON_UNILATER
! IN  NBGDCU : NOM JEVEUX DE LA SD INFOS POINTEURS GRANDEURS
!       ZI(JNBGD+IOCC-1): INDICE DEBUT DANS LISTE DES NOMS DES GRANDEURS
!                       POUR ZONE IOCC
!       ZI(JNBGD+IOCC) - ZI(JNBGD+IOCC-1): NOMBRE DE GRANDEURS DE LA
!                       ZONE IOCC
! IN  COEFCU : NOM JEVEUX DE LA SD CONTENANT LES COEFFICIENTS DES
!              GRANDEURS DE MEMBRE DE DROITE
!              VECTEUR TYPE ZR OU ZK8 SUIVANT FONREE
!       Z*(JCOEF+IOCC-1): VALEUR OU NOM FONCTION DU MEMBRE DE DROITE
! IN  COMPCU : NOM JEVEUX DE LA SD CONTENANT LES GRANDEURS DU MEMBRE
!              DE GAUCHE
!              LONGUEUR = ZI(JDUME+3)
!              INDEXE PAR NBGDCU:
!       ZI(JNBGD+IOCC-1): INDEX DEBUT POUR ZONE IOCC
!       ZI(JDUME+2*(IOCC-1)+5) = ZI(JNBGD+IOCC)-ZI(JNBGD+IOCC-1):
!                         NOMBRE GRANDEURS A GAUCHE POUR ZONE IOCC
!       ZK8(JCMPG-1+INDEX+ICMP-1): NOM ICMP-EME GRANDEUR
! IN  MULTCU : NOM JEVEUX DE LA SD CONTENANT LES COEF DU MEMBRE
!              DE GAUCHE
!              VECTEUR TYPE ZR OU ZK8 SUIVANT FONREE
!              MEME ACCES QUE COMPCU
! IN  PENACU : NOM JEVEUX DE LA SD CONTENANT LES COEF DE PENALITE
!              VECTEUR TYPE ZR, MEME ACCES QUE COEFCU
! OUT NTCMP  : NOMBRE TOTAL DE COMPOSANTES SUR TOUTES LES ZONES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: contForm = CONT_FORM_UNIL
    integer(kind=8), parameter :: ntCmpMaxi = 30, iterContMaxi = 10
    character(len=8) :: cmpName(ntCmpMaxi), k8bid, coefImpo, coefMult(ntCmpMaxi)
    integer(kind=8) :: noc, nbCmp, nbCoefMult
    integer(kind=8) :: iUnilZone, icmp
    character(len=16) :: algoCont
    real(kind=8) :: coefPena
    character(len=24) :: sdcontParaciJv
    integer(kind=8), pointer :: sdcontParaci(:) => null()
    real(kind=8), pointer :: penacu(:) => null()
    character(len=8), pointer :: coefcu(:) => null()
    integer(kind=8), pointer :: nbgdcu(:) => null()
    character(len=8), pointer :: compcu(:) => null()
    character(len=8), pointer :: multcu(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    ntCmp = 0

! - Datastructure for contact definition
    sdcontParaciJv = sdcont(1:8)//'.PARACI'
    call jeveuo(sdcontParaciJv, 'E', vi=sdcontParaci)

! - Set the formulation
    sdcontParaci(4) = contForm

! - Get algorithm
    call cazouu(zoneKeyword, nbUnilZone, 'ALGO_CONT', 'T')
    call getvtx(zoneKeyword, 'ALGO_CONT', iocc=1, scal=algoCont)
    if (algoCont .eq. 'CONTRAINTE') then
        sdcontParaci(30) = 1
    else if (algoCont .eq. 'PENALISATION') then
        sdcontParaci(30) = 4
        call wkvect(penacuJv, 'V V R', nbUnilZone, vr=penacu)
        do iUnilZone = 1, nbUnilZone
            call getvr8(zoneKeyword, 'COEF_PENA', iocc=iUnilZone, scal=coefPena)
            penacu(iUnilZone) = coefPena
        end do
    else
        ASSERT(ASTER_FALSE)
    end if

! - Parameter ITER_CONT_MAXI
    sdcontParaci(3) = iterContMaxi

! - Right-hand side
    call wkvect(coefcuJv, 'V V K8', nbUnilZone, vk8=coefcu)
    do iUnilZone = 1, nbUnilZone
        call getvid(zoneKeyword, 'COEF_IMPO', iocc=iUnilZone, scal=coefImpo, nbret=noc)
        coefcu(iUnilZone) = coefImpo
    end do

! - Left-hand side - Count DOF
    call wkvect(nbgdcuJv, 'V V I', nbUnilZone+1, vi=nbgdcu)
    nbgdcu(1) = 1
    ntCmp = 0
    do iUnilZone = 1, nbUnilZone
        call getvtx(zoneKeyword, 'NOM_CMP', iocc=iUnilZone, scal=k8bid, nbret=nbCmp)
        call getvid(zoneKeyword, 'COEF_MULT', iocc=iUnilZone, scal=k8bid, nbret=nbCoefMult)
        if (nbCmp .ne. nbCoefMult) then
            call utmess('F', 'UNILATER_42')
        end if
        nbCmp = abs(nbCmp)
        ntCmp = ntCmp+nbCmp
        if (ntCmp .gt. ntCmpMaxi) then
            call utmess('F', 'UNILATER_43')
        end if
        nbgdcu(iUnilZone+1) = nbgdcu(iUnilZone)+nbCmp
    end do

! - Left-hand side - List of parameters
    call wkvect(compcuJv, 'V V K8', ntCmp, vk8=compcu)
    call wkvect(multcuJv, 'V V K8', ntCmp, vk8=multcu)
    do iUnilZone = 1, nbUnilZone
        nbCmp = nbgdcu(iUnilZone+1)-nbgdcu(iUnilZone)
        call getvtx(zoneKeyword, 'NOM_CMP', iocc=iUnilZone, nbval=nbCmp, vect=cmpName, &
                    nbret=noc)
        call getvid(zoneKeyword, 'COEF_MULT', iocc=iUnilZone, nbval=nbCmp, vect=coefMult, &
                    nbret=noc)
        do icmp = 1, nbCmp
            compcu(nbgdcu(iUnilZone)+icmp-1) = cmpName(icmp)
            multcu(nbgdcu(iUnilZone)+icmp-1) = coefMult(icmp)
        end do
    end do
!
    call jedema()
!
end subroutine
