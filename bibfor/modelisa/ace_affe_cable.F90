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

subroutine ace_affe_cable(nbocc, infoconcept, infocarte, grplmax, grpnbma, lesmailles)
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ace_affe_verif_elem.h"
#include "asterfort/getvtx.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nocart.h"
!
    integer(kind=8)         :: nbocc
    type(cara_elem_info)    :: infoconcept
    type(cara_elem_carte)   :: infocarte(*)
    character(len=24)       :: grplmax(*)
    integer(kind=8)         :: grpnbma(*)
    integer(kind=8)         :: lesmailles(*)
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM : AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT CABLE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)   :: GroupeMaxOccur
    character(len=8)    :: nomu, noma
    character(len=19)   :: cartelem, cartefcx
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)   :: ii, jj, ioc, jdcelems, jdccfcx, jdvelems, jdvcfcx, jdme
    integer(kind=8)   :: ng, ibid, nbcara1, coderet
    character(len=8)  :: k8_fcx
!
    integer(kind=8)   :: TForme, TVaria, TCodeFV
!
    integer(kind=8), parameter  :: nbMaxCar = 5
    integer(kind=8)             :: ncara, nvale, nbCara
    real(kind=8)                :: vale(nbMaxCar)
    character(len=8)            :: cara(nbMaxCar)
    character(len=8)            :: NomCara(nbMaxCar)
    character(len=8)            :: NomCata(nbMaxCar)
!
!-----------------------------------------------------------------------
!
    if (nbocc .eq. 0) goto 999
!
    call jemarq()
!   Récupère les informations
    nomu = infoconcept%nomu
    noma = infoconcept%maillage
    jdme = infoconcept%jmodmail
    GroupeMaxOccur = infoconcept%GroupeMaxOccur
!   La carte est déjà alouée : ace_crea_carte
    cartelem = infocarte(ACE_CAR_CABLE)%nom_carte
    jdcelems = infocarte(ACE_CAR_CABLE)%adr_cmp
    jdvelems = infocarte(ACE_CAR_CABLE)%adr_val
!   La carte est déjà alouée : ace_crea_carte
    cartefcx = infocarte(ACE_CAR_CVCXF)%nom_carte
    jdccfcx = infocarte(ACE_CAR_CVCXF)%adr_cmp
    jdvcfcx = infocarte(ACE_CAR_CVCXF)%adr_val
!
    nbCara = 2
    NomCara(1:nbCara) = [character(len=6)::'AIRE', 'N_INIT']
    NomCata(1:nbCara) = [character(len=4)::'SECT', 'TENS']
!
    TForme = ACE_SECTION_CERCLE
    TVaria = ACE_SECTION_CONSTANTE
    TCodeFV = ACE_CodeFormeVaria(TForme, TVaria)
!
!   Lecture des valeurs et affectation dans la carte
    bioc: do ioc = 1, nbocc
!       Les GROUP_MA concernés
        call getvtx('CABLE', 'GROUP_MA', ioc, nbval=GroupeMaxOccur, vect=grplmax, nbret=ng)
!       Vérification que les affectations ne concernent que les CABLE
        call ace_affe_verif_elem(noma, jdme, lesmailles, ng, grplmax, grpnbma, &
                                 ACE_NU_CABLE, 'CABLE', coderet)
!       En // ou pas : si tous les groupes de mailles sont vides ==> on cycle
        if (coderet .eq. 0) cycle bioc
!       Les cartes vont être remplies
        infocarte(ACE_CAR_CABLE)%utilise = ASTER_TRUE
        infocarte(ACE_CAR_CVCXF)%utilise = ASTER_TRUE
!
        call getvtx('CABLE', 'CARA', iocc=ioc, nbval=nbMaxCar, vect=cara, nbret=ncara)
        call getvr8('CABLE', 'VALE', iocc=ioc, nbval=nbMaxCar, vect=vale, nbret=nvale)
!
!       On commence à Zéro
        nbcara1 = 0
        do ii = 1, nbCara
            do jj = 1, ncara
                if (NomCara(ii) .eq. cara(jj)) then
                    zk8(jdcelems+nbcara1) = NomCata(ii)
                    zr(jdvelems+nbcara1) = vale(jj)
                    nbcara1 = nbcara1+1
                end if
            end do
        end do
        ASSERT(nbcara1 .eq. nbCara)
!       Optionnels, avec une valeur par défaut sur tous les éléments
        k8_fcx = '.'
        call getvid('CABLE', 'FCX', iocc=ioc, scal=k8_fcx, nbret=ibid)
        zk8(jdccfcx) = 'FCXP'
        zk8(jdvcfcx) = k8_fcx
!       "GROUP_MA" = toutes les mailles de la liste de groupes mailles
        do ii = 1, ng
            if (grpnbma(ii) .gt. 0) then
                call nocart(cartelem, 2, nbcara1, groupma=grplmax(ii))
                call nocart(cartefcx, 2, 1, groupma=grplmax(ii))
            end if
        end do
!
    end do bioc
    call jedema()
!
999 continue
end subroutine
