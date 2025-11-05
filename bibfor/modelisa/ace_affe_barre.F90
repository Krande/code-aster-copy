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

subroutine ace_affe_barre(nbocc, infoconcept, infocarte, grplmax, grpnbma, lesmailles)
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/ace_affe_verif_elem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nocart.h"
#include "asterfort/tbcarapou.h"
!
    integer(kind=8)         :: nbocc
    type(cara_elem_info)    :: infoconcept
    type(cara_elem_carte)   :: infocarte(*)
    character(len=24)       :: grplmax(*)
    integer(kind=8)         :: grpnbma(*)
    integer(kind=8)         :: lesmailles(*)
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM : AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT BARRE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)   :: GroupeMaxOccur, jdme
    character(len=8)  :: nomu, noma
    character(len=19) :: cartelem, cartefcx
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)   :: ii, jj, ioc, jdcelems, jdccfcx, jdvelems, jdvcfcx
    integer(kind=8)   :: ng, ibid, nbcara1, nbret, nbr_cmp, coderet
    character(len=8)  :: k8bid
    character(len=16) :: typsec
!
    character(len=19) :: tabcar
    character(len=8)  :: nomsec
!
    integer(kind=8)   :: TForme, TVaria, TCodeFV
    integer(kind=8), parameter  :: nbval = 10
    integer(kind=8)   :: valok(nbval)
    character(len=8)  :: valp(nbval)
    real(kind=8)      :: valr(nbval)
!
    integer(kind=8), parameter  :: nbMaxCar = 10
    integer(kind=8)             :: ncara, nvale, nbCara
    real(kind=8)                :: vale(nbMaxCar)
    character(len=8)            :: cara(nbMaxCar)
    character(len=8)            :: NomCara(nbMaxCar)
!-----------------------------------------------------------------------
!
    if (nbocc .eq. 0) goto 999
!
    call jemarq()
!   Récupère les informations
    nomu = infoconcept%nomu
    noma = infoconcept%maillage
    GroupeMaxOccur = infoconcept%GroupeMaxOccur
    jdme = infoconcept%jmodmail
!   La carte est déjà alouée : ace_crea_carte
    cartelem = infocarte(ACE_CAR_BARRE)%nom_carte
    jdcelems = infocarte(ACE_CAR_BARRE)%adr_cmp
    jdvelems = infocarte(ACE_CAR_BARRE)%adr_val
    nbr_cmp = infocarte(ACE_CAR_BARRE)%nbr_cmp
!   La carte est déjà alouée : ace_crea_carte
    cartefcx = infocarte(ACE_CAR_CVCXF)%nom_carte
    jdccfcx = infocarte(ACE_CAR_CVCXF)%adr_cmp
    jdvcfcx = infocarte(ACE_CAR_CVCXF)%adr_val
!
!   lecture des valeurs et affectation dans la carte
    bioc: do ioc = 1, nbocc
!       Les GROUP_MA concernés
        call getvtx('BARRE', 'GROUP_MA', ioc, nbval=GroupeMaxOccur, vect=grplmax, nbret=ng)
!       Le type de section
        call getvtx('BARRE', 'SECTION', iocc=ioc, scal=typsec, nbret=nbret)
        TVaria = ACE_SECTION_CONSTANTE
        if (typsec .eq. 'GENERALE') then
            TForme = ACE_SECTION_GENERALE
            TVaria = ACE_SECTION_CONSTANTE
        else if (typsec .eq. 'RECTANGLE') then
            TForme = ACE_SECTION_RECTANGLE
            TVaria = ACE_SECTION_CONSTANTE
        else if (typsec .eq. 'CERCLE') then
            TForme = ACE_SECTION_CERCLE
            TVaria = ACE_SECTION_CONSTANTE
        else
            ASSERT(.false.)
        end if
        TCodeFV = ACE_CodeFormeVaria(TForme, TVaria)
!       Vérification que les affectations ne concernent que les BARRE
        call ace_affe_verif_elem(noma, jdme, lesmailles, ng, grplmax, grpnbma, &
                                 ACE_NU_BARRE, 'BARRE', coderet, TCodeFV)
!       En // ou pas : si tous les groupes de mailles sont vides ==> on cycle
        if (coderet .eq. 0) cycle bioc
!       Les cartes vont être remplies
        infocarte(ACE_CAR_BARRE)%utilise = ASTER_TRUE
        infocarte(ACE_CAR_CVCXF)%utilise = ASTER_TRUE
!
        if (TForme .eq. ACE_SECTION_GENERALE) then
            valr(:) = 0.0; valp(:) = '...'
            nbCara = 1
            NomCara(1:nbCara) = [character(len=3)::'A']
            call getvid('BARRE', 'TABLE_CARA', iocc=ioc, scal=tabcar, nbret=nbret)
            if (nbret .eq. 1) then
                call getvtx('BARRE', 'NOM_SEC', iocc=ioc, scal=nomsec, nbret=nbret)
                ASSERT(nbret .eq. 1)
                cara(1) = 'A'
                call tbcarapou(tabcar, nomsec, 1, cara, vale, valok)
                ncara = 1
            else
                call getvtx('BARRE', 'CARA', iocc=ioc, nbval=nbMaxCar, vect=cara, nbret=ncara)
                call getvr8('BARRE', 'VALE', iocc=ioc, nbval=nbMaxCar, vect=vale, nbret=nvale)
            end if
!
            ASSERT(ncara .eq. nbCara)
            do ii = 1, nbCara
                do jj = 1, ncara
                    if (NomCara(ii) .eq. cara(jj)) then
                        valr(ii) = vale(jj)
                        valp(ii) = cara(jj)
                    end if
                end do
            end do
!           On commence à Zéro
            nbcara1 = 0
            zk8(jdcelems+nbcara1) = 'A1'
            zr(jdvelems+nbcara1) = valr(1)
            nbcara1 = nbcara1+1
!
!             zk8(jdcelems+nbcara1) = 'TSEC'
!             zr(jdvelems+nbcara1) = ACE_SECTION_GENERALE
!             nbcara1 = nbcara1+1
!
        else if (TForme .eq. ACE_SECTION_RECTANGLE) then
            valr(:) = 0.0; valp(:) = '...'
            nbCara = 4
            NomCara(1:nbCara) = [character(len=3)::'HY', 'HZ', 'EPY', 'EPZ']
            call getvtx('BARRE', 'CARA', iocc=ioc, nbval=nbMaxCar, vect=cara, nbret=ncara)
            call getvr8('BARRE', 'VALE', iocc=ioc, nbval=nbMaxCar, vect=vale, nbret=nvale)
!
            ASSERT(ncara .eq. nbCara)
            do ii = 1, nbCara
                do jj = 1, ncara
                    if (NomCara(ii) .eq. cara(jj)) then
                        valr(ii) = vale(jj)
                        valp(ii) = cara(jj)
                    end if
                end do
            end do
!           On commence à Zéro
            nbcara1 = 0
            zk8(jdcelems+nbcara1) = 'A1'
            zr(jdvelems+nbcara1) = valr(1)*valr(2)-(valr(1)-valr(3)*2.0)*(valr(2)-valr(4)*2.0)
            nbcara1 = nbcara1+1
!
        else if (TForme .eq. ACE_SECTION_CERCLE) then
            valr(:) = 0.0; valp(:) = '...'
            nbCara = 2
            NomCara(1:nbCara) = [character(len=2)::'R', 'EP']
            call getvtx('BARRE', 'CARA', iocc=ioc, nbval=nbMaxCar, vect=cara, nbret=ncara)
            call getvr8('BARRE', 'VALE', iocc=ioc, nbval=nbMaxCar, vect=vale, nbret=nvale)
!
            ASSERT(ncara .eq. nbCara)
            do ii = 1, nbCara
                do jj = 1, ncara
                    if (NomCara(ii) .eq. cara(jj)) then
                        valr(ii) = vale(jj)
                        valp(ii) = cara(jj)
                    end if
                end do
            end do
!           On commence à Zéro
            nbcara1 = 0
            zk8(jdcelems+nbcara1) = 'A1'
            zr(jdvelems+nbcara1) = r8pi()*(valr(1)*valr(1)-(valr(1)-valr(2))*(valr(1)-valr(2)))
            nbcara1 = nbcara1+1
        else
            ASSERT(.false.)
        end if
!       Optionnels, avec une valeur par défaut sur tous les éléments
        k8bid = '.'
        call getvid('BARRE', 'FCX', iocc=ioc, scal=k8bid, nbret=ibid)
        zk8(jdccfcx) = 'FCXP'
        zk8(jdvcfcx) = k8bid
!       "GROUP_MA" = toutes les mailles de la liste de groupes mailles présente sur le proc
        do ii = 1, ng
            if (grpnbma(ii) .gt. 0) then
                call nocart(cartelem, 2, nbcara1, groupma=grplmax(ii))
                call nocart(cartefcx, 2, 1, groupma=grplmax(ii))
            end if
        end do
        !
    end do bioc
!
    call jedema()
!
999 continue
end subroutine
