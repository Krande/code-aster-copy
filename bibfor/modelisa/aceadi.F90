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

subroutine aceadi(nbocc, infoconcept, infcarte, mcf)
!
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!
!     Affectation des caractéristiques pour les éléments discrets
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
#include "jeveux.h"
#include "asterfort/affdis.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
#include "asterfort/verima.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbocc
    type(cara_elem_info) :: infoconcept
    type(cara_elem_carte) :: infcarte(*)
    character(len=*) :: mcf
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbval, nrd, lmax, ifm
!   nbval : nombre de valeurs maximum à lire pour *_TR_L non-symétrique 12*12 = 144
    parameter(nbval=150, nrd=2)
    integer(kind=8) :: jdc(3), jdv(3), dimcar, ii, jj, ikma, iv, ndim, jmail, nbmail, numa
    integer(kind=8) :: jdcinf, jdvinf, ncmp, nutyma
    integer(kind=8) :: nsym, neta, nrep, ivr(4)
    integer(kind=8) :: ioc, irep, isym, ng
    integer(kind=8) :: nval
    real(kind=8) :: val(nbval), eta
    character(len=1) :: kma(3)
    character(len=4) :: letype
    character(len=8) :: noma, nomo, cara, typel, nomail
    character(len=16) :: rep, repdis(nrd), sym, symdis(nrd)
    character(len=19) :: cart(3), cartdi
    character(len=24) :: grmama
!
    character(len=24) :: valk(4)
! --------------------------------------------------------------------------------------------------
    integer(kind=8), pointer :: typmail(:) => null()
    character(len=24), pointer  :: group_ma(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    data repdis/'GLOBAL          ', 'LOCAL           '/
    data symdis/'OUI             ', 'NON             '/
    data kma/'K', 'M', 'A'/
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    nomo = infoconcept%modele
    noma = infoconcept%maillage
    ndim = infoconcept%dimmod
    lmax = infoconcept%GroupeMaxOccur
    ivr(:) = infoconcept%ivr(:)
!
!   Pour controle
    grmama = noma//'.GROUPEMA'
!   Vecteur du type des mailles du maillage
    if (infoconcept%VerifMaille .or. infoconcept%IsParaMesh) then
        call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
    end if
!
    ASSERT((mcf .eq. 'DISCRET_2D') .or. (mcf .eq. 'DISCRET'))
!
    AS_ALLOCATE(vk24=group_ma, size=lmax)
!
!   Les cartes sont déjà construites : ace_crea_carte
    cartdi = infcarte(ACE_CAR_DINFO)%nom_carte
    jdcinf = infcarte(ACE_CAR_DINFO)%adr_cmp
    jdvinf = infcarte(ACE_CAR_DINFO)%adr_val
    dimcar = infcarte(ACE_CAR_DINFO)%nbr_cmp
!
    cart(1) = infcarte(ACE_CAR_DISCK)%nom_carte
    jdc(1) = infcarte(ACE_CAR_DISCK)%adr_cmp
    jdv(1) = infcarte(ACE_CAR_DISCK)%adr_val
!
    cart(2) = infcarte(ACE_CAR_DISCM)%nom_carte
    jdc(2) = infcarte(ACE_CAR_DISCM)%adr_cmp
    jdv(2) = infcarte(ACE_CAR_DISCM)%adr_val
!
    cart(3) = infcarte(ACE_CAR_DISCA)%nom_carte
    jdc(3) = infcarte(ACE_CAR_DISCA)%adr_cmp
    jdv(3) = infcarte(ACE_CAR_DISCA)%adr_val
!
    ifm = ivr(4)
!   Boucle sur les occurences de discret
    do ioc = 1, nbocc
        eta = 0.0d0
        val(:) = 0.0d0
        call getvtx(mcf, 'GROUP_MA', iocc=ioc, nbval=lmax, vect=group_ma, nbret=ng)
        ! En //, il faut vérifier que les groupes existent sur le proc
        !   en sortie de verima ==> les groupes présent sur le proc
        if (infoconcept%VerifMaille .or. infoconcept%IsParaMesh) then
            call verima(noma, group_ma, ng, 'GROUP_MA')
        end if
        if (ng .eq. 0) cycle
        !
        call getvr8(mcf, 'VALE', iocc=ioc, nbval=nbval, vect=val, nbret=nval)
        ASSERT(abs(nval) .le. nbval)
        call getvtx(mcf, 'CARA', iocc=ioc, scal=cara)
        call getvr8(mcf, 'AMOR_HYST', iocc=ioc, scal=eta, nbret=neta)
!
!       Repère : par défaut GLOBAL
        irep = 1
        call getvtx(mcf, 'REPERE', iocc=ioc, scal=rep, nbret=nrep)
        if (nrep .ge. 1) then
            do ii = 1, nrd
                if (rep .eq. repdis(ii)) irep = ii
            end do
        end if
!       Matrice symétrique ou non-symétrique : par défaut symétrique
        isym = 1
        call getvtx(mcf, 'SYME', iocc=ioc, scal=sym, nbret=nsym)
        if (nsym .ge. 1) then
            do ii = 1, nrd
                if (sym .eq. symdis(ii)) isym = ii
            end do
        end if
!
!       Vérification du bon type de maille en fonction de cara
        if (infoconcept%VerifMaille .or. infoconcept%IsParaMesh) then
            if ((cara(2:7) .eq. '_T_D_N') .or. (cara(2:8) .eq. '_TR_D_N') .or. &
                (cara(2:5) .eq. '_T_N') .or. (cara(2:6) .eq. '_TR_N')) then
                letype = 'POI1'
            else
                letype = 'SEG2'
            end if
            !
            do ii = 1, ng
                call jelira(jexnom(grmama, group_ma(ii)), 'LONUTI', nbmail)
                call jeveuo(jexnom(grmama, group_ma(ii)), 'L', jmail)
                do jj = 1, nbmail
                    numa = zi(jmail+jj-1)
                    nutyma = typmail(numa)
                    call jenuno(jexnum('&CATA.TM.NOMTM', nutyma), typel)
                    if (typel(1:4) .ne. letype) then
                        nomail = int_to_char8(numa)
                        valk(1) = nomail
                        valk(2) = letype
                        valk(3) = typel
                        valk(4) = cara
                        call utmess('F', 'MODELISA_56', nk=4, valk=valk)
                    end if
                end do
            end do
        end if
!
        if (ivr(3) .eq. 2) then
            if (isym .eq. 1) then
                write (ifm, 100) rep, 'SYMETRIQUE', ioc
            else
                write (ifm, 100) rep, 'NON-SYMETRIQUE', ioc
            end if
        end if
!
!       GROUP_MA = toutes les mailles de tous les groupes de mailles
        iv = 1
        call affdis(ndim, irep, eta, cara, val, jdc, jdv, ivr, iv, kma, &
                    ncmp, ikma, jdcinf, jdvinf, isym)
        do ii = 1, ng
            call nocart(cartdi, 2, dimcar, groupma=group_ma(ii))
            call nocart(cart(ikma), 2, ncmp, groupma=group_ma(ii))
        end do
        if (cara(1:1) .eq. 'K') then
            infcarte(ACE_CAR_DISCK)%utilise = .true.
        else if (cara(1:1) .eq. 'M') then
            infcarte(ACE_CAR_DISCM)%utilise = .true.
        else if (cara(1:1) .eq. 'A') then
            infcarte(ACE_CAR_DISCA)%utilise = .true.
        end if
    end do
!
    AS_DEALLOCATE(vk24=group_ma)
!
    call jedema()
!
100 format(/, 3x, '<DISCRET> MATRICES (REPERE ', a6, ') AFFECTEES AUX ELEMENTS DISCRETS ', &
            '(TYPE ', a, '), OCCURENCE ', i4)
end subroutine
