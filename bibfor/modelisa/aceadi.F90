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

subroutine aceadi(noma, nomo, mcf, lmax, nbocc, infcarte, ivr)
!
!
! --------------------------------------------------------------------------------------------------
!
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES ELEMENTS DISCRET
!
! --------------------------------------------------------------------------------------------------
!
!  IN
!     NOMA   : NOM DU MAILLAGE
!     NOMO   : NOM DU MODELE
!     LMAX   : NOMBRE MAX DE MAILLE OU GROUPE DE MAILLE
!     NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE DISCRET
!     IVR    : TABLEAU DES INDICES DE VERIFICATION
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    use cara_elem_parameter_module
    use cara_elem_carte_type
    implicit none
    character(len=8) :: noma, nomo
    integer(kind=8) :: lmax, nbocc, ivr(*), ifm
    type(cara_elem_carte) :: infcarte(*)
    character(len=*) :: mcf
!
#include "jeveux.h"
!
#include "asterc/getres.h"
!
#include "asterfort/affdis.h"
#include "asterfort/assert.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/nocart.h"
#include "asterfort/verdis.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: nbcar, nbval, nrd
    parameter(nbcar=100, nbval=1000, nrd=2)
    integer(kind=8) :: jdc(3), jdv(3), dimcar, nm, ii, l, iv, ndim
    integer(kind=8) :: jdcinf, jdvinf, ncmp, nn
    integer(kind=8) :: nsym, neta, nrep, i3d, i2d, ier
    integer(kind=8) :: jdls, i, ioc, irep, isym, ng, nj, ncar
    integer(kind=8) :: nval, jdls2
    real(kind=8) :: val(nbval), eta
    character(len=1) :: kma(3)
    character(len=8) :: nomu
    character(len=9) :: car(nbcar)
    character(len=16) :: rep, repdis(nrd), concep, cmd, sym, symdis(nrd)
    character(len=19) :: cart(3), ligmo, cartdi
    character(len=24) :: tmpdis, mlggno
!
! --------------------------------------------------------------------------------------------------
    data repdis/'GLOBAL          ', 'LOCAL           '/
    data symdis/'OUI             ', 'NON             '/
    data kma/'K', 'M', 'A'/
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call getres(nomu, concep, cmd)
    tmpdis = nomu//'.DISCRET'
    mlggno = noma//'.GROUPENO'
    ligmo = nomo//'.MODELE    '
!
!   Vérification des dimensions / modélisations
    ier = 0
    call verdis(nomo, noma, 'F', i3d, i2d, ndim, ier)
    ASSERT((mcf .eq. 'DISCRET_2D') .or. (mcf .eq. 'DISCRET'))
!
    call wkvect('&&TMPDISCRET', 'V V K24', lmax, jdls)
    call wkvect('&&TMPDISCRET2', 'V V K8', lmax, jdls2)
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
    nj = 0
    nn = 0
    do ioc = 1, nbocc
        eta = 0.0d0
        irep = 1
        isym = 1
        val(:) = 0.0d0
        call getvem(noma, 'GROUP_MA', mcf, 'GROUP_MA', ioc, lmax, zk24(jdls), ng)
        call getvem(noma, 'MAILLE', mcf, 'MAILLE', ioc, lmax, zk8(jdls2), nm)
        call getvr8(mcf, 'VALE', iocc=ioc, nbval=nbval, vect=val, nbret=nval)
        ASSERT(nbval .ge. 1)
        call getvtx(mcf, 'CARA', iocc=ioc, nbval=nbcar, vect=car, nbret=ncar)
        ASSERT(ncar .eq. 1)
!
        call getvtx(mcf, 'REPERE', iocc=ioc, scal=rep, nbret=nrep)
        call getvr8(mcf, 'AMOR_HYST', iocc=ioc, scal=eta, nbret=neta)
        if (ioc .eq. 1 .and. nrep .eq. 0) rep = repdis(1)
        do i = 1, nrd
            if (rep .eq. repdis(i)) irep = i
        end do
!
!       Matrice symétrique ou non-symétrique : par défaut symétrique
        call getvtx(mcf, 'SYME', iocc=ioc, scal=sym, nbret=nsym)
        if (nsym .eq. 0) sym = symdis(1)
        do i = 1, nrd
            if (sym .eq. symdis(i)) isym = i
        end do
!
        if (ivr(3) .eq. 2) then
            if (isym .eq. 1) then
                write (ifm, 100) rep, 'SYMETRIQUE', ioc
            else
                write (ifm, 100) rep, 'NON-SYMETRIQUE', ioc
            end if
        end if
!       GROUP_MA = toutes les mailles de tous les groupes de mailles
        if (ng .gt. 0) then
            iv = 1
            do i = 1, ncar
                call affdis(ndim, irep, eta, car(i), val, jdc, jdv, ivr, iv, kma, &
                            ncmp, l, jdcinf, jdvinf, isym)
                do ii = 1, ng
                    call nocart(cartdi, 2, dimcar, groupma=zk24(jdls+ii-1))
                    call nocart(cart(l), 2, ncmp, groupma=zk24(jdls+ii-1))
                end do
            end do
        end if
!       MAILLE = toutes les mailles de la liste de mailles
        if (nm .gt. 0) then
            iv = 1
            do i = 1, ncar
                call affdis(ndim, irep, eta, car(i), val, jdc, jdv, ivr, iv, kma, &
                            ncmp, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=nm, limano=zk8(jdls2))
                call nocart(cart(l), 3, ncmp, mode='NOM', nma=nm, limano=zk8(jdls2))
            end do
        end if
    end do
!
    call jedetr('&&TMPDISCRET')
    call jedetr('&&TMPDISCRET2')
!
    call jedema()
!
100 format(/, 3x, '<DISCRET> MATRICES (REPERE ', a6, ') AFFECTEES AUX ELEMENTS DISCRETS ', &
            '(TYPE ', a, '), OCCURENCE ', i4)
end subroutine
