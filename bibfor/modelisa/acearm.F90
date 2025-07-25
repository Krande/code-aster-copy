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

subroutine acearm(infdonn, lmax, nbocc, infcarte, ivr)
!
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
#include "jeveux.h"
#include "asterfort/affdis.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/irmifr.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/nocart.h"
#include "asterfort/r8inir.h"
#include "asterfort/rigmi1.h"
#include "asterfort/rigmi2.h"
#include "asterfort/ulisop.h"
#include "asterfort/ulopen.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    type(cara_elem_info) :: infdonn
    type(cara_elem_carte) :: infcarte(*)

    integer(kind=8) :: lmax, nbocc, ivr(*), ifm

! ----------------------------------------------------------------------
!     AFFE_CARA_ELEM
!     AFFECTATION DES CARACTERISTIQUES POUR LES ELEMENTS DISCRET
! ----------------------------------------------------------------------
! IN  : NOMA   : NOM DU MAILLAGE
! IN  : LMAX   : NOMBRE MAX DE MAILLE OU GROUPE DE MAILLE
! IN  : NBOCC  : NOMBRE D'OCCURENCES DU MOT CLE DISCRET
! IN  : IVR    : TABLEAU DES INDICES DE VERIFICATION
! ----------------------------------------------------------------------
    integer(kind=8) :: nrd
    parameter(nrd=2)
    integer(kind=8) :: jdc(3), jdv(3), dimcar, irgma, irgm2, irgm3
    integer(kind=8) :: itbmp, ndim, jdcinf, jdvinf, i, ioc
    integer(kind=8) :: irep, isym, impris, nu, nfr, ngp, ngl, ifreq, nma, ldgm, nbpo
    integer(kind=8) :: in, nfreq, iv, jd, ncmp, l, nbli, ncmp2, icf
    real(kind=8) :: eta, vale(3), freq, coef, zero(5)
    character(len=1) :: kma(3)
    character(len=7) :: ledisc
    character(len=8) :: nommai, noma
    character(len=9) :: cara
    character(len=16) :: rep, repdis(nrd), k16nom
    character(len=19) :: cart(3), cartdi
    character(len=24) :: nogp, nogl
!
!     integer :: ixckma, ixci
!     real(kind=8) :: r8bid
!     character(len=8) :: k8bid
!     character(len=24) :: tmpnd(3), tmpvd(3), tmcinf, tmvinf
!
    data repdis/'GLOBAL          ', 'LOCAL           '/
    data kma/'K', 'M', 'A'/
!     ------------------------------------------------------------------
!
    call jemarq()
!
    noma = infdonn%maillage
    ndim = infdonn%dimmod
!   Pour miss3d c'est obligatoirement du 3d
    ASSERT(ndim .eq. 3)
!
    call wkvect('&&TMPRIGMA', 'V V R', 3*lmax, irgma)
    call wkvect('&&TMPRIGM2', 'V V R', 3*lmax, irgm2)
    call wkvect('&&TMPRIGM3', 'V V R', 3*lmax, irgm3)
    call wkvect('&&TMPTABMP', 'V V K8', lmax, itbmp)
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
!   On ne peut faire qu'une occurrence de RIGI_MISS_3D : catalogue
!   On garde la boucle, au cas où l'on souhaiterait faire des évolutions
    ASSERT(nbocc .eq. 1)
    do ioc = 1, nbocc
        eta = 0.0d0
!       PAR DEFAUT ON EST DANS LE REPERE GLOBAL, MATRICES SYMETRIQUES
        irep = 1
        isym = 1
        rep = repdis(1)
        call getvis('RIGI_MISS_3D', 'UNITE_RESU_IMPE', iocc=ioc, scal=impris, nbret=nu)
        k16nom = ' '
        if (ulisop(impris, k16nom) .eq. 0) then
            call ulopen(impris, ' ', ' ', 'NEW', 'O')
        end if
        call getvr8('RIGI_MISS_3D', 'FREQ_EXTR', iocc=ioc, scal=freq, nbret=nfr)
        call getvtx('RIGI_MISS_3D', 'GROUP_MA_POI1', iocc=ioc, scal=nogp, nbret=ngp)
        call getvtx('RIGI_MISS_3D', 'GROUP_MA_SEG2', iocc=ioc, scal=nogl, nbret=ngl)
        do i = 1, nrd
            if (rep .eq. repdis(i)) irep = i
        end do
        if (ivr(3) .eq. 2) then
            write (ifm, 100) rep, ioc
        end if
100     format(/, 3x, &
                '<DISCRET> MATRICES AFFECTEES AUX ELEMENTS DISCRET ', &
                '(REPERE ', a6, '), OCCURENCE ', i4)
        call irmifr(impris, freq, ifreq, nfreq, icf)
!
! ---    "GROUP_MA" = TOUTES LES MAILLES DE TOUS LES GROUPES DE MAILLES
        if (ngl .ne. 0) then
            call rigmi2(noma, nogl, ifreq, nfreq, impris, &
                        zr(irgm2), zr(irgm3))
        end if
!
        call r8inir(5, 0.0d0, zero, 1)
!
        cara = 'K_T_D_N'
        if (ngp .ne. 0) then
            call jelira(jexnom(noma//'.GROUPEMA', nogp), 'LONMAX', nma)
            call jeveuo(jexnom(noma//'.GROUPEMA', nogp), 'L', ldgm)
            nbpo = nma
            call rigmi1(noma, nogp, ifreq, nfreq, impris, &
                        zr(irgma), zr(irgm3))
            do in = 0, nma-1
                nommai = int_to_char8(zi(ldgm+in))
                zk8(itbmp+in) = nommai
            end do
            do i = 1, nbpo
                iv = 1
                jd = itbmp+i-1
                call affdis(ndim, irep, eta, cara, zr(irgma+3*i-3), &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
!              AFFECTATION MATRICE MASSE NULLE
                iv = 1
                ledisc = 'M_T_D_N'
                call affdis(ndim, irep, eta, ledisc, zero, &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
!              AFFECTATION MATRICE AMORTISSEMENT NULLE
                iv = 1
                ledisc = 'A_T_D_N'
                call affdis(ndim, irep, eta, ledisc, zero, &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
            end do
        end if
!
        cara = 'K_T_D_L'
        if (ngl .ne. 0) then
            coef = 20.d0
            call jelira(jexnom(noma//'.GROUPEMA', nogl), 'LONMAX', nma)
            call jeveuo(jexnom(noma//'.GROUPEMA', nogl), 'L', ldgm)
            nbli = nma
            do in = 0, nma-1
                nommai = int_to_char8(zi(ldgm+in))
                zk8(itbmp+in) = nommai
            end do
            call r8inir(3, 0.d0, vale, 1)
            do i = 1, nbli
                iv = 1
                jd = itbmp+i-1
                vale(1) = -zr(irgm2+3*i-3)*coef
                vale(2) = -zr(irgm2+3*i-2)*coef
                vale(3) = -zr(irgm2+3*i-1)*coef
                call affdis(ndim, irep, eta, cara, vale, &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp2, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp2, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
!               AFFECTATION MATRICE MASSE NULLE
                iv = 1
                ledisc = 'M_T_D_L'
                call affdis(ndim, irep, eta, ledisc, zero, &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp2, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp2, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
!               AFFECTATION MATRICE AMORTISSEMENT NULLE
                iv = 1
                ledisc = 'A_T_D_L'
                call affdis(ndim, irep, eta, ledisc, zero, &
                            jdc, jdv, ivr, iv, kma, &
                            ncmp2, l, jdcinf, jdvinf, isym)
                call nocart(cartdi, 3, dimcar, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
                call nocart(cart(l), 3, ncmp2, mode='NOM', nma=1, &
                            limano=[zk8(jd)])
            end do
        end if
    end do
!
    call jedetr('&&TMPRIGMA')
    call jedetr('&&TMPRIGM2')
    call jedetr('&&TMPRIGM3')
    call jedetr('&&TMPTABMP')
    call jedetr('&&ACEARM.RIGM')
!
    call jedema()
end subroutine
