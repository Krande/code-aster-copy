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

subroutine pgpcrt(sd_pgp)
    implicit none
! Initialzes a table data structure as a result of the POST_GENE_PHYS
! command
! ----------------------------------------------------------------------
! person_in_charge: hassan.berro at edf.fr
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/pgpget.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"

!   ====================================================================
!   = 0 =   Variable declarations and initialization
!   ====================================================================
!   -0.1- Input/output arguments
    character(len=8)  :: sd_pgp
!   -0.2- Local variables
    aster_logical :: exist_complex
    integer(kind=8) :: nbparams
    parameter(nbparams=10)
    integer(kind=8) :: iobs, ipar, nbobs, nblines, nord
    integer(kind=8) :: i, iord, physlen, dec1, lc, dec
    integer(kind=8) :: jtab, jvec, jlog, jlogr, jlogc, jdsc
    integer(kind=8) :: lgnoeu, lgmail, lgpoin
    character(len=4)  :: typcha, typresin, typsc
    character(len=8)  :: result, partyp(nbparams)
    character(len=16) :: params(nbparams), champ
    character(len=24) :: nomjv, discjv, nomlgs(nbparams), nomjvs(nbparams)

    integer(kind=8), pointer :: indic(:) => null()
    character(len=8), pointer :: lcmp(:) => null()
    character(len=8), pointer :: nomnoeu(:) => null()
    character(len=8), pointer :: nommail(:) => null()
    character(len=8), pointer :: nompoin(:) => null()

!
!   -0.3- Initialization
    data params/'NUME_OBSE', 'INST', 'FREQ', 'NOM_CHAM', 'NOM_CMP', &
        'NOEUD', 'MAILLE', 'POINT', 'VALE_R', 'VALE_C'/
    data partyp/'I', 'R', 'R', 'K16', 'K8', &
        'K8', 'K8', 'K8', 'R', 'C'/

!   ------------------------------------------------------------------------------------
!   Definition of statement functions
#define disc(i) zr(jdsc+i-1)

    call jemarq()

!   ====================================================================
!   = 1 = Initialization of the table's data structure
!   ====================================================================
    call pgpget(sd_pgp, 'RESU_OUT', kscal=result)
    call pgpget(sd_pgp, 'TYP_RESU', kscal=typresin)
    call pgpget(sd_pgp, 'NB_OBSER', iscal=nbobs)

    call tbcrsd(result, 'G')
    call tbajpa(result, nbparams, params, partyp)

    nblines = 0
    do iobs = 1, nbobs
        call pgpget(sd_pgp, 'REF_COMP', iobs=iobs, lonvec=physlen)
        call pgpget(sd_pgp, 'DISC', iobs=iobs, lonvec=nord)

        nblines = nblines+nord*physlen
    end do

    call jeveuo(result//'           .TBLP', 'E', jtab)

    do ipar = 1, nbparams
        nomjv = zk24(jtab+4*(ipar-1)+2)
        call juveca(nomjv, nblines)
        call jeecra(nomjv, 'LONUTI', nblines)
        call jelibe(nomjv)
        nomjvs(ipar) = nomjv

        nomjv = zk24(jtab+4*(ipar-1)+3)
        call juveca(nomjv, nblines)
        call jeecra(nomjv, 'LONUTI', nblines)
        call jelibe(nomjv)
        nomlgs(ipar) = nomjv
    end do

!   Memory optimisation, for params 4 and 5 (CHAMP, COMP) refer to the logicals of
!   NUME_OBSE, i.e. ipar=1 in the table. Delete the logical objects corresponding
!   to ipar=4 and ipar=5
    call jedetr(nomlgs(4))
    call jedetr(nomlgs(5))
    zk24(jtab+4*(4-1)+3) = zk24(jtab+4*(1-1)+3)
    zk24(jtab+4*(5-1)+3) = zk24(jtab+4*(1-1)+3)

!   Memory optimisation, prefill the INST (2), FREQ (3), VALE_R (9), and VALE_C (10)
!   logicals, based on the type of result (transient or harmonic)
    if (typresin .eq. 'TRAN') then
        call jedetr(nomlgs(2))
        zk24(jtab+4*(2-1)+3) = zk24(jtab+4*(1-1)+3)

!       set FREQ logicals to 0
        call jeveuo(nomlgs(3), 'E', jlog)
        do i = 1, nblines
            zi(jlog+i-1) = 0
        end do
        call jelibe(nomlgs(3))

!       In the transient case, it might still be possible that the fields of the basis
!       themselves are complex, the result of the modal combination would thus be
!       complex, it is thus important to do a further check of the TYP_SCAL saved for
!       each observation
        exist_complex = .false.

!       Initialize logicals for VALE_R and VALE_C to false (0)
        call jeveuo(nomlgs(9), 'E', jlogr)
        call jeveuo(nomlgs(10), 'E', jlogc)
        do i = 1, nblines
            zi(jlogr+i-1) = 0
            zi(jlogc+i-1) = 0
        end do

!       Properly fill the logicals for for VALE_R and VALE_C per observation
        dec = 0
        do iobs = 1, nbobs
            call pgpget(sd_pgp, 'TYP_SCAL ', iobs=iobs, kscal=typsc)
            call pgpget(sd_pgp, 'REF_COMP', iobs=iobs, lonvec=physlen)
            call pgpget(sd_pgp, 'DISC', iobs=iobs, lonvec=nord)
            if (typsc(1:1) .eq. 'R') then
                do i = 1, physlen*nord
                    zi(jlogr+dec+i-1) = 1
                end do
            else
                exist_complex = .true.
                do i = 1, physlen*nord
                    zi(jlogc+dec+i-1) = 1
                end do
            end if
            dec = dec+physlen*nord
        end do
        call jelibe(nomlgs(9))
        call jelibe(nomlgs(10))

!       If no complex field is found, optimize memory by deleting the logicals for VALE_*
!       and refer to those of NUME_OBSE(1's) and FREQ (0's)
        if (.not. (exist_complex)) then
            call jedetr(nomlgs(9))
            call jedetr(nomlgs(10))
            zk24(jtab+4*(9-1)+3) = zk24(jtab+4*(1-1)+3)
            zk24(jtab+4*(10-1)+3) = zk24(jtab+4*(3-1)+3)
        end if

    else
        call jedetr(nomlgs(3))
        call jedetr(nomlgs(10))
        zk24(jtab+4*(3-1)+3) = zk24(jtab+4*(1-1)+3)
        zk24(jtab+4*(10-1)+3) = zk24(jtab+4*(1-1)+3)

!       set INST logicals to 0
        call jeveuo(nomlgs(2), 'L', jlog)
        do i = 1, nblines
            zi(jlog+i-1) = 0
        end do
        call jelibe(nomlgs(2))

        call jedetr(nomlgs(9))
        zk24(jtab+4*(9-1)+3) = zk24(jtab+4*(2-1)+3)
    end if

!   Line counter, across the whole table (for different observations)
    lc = 0
    do iobs = 1, nbobs

        call pgpget(sd_pgp, 'DISC ', iobs=iobs, lonvec=nord)
        call pgpget(sd_pgp, 'DISC ', iobs=iobs, savejv=discjv)

        call jeveuo(discjv, 'L', jdsc)

        call pgpget(sd_pgp, 'NOM_CHAM ', iobs=iobs, kscal=champ)
        call pgpget(sd_pgp, 'TYP_CHAM ', iobs=iobs, kscal=typcha)

        call pgpget(sd_pgp, 'REF_COMP', iobs=iobs, lonvec=physlen)
        AS_ALLOCATE(vk8=lcmp, size=physlen)
        call pgpget(sd_pgp, 'REF_COMP', iobs=iobs, kvect=lcmp)

        AS_ALLOCATE(vi=indic, size=physlen)
        call pgpget(sd_pgp, 'REF_INDI', iobs=iobs, ivect=indic)

!
        if (typcha .eq. 'NOEU') then
            AS_ALLOCATE(vk8=nomnoeu, size=physlen)
            call pgpget(sd_pgp, 'REF_SUP1', iobs=iobs, kvect=nomnoeu)

            AS_ALLOCATE(vk8=nommail, size=physlen)
            AS_ALLOCATE(vk8=nompoin, size=physlen)
            do i = 1, physlen
                nommail(i) = ' '
                nompoin(i) = ' '
            end do

            lgnoeu = 1
            lgmail = 0
            lgpoin = 0
        else if (typcha .eq. 'ELNO') then
            AS_ALLOCATE(vk8=nommail, size=physlen)
            call pgpget(sd_pgp, 'REF_SUP1', iobs=iobs, kvect=nommail)
            AS_ALLOCATE(vk8=nomnoeu, size=physlen)
            call pgpget(sd_pgp, 'REF_SUP2', iobs=iobs, kvect=nomnoeu)

            AS_ALLOCATE(vk8=nompoin, size=physlen)
            do i = 1, physlen
                nompoin(i) = ' '
            end do

            lgnoeu = 1
            lgmail = 1
            lgpoin = 0
        else
            AS_ALLOCATE(vk8=nommail, size=physlen)
            call pgpget(sd_pgp, 'REF_SUP1', iobs=iobs, kvect=nommail)
            AS_ALLOCATE(vk8=nompoin, size=physlen)
            call pgpget(sd_pgp, 'REF_SUP2', iobs=iobs, kvect=nompoin)

            AS_ALLOCATE(vk8=nomnoeu, size=physlen)
            do i = 1, physlen
                nomnoeu(i) = ' '
            end do

            lgnoeu = 0
            lgmail = 1
            lgpoin = 1
        end if
!
!       Memory optimisation, treatment of each parameter separately
!
!       Parameter 1 : NUME_OBS
        call jeveuo(nomjvs(1), 'E', jvec)
        call jeveuo(nomlgs(1), 'E', jlog)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zi(jvec+dec1+i-1) = iobs
                zi(jlog+dec1+i-1) = indic(i)
            end do
        end do
        call jelibe(nomjvs(1))
        call jelibe(nomlgs(1))

!       Parameter 2 : INSTANT
        call jeveuo(nomjvs(2), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zr(jvec+dec1+i-1) = disc(iord)
            end do
        end do
        call jelibe(nomjvs(2))

!       Parameter 3 : FREQ
        call jeveuo(nomjvs(3), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zr(jvec+dec1+i-1) = disc(iord)
            end do
        end do
        call jelibe(nomjvs(3))

        call jelibe(discjv)

!       Parameter 4 : CHAMP
        call jeveuo(nomjvs(4), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zk16(jvec+dec1+i-1) = champ
            end do
        end do
        call jelibe(nomjvs(4))

!       Parameter 5 : COMPOSANT
        call jeveuo(nomjvs(5), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zk8(jvec+dec1+i-1) = lcmp(i)
            end do
        end do
        call jelibe(nomjvs(5))

!       Parameter 6 : NOEUD
        call jeveuo(nomjvs(6), 'E', jvec)
        call jeveuo(nomlgs(6), 'E', jlog)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zi(jlog+dec1+i-1) = lgnoeu*indic(i)
                zk8(jvec+dec1+i-1) = nomnoeu(i)
            end do
        end do
        call jelibe(nomjvs(6))
        call jelibe(nomlgs(6))

!       Parameter 7 : MAILLE
        call jeveuo(nomjvs(7), 'E', jvec)
        call jeveuo(nomlgs(7), 'E', jlog)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zi(jlog+dec1+i-1) = lgmail*indic(i)
                zk8(jvec+dec1+i-1) = nommail(i)
            end do
        end do
        call jelibe(nomjvs(7))
        call jelibe(nomlgs(7))

!       Parameter 8 : POINT
        call jeveuo(nomjvs(8), 'E', jvec)
        call jeveuo(nomlgs(8), 'E', jlog)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zi(jlog+dec1+i-1) = lgpoin*indic(i)
                zk8(jvec+dec1+i-1) = nompoin(i)
            end do
        end do
        call jelibe(nomjvs(8))
        call jelibe(nomlgs(8))

!       Parameter 9 : VALE_R
        call jeveuo(nomjvs(9), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zr(jvec+dec1+i-1) = 0.d0
            end do
        end do
        call jelibe(nomjvs(9))

!       Parameter 10 : VALE_C
        call jeveuo(nomjvs(10), 'E', jvec)
        do iord = 1, nord
            dec1 = lc+(iord-1)*physlen
            do i = 1, physlen
                zc(jvec+dec1+i-1) = dcmplx(0.d0, 0.d0)
            end do
        end do
        call jelibe(nomjvs(10))
!
        lc = lc+nord*physlen
!

        AS_DEALLOCATE(vi=indic)
        AS_DEALLOCATE(vk8=lcmp)
        AS_DEALLOCATE(vk8=nomnoeu)
        AS_DEALLOCATE(vk8=nommail)
        AS_DEALLOCATE(vk8=nompoin)

    end do

    call jeveuo(result//'           .TBNP', 'E', jtab)

!   Finally update the number of lines saved
    zi(jtab+1) = nblines

    call jedema()

end subroutine
