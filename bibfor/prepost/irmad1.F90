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
!
subroutine irmad1(ifi, versio, nbno, prno, nueq, &
                  nec, dg, ncmpmx, itype, nstat, &
                  chamno, nomcmp, nomsym, numnoe)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/exisdg.h"
#include "asterfort/irgags.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ifi, nbno, prno(*), nueq(*), nec, dg(*), ncmpmx, numnoe(*)
    integer(kind=8) :: versio, itype, nstat
    character(len=*) :: nomcmp(*), nomsym, chamno(*)
!        ECRITURE D'UN CHAM_NO SUR FICHIER UNIVERSEL, DATASET TYPE 252
!        A VALEURS REELLES OU COMPLEXES
!--------------------------------------------------------------------
!     ------------------------------------------------------------------
    character(len=4) :: nomgds(50), nomchs(50)
    character(len=19) :: chamn
    integer(kind=8) :: nbchs
    integer(kind=8) :: nbcmps(50), ipcmps(50, 50), impre
    aster_logical :: ltabl(50), afaire
!
!  --- INITIALISATIONS ----
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iavale, ichs, icmp, icms, icmsup, icompt
    integer(kind=8) :: icp, iec, imat, inno, ino, irval, ist
    integer(kind=8) :: ival, k1, k2, mfor, mkey, mtyp, nbcmpt
    integer(kind=8) :: ncmp, ncol, ndim, nrow
!-----------------------------------------------------------------------
    call jemarq()
    do i = 1, ncmpmx
        ltabl(i) = .false.
    end do
!
    nbcmpt = 0
    do inno = 1, nbno
        ino = numnoe(inno)
        do iec = 1, nec
            dg(iec) = prno((ino-1)*(nec+2)+2+iec)
        end do
        ncmp = prno((ino-1)*(nec+2)+2)
        if (ncmp .eq. 0) goto 100
        icompt = 0
        do icmp = 1, ncmpmx
            if (exisdg(dg, icmp)) icompt = icompt+1
        end do
        nbcmpt = max(nbcmpt, icompt)
100     continue
    end do
    nrow = nbcmpt
    ncol = nbno*nstat
    ndim = ncol*nrow
!
! --- ALLOCATION DES TABLEAUX DE TRAVAIL ---
!
    ASSERT((itype .eq. 1) .or. (itype .eq. 2))
    if (itype .eq. 1) then
        call wkvect('&&IRMAD1.VAL', 'V V R', ndim, irval)
    else if (itype .eq. 2) then
        call wkvect('&&IRMAD1.VAL', 'V V C', ndim, irval)
    end if
!
! ---- RECHERCHE DES GRANDEURS SUPERTAB -----
!
    call irgags(ncmpmx, nomcmp, nomsym, nbchs, nomchs, &
                nbcmps, nomgds, ipcmps)
    do ichs = 1, 50
        do ist = 1, 50
            ipcmps(ichs, ist) = -1
        end do
    end do
!
! ---- BOUCLE SUR LES DIVERSES GRANDEURS SUPERTAB ----
    impre = 0
    do ichs = 1, nbchs
        if (ichs .gt. 1) then
            afaire = .false.
            do icp = 1, nbcmps(ichs)
                afaire = (afaire .or. ltabl(ipcmps(ichs, icp)))
            end do
            if (.not. afaire) goto 20
        end if
        impre = impre+1
        do ist = 1, nstat
            chamn = chamno(ist)
            call jeveuo(chamn//'.VALE', 'L', iavale)
            do inno = 1, nbno
                ino = numnoe(inno)
                do iec = 1, nec
                    dg(iec) = prno((ino-1)*(nec+2)+2+iec)
                end do
                ival = prno((ino-1)*(nec+2)+1)
                ncmp = prno((ino-1)*(nec+2)+2)
                if (ncmp .eq. 0) goto 40
                icompt = 0
                do icmp = 1, ncmpmx
                    if (exisdg(dg, icmp)) then
                        if (ichs .eq. 1) ltabl(icmp) = .true.
                        icompt = icompt+1
                        k1 = nueq(ival-1+icompt)
                        do icms = 1, nbcmps(ichs)
                            icmsup = ipcmps(ichs, icms)
                            if (icmp .eq. icmsup) then
                                k2 = icms+(inno-1)*nbcmpt+(ist-1)*nbcmpt*nbno
                                if (itype .eq. 1) then
                                    zr(irval-1+k2) = zr(iavale-1+k1)
                                else
                                    zc(irval-1+k2) = zc(iavale-1+k1)
                                end if
                                goto 44
                            end if
                        end do
                    end if
44                  continue
                end do
40              continue
            end do
        end do
20      continue
    end do
!
    ASSERT(impre .le. 1)
!
    if (versio .eq. 5 .and. impre .eq. 1) then
        imat = 147
        if (itype .eq. 1) then
            mtyp = 4
        else
            mtyp = 6
        end if
        mfor = 3
        mkey = 2
        write (ifi, '(A)') '    -1'
        write (ifi, '(A)') '   252'
        write (ifi, '(I10)') imat
        write (ifi, '(5I10)') mtyp, mfor, nrow, ncol, mkey
        if (itype .eq. 1) then
            write (ifi, '(1P,4D20.12)') (zr(irval+i), i=0, ndim-1)
        else
            write (ifi, '(1P,2(2D20.12))') (zc(irval+i), i=0, ndim-1 &
                                            )
        end if
        mfor = 3
        write (ifi, '(A)') '    -1'
    end if
!
    call jedetr('&&IRMAD1.VAL')
    call jedema()
end subroutine
