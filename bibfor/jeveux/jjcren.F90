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
subroutine jjcren(nomlu, icre, iret)
!
! In case of failure in this subroutine, check the string `nomlu`.
!
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjarep.h"
#include "asterfort/jxhcod.h"
#include "asterfort/utmess.h"
    character(len=*) :: nomlu
    integer(kind=8) :: icre, iret
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icla, iclasi, in, iref, isg, j
    integer(kind=8) :: jcara, jdate, jdocu, jgenr, jhcod, jiadd, jiadm
    integer(kind=8) :: jlong, jlono, jltyp, jluti, jmarq, jorig, jrnom
    integer(kind=8) :: jtype, lorep, lorepa, n, ne, nfic
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!     ------------------------------------------------------------------
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    integer(kind=8) :: nbcla
    common/nficje/nbcla
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
!     ------------------------------------------------------------------
    character(len=32) :: clel, cle, d32, valk(3)
    aster_logical :: linser, rinser
    integer(kind=8) :: iclain, idatin, iin
    data d32/'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/
! DEB ------------------------------------------------------------------
    if (icre .ne. 0) then
        if (nomlu(1:1) .eq. ' ') then
            call utmess('F', 'JEVEUX_19', sk=nomlu)
        end if
    end if
500 continue
    iclasi = iclas
    iret = 0
    linser = .false.
    rinser = .false.
    nfic = nbcla
    lorepa = 0
    do icla = 1, nfic
        if (classe(icla:icla) .ne. ' ') then
            lorep = nrhcod(icla)
            clel = nomlu
            if (lorep .ne. lorepa) then
                iref = jxhcod(clel, lorep)
                lorepa = lorep
            end if
            ne = 1
            i = iref
5           continue
            if (hcod(jhcod(icla)+i) .eq. 0 .and. .not. rinser) then
                if (icre .eq. 1 .or. icre .eq. 2) then
                    if (icla .eq. iclas) then
                        if (nreuti(icla) .ge. nremax(icla)) then
                            call jjarep(icla, 2*nremax(icla))
                            goto 500
                        end if
                        linser = .true.
                        j = nreuti(icla)+1
                        idatin = j
                        iclain = icla
                        iin = i
                        isg = 1
                    end if
                else
                    if (icla .eq. nfic) then
                        iret = 0
                    end if
                end if
            else
                j = hcod(jhcod(icla)+i)
                cle = '!'
                if (j .ne. 0) cle = rnom(jrnom(icla)+abs(j))
                if (cle .eq. clel) then
                    if (icre .eq. 1 .or. icre .eq. 2) then
                        valk(1) = clel
                        valk(2) = nomfic(icla)
                        call utmess('F', 'JEVEUX_10', nk=2, valk=valk)
                    else
                        if (icre .eq. -1 .or. icre .eq. -2) then
                            hcod(jhcod(icla)+i) = -j
                            rnom(jrnom(icla)+j) = '?'
                        end if
                        if (genr(jgenr(icla)+j) .ne. 'X') then
                            iclaos = icla
                            idatos = j
                            iret = 1
                        else
                            iclaco = icla
                            idatco = j
                            iret = 2
                        end if
                        goto 15
                    end if
                else
                    if (j .lt. 0 .and. .not. rinser) then
                        if (icre .eq. 1 .or. icre .eq. 2) then
                            if (icla .eq. iclas) then
                                linser = .true.
                                rinser = .true.
                                idatin = -j
                                iclain = icla
                                iin = i
                                isg = 0
                            end if
                        end if
                    end if
                    if (ne .eq. 1) in = jxhcod(clel, lorep-2)
                    ne = ne+1
                    i = 1+mod(i+in, lorep)
                    if (ne .le. lorep) then
                        j = hcod(jhcod(icla)+i)
                        if (j .eq. 0 .and. rinser) goto 10
                        goto 5
                    else
                        if (icre .eq. 1 .or. icre .eq. 2) then
                            if (icla .eq. iclas) then
                                call jjarep(icla, 2*nremax(icla))
                                lorep = nrhcod(icla)
                                goto 500
                            end if
                        else if (icla .eq. nfic) then
                            iret = 0
                        end if
                    end if
                end if
            end if
        end if
10      continue
    end do
    if (linser) then
        iclas = iclain
        if (icre .eq. 1) then
            iclaos = iclain
            idatos = idatin
            iret = 1
        else if (icre .eq. 2) then
            iclaco = iclain
            idatco = idatin
            iret = 2
        end if
        nreuti(iclas) = nreuti(iclas)+isg
        hcod(jhcod(iclas)+iin) = idatin
        rnom(jrnom(iclas)+idatin) = nomlu
    end if
15  continue
    if (iret .eq. 1) then
        nomos = nomlu
        if (iclasi .ne. iclaos) then
            nomco = d32
            nomoc = d32
        end if
    else if (iret .eq. 2) then
        nomco = nomlu
        nomoc = d32
        if (iclasi .ne. iclaco) then
            nomos = d32
        end if
    end if
! FIN ------------------------------------------------------------------
end subroutine
