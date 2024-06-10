! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine jedema()
    implicit none
    ! ----------------------------------------------------------------------
    ! DECREMENTE LA MARQUE N ET LIBERE LES OBJETS MARQUES PAR N
    ! ON OPERE EN DEUX TEMPS : -1 LIBERATION DES COLLECTION PAR JJLIDE
    !                          -2 LIBERATION DES OBJETS SIMPLES
    ! ----------------------------------------------------------------------
#include "jeveux_private.h"
#include "asterfort/jjlide.h"
#include "asterfort/utmess.h"
    integer :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    !-----------------------------------------------------------------------
    integer :: jcara, jdate, jdocu, jgenr, jhcod, jiadd
    integer :: jiadm, jlong, jlono, jltyp, jluti, jmarq, jorig
    integer :: jrnom, jtype, n
    !-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
            &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    !
    integer :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer :: istat
    common/istaje/istat(4)
    integer :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
    integer :: lundef, idebug
    common/undfje/lundef, idebug
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
    ! ----------------------------------------------------------------------
    integer :: k, iadmi, ideb, ifin, idos, idco, ic, is
    character(len=8) :: ksuf
    character(len=24) :: d24
    data d24/'$$$$$$$$$$$$$$$$$$$$$$$$'/
    ! Initialize arrays to catch uninitialized memory issues
    iszon = 0
    print *, 'Entering jedema'
    print *, 'ipgc:', ipgc
    if (ipgc .eq. 0) then
        call utmess('F', 'JEVEUX_06')
    else
        ideb = iszon(jiszon+kposma(1)+ipgc-1)
        print *, 'ideb:', ideb
    end if
    ifin = lgduti-1
    print *, 'ifin:', ifin
    !
    ! Print array contents to check for valid values
    print *, 'Array contents before processing:'
    do k = 1, size(iszon)
        print *, 'iszon(', k, '):', iszon(k)
    end do
    !
    ! --- ON TRAITE D'ABORD LES COLLECTIONS
    !
    do k = ideb, ifin
        iadmi = iszon(jiszon+kdesma(1)+k)
        print *, 'iadmi:', iadmi, 'at k:', k
        if (iadmi .ne. 0) then
            idos = iszon(jiszon+iadmi-2)
            is = iszon(jiszon+iadmi-4)
            idco = iszon(jiszon+is-3)
            ic = iszon(jiszon+is-2)
            print *, 'idos:', idos, 'is:', is, 'idco:', idco, 'ic:', ic
            if (idco .gt. 0) then
                iclaco = ic
                idatco = idco
                nomco = d24
                call jjlide('JELIBE', rnom(jrnom(ic)+idco), 2)
            else if (idos .gt. 0) then
                if (genr(jgenr(ic)+idos) .eq. 'X') then
                    iclaco = ic
                    idatco = idos
                    nomco = d24
                    call jjlide('JELIBE', rnom(jrnom(ic)+idos) (1:24), 2)
                end if
            end if
        end if
    end do
    !
    ! --- ON TRAITE MAINTENANT LES OBJETS SIMPLES
    ! --- ET LE $$DESO DES COLLECTIONS CONTIGUES
    !
    do k = ideb, ifin
        iadmi = iszon(jiszon+kdesma(1)+k)
        print *, 'Processing simple objects, iadmi:', iadmi, 'at k:', k
        if (iadmi .ne. 0) then
            idos = iszon(jiszon+iadmi-2)
            is = iszon(jiszon+iadmi-4)
            idco = iszon(jiszon+is-3)
            ic = iszon(jiszon+is-2)
            print *, 'idos:', idos, 'is:', is, 'idco:', idco, 'ic:', ic
            if (idco .gt. 0) goto 200
            if (idos .gt. 0) then
                ksuf = rnom(jrnom(ic)+idos) (25:32)
                print *, 'ksuf:', ksuf
                if ((ksuf(1:2) .eq. '$$' .and. ksuf(3:6) .ne. 'DESO') .or. ksuf(1:2) .eq. '&&' &
                        .or. genr(jgenr(ic)+idos) .eq. 'X') goto 200
                if (idebug .eq. 1) then
                    iclaos = ic
                    idatos = idos
                    nomos = d24
                    call jjlide('JELIBE', rnom(jrnom(ic)+idos) (1:24), 1)
                else
                    imarq(jmarq(ic)+2*idos-1) = 0
                    imarq(jmarq(ic)+2*idos) = 0
                    iszon(jiszon+iadmi-1) = istat(1)
                    iszon(jiszon+kdesma(1)+k) = 0
                    svuse = svuse-(iszon(jiszon+iadmi-4)-iadmi+4)
                    smxuse = max(smxuse, svuse)
                end if
                end if
            end if
            200     continue
    end do
    lgputi = lgputi-1
    lgduti = ideb
    iszon(jiszon+kposma(1)+ipgc-1) = 0
    ipgc = ipgc-1
    print *, 'Exiting jedema'
    ! FIN ------------------------------------------------------------------
end subroutine
