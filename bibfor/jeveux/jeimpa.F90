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
subroutine jeimpa(unit, nomlu, com)
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jelira.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: unit
    character(len=*) :: nomlu, com
!
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
!
    integer(kind=8) :: numec
    common/inumje/numec
    character(len=24) :: nomec
    common/knomje/nomec
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibacol, ic, icol, id, ilon, ipgcex
    integer(kind=8) :: ixiadd, ixlong, j, jcol, jdocu, jgenr, jlon
    integer(kind=8) :: jorig, jrnom, jtype, k, n, nnac, nnaci
    integer(kind=8) :: nnao
!-----------------------------------------------------------------------
    parameter(n=5)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    character(len=8) :: nume, nome
!
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!     -----------------------------------------------------------------
    integer(kind=8) :: iddeso, idiadd, idlong
    parameter(iddeso=1, idiadd=2,&
     &                            idlong=7)
!
    character(len=72) :: coml
    character(len=32) :: noml32
    character(len=33) :: cval
    integer(kind=8) :: icre, iret, ival
!
    parameter(nnao=15)
    parameter(nnac=6)
    character(len=8) :: nac(nnac), nao(nnao)
    character(len=1) :: tac(nnac), tao(nnao), genri
    integer(kind=8) :: lac(nnac), lao(nnao)
    aster_logical :: tab1(5, 4), tab2(5, 4), tab3(2, 3)
    aster_logical :: lconst, lconti, lcol
!
    data nume, nome&
     &               /'$$XNUM  ', '$$XNOM  '/
    data nao/&
     &     'CLAS    ', 'GENR    ', 'TYPE    ', 'LTYP    ',&
     &     'DOCU    ', 'DATE    ', 'LONMAX  ',&
     &     'NOMMAX  ', 'LONUTI  ', 'NOMUTI  ', 'LONO    ',&
     &     'IADM    ', 'IADD    ', 'LADD    ', 'USAGE   '/
    data nac/'ACCES   ', 'STOCKAGE', 'MODELONG',&
     &                 'NMAXOC  ', 'NUTIOC  ', 'LONT    '/
    data tao/&
     &     'K', 'K', 'K', 'I',&
     &     'K', 'I', 'I',&
     &     'I', 'I', 'I', 'I',&
     &     'I', 'I', 'I', 'K'/
    data tac/&
     &  'K', 'K', 'K', 'I', 'I', 'I'/
    data lao/&
     &     1, 1, 1, 0,&
     &     4, 0, 0,&
     &     0, 0, 0, 0,&
     &     0, 0, 0, 3/
    data lac/&
     &    33, 8, 33, 0, 0, 0/
! 1 : CONT CSTE - 2 : DISP CSTE - 3 : CONT VARI - 4 : DISP VARI
!     - IRET = 3 - CONDITION D'ACCES A LONO / IADM / IADD / LADD / USAGE
    data((tab1(i, j), i=1, 5), j=1, 4)/&
     &     .false., .false., .false., .false., .true.,&
     &     .true., .true., .true., .true., .true.,&
     &     .false., .false., .false., .false., .true.,&
     &     .true., .true., .true., .true., .true./
!     - IRET = 2 - CONDITION D'ACCES A LONO / IADM / IADD / LADD / USAGE
    data((tab2(i, j), i=1, 5), j=1, 4)/&
     &     .true., .true., .true., .true., .true.,&
     &     .false., .false., .false., .false., .false.,&
     &     .true., .true., .true., .true., .true.,&
     &     .false., .false., .false., .false., .false./
!     ------------------- CONDITION D'ACCES A LON... / NOM...
    data((tab3(i, j), i=1, 2), j=1, 3)/&
     &     .true., .false.,&
     &     .true., .false.,&
     &     .false., .true./
! DEB -----------------------------------------------------------------
    ipgcex = ipgc
    ipgc = -2
    noml32 = nomlu
    coml = com
    icre = 0
    iret = 0
    jcol = 1
    ilon = 1
    lconst = .false.
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_26', sk=noml32(1:24))
    else if (iret .eq. 1) then
        lcol = .false.
        ic = iclaos
        id = idatos
    else if (iret .eq. 2) then
        ic = iclaco
        lcol = .true.
        call jjallc(iclaco, idatco, 'L', ibacol)
        id = iszon(jiszon+ibacol+iddeso)
        ixlong = iszon(jiszon+ibacol+idlong)
        lconst = ixlong .eq. 0
        ixiadd = iszon(jiszon+ibacol+idiadd)
        lconti = ixiadd .eq. 0
        if (.not. lconti .and. lconst) jcol = 2
        if (lconti .and. .not. lconst) jcol = 3
        if (.not. lconti .and. .not. lconst) jcol = 4
        if (noml32(25:32) .ne. '        ') then
            call jjcroc(noml32(25:32), icre)
            iret = 3
        end if
    end if
    genri = genr(jgenr(ic)+id)
    jlon = 1
    if (genri .eq. 'V') jlon = 2
    if (genri .eq. 'N') jlon = 3
!
    write (unit, '(A)') 'JEIMPA  IMPRESSION DES ATTRIBUTS DE >'&
     &                  //noml32(1:24)//'<'
    write (unit, '(A1,A72)') ' ', coml
    if (iret .eq. 3) then
        if (noml32(25:32) .eq. nome) then
            write (unit, '(A,A8)') 'NOM OC', nomec
        else if (noml32(25:32) .eq. nume) then
            write (unit, '(A,I12)') 'NUM OC', numec
        end if
    end if
    if (iret .eq. 2) then
        nnaci = nnac
        if (.not. lconti) nnaci = nnac-1
        do k = 1, nnaci
            call jelira(noml32, nac(k), ival, cval)
            if (tac(k) .eq. 'I') then
                write (unit, '(A8,I12)') nac(k), ival
            else
                write (unit, '(A8,A)') nac(k), cval(1:lac(k))
            end if
        end do
    end if
    do k = 1, nnao
        icol = k-10
        if (nao(k) (1:3) .eq. 'LON') ilon = 1
        if (nao(k) (1:3) .eq. 'NOM') ilon = 2
        if ((k .le. 6) .or. (k .gt. 6 .and. k .le. 10 .and. lconst .and. tab3(ilon, jlon)) .or. &
          ((iret .eq. 1 .or. iret .eq. 3) .and. (k .gt. 6 .and. k .le. 10) .and. tab3(ilon, jlon)) &
            .or. (iret .eq. 2 .and. (k .gt. 10 .and. tab2(icol, jcol))) .or. &
            (iret .eq. 3 .and. (k .gt. 10 .and. tab1(icol, jcol))) .or. &
            (iret .eq. 1 .and. (k .gt. 10))) then
            call jelira(noml32, nao(k), ival, cval)
            if (tao(k) .eq. 'I') then
                write (unit, '(A8,I12)') nao(k), ival
            else
                write (unit, '(A8,A)') nao(k), cval(1:lao(k))
            end if
        end if
    end do
    if (lcol) then
        call jjlide('JEIMPA', noml32(1:24), 2)
    end if
    ipgc = ipgcex
end subroutine
