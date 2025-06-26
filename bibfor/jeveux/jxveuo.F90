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

subroutine jxveuo(cel, itab, inat, jitab)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjlirs.h"
#include "asterfort/jjprem.h"
#include "asterfort/jxliro.h"
#include "asterfort/jxlocs.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: itab(*), inat, jitab
    character(len=*) :: cel
! ----------------------------------------------------------------------
! MISE EN MEMOIRE D'UN SEGMENT DE VALEUR
! ROUTINE AVEC ADHERENCE SYSTEME CRAY : SHIFTR AND
!
! IN  CEL    : ACCES : 'L' OU 'E'
! IN  ITAB   : TABLEAU PAR RAPPORT AUQUEL L'ADRESSE EST CALCULEE
! IN  INAT   : TYPE D'OBJET 1:OS, 2:CO, 3:OC
! OUT JITAB  : ADRESSE PAR RAPPORT A ITAB
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadyn, ibacol, ibiadd, ibiadm, iblong, iblono, ibluti
    integer(kind=8) :: ibmarq, ic, idco, idos, ista, iusa, ixdeso
    integer(kind=8) :: ixiadd, ixiadm, ixlong, ixlono, ixluti, ixmarq, jcara
    integer(kind=8) :: jdate, jdocu, jgenr, jhcod, jiadd, jiadm, jlong
    integer(kind=8) :: jlono, jltyp, jluti, jmarq, jorig, jrnom, jtype
    integer(kind=8) :: k, longj, lonoj, lonok, lonti, lutilo, n
!
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: idehc
    parameter(idehc=6)
! ----------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    integer(kind=8) :: istat
    common/istaje/istat(4)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
! ----------------------------------------------------------------------
    character(len=1) :: typei, genri
    integer(kind=8) :: ltypi, iaddi(2), iadmi, lonoi, irt
    aster_logical :: ldeps, lconst
! ----------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idlong, idlono, idluti
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4, idlong=7,&
     &               idlono=8, idluti=9)
! DEB ------------------------------------------------------------------
    jitab = 0
    irt = 0
    select case (inat)
!
! ----- INAT =  1 : OBJET SIMPLE
!
    case (1)
        ic = iclaos
        idos = idatos
        idco = 0
        ixdeso = idatos
        genri = genr(jgenr(ic)+idos)
        typei = type(jtype(ic)+idos)
        ltypi = ltyp(jltyp(ic)+idos)
        lonoi = lono(jlono(ic)+idos)*ltypi
        iadmi = iadm(jiadm(ic)+2*idos-1)
        iadyn = iadm(jiadm(ic)+2*idos)
        iaddi(1) = iadd(jiadd(ic)+2*idos-1)
        iaddi(2) = iadd(jiadd(ic)+2*idos)
!
! ----- INAT = 2 : COLLECTION
!
    case (2)
        ic = iclaco
        ibacol = iadm(jiadm(ic)+2*idatco-1)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        idos = ixdeso
        idco = 0
        genri = genr(jgenr(ic)+ixdeso)
        typei = type(jtype(ic)+ixdeso)
        ltypi = ltyp(jltyp(ic)+ixdeso)
        ixlong = iszon(jiszon+ibacol+idlong)
        lconst = (ixlong .eq. 0)
        ixlono = iszon(jiszon+ibacol+idlono)
        ixluti = iszon(jiszon+ibacol+idluti)
        if (lconst) then
            if (long(jlong(ic)+ixdeso) .ne. 0) then
                lonoi = lono(jlono(ic)+ixdeso)*ltypi
            else
                call utmess('F', 'JEVEUX1_62')
            end if
        else
            iblono = iadm(jiadm(ic)+2*ixlono-1)
            iblong = iadm(jiadm(ic)+2*ixlong-1)
            lonti = lono(jlono(ic)+ixdeso)
            lutilo = luti(jluti(ic)+ixlono)
            if (lutilo .eq. 0) then
                k = 1
                iszon(jiszon+iblono-1+k) = 1
5               continue
                if (k .le. iszon(jiszon+ibacol+ivnmax)) then
                    longj = iszon(jiszon+iblong-1+k)
                    if (longj .gt. 0) then
                        if (genri .eq. 'V') then
                            lonoj = longj
                        else if (genri .eq. 'N') then
                            lonok = (idehc+jjprem(longj, irt))*lois+(longj+1)*ltypi
                            if (mod(lonok, ltypi) .gt. 0) then
                                lonok = (lonok/ltypi+1)
                            else
                                lonok = lonok/ltypi
                            end if
                            lonoj = lonok
                            ibluti = iadm(jiadm(ic)+2*ixluti-1)
                            iszon(jiszon+ibluti-1+k) = 0
                        end if
                    else
                        lonoj = 0
                    end if
                    iszon(jiszon+iblono-1+k+1) = lonoj+iszon(jiszon+iblono- &
                                                             1+k)
                    k = k+1
                    lutilo = lutilo+1
                    goto 5
                end if
                luti(jluti(ic)+ixlono) = lutilo
            end if
            if (lonti .ne. 0) then
                lonoi = lonti*ltypi
            else
                lonti = iszon(jiszon+iblono-1+lutilo+1)-1
                lonoi = ltypi*lonti
                lono(jlono(ic)+ixdeso) = lonoi/ltypi
                luti(jluti(ic)+ixdeso) = lutilo
            end if
        end if
        iadmi = iadm(jiadm(ic)+2*ixdeso-1)
        iadyn = iadm(jiadm(ic)+2*ixdeso)
        iaddi(1) = iadd(jiadd(ic)+2*ixdeso-1)
        iaddi(2) = iadd(jiadd(ic)+2*ixdeso)
        if (iadmi .eq. 0) then
            if (iaddi(1) .ne. 0) then
                call jjalls(lonoi, ic, genri, typei, ltypi, &
                            'NOINIT', itab, jitab, iadmi, iadyn)
                call jxliro(ic, iadmi, iaddi, lonoi)
            else
                call jjalls(lonoi, ic, genri, typei, ltypi, &
                            'INIT  ', itab, jitab, iadmi, iadyn)
            end if
            iadm(jiadm(ic)+2*ixdeso-1) = iadmi
            iadm(jiadm(ic)+2*ixdeso) = iadyn
            call jjecrs(iadmi, ic, ixdeso, 0, cel, &
                        imarq(jmarq(ic)+2*ixdeso-1))
        end if
!
! ----- INAT = 3 : OBJET DE COLLECTION
!
    case (3)
        ic = iclaco
        idco = idatco
        idos = idatoc
        ibacol = iadm(jiadm(ic)+2*idatco-1)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        genri = genr(jgenr(ic)+ixdeso)
        typei = type(jtype(ic)+ixdeso)
        ltypi = ltyp(jltyp(ic)+ixdeso)
        ixlono = iszon(jiszon+ibacol+idlono)
        if (ixlono .eq. 0) then
            lonoi = lono(jlono(ic)+ixdeso)*ltypi
        else
            iblono = iadm(jiadm(ic)+2*ixlono-1)
            lonoi = iszon(jiszon+iblono-1+idatoc)*ltypi
        end if
        ixiadm = iszon(jiszon+ibacol+idiadm)
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixmarq = iszon(jiszon+ibacol+idmarq)
        ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
        ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
        ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
        iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
        iadyn = iszon(jiszon+ibiadm-1+2*idatoc)
        iaddi(1) = iszon(jiszon+ibiadd-1+2*idatoc-1)
        iaddi(2) = iszon(jiszon+ibiadd-1+2*idatoc)
!
    end select
!
    if (iadmi .eq. 0) then
!
! ----- PAS DE SEGMENT EN MEMOIRE
!
        if (iaddi(1) .eq. 0) then
!
! ------- PAS D'IMAGE DISQUE
!
            if (cel .eq. 'E') then
                call jjalls(lonoi, ic, genri, typei, ltypi, &
                            'INIT', itab, jitab, iadmi, iadyn)
            else
                call utmess('F', 'JEVEUX1_61')
            end if
        else
!
! ------- AVEC  IMAGE DISQUE
!
            call jjalls(lonoi, ic, genri, typei, ltypi, &
                        'NOINIT', itab, jitab, iadmi, iadyn)
            call jxliro(ic, iadmi, iaddi, lonoi)
        end if
    else
!
! ----- SEGMENT EN MEMOIRE
!
        call jjlirs(iadmi, ic, idos, iusa, ista)
        ldeps = .false.
        if (iusa .ne. istat(2)) ldeps = .true.
        call jxlocs(itab, genri, ltypi, lonoi, iadmi, &
                    ldeps, jitab)
    end if
!
    if (inat .eq. 3) then
        iszon(jiszon+ibiadm-1+2*idatoc-1) = iadmi
        iszon(jiszon+ibiadm-1+2*idatoc) = iadyn
        iszon(jiszon+ibiadd-1+2*idatoc-1) = iaddi(1)
        iszon(jiszon+ibiadd-1+2*idatoc) = iaddi(2)
        call jjecrs(iadmi, ic, idos, idco, cel, &
                    iszon(jiszon+ibmarq-1+2*idatoc-1))
    else
        iadm(jiadm(ic)+2*ixdeso-1) = iadmi
        iadm(jiadm(ic)+2*ixdeso) = iadyn
        iadd(jiadd(ic)+2*ixdeso-1) = iaddi(1)
        iadd(jiadd(ic)+2*ixdeso) = iaddi(2)
        call jjecrs(iadmi, ic, idos, 0, cel, &
                    imarq(jmarq(ic)+2*idos-1))
    end if
! FIN ------------------------------------------------------------------
end subroutine
