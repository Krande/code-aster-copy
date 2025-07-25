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
subroutine jeprat(unit, nomlu, cidatr, mess)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjalty.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjimpo.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: unit
    character(len=*) :: nomlu, cidatr, mess
! ----------------------------------------------------------------------
! ROUTINE D'IMPRESSION DES OBJETS SYSTEME OU DES OBJETS ATTRIBUT DE
! COLLECTION
!
! IN  UNIT  : UNITE LOGIQUE D'IMPRESSION
! IN  NOMLU : NOM DE L'OBJET A IMPRIMER OU NOM DE CLASSE
! IN  CIDATR: NOM DE L'ATTRIBUT
! IN  MESS  : MESSAGE UTILISATEUR
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     -----------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmex, iadmi, ibatr, idatr, ideci, ipgcex, iret2
    integer(kind=8) :: ixatr, jcara, jdate, jdocu, jgenr, jhcod, jiadd
    integer(kind=8) :: jiadm, jlong, jlono, jltyp, jluti, jmarq, jorig
    integer(kind=8) :: jrnom, jtype, k, n
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
!     ------------------------------------------------------------------
    character(len=32) :: noml32, valk(2)
    character(len=1) :: genri, typei
    character(len=8) :: nom
    integer(kind=8) :: icre, iret, jctab, ltypi, lonoi
    integer(kind=8) :: ibacol
    aster_logical :: lcol
!     ------------------------------------------------------------------
    integer(kind=8) :: idnum
    parameter(idnum=10)
    character(len=8) :: cidnom(idnum)
    integer(kind=8) :: idpar
    parameter(idpar=3)
    character(len=8) :: cidpar(idpar)
    integer(kind=8) :: lidbas
    parameter(lidbas=20)
    character(len=8) :: cidbas(lidbas)
    data cidnom/'$$DESO  ', '$$IADD  ', '$$IADM  ', '$$MARQ  ',&
     &               '$$NOM   ', '$$XXXX  ', '$$LONG  ', '$$LONO  ',&
     &               '$$LUTI  ', '$$NUM   '/
    data cidpar/'&&LONO  ', '&&LUTI  ', '&&PART  '/
    data cidbas/'$$CARA  ', '$$IADD  ', '$$GENR  ', '$$TYPE  ',&
     &               '$$DOCU  ', '$$ORIG  ', '$$RNOM  ', '$$LTYP  ',&
     &               '$$LONG  ', '$$LONO  ', '$$DATE  ', '$$LUTI  ',&
     &               '$$HCOD  ', '$$USADI ', '$$ACCE  ', '$$MARQ  ',&
     &               '$$XXXX  ', '$$TLEC  ', '$$TECR  ', '$$IADM  '/
! DEB ------------------------------------------------------------------
    ipgcex = ipgc
    ipgc = -2
    noml32 = nomlu(1:min(24, len(nomlu)))
    nom = cidatr
!
    if (noml32(1:1) .eq. '$') then
        iclas = index(classe, noml32(2:2))
        if (iclas .eq. 0) then
            call utmess('F', 'JEVEUX1_15', sk=noml32(2:2))
        end if
        do k = 1, lidbas
            if (nom .eq. cidbas(k)) then
                idatr = k
                ideci = 0
                iadmi = iadm(jiadm(iclas)+2*idatr-1)
                genri = genr(jgenr(iclas)+idatr)
                typei = type(jtype(iclas)+idatr)
                ltypi = ltyp(jltyp(iclas)+idatr)
                lonoi = lono(jlono(iclas)+idatr)*ltypi
                call jjimpo(unit, iadmi, ideci, 0, genri, &
                            typei, ltypi, lonoi, mess)
                goto 10
            end if
        end do
        call utmess('F', 'JEVEUX1_16', sk=nom)
10      continue
!
    else
        lcol = .false.
        icre = 0
        call jjvern(noml32, icre, iret)
!
        if (iret .eq. 1 .and. nom(1:2) .eq. '&&') then
            do k = 1, idpar
                if (nom .eq. cidpar(k)) then
                    idatr = k
                    goto 20
                end if
            end do
            call utmess('F', 'JEVEUX1_17', sk=nom)
20          continue
            call jjcren(noml32(1:24)//nom, 0, iret2)
            if (iret2 .eq. 0) then
                call utmess('F', 'JEVEUX_26', sk=noml32(1:24))
            end if
            iadmi = iadm(jiadm(iclaos)+2*idatos-1)
            iadmex = iadmi
            genri = genr(jgenr(iclaos)+idatos)
            typei = type(jtype(iclaos)+idatos)
            ltypi = ltyp(jltyp(iclaos)+idatos)
            lonoi = lono(jlono(iclaos)+idatos)*ltypi
            if (iadmex .eq. 0) then
                call jjalty(typei, ltypi, 'L', 1, jctab)
                iadmi = iadm(jiadm(iclaos)+2*idatos-1)
            end if
            ideci = 0
            call jjimpo(unit, iadmi, ideci, 0, genri, &
                        typei, ltypi, lonoi, mess)
            if (iadmex .eq. 0) then
                call jjlide('JEIMPO', noml32(1:24)//nom, 1)
                ipgc = ipgcex
            end if
        else if (iret .ne. 2) then
            valk(1) = nom
            valk(2) = noml32
            call utmess('F', 'JEVEUX1_18', nk=2, valk=valk)
        else
            lcol = .true.
            call jjallc(iclaco, idatco, 'L', ibacol)
            do k = 1, idnum
                if (nom .eq. cidnom(k)) then
                    idatr = k
                    goto 30
                end if
            end do
            call utmess('F', 'JEVEUX1_17', sk=nom)
30          continue
            ixatr = iszon(jiszon+ibacol+idatr)
            if (ixatr .gt. 0) then
                ibatr = iadm(jiadm(iclaco)+2*ixatr-1)
                if (ibatr .eq. 0) then
                    call utmess('F', 'JEVEUX1_19', sk=nom)
                end if
                ideci = 0
                genri = genr(jgenr(iclaco)+ixatr)
                typei = type(jtype(iclaco)+ixatr)
                ltypi = ltyp(jltyp(iclaco)+ixatr)
                lonoi = lono(jlono(iclaco)+ixatr)*ltypi
                call jjimpo(unit, ibatr, ideci, 0, genri, &
                            typei, ltypi, lonoi, mess)
            end if
        end if
        if (lcol) then
            call jjlide('JEIMPO', noml32, 2)
        end if
    end if
    ipgc = ipgcex
! FIN -----------------------------------------------------------------
end subroutine
