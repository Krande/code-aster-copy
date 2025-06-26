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
subroutine jjlide(nomap, nomlu, itype)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jxecro.h"
    character(len=*) :: nomap, nomlu
    integer(kind=8) :: itype
! ----------------------------------------------------------------------
! LIBERATION D'UN OBJET JEVEUX
!
! IN  NOMAP : NOM DE LA ROUTINE APPELANTE (JELIBE,JETASS,JELIBF)
! IN  NOMLU : NOM DE L'OBJET A LIBERER
! IN  ITYPE  : TYPE D'OBJET: 1, 2 OU 3
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: iadit, iadmar, iady1, iady2, iasig, ibacol, ibiadd
    integer(kind=8) :: ibiadm, iblono, ibmarq, ic, icre, idco, idos
    integer(kind=8) :: idyni, ijit, ikit, is, itrold, ix, ixdeso
    integer(kind=8) :: ixiadd, ixiadm, ixlono, ixmarq, jcara, jdate, jdocu
    integer(kind=8) :: jgenr, jhcod, jiadd, jiadm, jit, jlong, jlono
    integer(kind=8) :: jltyp, jluti, jmarq, jorig, jrnom, jtype, k
    integer(kind=8) :: kit, kk, ldynol, lonoi, marqi, n, nadm
    integer(kind=8) :: nalloc, nldo, nmax, nnn
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
! ----------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: datei
    common/iheuje/datei
    integer(kind=8) :: istat
    common/istaje/istat(4)
    integer(kind=8) :: lundef, idebug
    common/undfje/lundef, idebug
    integer(kind=8) :: idinit, idxaxd, itrech, itiad, itcol, lmots, idfr
    common/ixadje/idinit(2), idxaxd(2), itrech, itiad, itcol, lmots, idfr
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
! ----------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idlono, idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4,&
     &               idlono=8, idnum=10)
! ----------------------------------------------------------------------
    integer(kind=8) :: iiadm, iiadd, iidos, iidco, idate, iorig, ilono, iltyp, iimar
    integer(kind=8) :: isauv, nparm
    parameter(iiadm=0, iiadd=1, iidos=2, iidco=3, idate=4,&
     &                 iorig=5, ilono=6, iltyp=7, iimar=8, isauv=9,&
     &                 nparm=10)
! ----------------------------------------------------------------------
    integer(kind=8) :: ista1, ista2, ipgcl
    character(len=4) :: fonc
    character(len=32) :: noml32
    integer(kind=8) :: iadmi, iaddi(2), nld, it(nparm), kt(nparm)
    aster_logical :: lsauv, ldate, lmarq, llibp, ltout, lattr, lxu, lad
! DEB ------------------------------------------------------------------
    do k = 1, nparm
        kt(k) = 0
    end do
    noml32 = nomlu
    iadit = 0
    iasig = 0
    iady1 = 0
    iady2 = 0
    ipgcl = ipgc
    if (nomap .eq. 'JELIBE') then
        fonc = 'LIBE'
        if (idebug .eq. 1) then
            fonc = 'DEBG'
        end if
    else if (nomap .eq. 'JELIBZ') then
        fonc = 'LIBE'
        ipgcl = -1
    else if (nomap .eq. 'JETASS') then
        fonc = 'TASS'
!
!       APPEL A JJLIDE PAR JETASS INVALIDE POUR TYPE =/= 1
        ASSERT(itype .eq. 1)
    else if (nomap .eq. 'JELIBF') then
        fonc = 'LIBF'
    else if (nomap .eq. 'SYSTEM') then
        fonc = 'LIBS'
    else if (nomap .eq. 'JEIMPO' .or. nomap .eq. 'JEIMPA' .or. &
             nomap .eq. 'JENUNO' .or. nomap .eq. 'JENONU') then
        fonc = 'LIBE'
        ipgcl = -2
    else
        goto 101
    end if
!
    nld = 0
    nldo = 1
    jit = 1
    kit = 1
!
! --- CAS D'UN OBJET SIMPLE
!
    if (itype .eq. 1) then
        ic = iclaos
        nld = 1
        it(jit+iiadm) = jiadm(ic)+2*idatos-1
        iadmi = iadm(jiadm(ic)+2*idatos-1)
        if (iadmi .eq. 0) then
            goto 101
        else
            ista1 = iszon(jiszon+iadmi-1)
            is = jiszon+iszon(jiszon+iadmi-4)
            ista2 = iszon(is-4)
            if (ista1 .eq. istat(1) .and. ista2 .eq. istat(3)) goto 101
        end if
        it(jit+iiadd) = jiadd(ic)+2*idatos-1
        it(jit+iidos) = idatos
        it(jit+iidco) = 0
        it(jit+idate) = jdate(ic)+idatos
        it(jit+iorig) = jorig(ic)+idatos
        it(jit+iimar) = jmarq(ic)+2*idatos-1
        it(jit+ilono) = jlono(ic)+idatos
        it(jit+iltyp) = jltyp(ic)+idatos
        it(jit+isauv) = 1
!
! --- CAS D'UNE COLLECTION
!
    else if (itype .eq. 2) then
        ic = iclaco
!
! ----- CAS D'UNE COLLECTION ENTIERE
!
        if (noml32(25:32) .eq. '        ') then
            ibacol = iadm(jiadm(ic)+2*idatco-1)
            if (ibacol .eq. 0) then
                goto 101
            else
                ista1 = iszon(jiszon+ibacol-1)
                is = jiszon+iszon(jiszon+ibacol-4)
                ista2 = iszon(is-4)
                if (ista1 .eq. istat(1) .and. ista2 .eq. istat(3)) goto 101
            end if
            ixiadd = iszon(jiszon+ibacol+idiadd)
            ixdeso = iszon(jiszon+ibacol+iddeso)
            ixiadm = iszon(jiszon+ibacol+idiadm)
            ixmarq = iszon(jiszon+ibacol+idmarq)
            nalloc = idnum
            nmax = 0
!
! ------- CAS D'UNE COLLECTION DISPERSEE
!
            if (ixiadd .ne. 0) then
                nmax = iszon(jiszon+ibacol+ivnmax)
                nalloc = nalloc+nmax
            end if
!
! ------- TOUTE FORME DE COLLECTION
! ---     ALLOCATION EN EVITANT L'APPEL A JJLDYN
! ---          (ET LES APPELS RECURSIFS)
!
            ldynol = ldyn
            if (ldyn .eq. 1) then
                ldyn = 2
            end if
            itrold = itrech
            itrech = 2
            nnn = nalloc*nparm*lois
            call jjalls(nnn, 0, 'V', 'I', lois, &
                        'INIT', it, jit, iadit, iady1)
            iszon(jiszon+iadit-1) = istat(2)
            iszon(jiszon+iszon(jiszon+iadit-4)-4) = istat(4)
            svuse = svuse+(iszon(jiszon+iadit-4)-iadit+4)
            call jjalls(nnn, 0, 'V', 'I', lois, &
                        'INIT', kt, kit, iasig, iady2)
            iszon(jiszon+iasig-1) = istat(2)
            iszon(jiszon+iszon(jiszon+iasig-4)-4) = istat(4)
            svuse = svuse+(iszon(jiszon+iasig-4)-iasig+4)
            smxuse = max(smxuse, svuse)
            itrech = itrold
            ldyn = ldynol
! --------OBJETS DE COLLECTION
!
            if (nmax .ne. 0) then
                ixlono = iszon(jiszon+ibacol+idlono)
                ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
                ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
                ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                do k = 1, nmax
                    nadm = jiszon+ibiadm-1+2*k-1
                    if (iszon(nadm) .ne. 0) then
                        ijit = jit+nld*nparm
                        ikit = kit+nld*nparm
                        nld = nld+1
                        it(ijit+iiadm) = nadm
                        kt(ikit+iiadm) = -1
                        it(ijit+iimar) = jiszon+ibmarq-1+2*k-1
                        kt(ikit+iimar) = -1
                        it(ijit+iiadd) = jiszon+ibiadd-1+2*k-1
                        kt(ikit+iiadd) = -1
                        if (ixlono .ne. 0) then
                            iblono = iadm(jiadm(ic)+2*ixlono-1)
                            it(ijit+ilono) = jiszon+iblono-1+k
                            kt(ikit+ilono) = -1
                        else
                            it(ijit+ilono) = jlono(ic)+ixdeso
                        end if
                        it(ijit+iidos) = k
                        it(ijit+iidco) = idatco
                        it(ijit+idate) = jdate(ic)+ixdeso
                        it(ijit+iorig) = jorig(ic)+ixdeso
                        it(ijit+iltyp) = jltyp(ic)+ixdeso
                        it(ijit+isauv) = 1
                    end if
                end do
            end if
            nldo = nld
!
! --------OBJETS ATTRIBUTS DE COLLECTION
!
            do k = 1, idnum
                ix = iszon(jiszon+ibacol+k)
                if (ix .gt. 0) then
                    if (rnom(jrnom(ic)+ix) (25:26) .eq. '$$' .or. ipgcl .eq. -2) then
!
! ----------- UNIQUEMENT LES OBJETS $$ PRESENTS EN MEMOIRE
! ----------- LES POINTEURS PARTAGES DOIVENT ETRE LIBERES EXPLICITEMENT
!
                        if (iadm(jiadm(ic)+2*ix-1) .ne. 0) then
                            ijit = jit+nld*nparm
                            ikit = kit+nld*nparm
                            nld = nld+1
                            it(ijit+iiadm) = jiadm(ic)+2*ix-1
                            it(ijit+iiadd) = jiadd(ic)+2*ix-1
                            it(ijit+iidos) = ix
                            it(ijit+iidco) = 0
                            it(ijit+idate) = jdate(ic)+ix
                            it(ijit+iorig) = jorig(ic)+ix
                            it(ijit+iimar) = jmarq(ic)+2*ix-1
                            it(ijit+ilono) = jlono(ic)+ix
                            it(ijit+iltyp) = jltyp(ic)+ix
                            if (ix .ne. ixiadm .and. ix .ne. ixmarq) then
                                it(ijit+isauv) = 1
                            else
                                it(ijit+isauv) = 0
                            end if
                        end if
                    end if
                end if
            end do
            ijit = jit+nld*nparm
            ikit = kit+nld*nparm
            nld = nld+1
            it(ijit+iiadm) = jiadm(ic)+2*idatco-1
            it(ijit+iiadd) = jiadd(ic)+2*idatco-1
            it(ijit+iidos) = idatco
            it(ijit+iidco) = 0
            it(ijit+idate) = jdate(ic)+idatco
            it(ijit+iorig) = jorig(ic)+idatco
            it(ijit+iimar) = jmarq(ic)+2*idatco-1
            it(ijit+ilono) = jlono(ic)+idatco
            it(ijit+iltyp) = jltyp(ic)+idatco
            it(ijit+isauv) = 1
!
! ----- CAS D'UN OBJET DE COLLECTION ------
!
        else
            icre = 0
            call jjcroc(noml32(25:32), icre)
            ibacol = iadm(jiadm(ic)+2*idatco-1)
            if (ibacol .eq. 0) goto 101
            ixiadd = iszon(jiszon+ibacol+idiadd)
            ixiadm = iszon(jiszon+ibacol+idiadm)
            ixmarq = iszon(jiszon+ibacol+idmarq)
            ixdeso = iszon(jiszon+ibacol+iddeso)
            ixlono = iszon(jiszon+ibacol+idlono)
!
!         LIBERATION D''UN OBJET DE COLLECTION CONTIGUE REFUSEE
            ASSERT(ixiadd .ne. 0)
            ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
            ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
            ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
            nld = 1
            it(jit+iiadm) = jiszon+ibiadm-1+2*idatoc-1
            kt(kit+iiadm) = -1
            iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
            if (iadmi .ne. 0) then
                ista1 = iszon(jiszon+iadmi-1)
                is = jiszon+iszon(jiszon+iadmi-4)
                ista2 = iszon(is-4)
                if (ista1 .eq. istat(1) .and. ista2 .eq. istat(3)) goto 101
            else
                goto 101
            end if
            it(jit+iiadd) = jiszon+ibiadd-1+2*idatoc-1
            kt(kit+iiadd) = -1
            it(jit+iimar) = jiszon+ibmarq-1+2*idatoc-1
            kt(kit+iimar) = -1
            it(jit+iidos) = idatoc
            it(jit+iidco) = idatco
            it(jit+idate) = jdate(ic)+ixdeso
            it(jit+iorig) = jorig(ic)+ixdeso
            if (ixlono .ne. 0) then
                iblono = iadm(jiadm(ic)+2*ixlono-1)
                it(jit+ilono) = jiszon+iblono-1+idatoc
                kt(kit+ilono) = -1
            else
                it(jit+ilono) = jlono(ic)+ixdeso
            end if
            it(jit+iltyp) = jltyp(ic)+ixdeso
            it(jit+isauv) = 1
        end if
    end if
!
    if (fonc .eq. 'LIBE') then
        if (itiad .eq. 3 .and. itype .ne. 2) then
            if (iadmi .lt. idfr) then
                idxaxd(1) = iadmi-4
            else
                idxaxd(2) = iadmi-4
            end if
        end if
    end if
!
    ltout = .true.
    do k = 1, nld
        lattr = (k .gt. nldo .and. nldo .gt. 0)
        llibp = .false.
        lsauv = .false.
        lmarq = .false.
        ldate = .false.
        lxu = .false.
        lad = .false.
        kk = (k-1)*nparm
!
! ----- PREPARATION AUX DIVERSES OPERATIONS
!
        if (kt(kit+kk+iiadm) .eq. 0) then
            iadmi = iadm(it(jit+kk+iiadm))
            idyni = iadm(it(jit+kk+iiadm)+1)
        else
            iadmi = iszon(it(jit+kk+iiadm))
            idyni = iszon(it(jit+kk+iiadm)+1)
        end if
        if (kt(kit+kk+iimar) .eq. 0) then
            marqi = imarq(it(jit+kk+iimar))
        else
            marqi = iszon(it(jit+kk+iimar))
        end if
!
        if (kt(kit+kk+ilono) .eq. 0) then
            lonoi = lono(it(jit+kk+ilono))*ltyp(it(jit+kk+iltyp))
        else
            lonoi = iszon(it(jit+kk+ilono))*ltyp(it(jit+kk+iltyp))
        end if
!
        idos = it(jit+kk+iidos)
        idco = it(jit+kk+iidco)
        ista1 = iszon(jiszon+iadmi-1)
        is = jiszon+iszon(jiszon+iadmi-4)
        ista2 = iszon(is-4)
        if (ista1 .eq. istat(1) .and. ista2 .eq. istat(1)) then
            goto 100
        end if
!
! ----- OPERATION LIBE
!       ==============
! ----- LES SEGMENTS XA OU XD RESTENT EN L'ETAT
! ----- SI LA MARQUE COURANTE CORRESPOND A CELLE DU SEGMENT DE
! ----- VALEUR :
! -----    POUR LES OBJETS SIMPLES :
! -----    - LES SEGMENTS UA PASSENT EN XA
! -----    - LES SEGMENTS UD PASSENT EN XD
! -----    POUR LES OBJETS ATTRIBUTS DE COLLECTION :
! -----    LES OBJETS ATTRIBUTS DE COLLECTION RESTENT U (UA OU UD)
! -----    SAUF SI IL N'Y A AUCUN OBJET DE COLLECTION EN  MEMOIRE
!
        if (fonc .eq. 'LIBE') then
            if (ista1 .eq. istat(2)) then
                if (ipgcl .eq. marqi) then
                    if (lattr) then
                        if (ltout .and. ista2 .eq. istat(3)) lxu = .true.
                    else
                        ltout = .false.
                        lxu = .true.
                    end if
                    if (ista2 .eq. istat(4) .and. lxu) then
                        ldate = .true.
                    end if
                    lmarq = .true.
                else
                    ltout = .false.
                    goto 100
                end if
            else
                if (.not. lattr) ltout = .false.
                goto 100
            end if
!
! ----- OPERATION TASS
!       ==============
! -----    UNIQUEMENT UTILISEE POUR ECRITURE DANS JETASS (POUR DES
! -----    OBJETS SIMPLES ADRESSE DISQUE DES OBJETS DE COLLECTION)
!
        else if (fonc .eq. 'TASS') then
            lmarq = .false.
            llibp = .false.
            lsauv = .true.
            ldate = .true.
!
! ----- OPERATION DEBG
!       ==============
! -----    UTILISEE POUR ECRITURE IMMEDIATE AVEC L'OPTION DE DEBUG
!
        else if (fonc .eq. 'DEBG') then
            if (ista1 .eq. istat(2)) then
                if (ipgcl .eq. marqi) then
                    if (lattr) then
                        if (ltout .and. ista2 .eq. istat(3)) lxu = .true.
                    else
                        lxu = .true.
                    end if
                    if (ista2 .eq. istat(4) .and. lxu) then
                        lad = .true.
                        lsauv = .true.
                        ldate = .true.
                        llibp = .true.
                    else if (lxu) then
                        llibp = .true.
                    end if
                    lmarq = .true.
                else
                    ltout = .false.
                    goto 100
                end if
            else if (ista2 .eq. istat(4)) then
                if (.not. lattr .or. (lattr .and. ltout)) then
                    lad = .true.
                    lsauv = .true.
                    ldate = .true.
                    llibp = .true.
                end if
            else
                if (.not. lattr) ltout = .false.
                goto 100
            end if
!
! ----- OPERATIONS LIBF
!       ===============
! ----- ON NE S'OCCUPE PAS DE LA MARQUE ASSOCIEE, ON FORCE LA
! ----- LIBERATION
!
        else if (fonc .eq. 'LIBF' .or. fonc .eq. 'LIBS') then
            llibp = .true.
            if (ista1 .eq. istat(2)) then
                lxu = .true.
                if (ista2 .eq. istat(4)) then
                    lad = .true.
                    lsauv = .true.
                    ldate = .true.
                end if
                lmarq = .true.
            else
                if (ista2 .eq. istat(4)) then
                    lad = .true.
                    lsauv = .true.
                else
                    llibp = .true.
                end if
            end if
            if (fonc .eq. 'LIBS') llibp = .false.
        end if
!
! ----- ACTUALISATION DES ATTRIBUTS
!
        if (lmarq) then
            if (kt(kit+kk+iimar) .eq. 0) then
                imarq(it(jit+kk+iimar)) = 0
                iadmar = imarq(it(jit+kk+iimar)+1)
                if (iadmar .ne. 0) then
                    iszon(jiszon+kdesma(1)+iadmar-1) = 0
                    imarq(it(jit+kk+iimar)+1) = 0
                end if
            else
                iszon(it(jit+kk+iimar)) = 0
                iadmar = iszon(it(jit+kk+iimar)+1)
                if (iadmar .ne. 0) then
                    iszon(jiszon+kdesma(1)+iadmar-1) = 0
                    iszon(it(jit+kk+iimar)+1) = 0
                end if
            end if
            if (lxu) then
                iszon(jiszon+iadmi-1) = istat(1)
                svuse = svuse-(iszon(jiszon+iadmi-4)-iadmi+4)
                smxuse = max(smxuse, svuse)
            end if
            if (lad) then
                is = iszon(jiszon+iadmi-4)
                iszon(jiszon+is-4) = istat(3)
            end if
        end if
        if (lsauv .and. it(jit+kk+isauv) .eq. 1) then
            if (kt(kit+kk+iiadd) .eq. 0) then
                iaddi(1) = iadd(it(jit+kk+iiadd))
                iaddi(2) = iadd(it(jit+kk+iiadd)+1)
            else
                iaddi(1) = iszon(it(jit+kk+iiadd))
                iaddi(2) = iszon(it(jit+kk+iiadd)+1)
            end if
            idco = it(jit+kk+iidco)
            idos = it(jit+kk+iidos)
            call jxecro(ic, iadmi, iaddi, lonoi, idco, &
                        idos)
            if (kt(kit+kk+iiadd) .eq. 0) then
                iadd(it(jit+kk+iiadd)) = iaddi(1)
                iadd(it(jit+kk+iiadd)+1) = iaddi(2)
            else
                iszon(it(jit+kk+iiadd)) = iaddi(1)
                iszon(it(jit+kk+iiadd)+1) = iaddi(2)
            end if
        end if
        if (ldate) date(it(jit+kk+idate)) = datei
        if (llibp) then
            call jjlidy(idyni, iadmi)
            if (kt(kit+kk+iiadm) .eq. 0) then
                iadm(it(jit+kk+iiadm)) = 0
                iadm(it(jit+kk+iiadm)+1) = 0
            else
                iszon(it(jit+kk+iiadm)) = 0
                iszon(it(jit+kk+iiadm)+1) = 0
            end if
        end if
100     continue
    end do
101 continue
    call jjlidy(iady1, iadit)
    call jjlidy(iady2, iasig)
! FIN ------------------------------------------------------------------
end subroutine
