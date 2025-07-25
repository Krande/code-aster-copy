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
subroutine jetass(clas)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jesvos.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jxecrb.h"
#include "asterfort/jxecro.h"
#include "asterfort/jxlibd.h"
#include "asterfort/jxliro.h"
!
!
    character(len=1), intent(in) :: clas
! ----------------------------------------------------------------------
! COMPRESSION D'UNE BASE DE DONNEES PAR RECUPERATION DES ENREGISTREMENTS
! LIBERES LORS D'UNE DESTRUCTION
!
! IN  CLAS   : NOM DE CLASSE ASSOCIEE
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     ------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
    integer(kind=8) :: istat
    common/istaje/istat(4)
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ibacol, ibiadd, ic, idco, idcol, idcop, idec
    integer(kind=8) :: idos, idosl, idosp, ixiadd, jcara, jdate
    integer(kind=8) :: jdocu, jgenr, jhcod, jiadd, jiadm, jlong, jlono
    integer(kind=8) :: jltyp, jluti, jmarq, jorig, jrnom, jtype, jusadi
    integer(kind=8) :: k, kadd, klib, ladd, ld, lgl, n
    integer(kind=8) :: ncla1, ncla2
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    aster_logical :: litlec
    common/lficje/litlec(n)
    common/jusadi/jusadi(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    real(kind=8) :: svuse, smxuse
    common/statje/svuse, smxuse
!     ------------------------------------------------------------------
    integer(kind=8) :: idiadd
    parameter(idiadd=2)
!     ------------------------------------------------------------------
    aster_logical :: libre, actu
    character(len=1) :: kclas
    integer(kind=8) :: itp(1), jitp, iaditp, iaddi(2), iaddib(2), lgbl, iadyn
! DEB ------------------------------------------------------------------
    iaddi(2) = 0
    iaddib(2) = 0
    kclas = clas
!
!   ON COMMENCE PAR S'ASSURER QUE TOUS LES OBJETS POSSEDENT UNE IMAGE DISQUE
!
    call jesvos(kclas)
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    do ic = ncla1, ncla2
        lgbl = 1024*longbl(ic)*lois
        call jjalls(lgbl, 0, 'V', 'I', lois, &
                    'INIT', itp, jitp, iaditp, iadyn)
        iszon(jiszon+iaditp-1) = istat(2)
        iszon(jiszon+iszon(jiszon+iaditp-4)-4) = istat(4)
        svuse = svuse+(iszon(jiszon+iaditp-4)-iaditp+4)
!
! ----- BOUCLE "TANT QUE" SUR LES ENREGISTREMENTS UTILISES
!
        k = 1
        klib = 0
        idosp = 0
        idcop = 0
200     continue
!
! ----- DECHARGEMENT DES TAMPONS DE LECTURE ET D'ECRITURE
! ----- AFIN D'ACTUALISER LES ADRESSES DISQUES DES COLLECTIONS
! ----- STOCKEES DANS DES PETITS OBJETS
!
        if (iitecr(ic) .gt. 0) then
            call jxecrb(ic, iitecr(ic), kitecr(ic)+1, lgbl, 0, &
                        0)
            iitecr(ic) = 0
            nitecr(ic) = 0
        end if
        if (litlec(ic)) then
            call jxecrb(ic, iitlec(ic), kitlec(ic)+1, lgbl, 0, &
                        0)
            litlec(ic) = .false.
            iitlec(ic) = 0
        end if
        k = k+1
        if (k .le. nbluti(ic)) then
            libre = iusadi(jusadi(ic)+3*k-2) .lt. -1 .or. iusadi(jusadi(ic)+3*k-1) .lt. -1
!
! ------- "GROS" OBJET DETRUIT
            if (libre) then
                if (klib .eq. 0) then
                    klib = k
!
! ------- ON POSITIONNE LES INDICATEURS POUR EVITER D'UTILISER
!-------- L'ENREGISTREMENT DANS JXECRO
                    iusadi(jusadi(ic)+3*klib-2) = 0
                    iusadi(jusadi(ic)+3*klib-1) = 0
                end if
            else
!
! --------- ENREGISTREMENT A DEPLACER
                if (klib .ne. 0) then
                    idco = iusadi(jusadi(ic)+3*k-2)
                    idos = iusadi(jusadi(ic)+3*k-1)
                    if (idos .gt. 0 .or. idco .gt. 0) then
!
! ----------- L'ENREGISTREMENT CONTIENT UN OU UNE PARTIE D'UN "GROS"
! ----------- OBJET ON ACTUALISE L'ADRESSE DISQUE ET ON REECRIT
!
                        iaddi(1) = k
                        call jxliro(ic, iaditp, iaddi, lgbl)
                        idosl = idos
                        idcol = idco
                        iaddib(1) = klib
                        call jxecro(ic, iaditp, iaddib, lgbl, idco, &
                                    idos)
                        if (idosl .ne. idosp .or. idcol .ne. idcop) then
!
! ------------- ON ACTUALISE L'ADRESSE DISQUE QUE SI L'UN DES
! ------------- IDENTIFICATEURS A ETE MODIFIE
!
                            idosp = idosl
                            idcop = idcol
                            if (idco .eq. 0) then
                                iadd(jiadd(ic)+2*idos-1) = klib
                            else
                                call jjallc(ic, idco, 'E', ibacol)
                                ixiadd = iszon(jiszon+ibacol+idiadd)
                                if (ixiadd .gt. 0) then
                                    ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                                    iszon(jiszon+ibiadd-1+2*idos-1) = &
                                        klib
                                end if
!
! ----------------- MISE A JOUR DU COMMON /IATCJE/ POUR APPEL JJLIDE
                                iclas = ic
                                iclaos = ic
                                idatos = ixiadd
                                nomos = rnom(jrnom(ic)+ixiadd)
                                call jjlide('JETASS', rnom(jrnom(ic)+ixiadd), 1)
                            end if
                        end if
                        call jxlibd(idcol, idosl, ic, iaddi, lois)
                        klib = min(klib+1, k)
                        iusadi(jusadi(ic)+3*klib-2) = 0
                        iusadi(jusadi(ic)+3*klib-1) = 0
!
! ----------- L'ENREGISTREMENT CONTIENT DES PETITS OBJETS
! ----------- ON ACTUALISE LES ADRESSES DISQUE ET ON REECRIT
! ----------- IL N'Y A PAS DE RETASSAGE AU SEIN DE L'ENREGISTREMENT
! ----------- ON EXPLORE L'ENREGISTREMENT
!
                    else if (idco .eq. 0 .and. idos .eq. 0) then
                        iaddi(1) = k
                        call jxliro(ic, iaditp, iaddi, lgbl)
                        actu = .false.
                        idec = 0
300                     continue
                        idcol = iszon(jiszon+iaditp+idec)
                        idosl = iszon(jiszon+iaditp+idec+1)
                        lgl = iszon(jiszon+iaditp+idec+2)
                        if (idcol .eq. 0 .and. idosl .eq. 0) then
                            goto 350
                        else if (idcol .lt. 0 .or. idosl .lt. 0) &
                            then
                            goto 320
                        end if
                        actu = .true.
                        if (idcol .eq. 0) then
                            iadd(jiadd(ic)+2*idosl-1) = klib
                        else
                            call jjallc(ic, idcol, 'E', ibacol)
                            ixiadd = iszon(jiszon+ibacol+idiadd)
                            kadd = iadd(jiadd(ic)+2*ixiadd-1)
                            if (ixiadd .gt. 0) then
                                ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                                ladd = iszon(jiszon+ibiadd-1+2*idosl)
                                iszon(jiszon+ibiadd-1+2*idosl-1) = &
                                    klib
                                if (kadd .ne. k) then
                                    iclas = ic
                                    iclaos = ic
                                    idatos = ixiadd
                                    nomos = rnom(jrnom(ic)+ixiadd)
                                    call jjlide('JETASS', rnom(jrnom(ic)+ixiadd), 1)
                                else
                                    ld = iadd(jiadd(ic)+2*ixiadd)
                                    iszon(jiszon+iaditp+(ld/lois)-1+2* &
                                          idosl-1) = klib
                                    iszon(jiszon+iaditp+(ld/lois)-1+2* &
                                          idosl) = ladd
                                end if
                            end if
                        end if
320                     continue
                        idec = idec+lgl+3
                        goto 300
350                     continue
                        if (actu) then
                            call jxlibd(idcol, idosl, ic, iaddi, lois)
                            iaddib(1) = klib
                            call jxecro(ic, iaditp, iaddib, lgbl, idco, &
                                        idos)
                            iusadi(jusadi(ic)+3*klib) = iusadi(jusadi(ic)+3*k)
                            klib = min(klib+1, k)
                            iusadi(jusadi(ic)+3*klib-2) = 0
                            iusadi(jusadi(ic)+3*klib-1) = 0
                        end if
                    end if
                end if
            end if
            goto 200
        end if
        if (klib .gt. 0) then
            do k = klib, nbluti(ic)
                iusadi(jusadi(ic)+3*k-2) = -1
                iusadi(jusadi(ic)+3*k-1) = -1
                iusadi(jusadi(ic)+3*k) = 0
            end do
            nbluti(ic) = klib-1
        end if
        call jjlidy(iadyn, iaditp)
    end do
! FIN ------------------------------------------------------------------
end subroutine
