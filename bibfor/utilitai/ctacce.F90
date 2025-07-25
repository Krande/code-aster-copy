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
subroutine ctacce(nsymb, typac, nbval, nival, nrval, &
                  niord, nkcha, resu)
    implicit none
#include "jeveux.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbval
    character(len=24) :: nival, nrval, niord, nkcha
    character(len=8) :: typac, resu
    character(len=16) :: nsymb
!     ----- OPERATEUR CREA_TABLE , MOT-CLE FACTEUR RESU   --------------
!
!        BUT : RECUPERER LES ACCES DE LA SD RESULTATS
!
!        IN/OUT : NIVAL (K24): OBJET DES VALEURS D'ACCES (ENTIERS)
!                 NRVAL (K24): OBJET DES VALEURS D'ACCES (REELS)
!                 NIORD (K24): OBJET DES NUMEROS D'ORDRE
!                 NKCHA (K24): OBJET DES NOMS DE CHAMP
!           OUT : RESU  (K8) : NOM DU RESULTAT (SI RESULTAT, SINON ' ')
!                 TYPAC (K8) : ACCES (ORDRE,MODE,FREQ,INST)
!                 NSYMB (K16): NOM SYMBOLIQUE DU CHAMP
!                 NBVAL (I)  : NOMBRE DE VALEURS D'ACCES
!
! ----------------------------------------------------------------------
    character(len=16) :: concep
    character(len=8) :: k8b
    integer(kind=8) :: n0, n1, jkcha, jrval, jival, jniord, nbto, nbno, nblo, nbni, nbli
    integer(kind=8) :: nbnm, nblm, nbnf, nblf, nbis, nbr8, jlist, ibid, nbtrou, i
    integer(kind=8) :: vali, n2, jinst, kk, nuord, tord(1)
    real(kind=8) :: r8b, epsi, valr, rinst
    complex(kind=8) :: cbid
    character(len=8) :: crit
    character(len=16) :: valk
    character(len=24) :: nlist
!     ------------------------------------------------------------------
    call jemarq()
!
    call getvid('RESU', 'CHAM_GD', iocc=1, nbval=0, nbret=n0)
    call getvid('RESU', 'RESULTAT', iocc=1, nbval=0, nbret=n1)
!
! =============================================================
! -1- CAS: CHAM_GD
! =============================================================
!
    if (n0 .ne. 0) then
!
        call wkvect(nkcha, 'V V K24', 1, jkcha)
        call getvid('RESU', 'CHAM_GD', iocc=1, scal=zk24(jkcha), nbret=n0)
        call wkvect(nrval, 'V V R', 1, jrval)
        zr(jrval) = 0.0d0
        call wkvect(nival, 'V V I', 1, jival)
        zi(jival) = 0
        call wkvect(niord, 'V V I', 1, jniord)
        zi(jniord) = 0
        nbval = 1
        typac = ' '
        nsymb = ' '
        resu = ' '
!
! =============================================================
! -2- CAS: RESULTAT/NOM_CHAM
! =============================================================
    else if (n1 .ne. 0) then
!
! --- 2.1- ON DETERMINE :
!     ------------------
!         NBVAL = NOMBRE DE VALEURS D'ACCES
!         TYPAC = TYPE D'ACCES (ORDRE,INST,FREQ,MODE,TOUT)
!         NIVAL = TABLEAU DES VALEURS D'ACCES (ENTIERES)
!         NRVAL = TABLEAU DES VALEURS D'ACCES (REELLES)
!         NKCHA = TABLEAU DES NOMS DE CHAMP
!
        call getvid('RESU', 'RESULTAT', iocc=1, scal=resu, nbret=n1)
!
        call getvr8('RESU', 'PRECISION', iocc=1, scal=epsi, nbret=n1)
        call getvtx('RESU', 'CRITERE', iocc=1, scal=crit, nbret=n1)
        call getvtx('RESU', 'TOUT_ORDRE', iocc=1, nbval=0, nbret=nbto)
        call getvis('RESU', 'NUME_ORDRE', iocc=1, nbval=0, nbret=nbno)
        call getvid('RESU', 'LIST_ORDRE', iocc=1, nbval=0, nbret=nblo)
        call getvr8('RESU', 'INST', iocc=1, nbval=0, nbret=nbni)
        call getvid('RESU', 'LIST_INST', iocc=1, nbval=0, nbret=nbli)
        call getvis('RESU', 'MODE', iocc=1, nbval=0, nbret=nbnm)
        call getvid('RESU', 'LIST_MODE', iocc=1, nbval=0, nbret=nblm)
        call getvr8('RESU', 'FREQ', iocc=1, nbval=0, nbret=nbnf)
        call getvid('RESU', 'LIST_FREQ', iocc=1, nbval=0, nbret=nblf)
!
        nbto = -nbto
        nbno = -nbno
        nblo = -nblo
        nbni = -nbni
        nbli = -nbli
        nbnm = -nbnm
        nblm = -nblm
        nbnf = -nbnf
        nblf = -nblf
        nbis = nbno+nblo+nbnm+nblm
        nbr8 = nbni+nbli+nbnf+nblf
!
!    -- ACCES PAR ORDRE, MODE, LIST_ORDRE, LIST_MODE :
        if (nbis .ne. 0) then
            call wkvect(nrval, 'V V R', 1, jrval)
            zr(jrval) = 0.0d0
!        -- NUME_ORDRE
            if (nbno .ne. 0) then
                typac = 'ORDRE'
                nbval = nbno
                call wkvect(nival, 'V V I', nbno, jival)
                call getvis('RESU', 'NUME_ORDRE', iocc=1, nbval=nbno, vect=zi(jival), &
                            nbret=n1)
!        -- MODE
            else if (nbnm .ne. 0) then
                typac = 'MODE'
                nbval = nbnm
                call wkvect(nival, 'V V I', nbnm, jival)
                call getvis('RESU', 'MODE', iocc=1, nbval=nbnm, vect=zi(jival), &
                            nbret=n1)
!        -- LIST_ORDRE, LIST_MODE
            else if (nblo .ne. 0) then
                if (nblo .ne. 0) then
                    typac = 'ORDRE'
                    call getvid('RESU', 'LIST_ORDRE', iocc=1, scal=nlist, nbret=n1)
                else
                    typac = 'MODE'
                    call getvid('RESU', 'LIST_MODE', iocc=1, scal=nlist, nbret=n1)
                end if
                nlist(20:24) = '.VALE'
                call jelira(nlist, 'LONMAX', nbval)
                call jeveuo(nlist, 'L', jlist)
                call wkvect(nival, 'V V I', nbval, jival)
                do i = 1, nbval
                    zi(jival+i-1) = zi(jlist+i-1)
                end do
            end if
!
!    -- ACCES PAR INST, FREQ, LIST_INST, LIST_FREQ :
        else if (nbr8 .ne. 0) then
            call wkvect(nival, 'V V I', 1, jival)
            zi(jival) = 0
!        -- INST
            if (nbni .ne. 0) then
                typac = 'INST'
                nbval = nbni
                call wkvect(nrval, 'V V R', nbni, jrval)
                call getvr8('RESU', 'INST', iocc=1, nbval=nbni, vect=zr(jrval), &
                            nbret=n1)
!        -- FREQ
            else if (nbnf .ne. 0) then
                typac = 'FREQ'
                nbval = nbnf
                call wkvect(nrval, 'V V R', nbnf, jrval)
                call getvr8('RESU', 'FREQ', iocc=1, nbval=nbnf, vect=zr(jrval), &
                            nbret=n1)
            else
                if (nbli .ne. 0) then
                    typac = 'INST'
                    call getvid('RESU', 'LIST_INST', iocc=1, scal=nlist, nbret=n1)
                else
                    typac = 'FREQ'
                    call getvid('RESU', 'LIST_FREQ', iocc=1, scal=nlist, nbret=n1)
                end if
                nlist(20:24) = '.VALE'
                call jelira(nlist, 'LONMAX', nbval)
                call jeveuo(nlist, 'L', jlist)
                call wkvect(nrval, 'V V R', nbval, jrval)
                do i = 1, nbval
                    zr(jrval+i-1) = zr(jlist+i-1)
                end do
            end if
!
!    -- ACCES TOUT_ORDRE :
        else
!
            cbid = dcmplx(0, 0)
            call rsorac(resu, 'LONUTI', ibid, r8b, k8b, &
                        cbid, r8b, k8b, tord, 1, &
                        nbtrou)
            nbval = tord(1)
            call wkvect(nival, 'V V I', nbval, jival)
            call rsorac(resu, 'TOUT_ORDRE', ibid, r8b, k8b, &
                        cbid, r8b, k8b, zi(jival), nbval, &
                        nbtrou)
!
!         -- SI LE RESULTAT EST UN EVOL_XXX, ON FORCE TYPAC='INST'
!            POUR QUE LES INSTANTS SOIENT IMPRIMES DANS LA TABLE :
            call gettco(resu, concep)
            if ((concep(1:5) .eq. 'EVOL_') .or. (concep(1:10) .eq. 'DYNA_TRANS')) then
                typac = 'INST'
                call wkvect(nrval, 'V V R', nbval, jrval)
                do kk = 1, nbval
                    nuord = zi(jival-1+kk)
                    call rsadpa(resu, 'L', 1, 'INST', nuord, &
                                0, sjv=jinst, styp=k8b)
                    rinst = zr(jinst)
                    zr(jrval-1+kk) = rinst
                end do
            else
                typac = 'ORDRE'
                call wkvect(nrval, 'V V R', 1, jrval)
                zr(jrval) = 0.0d0
            end if
!
        end if
!
!
! --- 2.2- ON DETERMINE :
!     ------------------
!         NKCHA = TABLEAU DES NOMS DE CHAMP
!
        call getvtx('RESU', 'NOM_CHAM', iocc=1, scal=nsymb, nbret=n1)
!
        call wkvect(niord, 'V V I', nbval, jniord)
        zi(jniord) = -1
        call wkvect(nkcha, 'V V K24', nbval, jkcha)
!
        if (typac .eq. 'ORDRE') then
            call jeveuo(nival, 'L', jival)
            do i = 1, nbval
                zi(jniord+i-1) = zi(jival+i-1)
                call rsexch(' ', resu, nsymb, zi(jival+i-1), zk24(jkcha+i-1), &
                            n1)
                if (n1 .ne. 0) then
                    valk = nsymb
                    vali = zi(jival+i-1)
                    call utmess('I', 'TABLE0_38', sk=valk, si=vali)
                    zk24(jkcha+i-1) = '&&CHAMP_INEXISTANT'
                end if
            end do
!
        else if (typac .eq. 'MODE') then
            call jeveuo(nival, 'L', jival)
            do i = 1, nbval
                call rsorac(resu, 'NUME_MODE', zi(jival+i-1), 0.0d0, k8b, &
                            cbid, epsi, crit, tord, 0, &
                            n1)
                if (n1 .eq. 0) then
                    zk24(jkcha+i-1) = '&&CHAMP_INEXISTANT'
                    valk = typac
                    vali = zi(jival+i-1)
                    zi(jniord+i-1) = -1
                    call utmess('I', 'TABLE0_39', sk=valk, si=vali)
                    goto 40
                end if
                n1 = -n1
                call rsorac(resu, 'NUME_MODE', zi(jival+i-1), 0.0d0, k8b, &
                            cbid, epsi, crit, zi(jniord+i-1), n1, &
                            n2)
                call rsexch(' ', resu, nsymb, zi(jniord+i-1), zk24(jkcha+i-1), &
                            n2)
                if (n2 .ne. 0) then
                    valk = nsymb
                    vali = zi(jival+i-1)
                    zi(jniord+i-1) = -1
                    call utmess('I', 'TABLE0_38', sk=valk, si=vali)
                    zk24(jkcha+i-1) = '&&CHAMP_INEXISTANT'
                end if
40              continue
            end do
!
        else if (typac .eq. 'INST' .or. typac .eq. 'FREQ') then
            call jeveuo(nrval, 'L', jrval)
            do i = 1, nbval
                call rsorac(resu, typac, 0, zr(jrval+i-1), k8b, &
                            cbid, epsi, crit, tord, 0, &
                            n1)
                if (n1 .eq. 0) then
                    zk24(jkcha+i-1) = '&&CHAMP_INEXISTANT'
                    valk = typac
                    valr = zr(jrval+i-1)
                    zi(jniord+i-1) = -1
                    call utmess('I', 'TABLE0_40', sk=valk, sr=valr)
                    goto 50
                end if
                n1 = -n1
                if (n1 .gt. 1) then
                    valk = typac
                    call utmess('F', 'TABLE0_46', sk=valk, sr=zr(jrval+i-1))
                end if
                call rsorac(resu, typac, 0, zr(jrval+i-1), k8b, &
                            cbid, epsi, crit, zi(jniord+i-1), n1, &
                            n2)
                call rsexch(' ', resu, nsymb, zi(jniord+i-1), zk24(jkcha+i-1), &
                            n2)
                if (n2 .ne. 0) then
                    valk = nsymb
                    vali = zi(jniord+i-1)
                    call utmess('I', 'TABLE0_38', sk=valk, si=vali)
                    zk24(jkcha+i-1) = '&&CHAMP_INEXISTANT'
                end if
50              continue
            end do
        end if
    end if
!
    call jedema()
!
end subroutine
