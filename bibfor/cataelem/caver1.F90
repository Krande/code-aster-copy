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
subroutine caver1()
    implicit none
!
! ----------------------------------------------------------------------
!     BUT: VERIFIER LA COHERENCE DES OBJETS DES CATALOGUES.
!
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/kndoub.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: opt, te
    aster_logical :: error
    character(len=8) :: para, typmai
    character(len=16) :: nomopt, nomte
    character(len=24) :: valk(4)
    character(len=8) :: gd1, gd2, tgd1(10), tgd2(10), typout, typou2
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iadesc, iamolo, iaopmo, iaopno, iapara
    integer(kind=8) :: icode, ier, igd, igdop, imolo, ioptte, ipara
    integer(kind=8) :: iret, itrou, jnbno, jnocm1, jnocm2, k
    integer(kind=8) :: kk, lgco, n1, n2, nbgd, nbin, nbinte
    integer(kind=8) :: nbno, nbopt, nbout, nboute, nbpt1, nbpt2, nbte
    integer(kind=8) :: nbvol, nucalc
    integer(kind=8), pointer :: nbligcol(:) => null()
    integer(kind=8), pointer :: optte(:) => null()
    character(len=8), pointer :: typema(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    call jelira('&CATA.OP.NOMOPT', 'NOMMAX', nbopt)
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nbte)
    call jeveuo('&CATA.TE.OPTTE', 'L', vi=optte)
    call jeveuo('&CATA.TE.NBLIGCOL', 'L', vi=nbligcol)
    lgco = nbligcol(1)
!
    ier = 0
!
    call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typema)
!
!
    do opt = 1, nbopt
        call jenuno(jexnum('&CATA.OP.NOMOPT', opt), nomopt)
        call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iadesc)
        call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iapara)
        nbin = zi(iadesc-1+2)
        nbout = zi(iadesc-1+3)
        nbvol = zi(iadesc-1+4)
        if (nbvol .ne. 0) then
            call utmess('E', 'CATAELEM_1', sk=nomopt)
            ier = ier+1
        end if
!
!
!       -- ON VERIFIE QU'UNE OPTION N'A JAMAIS 2 PARAMETRES
!          DE MEME NOM
!          -------------------------------------------------
        call kndoub(8, zk8(iapara), nbin+nbout, iret)
        if (iret .gt. 0) then
            call utmess('E', 'CATAELEM_2', sk=nomopt)
        end if
!
!
        do te = 1, nbte
            call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
            ioptte = optte((te-1)*lgco+opt)
            if (ioptte .eq. 0) goto 30
            call jeveuo(jexnum('&CATA.TE.OPTMOD', ioptte), 'L', iaopmo)
            nucalc = zi(iaopmo-1+1)
            nbinte = zi(iaopmo-1+2)
!
            typmai = typema(te)
            call jeveuo(jexnom('&CATA.TM.NBNO', typmai), 'L', jnbno)
            nbno = zi(jnbno)
!
!           -- on ne traque pas les erreurs si nucalc = 0, -1 ou -2
            if ((nucalc .le. 0) .and. (nucalc .ge. -2)) goto 30
!
            call jeveuo(jexnum('&CATA.TE.OPTNOM', ioptte), 'L', iaopno)
            do ipara = 1, nbinte
                para = zk8(iaopno-1+ipara)
                imolo = zi(iaopmo-1+3+ipara)
                if (imolo .eq. 0) then
                    valk(1) = para
                    valk(2) = nomopt
                    valk(3) = nomte
                    call utmess('E', 'CATAELEM_3', nk=3, valk=valk)
                    ier = ier+1
                    goto 10
                end if
!
                call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', iamolo)
                igd = zi(iamolo-1+2)
                itrou = indik8(zk8(iapara-1+1), para, 1, nbin)
                igdop = zi(iadesc-1+4+itrou)
                if ((itrou .eq. 0) .or. (igdop .ne. igd)) then
                    if (itrou .eq. 0) then
                        valk(1) = para
                        valk(2) = nomopt
                        valk(3) = nomte
                        call utmess('E', 'CATAELEM_4', nk=3, valk=valk)
                        ier = ier+1
                    end if
                    if (igdop .ne. igd) then
                        valk(1) = para
                        valk(2) = nomopt
                        valk(3) = nomte
                        call utmess('E', 'CATAELEM_5', nk=3, valk=valk)
                        ier = ier+1
                    end if
                end if
!
!              -- ON VERIFIE QUE POUR LES MODE LOCAUX AUX NOEUDS
!                 LE NOMBRE DE NOEUDS EST LE NOMBRE DE NOEUDS DE
!                 LA MAILLE SUPPORT :
!              ---------------------------------------------------
                icode = zi(iamolo-1+1)
                nbpt2 = -1
                if (icode .eq. 2) then
                    nbpt2 = mod(zi(iamolo-1+4), 10000)
                else if (icode .eq. 3) then
                    nbpt1 = zi(iamolo-1+4)
                    if (nbpt1 .lt. 0) then
                        nbpt2 = mod(abs(nbpt1), 10000)
                    end if
                end if
                if (nbpt2 .ge. 0) then
                    if (nbpt2 .ne. nbno) then
                        valk(1) = para
                        valk(2) = nomopt
                        valk(3) = nomte
                        call utmess('E', 'CATAELEM_6', nk=3, valk=valk)
                        ier = ier+1
                    end if
                end if
!
!
10              continue
            end do
!
!
!         -- VERIFICATION DES MODES LOCAUX "OUT" DES TE/OPTIONS:
!         ------------------------------------------------------
            nboute = zi(iaopmo-1+3)
            do ipara = 1, nboute
                para = zk8(iaopno-1+nbinte+ipara)
                imolo = zi(iaopmo-1+3+nbinte+ipara)
                if (imolo .eq. 0) then
                    valk(1) = para
                    valk(2) = nomopt
                    valk(3) = nomte
                    call utmess('E', 'CATAELEM_3', nk=3, valk=valk)
                    ier = ier+1
                    goto 20
                end if
                call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', iamolo)
                igd = zi(iamolo-1+2)
                itrou = indik8(zk8(iapara-1+nbin+1), para, 1, nbout)
                igdop = zi(iadesc-1+4+nbin+itrou)
                if ((itrou .eq. 0) .or. (igdop .ne. igd)) then
                    if (itrou .eq. 0) then
                        valk(1) = para
                        valk(2) = nomopt
                        valk(3) = nomte
                        call utmess('E', 'CATAELEM_4', nk=3, valk=valk)
                        ier = ier+1
                    end if
                    if (igdop .ne. igd) then
                        valk(1) = para
                        valk(2) = nomopt
                        valk(3) = nomte
                        call utmess('E', 'CATAELEM_5', nk=3, valk=valk)
                        ier = ier+1
                    end if
                end if
!
!           -- ON VERIFIE QUE  LE TYPE DU CHAMP LOCAL EST COHERENT
!               AVEC CELUI DECLARE DANS L'OPTION :
!           ---------------------------------------------------
                icode = zi(iamolo-1+1)
                typou2 = zk8(iapara-1+nbin+nbout+itrou)
                typout = '????'
                if (icode .ge. 4) then
                    typout = 'RESL'
                else if (icode .eq. 3) then
                    typout = 'ELGA'
                else if (icode .eq. 2) then
                    typout = 'ELNO'
                else if (icode .eq. 1) then
                    typout = 'ELEM'
                end if
!
                if (typout .ne. typou2) then
                    valk(1) = para
                    valk(2) = nomopt
                    valk(3) = nomte
                    valk(4) = typou2
                    call utmess('E', 'CATAELEM_7', nk=4, valk=valk)
!             IER = IER + 1
                end if
20              continue
            end do
!
30          continue
        end do
!
    end do
!
!    -- ON VERIFIE QUE CERTAINES GRANDEURS SONT "PARALLELES" :
!       ELLES DOIVENT AVOIR EXACTEMENT LES MEMES CMPS :
!    ----------------------------------------------------------
    nbgd = 2
    tgd1(1) = 'DEPL_R'
    tgd2(1) = 'DEPL_C'
    tgd1(2) = 'DEPL_R'
    tgd2(2) = 'DEPL_F'
    error = .false.
    do k = 1, nbgd
        gd1 = tgd1(k)
        gd2 = tgd2(k)
        call jelira(jexnom('&CATA.GD.NOMCMP', gd1), 'LONMAX', n1)
        call jelira(jexnom('&CATA.GD.NOMCMP', gd2), 'LONMAX', n2)
        if (n1 .ne. n2) then
            error = .true.
        else
            call jeveuo(jexnom('&CATA.GD.NOMCMP', gd1), 'L', jnocm1)
            call jeveuo(jexnom('&CATA.GD.NOMCMP', gd2), 'L', jnocm2)
            do kk = 1, n1
                if (zk8(jnocm1-1+kk) .ne. zk8(jnocm1-1+kk)) error = .true.
            end do
        end if
        if (error) then
            valk(1) = gd1
            valk(2) = gd2
            call utmess('E', 'CATAELEM_8', nk=2, valk=valk)
            ier = ier+1
        end if
    end do
!
!
    if (ier .gt. 0) then
        call utmess('F', 'CATAELEM_9')
    end if
!
    call jedema()
end subroutine
