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
subroutine trmult(modsta, numexi, mailla, neq, iddeeq, &
                  pside, numddl)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/compno.h"
#include "asterfort/copmod.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/r8inir.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsorac.h"
#include "asterfort/rsvpar.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: modsta, mailla, numddl
    integer(kind=8) :: numexi, neq, iddeeq
    real(kind=8) :: pside(neq)
!     OPERATEUR :   DYNA_TRAN_MODAL
!
!     CREE ET CALCULE LE VECTEUR PSI*DIRECTION DANS LE CAS D'UN CALCUL
!     SISMIQUE D UNE STRUCTURE MULTI-SUPPORTEE
!     ------------------------------------------------------------------
! IN  : MODSTA : NOM DU CONCEPT MODES STATIQUES
! IN  : NUMEXI : NUMERO D'OCCURENCE DU MOT CLE EXCIT
! IN  : MAILLA : NOM DU MAILLAGE
! IN  : NEQ    : NOMBRE D'EQUATIONS
! IN  : IDDEEQ : INDICE DE L'EQUATION
! OUT : PSIDE  : VALEURS DU VECTEUR PSI*DELTA
    real(kind=8) :: xnorm, depl(6)
    character(len=8) :: nomnoe
    character(len=24) :: magrno
    character(len=24) :: valk(3)
    character(len=8) :: kbid
    integer(kind=8) :: ibid, iordr(1), ier
    real(kind=8) :: r8b, epsi
    character(len=8) :: cmp(6), crit
    character(len=16) :: acces
    character(len=19) :: chamno
    complex(kind=8) :: c16b
    integer(kind=8) :: imode
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, id, idno, ii, in
    integer(kind=8) :: iret, ldgn, nb, nbd, nbdir, nbgr, nbno
    integer(kind=8) :: nbtrou, nbv, tmod(1)
    real(kind=8) :: xd
    character(len=24), pointer :: group_no(:) => null()
    real(kind=8), pointer :: base(:) => null()
!-----------------------------------------------------------------------
    data cmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!
!     --- RECUPERATION DES ARGUMENTS DE LA COMMANDE ---
!
    call jemarq()
    epsi = 1.d-4
    magrno = ' '
    ier = 0
!
!     --- RECUPERATION DE LA DIRECTION SISMIQUE  ---
!
    call getvr8('EXCIT', 'DIRECTION', iocc=numexi, nbval=0, nbret=nbd)
    nbdir = -nbd
    call getvr8('EXCIT', 'DIRECTION', iocc=numexi, nbval=nbdir, vect=depl, &
                nbret=nbd)
!     --- ON NORMALISE LE VECTEUR ---
    xnorm = 0.d0
    do i = 1, nbdir
        xnorm = xnorm+depl(i)*depl(i)
    end do
    xnorm = sqrt(xnorm)
    if (xnorm .lt. 0.d0) then
        call utmess('F', 'ALGORITH9_81')
    end if
    do i = 1, nbdir
        depl(i) = depl(i)/xnorm
    end do
!
!     --- RECUPERATION DES POINTS D'ANCRAGE ---
!
    call getvem(mailla, 'NOEUD', 'EXCIT', 'NOEUD', numexi, &
                0, kbid, nbno)
    if (nbno .ne. 0) then
!        --- ON RECUPERE UNE LISTE DE NOEUD ---
        nbno = -nbno
        call wkvect('&&TRMULT.NOEUD', 'V V K8', nbno, idno)
        call getvem(mailla, 'NOEUD', 'EXCIT', 'NOEUD', numexi, &
                    nbno, zk8(idno), nbv)
    else
!        --- ON RECUPERE UNE LISTE DE GROUP_NO ---
        call getvem(mailla, 'GROUP_NO', 'EXCIT', 'GROUP_NO', numexi, &
                    0, kbid, nbgr)
        nbgr = -nbgr
        if (nbgr .ne. 0) then
            AS_ALLOCATE(vk24=group_no, size=nbgr)
            call getvem(mailla, 'GROUP_NO', 'EXCIT', 'GROUP_NO', numexi, &
                        nbgr, group_no, nbv)
!           --- ECLATE LE GROUP_NO EN NOEUD ---
            call compno(mailla, nbgr, group_no, nbno)
            call wkvect('&&TRMULT.NOEUD', 'V V K8', nbno, idno)
            magrno = mailla//'.GROUPENO'
            ii = -1
            do i = 1, nbgr
                call jelira(jexnom(magrno, group_no(i)), 'LONUTI', nb)
                call jeveuo(jexnom(magrno, group_no(i)), 'L', ldgn)
                do in = 0, nb-1
                    nomnoe = int_to_char8(zi(ldgn+in))
                    ii = ii+1
                    zk8(idno+ii) = nomnoe
                end do
            end do
        end if
    end if
!
    call rsorac(modsta, 'LONUTI', 0, r8b, kbid, &
                c16b, r8b, kbid, tmod, 1, &
                ibid)
    AS_ALLOCATE(vr=base, size=neq*tmod(1))
    call copmod(modsta, bmodr=base, numer=numddl)
!
!     --- CALCUL DU VECTEUR PSI*DELTA ---
!
    call r8inir(neq, 0.d0, pside, 1)
    do id = 1, nbdir
        xd = depl(id)
        if (abs(xd) .gt. epsi) then
            do in = 1, nbno
                acces(1:8) = zk8(idno+in-1)
                acces(9:16) = cmp(id)
!
!              --- ON RECUPERE LE MODE STATIQUE ASSOCIE AU NOEUD ---
                call rsorac(modsta, 'NOEUD_CMP', ibid, r8b, acces, &
                            c16b, epsi, crit, iordr, 1, &
                            nbtrou)
                if (nbtrou .ne. 1) then
                    ier = ier+1
                    valk(1) = acces(1:8)
                    valk(2) = acces(9:16)
                    call utmess('F', 'ALGELINE4_61', nk=2, valk=valk)
                    goto 40
                end if
                imode = iordr(1)
!
                call rsvpar(modsta, imode, 'TYPE_DEFO', ibid, r8b, &
                            'DEPL_IMPO', iret)
                if (iret .ne. 100) then
                    ier = ier+1
                    valk(1) = 'MODE_MECA'
                    valk(2) = acces(1:8)
                    valk(3) = acces(9:16)
                    call utmess('F', 'ALGELINE4_62', nk=3, valk=valk)
                    goto 40
                end if
                call rsexch('F', modsta, 'DEPL', imode, chamno, &
                            iret)
!
!              --- ON EFFECTUE LE PRODUIT  MODE_STAT * DIR ---
                do i = 1, neq
                    pside(i) = pside(i)+xd*base((imode-1)*neq+i)
                end do
40              continue
            end do
        end if
    end do
!
!     --- MISE A ZERO DES DDL DE LAGRANGE
    call zerlag(neq, zi(iddeeq), vectr=pside(1))
!
    call jedetr('&&TRMULT.NOEUD')
    AS_DEALLOCATE(vk24=group_no)
    AS_DEALLOCATE(vr=base)
!
    call jedema()
end subroutine
