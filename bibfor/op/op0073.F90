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
subroutine op0073()
    implicit none
!
!     DEFINITION D UN OBSTACLE DE CHOC DISCRETISE PAR FACETTES
!
!-----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=8) :: nomres
    integer(kind=8) :: nbpara, nbinfo
    parameter(nbpara=3)
    character(len=8) :: typara(nbpara)
    character(len=16) :: typres, nomcom, nopara(nbpara)
    character(len=24) :: type, tabk(nbpara)
    character(len=19) :: nomfon
    integer(kind=8) :: lval, lpro, lfon, nbval, nbpair
    integer(kind=8) :: ibid, idtemp, i
    integer(kind=8) :: ifm, niv
    real(kind=8) :: r8bid, rad
    complex(kind=8) :: cbid
    aster_logical :: crprol
!     ------------------------------------------------------------------
    data nopara/'LIEU', 'TYPE', 'FONCTION'/
    data typara/'K8', 'K24', 'K24'/
!     ------------------------------------------------------------------
!
    call jemarq()
    cbid = (0.d0, 0.d0)
    r8bid = 0.d0
    call infmaj()
    call infniv(ifm, niv)
!
    call getres(nomres, typres, nomcom)
!
!     --- VERIFICATIONS DE PREMIER NIVEAU ---
    call getvr8(' ', 'VALE', nbval=0, nbret=nbval)
    nbval = -nbval
    if ((nbval/2)*2 .ne. nbval) then
        call utmess('F', 'ALGORITH9_43')
    end if
!
! --- CREATION DE LA TABLE
    call tbcrsd(nomres, 'G')
    call tbajpa(nomres, nbpara, nopara, typara)
!
! --- TYPE DE L'OBSTACLE
    call getvtx(' ', 'TYPE', scal=type, nbret=ibid)
!
! --- FONCTION R=F(THETA EN RADIAN) DECRIVANT LA GEOMETRIE
    nomfon = nomres//'_INITIAL'
    crprol = .true.
!
! --- LIGNE DESCRIPTIVE
    nbinfo = nbpara
!     ON LIMITERA AU 2 PREMIERS PARAMETRES S'IL N'Y A PAS DE FONCTION...
    tabk(1) = 'DEFIOBST'
    tabk(2) = type
    tabk(3) = nomfon
!
! ===================================================================
!
! --- DIMENSIONNEMENT DES OBJETS DE STOCKAGE ---
    rad = r8dgrd()
    nbpair = nbval/2
!
    if (type(1:7) .eq. 'DISCRET') then
        if (nbval .gt. 0) then
            call wkvect('&&OP0073.TEMP', 'V V R', nbval, idtemp)
            call getvr8(' ', 'VALE', nbval=nbval, vect=zr(idtemp), nbret=ibid)
!
            call wkvect(nomfon//'.VALE', 'G V R', nbval, lval)
            lfon = lval+nbpair
            do i = 1, nbpair
                zr(lval-1+i) = zr(idtemp+2*(i-1))*rad
                zr(lfon-1+i) = zr(idtemp+2*(i-1)+1)
            end do
        end if
!
! --- CAS CERCLE, PLAN... SEUL LE .REFO ETAIT PRODUIT DANS L'ANCIENNE SD
    else
        crprol = .false.
        nbinfo = 2
    end if
!
    if (crprol) then
        ASSERT(lxlgut(nomfon) .le. 24)
        call wkvect(nomfon//'.PROL', 'G V K24', 6, lpro)
        zk24(lpro) = 'FONCTION'
        zk24(lpro+1) = 'LIN LIN'
        zk24(lpro+2) = 'THETA'
        zk24(lpro+3) = 'R'
        zk24(lpro+4) = 'EE'
        zk24(lpro+5) = nomfon
    end if
!
! --- INSERTION EFFECTIVE DE LA LIGNE DANS LA TABLE
    call tbajli(nomres, nbinfo, nopara, [ibid], [r8bid], &
                [cbid], tabk, 0)
!
    call jedema()
end subroutine
