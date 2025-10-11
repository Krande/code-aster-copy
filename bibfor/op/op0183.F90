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
subroutine op0183()
!
    use listLoad_type
    implicit none
!
!-----------------------------------------------------------------------
!     COMMANDE :  CALC_FORC_NONL
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/asasve.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/comp_info.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infbav.h"
#include "asterfort/infmaj.h"
#include "asterfort/infmue.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medomm.h"
#include "asterfort/nmdocc.h"
#include "asterfort/nmdoch.h"
#include "asterfort/onerrf.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/utmess.h"
#include "asterfort/vefnme.h"
#include "asterfort/vrcins.h"
#include "asterfort/vtcreb.h"
    integer(kind=8) :: ibid
    integer(kind=8) :: i, iad
    integer(kind=8) :: iordr, iret, iret2, j
    integer(kind=8) :: jfo
    integer(kind=8) :: jordr, ifm, niv
    integer(kind=8) :: lonch, lvafon, n0
    integer(kind=8) :: nbddl, nbordr, nc, nh, np
    integer(kind=8) :: ltps, ltps2
    real(kind=8) :: time, prec
    character(len=1), parameter :: jvBase = "V"
    character(len=2) :: codret
    character(len=6) :: nompro
    character(len=8) :: ctyp, crit, materField
    character(len=8) :: kiord
    character(len=16) :: resultType, oper, k16bid
    character(len=16) :: compex
    character(len=19) :: resultIn, knum, ligrel, resultOut, chdep2
    character(len=19), parameter :: listLoad = '&&OP0183.LIST_LOAD'
    character(len=16), parameter :: option = 'FONL_NOEU'
    character(len=24) :: model, caraElem, chamno, mateco
    character(len=24) :: nume, vafono, sigma, chdepl, k24bid
    character(len=24) :: compor, chvive, chacve, raide
    character(len=24), parameter :: vefnod = '&&OP0183.VEFNOD'
    character(len=24) :: chvarc
    character(len=24) :: numref, valk(3)
    aster_logical :: l_etat_init
    type(ListLoad_Prep) :: listLoadPrep
!     ------------------------------------------------------------------
    parameter(nompro='OP0183')
!     ------------------------------------------------------------------
!
    aster_logical :: exitim
    real(kind=8), pointer :: fono(:) => null()
    real(kind=8), pointer :: noch(:) => null()
!
!     ------------------------------------------------------------------
    data k24bid/' '/
    data chvarc/'&&OP0183.CHVARC'/
!     ------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
! --- ON STOCKE LE COMPORTEMENT EN CAS D'ERREUR AVANT MNL : COMPEX
! --- PUIS ON PASSE DANS LE MODE "VALIDATION DU CONCEPT EN CAS D'ERREUR"
    call onerrf(' ', compex, ibid)
    call onerrf('EXCEPTION+VALID', k16bid, ibid)
!
    call infmue()

! - Get result datastructures
    call getres(resultOut, resultType, oper)
    ASSERT(resultType .eq. 'DYNA_TRANS')
    call getvid(' ', 'RESULTAT', scal=resultIn, nbret=n0)
    ASSERT(resultIn .ne. resultOut)

!
    knum = '&&'//nompro//'.NUME_ORDRE'

    call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
    call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
    call rsutnu(resultIn, ' ', 0, knum, nbordr, &
                prec, crit, iret)
    if (iret .eq. 10) then
        call utmess('F', 'CALCULEL4_8', sk=resultIn)
        goto 60
    end if
    if (iret .ne. 0) then
        call utmess('F', 'ALGORITH3_41')
        goto 60
    end if
    call jeveuo(knum, 'L', jordr)
!
    exitim = .true.
    call rscrsd('G', resultOut, resultType, nbordr)

! - Get main parameters from user
    call medomm(model, materField, mateco, caraElem)
    ASSERT(model .ne. ' ')
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)

! - Get loads/BC and create list of loads datastructure
    listLoadPrep%model = model(1:8)
    call nmdoch(listLoadPrep, listLoad, jvBase)
!
    time = 0.d0
    l_etat_init = .false.
!
    numref = ' '
    call refdcp(resultIn, resultOut)
    call dismoi('REF_RIGI_PREM', resultOut, 'RESU_DYNA', repk=raide, arret='C')
    if (raide .ne. ' ') then
        call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=numref)
    end if
!
    do i = 1, nbordr
        call jemarq()
        iordr = zi(jordr+i-1)
        vafono = ' '
        nh = 0
!
        call rsexch(' ', resultIn, 'SIEF_ELGA', iordr, sigma, &
                    iret)
        if (iret .ne. 0) then
            call rsexch(' ', resultIn, 'SIEF_ELGA', iordr, sigma, &
                        iret2)
            if (iret2 .ne. 0 .and. option .ne. 'FONL_NOEU') then
                call codent(iordr, 'G', kiord)
                valk(1) = kiord
                valk(2) = option
                call utmess('A', 'PREPOST5_2', nk=2, valk=valk)
                goto 40
!
            end if
            if (iret2 .ne. 0 .and. option .eq. 'FONL_NOEU') then
                sigma = ' '
            end if
        end if
!
        call rsexch(' ', resultIn, 'DEPL', iordr, chdepl, &
                    iret)
        if (iret .ne. 0) then
            call codent(iordr, 'G', kiord)
            valk(1) = kiord
            valk(2) = option
            call utmess('A', 'PREPOST5_3', nk=2, valk=valk)
            goto 40
        else
!         CREATION D'UN VECTEUR ACCROISSEMENT DE DEPLACEMENT NUL
!         POUR LE CALCUL DE FORC_NODA DANS LES POU_D_T_GD
!
            chdep2 = '&&'//nompro//'.CHDEP_NUL'
            call copisd('CHAMP_GD', 'V', chdepl, chdep2)
            call jelira(chdep2//'.VALE', 'LONMAX', nbddl)
            call jerazo(chdep2//'.VALE', nbddl, 1)
        end if
!
!       -- CALCUL D'UN NUME_DDL "MINIMUM" POUR ASASVE :
        nume = numref(1:14)//'.NUME'
!
        call rsexch(' ', resultIn, 'VITE', iordr, chvive, &
                    iret)
        if (iret .eq. 0) then
            chvive = '&&'//nompro//'.CHVIT_NUL'
            call copisd('CHAMP_GD', 'V', chdepl, chvive)
            call jelira(chvive(1:19)//'.VALE', 'LONMAX', nbddl)
            call jerazo(chvive(1:19)//'.VALE', nbddl, 1)
        end if
        call rsexch(' ', resultIn, 'ACCE', iordr, chacve, &
                    iret)
        if (iret .eq. 0) then
            chacve = '&&'//nompro//'.CHACC_NUL'
            call copisd('CHAMP_GD', 'V', chdepl, chacve)
            call jelira(chacve(1:19)//'.VALE', 'LONMAX', nbddl)
            call jerazo(chacve(1:19)//'.VALE', nbddl, 1)
        end if
!
        if (exitim) then
            call rsadpa(resultIn, 'L', 1, 'INST', iordr, &
                        0, sjv=iad, styp=ctyp)
            time = zr(iad)
        end if
!
        call vrcins(model, materField, caraElem, time, chvarc(1:19), &
                    codret)
!
!       --- CALCUL DES VECTEURS ELEMENTAIRES ---
        if (i .eq. 1) then
            compor = '&&OP0183.COMPOR'
            call nmdocc(model(1:8), materField, l_etat_init, compor, 'V')
            if (niv .ge. 2) then
                call comp_info(model(1:8), compor)
            end if
        end if
!
        call vefnme(option, model, mateco, caraElem, &
                    compor, nh, ligrel, chvarc, &
                    sigma, ' ', chdepl, 'V', &
                    vefnod)
!
!       --- ASSEMBLAGE DES VECTEURS ELEMENTAIRES ---
        call asasve(vefnod, nume, 'R', vafono)
!
        call rsexch(' ', resultOut, 'DEPL', iordr, chamno, &
                    iret)
        call rsadpa(resultOut, 'E', 1, 'INST', iordr, 0, sjv=ltps2)
        call rsadpa(resultIn, 'L', 1, 'INST', iordr, 0, sjv=ltps)
        zr(ltps2) = zr(ltps)
!
!
        call jeexin(chamno(1:19)//'.REFE', iret)
        if (iret .ne. 0) then
            call codent(iordr, 'G', kiord)
            valk(1) = option
            valk(2) = kiord
            call utmess('I', 'PREPOST5_1', nk=2, valk=valk)
            call detrsd('CHAM_NO', chamno(1:19))
        end if
        call vtcreb(chamno, 'G', 'R', nume_ddlz=nume)
        call jeveuo(chamno(1:19)//'.VALE', 'E', vr=noch)
!
        call jeveuo(vafono, 'L', jfo)
        call jeveuo(zk24(jfo) (1:19)//'.VALE', 'L', vr=fono)
        call jelira(zk24(jfo) (1:19)//'.VALE', 'LONMAX', lvafon)
        call jelira(chamno(1:19)//'.VALE', 'LONMAX', lonch)
!
        do j = 0, lonch-1
            noch(1+j) = fono(1+j)
        end do
!
        call rsnoch(resultOut, 'DEPL', iordr)
        call detrsd('CHAMP_GD', '&&'//nompro//'.SIEF')
        call detrsd('VECT_ELEM', vefnod)
40      continue
        call jedema()
    end do
!
    call jedetr(knum)
!
! --- ON REMET LE MECANISME D'EXCEPTION A SA VALEUR INITIALE
    call onerrf(compex, k16bid, ibid)
!
60  continue
    call infbav()
    call jedema()
!
end subroutine
