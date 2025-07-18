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

subroutine ccpara(option, modele, resuin, resuou, numord, &
                  nordm1, exitim, mater, carael)
    implicit none
!     --- ARGUMENTS ---
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mearcc.h"
#include "asterfort/mecact.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterc/isnnem.h"
!
    integer(kind=8) :: numord, nordm1
    character(len=8) :: modele, resuin, resuou, mater, carael
    character(len=16) :: option
    aster_logical :: exitim
! person_in_charge: nicolas.sellenet at edf.fr
    integer(kind=8) :: nparin, iret, ibid, jainst, jfreq, iaopds, iaoplo, ipara
    integer(kind=8) :: jnmo, opt, inume, jains2, iundef
!
    real(kind=8) :: tps(6), rundf, omega2, freq, time, zero
    parameter(zero=0.0d0)
!
    character(len=2) :: chdret
    character(len=8) :: mailla, k8b, nomcmp(6)
    character(len=16) :: vari
    character(len=24) :: curcha, chtime, chfreq, chome2, chharm, chvref
    character(len=24) :: chmass, chnova, chsigf, chsig, chvarc, chvac2
    parameter(chtime='&&CCPARA.CH_INST_R')
    parameter(chfreq='&&CCPARA.FREQ')
    parameter(chome2='&&CCPARA.OMEGA2')
    parameter(chharm='&&CCPARA.NUME_MODE')
    parameter(chvarc='&&CCPARA.VARI_INT_N')
    parameter(chvref='&&CCPARA.VARI_INT_REF')
    parameter(chvac2='&&CCPARA.VARI_INT_NM1')
    parameter(chmass='&&CCPARA.MASS_MECA_D')
    parameter(chnova='&&CCPARA.NOM_VARI')
    parameter(chsigf='&&CCPARA.CHAM_SI2D')
!
    data nomcmp/'INST    ', 'DELTAT  ', 'THETA   ', 'KHI     ',&
     &              'R       ', 'RHO     '/
! ======================================================================
!
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), opt)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopds)
    call jeveuo(jexnum('&CATA.OP.LOCALIS', opt), 'L', iaoplo)
    nparin = zi(iaopds-1+2)
!
    rundf = r8nnem()
    iundef = isnnem()
!
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mailla)
!
    if (exitim) then
        call rsadpa(resuin, 'L', 1, 'INST', numord, &
                    0, sjv=jainst, styp=k8b)
        time = zr(jainst)
    else
        time = zero
    end if
    call vrcref(modele, mater, carael, chvref(1:19))
    call vrcins(modele, mater, carael, time, chvarc(1:19), &
                chdret)
!
    do ipara = 1, nparin
        curcha = zk24(iaoplo+3*ipara-2)
        if (curcha .eq. chtime) then
            call jenonu(jexnom(resuin//'           .NOVA', 'INST'), iret)
            if (iret .ne. 0) then
                call rsadpa(resuin, 'L', 1, 'INST', numord, &
                            0, sjv=jainst, styp=k8b)
!
                tps(1) = zr(jainst)
                tps(2) = rundf
                tps(3) = rundf
                tps(4) = rundf
                tps(5) = rundf
                tps(6) = rundf
!
                call mecact('V', chtime, 'MAILLA', mailla, 'INST_R', &
                            ncmp=6, lnomcmp=nomcmp, vr=tps)
            end if
!
        else if (curcha .eq. chfreq) then
            call jenonu(jexnom(resuin//'           .NOVA', 'FREQ'), iret)
            if (iret .ne. 0) then
                call rsadpa(resuin, 'L', 1, 'FREQ', numord, &
                            0, sjv=jfreq, styp=k8b)
                freq = zr(jfreq)
            else
                freq = 1.d0
            end if
            call mecact('V', chfreq, 'MAILLA', mailla, 'FREQ_R', &
                        ncmp=1, nomcmp='FREQ', sr=freq)
!
        else if (curcha .eq. chome2) then
            call jenonu(jexnom(resuin//'           .NOVA', 'OMEGA2'), iret)
            if (iret .ne. 0) then
                call rsadpa(resuin, 'L', 1, 'OMEGA2', numord, &
                            0, sjv=jfreq, styp=k8b)
                omega2 = zr(jfreq)
            else
                omega2 = 1.0d0
            end if
            call mecact('V', chome2, 'MAILLA', mailla, 'OME2_R', &
                        ncmp=1, nomcmp='OMEG2', sr=omega2)
!
        else if (curcha .eq. chharm) then
            call jenonu(jexnom(resuin//'           .NOVA', 'NUME_MODE'), iret)
            if (iret .ne. 0) then
                call rsadpa(resuin, 'L', 1, 'NUME_MODE', numord, &
                            0, sjv=jnmo, styp=k8b, istop=0)
                if (zi(jnmo) .ne. iundef) then
                    call mecact('V', chharm, 'MAILLA', mailla, 'HARMON', &
                                ncmp=1, nomcmp='NH', si=zi(jnmo))

                end if
            end if
!
        else if (curcha .eq. chmass) then
            inume = 1
            call mecact('V', chmass, 'MAILLA', mailla, 'POSI', &
                        ncmp=1, nomcmp='POS', si=inume)
!
        else if (curcha .eq. chvac2) then
            if (exitim) then
                call rsadpa(resuin, 'L', 1, 'INST', nordm1, &
                            0, sjv=jains2, styp=k8b)
                time = zr(jains2)
            else
                time = zero
            end if
            call vrcins(modele, mater, carael, time, chvac2(1:19), &
                        chdret)
!
        else if (curcha .eq. chnova) then
            call getvtx(' ', 'NOM_VARI', scal=vari, nbret=ibid)
            call mecact('V', chnova, 'MAILLA', mailla, 'NEUT_K24', &
                        ncmp=1, nomcmp='Z1', sk=vari)
!
        else if (curcha .eq. chsigf) then
            call rsexch(' ', resuin, 'SIGM_ELNO', numord, chsig, &
                        iret)
            if (iret .ne. 0) then
                call rsexch('F', resuou, 'SIGM_ELNO', numord, chsig, &
                            iret)
            end if
            call mearcc(option, modele, chsig, chsigf)
!
        end if
    end do
!
end subroutine
