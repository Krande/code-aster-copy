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
subroutine op0038()
!
    use HHO_precalc_module, only: hhoAddInputField
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/sdmpic.h"
#include "asterfort/vrcins.h"
!
! --------------------------------------------------------------------------------------------------
!
! Command: CALC_CHAM_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ierd, iret, nh, nbRet, nbin
    real(kind=8) :: time, rundf
    character(len=1), parameter :: base = 'G'
    character(len=2) :: chdret
    character(len=8) :: model, caraElem, temp, mesh, kmpic, chmate
    character(len=8) :: lpain(10), lpaout(1)
    character(len=16) :: type, oper, option, phenom
    character(len=19) :: chelem, press, ligrel
    character(len=24) :: chgeom, chcara(18), chharm, mateco
    character(len=24) :: chtemp, chtime, chflug, chpres, chvarc
    character(len=24) :: lchin(10), lchout(1)
    aster_logical :: exitim, l_ther
    parameter(chvarc='&&OP0038.CHVARC')
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Initializations
    rundf = r8vide()

! - Output result
    call getres(chelem, type, oper)

! - Get main parameters
    model = ' '
    mateco = ' '
    caraElem = ' '
    chmate = ' '
    chtemp = ' '
    chpres = ' '
    call getvid(' ', 'MODELE', scal=model, nbret=nbRet)
    ASSERT(nbRet .eq. 1)
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenom)
    l_ther = phenom .eq. 'THERMIQUE'
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=nbRet)
    call getvid(' ', 'CHAM_MATER', scal=chmate, nbret=nbRet)
    mateco = ' '
    if (nbRet .ne. 0) then
        call rcmfmc(chmate, mateco, l_ther_=l_ther)
    end if
    call getvtx(' ', 'OPTION', scal=option, nbret=nbRet)
    temp = ' '
    call getvid(' ', 'TEMP', scal=temp, nbret=nbRet)
    if (nbRet .ne. 0) then
        chtemp = temp
        call chpver('F', chtemp, 'NOEU', 'TEMP_R', ierd)
    end if
    press = ' '
    call getvid(' ', 'PRES', scal=press, nbret=nbRet)
    if (nbRet .ne. 0) then
        chpres = press
        call chpver('F', chpres, 'NOEU', 'PRES_C', ierd)
    end if
    call getvr8(' ', 'INST', scal=time, nbret=nbRet)
    exitim = nbRet .ne. 0
    call getvis(' ', 'MODE_FOURIER', scal=nh, nbret=nbRet)
    if (nbRet .eq. 0) then
        nh = 0
    end if
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - List of cells for computation: all model
    call exlima(' ', 0, 'G', model, ligrel)

! - Prepare input field
    call mecham(option, model, caraElem, nh, chgeom, &
                chcara, chharm, iret)
    chtime = ' '
    if (exitim) then
        call mechti(mesh, time, rundf, rundf, chtime)
    end if

    if (iret .ne. 0) goto 10

! - Compute
    if (option(1:7) .eq. 'FLUX_EL') then
        chflug = '&&OP0038.FLUXGAUSS'
        lchin(1) = chgeom
        lpain(1) = 'PGEOMER'
        lchin(2) = mateco
        lpain(2) = 'PMATERC'
        lchin(3) = chcara(7)
        lpain(3) = 'PCACOQU'
        lchin(4) = chcara(12)
        lpain(4) = 'PCAMASS'
        lchin(5) = chtemp
        lpain(5) = 'PTEMPER'
        lchin(6) = chtime
        lpain(6) = 'PINSTR'
        lchin(7) = chharm
        lpain(7) = 'PHARMON'
        lchin(8) = ' '
        lpain(8) = 'PVARCPR'
        nbin = 8

        call hhoAddInputField(model, 10, lchin, lpain, nbin)

        lchout(1) = chflug
        lpaout(1) = 'PFLUXPG'
        call calcul('S', 'FLUX_ELGA', ligrel, nbin, lchin, &
                    lpain, 1, lchout, lpaout, 'V', &
                    'OUI')
        if (option .eq. 'FLUX_ELNO') then
            lchin(1) = chflug
            lpain(1) = 'PFLUXPG'
            lchout(1) = chelem
            lpaout(1) = 'PFLUXNO'
            call calcul('S', option, ligrel, 1, lchin, &
                        lpain, 1, lchout, lpaout, base, &
                        'OUI')

        else if (option .eq. 'FLUX_ELGA') then
            call copisd('CHAMP', 'G', chflug, chelem)

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (option .eq. 'COOR_ELGA') then
        lchin(1) = chgeom
        lpain(1) = 'PGEOMER'
        lchin(2) = chcara(1)
        lpain(2) = 'PCAORIE'
        lchin(3) = chcara(17)
        lpain(3) = 'PFIBRES'
        lchin(4) = chcara(16)
        lpain(4) = 'PNBSP_I'
        lchin(5) = chcara(7)
        lpain(5) = 'PCACOQU'
        lchin(6) = chcara(5)
        lpain(6) = 'PCAGEPO'
        lchout(1) = chelem
        lpaout(1) = 'PCOORPG'
        call cesvar(caraElem, ' ', ligrel, lchout(1))
        call calcul('S', option, ligrel, 6, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')

    else if (option .eq. 'ROCH_ELNO') then
        call vrcins(model, chmate, caraElem, time, chvarc(1:19), &
                    chdret)

        lchin(1) = mateco
        lpain(1) = 'PMATERC'
        lchin(2) = chcara(6)
        lpain(2) = 'PCAGNPO'
        lchin(3) = chcara(5)
        lpain(3) = 'PCAGEPO'
        lchin(4) = chvarc(1:19)
        lpain(4) = 'PVARCPR'
        lchout(1) = chelem
        lpaout(1) = 'PROCHRR'
        call calcul('S', option, ligrel, 4, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')
        call detrsd('CHAMP_GD', chvarc)

    else if (option .eq. 'PRAC_ELNO') then
        lpain(1) = 'PPRESSC'
        lchin(1) = chpres
        lchout(1) = chelem
        lpaout(1) = 'PPRAC_R'
        call calcul('S', option, ligrel, 1, lchin, &
                    lpain, 1, lchout, lpaout, 'G', &
                    'OUI')

    else
        ASSERT(ASTER_FALSE)

    end if
!
10  continue
!
!     -- SI CHELEM N'EST PAS MPI_COMPLET, ON LE COMPLETE :
    call dismoi('MPI_COMPLET', chelem, 'CHAM_ELEM', repk=kmpic)
    if (kmpic .eq. 'NON') call sdmpic('CHAM_ELEM', chelem)
!
    call jedema()
end subroutine
