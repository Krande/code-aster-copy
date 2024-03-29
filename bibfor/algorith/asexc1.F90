! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
subroutine asexc1(motfac, nbocc, nbmode, momec, amort, &
                  corfre, ndir, valspe, asyspe, nopara, &
                  nordr)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
!
    integer :: nbocc, nbmode, ndir(*), nordr(*)
    real(kind=8) :: amort(*), valspe(3, *), asyspe(*)
    character(len=*) :: motfac, momec
    character(len=24) :: nopara(*)
    aster_logical :: corfre
!     COMMANDE : COMB_SISM_MODAL
!                TRAITEMENT DU MOT-CLE "EXCIT" POUR LE MONO-APPUI
!     ------------------------------------------------------------------
! IN  : MOTFAC : MOT CLE FACTEUR
! IN  : NBOCC  : NOMBRE D'OCCURENCE DU MOT CLE FACTEUR
! IN  : NBMODE : NOMBRE DE MODES
! IN  : AMORT  : AMORTISSEMENTS MODAUX
! IN  : MOMEC  : MODES MECANIQUES
! IN  : CORFRE : CORRECTION FREQUENCE SI .TRUE.
! OUT : NDIR   : DIRECTION DU SEISME A ETUDIER
! OUT : VALSPE : VALEURS DU SPECTRE
! OUT : ASYSPE : VALEURS ASYMPTOTIQUES DU SPECTRE
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer :: nature(3), id, ier, ii, im, inat, ioc, n1
    integer :: nimpr, ival
    real(kind=8) :: amor, coef, deuxpi, echel, epsi, freq, dirspe(3), echspe(3)
    real(kind=8) :: valpu(2), omega, omega2, resu, un, uns2pi, xnorm, zero
    real(kind=8) :: fcoup
    character(len=1) :: dir(3)
    character(len=4) :: knat
    character(len=8) :: spect, nomspe(3), nompu(2)
    character(len=9) :: niveau
    real(kind=8) :: correc
!     ------------------------------------------------------------------
    data nompu/'AMOR', 'FREQ'/
    data dir/'X', 'Y', 'Z'/
!     ------------------------------------------------------------------
!
    call jemarq()
    ier = 0
    epsi = 1.d-03
    zero = 0.d0
    un = 1.d0
    deuxpi = r8depi()
    uns2pi = un/deuxpi
!
!     --- LECTURE MOT-CLE FACTEUR IMPRESSION ---
!
    call getvtx('IMPRESSION', 'NIVEAU', iocc=1, scal=niveau, nbret=nimpr)
    if (nimpr .eq. 0) niveau = 'TOUT     '
!
    call getvr8(' ', 'FREQ_COUP', iocc=1, scal=fcoup, nbret=n1)
    if (n1 .eq. 0) then
        call rsadpa(momec, 'L', 1, nopara(2), nordr(nbmode), &
                    0, sjv=ival, istop=0)
        fcoup = uns2pi*sqrt(zr(ival))
    end if
!
!
    do ioc = 1, nbocc
!
        echspe(1) = un
        echspe(2) = un
        echspe(3) = un
        dirspe(1) = un
        dirspe(2) = un
        dirspe(3) = un
        xnorm = un
!
!        --- RECUPERATION DE LA DIRECTION DU SPECTRE ---
        call getvr8(motfac, 'AXE', iocc=ioc, nbval=0, nbret=n1)
        if (n1 .ne. 0) then
            call getvr8(motfac, 'AXE', iocc=ioc, nbval=3, vect=dirspe, &
                        nbret=n1)
            xnorm = zero
            do id = 1, 3
                xnorm = xnorm+dirspe(id)*dirspe(id)
            end do
            if (xnorm .lt. epsi) then
                ier = ier+1
                call utmess('E', 'SEISME_4')
                goto 10
            end if
            xnorm = un/sqrt(xnorm)
            call getvid(motfac, 'SPEC_OSCI', iocc=ioc, scal=spect, nbret=n1)
            nomspe(1) = spect
            nomspe(2) = spect
            nomspe(3) = spect
            call getvr8(motfac, 'ECHELLE', iocc=ioc, scal=echel, nbret=n1)
            if (n1 .ne. 0) then
                echspe(1) = echel
                echspe(2) = echel
                echspe(3) = echel
            end if
!
        else
            call getvr8(motfac, 'TRI_AXE', iocc=ioc, nbval=0, nbret=n1)
            if (n1 .ne. 0) then
                call getvr8(motfac, 'TRI_AXE', iocc=ioc, nbval=3, vect=dirspe, &
                            nbret=n1)
                call getvid(motfac, 'SPEC_OSCI', iocc=ioc, scal=spect, nbret=n1)
                nomspe(1) = spect
                nomspe(2) = spect
                nomspe(3) = spect
                call getvr8(motfac, 'ECHELLE', iocc=ioc, scal=echel, nbret=n1)
                if (n1 .ne. 0) then
                    echspe(1) = echel
                    echspe(2) = echel
                    echspe(3) = echel
                end if
!
            else
!
                call getvid(motfac, 'SPEC_OSCI', iocc=ioc, nbval=3, vect=nomspe, &
                            nbret=n1)
                call getvr8(motfac, 'ECHELLE', iocc=ioc, nbval=3, vect=echspe, &
                            nbret=n1)
            end if
        end if
!
        call getvtx(motfac, 'NATURE', iocc=ioc, scal=knat, nbret=n1)
        if (knat .eq. 'ACCE') inat = 1
        if (knat .eq. 'VITE') inat = 2
        if (knat .eq. 'DEPL') inat = 3
!
        do id = 1, 3
            dirspe(id) = xnorm*dirspe(id)
            if (abs(dirspe(id)) .gt. epsi) then
                if (ndir(id) .ne. 0) then
                    ier = ier+1
                    call utmess('E', 'SEISME_5')
                    goto 10
                else
                    ndir(id) = 1
                end if
                nature(id) = inat
            end if
        end do
!
10      continue
    end do
!
    if (ier .ne. 0) then
        call utmess('F', 'SEISME_6')
    end if
!
!     --- INTERPOLATION DES SPECTRES ---
    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
        call utmess('I', 'SEISME_53')
    end if
    do im = 1, nbmode
        ii = 0
        amor = amort(im)
        call rsadpa(momec, 'L', 1, nopara(2), nordr(im), &
                    0, sjv=ival, istop=0)
        omega2 = zr(ival)
        omega = sqrt(omega2)
        freq = uns2pi*omega
        valpu(1) = amor
        valpu(2) = freq
        if (corfre) then
            correc = sqrt(un-amor*amor)
        else
            correc = 1.
        end if
        do id = 1, 3
            if (ndir(id) .eq. 1) then
                call fointe('F ', nomspe(id), 2, nompu, valpu, &
                            resu, ier)
                coef = dirspe(id)*echspe(id)
                if (nature(id) .eq. 1) then
                    valspe(id, im) = resu*coef*correc
                else if (nature(id) .eq. 2) then
                    valspe(id, im) = resu*coef*omega*correc
                else
                    valspe(id, im) = resu*coef*omega2*correc
                end if
                if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
                    if (ii .eq. 0) then
                        ii = 1
                        call utmess('I', 'SEISME_54', si=im, sk=dir(id), nr=3, &
                                    valr=[freq, amor, valspe(id, im)])
                    else
                        call utmess('I', 'SEISME_55', sk=dir(id), sr=valspe(id, im))
                    end if
                end if
            end if
        end do
    end do
!
!     --- VALEURS ASYMPTOTIQUES DES SPECTRES ---
    if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
        call utmess('I', 'SEISME_56')
    end if
    do id = 1, 3
        if (ndir(id) .eq. 1) then
            amor = amort(nbmode)
            valpu(1) = amor
            valpu(2) = fcoup
            omega = deuxpi*fcoup
            if (corfre) then
                correc = sqrt(un-amor*amor)
            else
                correc = 1.
            end if
            call fointe('F ', nomspe(id), 2, nompu, valpu, &
                        resu, ier)
            coef = dirspe(id)*echspe(id)
            if (nature(id) .eq. 1) then
                asyspe(id) = resu*coef*correc
            else if (nature(id) .eq. 2) then
                asyspe(id) = resu*coef*omega*correc
            else
                asyspe(id) = resu*coef*omega*omega*correc
            end if
            if (niveau .eq. 'TOUT     ' .or. niveau .eq. 'SPEC_OSCI') then
                call utmess('I', 'SEISME_57', sk=dir(id), sr=asyspe(id))
            end if
        end if
    end do
!
    call jedema()
end subroutine
