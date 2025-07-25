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
subroutine pdadom(xm0, xm2, xm4, dom)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/erfcam.h"
#include "asterc/r8pi.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/limend.h"
#include "asterfort/rccome.h"
#include "asterfort/rcpare.h"
#include "asterfort/rcvale.h"
#include "asterfort/utmess.h"
    real(kind=8) :: xm0, xm2, xm4, dom
!
!     CALCUL DOMMAGE EN FREQUENTIEL
!     ----------------------------------------------------------------
!
!
    integer(kind=8) :: icodwo, icodre(6)
    integer(kind=8) :: icodba, icodhs
    character(len=8) :: nommat, method, mecomp, nompar, kcorre, kbid
    character(len=16) :: nomres(6), cara
    character(len=32) :: pheno
    real(kind=8) :: delta, rvke, alpha, pi, salt, x, val(6), re(1)
    real(kind=8) :: valmin, valmax, pas, xireg, rundf, nrupt(1)
    integer(kind=8) :: ibask, ifonc, ihosin, nbval
    aster_logical :: endur
!
!     ----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: ipoint, nbpar, nbpoin
    real(kind=8) :: rbid, x1, x2, xnpoin
    real(kind=8) :: xp, y, y1, yd1, yd2, ypic1, ypic2
    real(kind=8), pointer :: dispics(:) => null()
    real(kind=8), pointer :: wohler2(:) => null()
!
!-----------------------------------------------------------------------
    rbid = 0.d0
    ifonc = 0
    ibask = 0
    ihosin = 0
    pi = r8pi()
    rundf = r8vide()
    call getvtx(' ', 'DOMMAGE', scal=method, nbret=nbval)
    call getvid(' ', 'MATER', scal=nommat, nbret=nbval)
    pheno = 'FATIGUE'
    call rccome(nommat, pheno, icodre(1))
    if (icodre(1) .eq. 1) then
        call utmess('F', 'FATIGUE1_24')
    end if
    cara = 'WOHLER'
    call rcpare(nommat, pheno, cara, icodwo)
    cara = 'A_BASQUIN'
    call rcpare(nommat, pheno, cara, icodba)
    cara = 'A0'
    call rcpare(nommat, pheno, cara, icodhs)
    if (icodwo .eq. 0) then
        ifonc = 1
    else if (icodba .eq. 0) then
        ibask = 1
    else if (icodhs .eq. 0) then
        ihosin = 1
    else
        call utmess('F', 'FATIGUE1_34')
    end if
!
!----  DEFINITION DES BORNES INTEGRATION
!
    call getvtx(' ', 'COMPTAGE', scal=mecomp, nbret=nbval)
    if (mecomp .eq. 'PIC     ' .and. xm4 .eq. rundf) then
        call utmess('F', 'FATIGUE1_35')
    end if
    if (ihosin .ne. 0) then
        nomres(6) = 'SL'
        nbpar = 0
        nompar = ' '
        call rcvale(nommat, 'FATIGUE', nbpar, nompar, [rbid], &
                    1, nomres(6), val(6), icodre(6), 2)
        valmin = val(6)
        valmax = 10*sqrt(xm0)
    else
        valmin = 0.d0
        valmax = 10*sqrt(xm0)
    end if
    pas = (valmax-valmin)/300.d0
    if (pas .eq. 0.0d0) then
        call utmess('F', 'FATIGUE1_36')
    end if
    xnpoin = (valmax-valmin)/pas
    nbpoin = int(xnpoin)+1
!
!------- CALCUL DES POINTS INTEGRATION
!
    if (xm2 .eq. 0.d0) then
        call utmess('F', 'FATIGUE1_37')
    end if
    if (mecomp .eq. 'PIC' .and. xm4 .eq. 0.d0) then
        call utmess('F', 'FATIGUE1_38')
    end if
    AS_ALLOCATE(vr=dispics, size=2*nbpoin)
    if (mecomp .eq. 'PIC     ') xireg = sqrt(xm2*xm2/xm0/xm4)
    do ipoint = 1, nbpoin
        x1 = valmin+(ipoint-1)*pas
        if (mecomp .eq. 'PIC     ') then
            alpha = xireg*x1/((sqrt(1.d0-xireg*xireg))*(sqrt(xm0)))
            alpha = (-1.d0/sqrt(2.d0))*alpha
            y1 = sqrt(1-xireg*xireg)*exp(-x1*x1/(2.d0*xm0*(1.d0-xireg* &
                                                           xireg)))
            xp = sqrt(pi/2.d0)*erfcam(alpha)
            y1 = y1+((xireg*x1/sqrt(xm0))*exp(-x1*x1/(2.d0*xm0)))*(xp)
            y1 = (sqrt(xm4)/(sqrt(xm2)*sqrt(xm0)))*y1
            y1 = (1.d0/(2.d0*pi))*(1.d0/sqrt(2.d0*pi))*y1
        else if (mecomp .eq. 'NIVEAU  ') then
            y1 = (1.d0/(2.d0*pi))*sqrt(xm2/(xm0*xm0*xm0))
            y1 = y1*x1*exp(-x1*x1/(2.d0*xm0))
        end if
        dispics(ipoint) = x1
        dispics(nbpoin+ipoint) = y1
    end do
!
!---------CORRECTION ELASTO-PLASTIQUE
!
    call getvtx(' ', 'CORR_KE', nbval=0, nbret=nbval)
    if (nbval .ne. 0) then
        call getvtx(' ', 'CORR_KE', scal=kcorre, nbret=nbval)
        call getvid(' ', 'MATER', scal=nommat, nbret=nbval)
        if (kcorre .eq. 'RCCM') then
            nomres(1) = 'N_KE'
            nomres(2) = 'M_KE'
            nomres(3) = 'SM'
            nbpar = 0
            nompar = ' '
            call rcvale(nommat, 'RCCM', nbpar, nompar, [rbid], &
                        3, nomres(1), val(1), icodre(1), 2)
            do ipoint = 1, nbpoin
                delta = dispics(ipoint)
                if (delta .le. 3.d0*val(3)) then
                    rvke = 1.d0
                elseif (delta .gt. 3.d0*val(3) .and. delta .lt. 3.d0*val(2)* &
                        val(3)) then
                    rvke = 1.d0+((1-val(1))/(val(1)*(val(2)-1)))*((delta/(3.d0*val(3)))-1.d0 &
                                                                  )
                else if (delta .ge. 3*val(2)*val(3)) then
                    rvke = 1.d0/val(1)
                end if
                dispics(ipoint) = rvke*dispics(ipoint)
            end do
        end if
    end if
!
! ----- INTERPOLATION
!
    if (method .eq. 'WOHLER') then
!
! --- INTERPOLATION SUR LA COURBE DE WOHLER ---
!
        AS_ALLOCATE(vr=wohler2, size=nbpoin)
        if (ifonc .ne. 0) then
            nomres(1) = 'WOHLER'
            nbpar = 1
            pheno = 'FATIGUE'
            nompar = 'SIGM'
            do ipoint = 1, nbpoin
                delta = dispics(ipoint)
                call limend(nommat, delta, 'WOHLER', kbid, endur)
                if (endur) then
                    wohler2(ipoint) = 0.d0
                else
                    call rcvale(nommat, pheno, nbpar, nompar, [delta], &
                                1, nomres(1), nrupt(1), icodre(1), 2)
                    wohler2(ipoint) = 1.d0/nrupt(1)
                end if
            end do
        else if (ibask .ne. 0) then
            nompar = ' '
            nbpar = 0
            nomres(1) = 'A_BASQUIN'
            nomres(2) = 'BETA_BASQUIN'
            call rcvale(nommat, 'FATIGUE', nbpar, nompar, [rbid], &
                        2, nomres, val, icodre, 2)
            do ipoint = 1, nbpoin
                wohler2(ipoint) = val(1)*dispics(ipoint)**val(2)
            end do
        else if (ihosin .ne. 0) then
            nomres(1) = 'E_REFE'
            nomres(2) = 'A0'
            nomres(3) = 'A1'
            nomres(4) = 'A2'
            nomres(5) = 'A3'
            nomres(6) = 'SL'
            nbpar = 0
            nompar = ' '
            call rcvale(nommat, 'FATIGUE', nbpar, nompar, [rbid], &
                        6, nomres, val, icodre, 2)
            nomres(1) = 'E'
            call rcvale(nommat, 'ELAS', nbpar, nompar, [rbid], &
                        1, nomres, re(1), icodre, 2)
            do ipoint = 1, nbpoin
                salt = (val(1)/re(1))*dispics(ipoint)
                if (salt .ge. val(6)) then
                    x = log10(salt)
                    y = val(2)+val(3)*x+val(4)*(x**2)+val(5)*(x**3)
                    wohler2(ipoint) = 1.d0/(10.d0**y)
                else
                    wohler2(ipoint) = 0.d0
                end if
            end do
        end if
    end if
!
! -- CALCUL INTEGRALE
!
    dom = 0.d0
    do ipoint = 2, nbpoin
        x2 = dispics(ipoint)
        x1 = dispics(ipoint-1)
        yd2 = wohler2(ipoint)
        yd1 = wohler2(ipoint-1)
        ypic2 = dispics(nbpoin+ipoint)
        ypic1 = dispics(nbpoin+ipoint-1)
        dom = dom+(yd2*ypic2+yd1*ypic1)*(x2-x1)/2.d0
    end do
!
    AS_DEALLOCATE(vr=dispics)
    AS_DEALLOCATE(vr=wohler2)
!
end subroutine
