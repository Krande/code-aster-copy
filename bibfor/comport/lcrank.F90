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
! aslint: disable=W0413
!
!
subroutine lcrank(ndim, typmod, imate, option, tmpm, &
                  tmpp, dstrai0, stresm0, stres, vim, &
                  vip, dsidep, codret)
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE RANKINE
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  IMATE   : NATURE DU MATERIAU
! IN  DSTRAI0 : INCREMENT DE DEFORMATION
! IN  STRESM0 : CONTRAINTE EN T-
! IN  VIM     : VARIABLES INTERNES EN T-
! IN  TMPP    : TEMPERATURE EN T+
! IN  TMPREF  : TEMPERATURE DE REFERENCE
! IN  OPTION  : OPTION DEMANDEE
!                 RIGI_MECA_TANG -> DSIDEP
!                 FULL_MECA      -> STRES DSIDEP VIP
!                 RAPH_MECA      -> STRES        VIP
!
! OUT STRES   : CONTRAINTE
! OUT VIP     : VARIABLES INTERNES:
!               1 -> DEFORMATION PLASTIQUE VOLUMIQUE
!               2 -> NORME DE LA DEFORMATION PLASTIQUE DEVIATORIQUE
!               3 -> INDICATEUR DE PLASTICITE = | 0 : elastique
!                                               | 1 : 1 mecanisme actif
!                                               | 2 : 2 mecanismes actifs
!                                               | 3 : 3 mecanismes actifs
!               4 -> COMPOSANTE MAJEURE DE L'INCREMENT
!                    DE DEFORMATION PLASTIQUE : DGAMA1
!               5 -> COMPOSANTE INTERMEDIAIRE DE L'INCREMENT
!                    DE DEFORMATION PLASTIQUE : DGAMA2
!               6 -> COMPOSANTE MINEURE DE L'INCREMENT
!                    DE DEFORMATION PLASTIQUE : DGAMA3
! OUT DSIDEP  : MATRICE TANGENTE
! OUT CODRET  : CODE RETOUR = | 0: OK
!                             | 1: NOOK
!
! ----------------------------------------------------------------------
    implicit none
! ======================================================================
!
    character(len=8) :: typmod(*)
    character(len=16) :: option
    integer(kind=8) :: ndim, imate, codret
    real(kind=8) :: tmpm, tmpp
    real(kind=8) :: dstrai(6), dstrai0(6)
    real(kind=8) :: stresm(6), stresm0(6), stres(6)
    real(kind=8) :: vim(*), vip(*)
    real(kind=8) :: dsidep(6, 6)
!
#include "asterf_types.h"
#include "asterfort/bptobg.h"
#include "asterfort/mcpstr.h"
#include "asterfort/mctgel.h"
#include "asterfort/ratg2d.h"
#include "asterfort/ratg3d.h"
#include "asterfort/rcvala.h"
!
! Declaration of constant parameters
    integer(kind=8) :: mmax, nmax
    parameter(mmax=3, nmax=6)
!
! Declaration of real type variables
    real(kind=8) :: rprops(nmax), strait(nmax), strain(nmax)
    real(kind=8) :: strest(nmax), young, poiss
    real(kind=8) :: cohe, eigprj(mmax, mmax), pstrs(mmax)
    real(kind=8) :: gmodu, bulk, r2g, r4g
    real(kind=8) :: r2bulk, r1d3, eetv
    real(kind=8) :: pt, eetvd3, pstrs1, pstrs2, pstrs3, smct
    real(kind=8) :: phia, phib, phic, res
    real(kind=8) :: consta, ddgama, dgama, dgamc
    real(kind=8) :: smcta, smctb, smctc, constb, drvaa
    real(kind=8) :: drvab, drvba, drvbb, r1ddet, ddgamb, dgamb
    real(kind=8) :: resnor, factor, aux1, aux2, aux3
    real(kind=8) :: ddepv, s1, s2, s3, sttv
    real(kind=8) :: r0, r1, r2, r3, r4, small, tol, sqr
    real(kind=8) :: apex, edge, ifplas, sufail
    real(kind=8) :: tr(nmax), p, r1nu, r1ny
!
! Declaration of integer type variables
    integer(kind=8) :: icode(nmax), ndt, ndi, i, ii, jj, mm
    integer(kind=8) :: mxiter, indic(mmax), nbmeca
!
! Declaration of integer type variables
    aster_logical :: epflag
    aster_logical :: tridim, outofp
!
! Declaration of character type variables
    character(len=8) :: mod
    character(len=16) :: nomres(3)
!
! Declaration of constant variables
    data r0, r1, r2, r3, r4, small, tol, sqr/ &
        0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 1.d-06, 1.d-6, &
        1.4142135623730951d0/
    data mxiter/50/
!
! Declaration of Common space variables
    common/debug/epflag
!
! Is the problem tridimensional?
    if (ndim .eq. 3) then
        tridim = .true.
    else
        tridim = .false.
    end if
!
! Modelisation type:
! 3D, D_PLAN, C_PLAN, AXIS
    mod = typmod(1)
! =================================================================
! --- NOMBRE DE COMPOSANTES ---------------------------------------
! =================================================================
    if (mod(1:2) .eq. '3D' .or. mod(1:6) .eq. 'C_PLAN') then
        ndt = 6
        ndi = 3
    else if (mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        ndt = 4
        ndi = 3
    end if
!
    if (mod(1:6) .eq. 'C_PLAN') then
        outofp = .false.
    else
        outofp = .true.
    end if
!
    epflag = .false.
!
!    write(6,*)
!     write(6,'(A,6(I2))')'! * NDI,NDT =',ndi,ndt
!
! Reading material linear elastic properties
! ALPHA must be constant
    nomres(1) = 'E       '
    nomres(2) = 'NU      '
    nomres(3) = 'ALPHA   '
    call rcvala(imate, ' ', 'ELAS', 0, '   ', &
                [tmpp], 3, nomres, rprops, icode, &
                2)
!
! Reading material Rankine properties
    nomres(1) = 'SIGMA_T'
    call rcvala(imate, ' ', 'RANKINE', 0, '   ', &
                [tmpp], 1, nomres, rprops(4), icode(4), &
                2)
!
! Remove SQRT(2) from extradiagonal terms of input strain and stress
! ------------------------------------------------------------------
    dstrai(1:nmax) = dstrai0(1:nmax)
    stresm(1:nmax) = stresm0(1:nmax)
!
! Thermal deformation rate
! -------------------------
!     write(6,'(4(A,E12.5))')'! * TEMP+=',tmpp,' TEMP-=',tmpm,' ALPHA=',rprops(3)
    if ((.not. isnan(tmpm)) .and. (.not. isnan(tmpp)) .and. (rprops(3) .gt. r0)) then
        do i = 1, ndi
            dstrai(i) = dstrai(i)-rprops(3)*(tmpp-tmpm)
        end do
    end if
!
    dstrai(4) = dstrai(4)/sqr
    stresm(4) = stresm(4)/sqr
    if (tridim) then
        dstrai(5) = dstrai(5)/sqr
        dstrai(6) = dstrai(6)/sqr
        stresm(5) = stresm(5)/sqr
        stresm(6) = stresm(6)/sqr
    end if
!
! Initialize some algorithmic and internal variables
    dgama = r0
    dgamb = r0
    dgamc = r0
    ifplas = r0
    sufail = r0
    edge = r0
    apex = r0
!
! Set some material properties
    young = rprops(1)
    poiss = rprops(2)
    cohe = rprops(4)
!
!      if (epflag) then
!         write(6,'(4(A,E12.5))')&
!         '! * YOUNG  =',young,' POISSON =',poiss,' COHESION =',cohe,&
!         ' ALPHA =',rprops(3)
!      endif
!
! Set some constants
    gmodu = young/(r2*(r1+poiss))
    bulk = young/(r3*(r1-r2*poiss))
    r2g = r2*gmodu
    r4g = r4*gmodu
    r2bulk = r2*bulk
    r1d3 = r1/r3
!
! Compute elastic trial state
! ---------------------------
    strest(:) = r0
    eetv = dstrai(1)+dstrai(2)+dstrai(3)
    pt = bulk*eetv
    eetvd3 = eetv*r1d3
    strest(1) = stresm(1)+r2g*(dstrai(1)-eetvd3)+pt
    strest(2) = stresm(2)+r2g*(dstrai(2)-eetvd3)+pt
    strest(3) = stresm(3)+r2g*(dstrai(3)-eetvd3)+pt
    strest(4) = (stresm(4)+r2g*dstrai(4))
    if (tridim) then
        strest(5) = (stresm(5)+r2g*dstrai(5))
        strest(6) = (stresm(6)+r2g*dstrai(6))
    end if
!
! Compute elastic trial strain
! ----------------------------
! Elastic trial volumetric strain and pressure stress
! Be careful that we have the following form of input strain:
!
!          (Sqrt(2)*EPXY Sqrt(2)*EPXZ Sqrt(2)*EPYZ)
!
    strait(:) = r0
    sttv = strest(1)+strest(2)+strest(3)
    r1nu = (r1+poiss)/young
    r1ny = poiss/young
    strait(1) = r1nu*strest(1)-r1ny*sttv
    strait(2) = r1nu*strest(2)-r1ny*sttv
    strait(3) = r1nu*strest(3)-r1ny*sttv
    strait(4) = r1nu*strest(4)
    if (tridim) then
        strait(5) = r1nu*strest(5)
        strait(6) = r1nu*strest(6)
    end if
!
!
! Compute eigen-stress for prediction: pstrs1, pstrs2, pstrs3
    call mcpstr(strest, tridim, pstrs, eigprj, ii, &
                jj, mm, codret)
!
    pstrs1 = pstrs(ii)
    pstrs2 = pstrs(mm)
    pstrs3 = pstrs(jj)
    indic(1) = ii
    indic(2) = mm
    indic(3) = jj
!
! ==================================================================
!
! 0. Is the behaviour plastic?
!
! ==================================================================
!
! Compute the number of active mechanisms: nbmeca
! --------------------------------------------------------------
    nbmeca = 0
    do i = 1, 3
        smct = pstrs(indic(i))
        phia = smct-cohe
        res = phia
        if (cohe .ne. r0) res = res/abs(cohe)
!
        if (res .gt. tol) then
            nbmeca = nbmeca+1
        end if
    end do
!
!     if (epflag) then
!         write(6,'(A,I2)')'! * NBMECA =',nbmeca
!     endif
!
    if (nbmeca .eq. 2) then
        goto 30
    else if (nbmeca .eq. 3) then
        goto 50
    else if (nbmeca .eq. 0) then
        goto 75
    end if
!
! Plastic step: Apply return mapping
! ==================================
    ifplas = r1
! ==================================================================
!
! 1. Apply one-vector return mapping first (return to MAIN PLANE)
!
! ==================================================================
    smct = pstrs1
    phia = smct-cohe
    res = phia
    if (cohe .ne. r0) res = res/abs(cohe)
!
    consta = r4g*r1d3+bulk
    constb = -r2g*r1d3+bulk
!
! Compute Plastic multiplier: DGAMA
    ddgama = phia/consta
    dgama = dgama+ddgama
!
! Compute new residual
    phia = smct-consta*dgama-cohe
!
! Check convergence
    resnor = abs(phia)
    factor = abs(smct)
    if (factor .ne. r0) resnor = resnor/factor
!
    if (resnor .le. tol) then
!
! Check validity of 1-vector return (check sextant of converged stress)
        s1 = pstrs1-consta*dgama
        s2 = pstrs2-constb*dgama
        s3 = pstrs3-constb*dgama
        p = (s1+s2+s3)*r1d3
!
!         if (epflag) then
!             write(6, '(A)') '!'
!             write(6, '(A)') '!!! RETURN TO MAIN PLANE !!!'
!             write(6, '(A)') '!'
!             write(6, '(4(A,E15.8))') &
!             '! * S1=', s1, ' S2=', s2, ' S3=', s3
!         endif
!
        goto 70
    end if
! failure of stress update procedure
    sufail = r1
    codret = 1
!
30  continue
! ==================================================================
!
! 2. Apply two-vector return mapping to appropriate EDGE
!
! ==================================================================
    ifplas = r2
!
    dgama = r0
    smcta = pstrs1
    smctb = pstrs2
    phia = smcta-cohe
    phib = smctb-cohe
    consta = r4g*r1d3+bulk
    constb = -r2g*r1d3+bulk
!
! Compute residual derivative matrix
    drvaa = -consta
    drvab = -constb
    drvba = -constb
    drvbb = -consta
!
! Compute Newton-Raphson increment and update variables DGAMA and DGAMB
    r1ddet = r1/(drvaa*drvbb-drvab*drvba)
    ddgama = (-drvbb*phia+drvab*phib)*r1ddet
    ddgamb = (drvba*phia-drvaa*phib)*r1ddet
    dgama = dgama+ddgama
    dgamb = dgamb+ddgamb
!
! Compute new residual
    phia = smcta-consta*dgama-constb*dgamb-cohe
    phib = smctb-constb*dgama-consta*dgamb-cohe
!
! Check convergence
    resnor = (abs(phia)+abs(phib))
    factor = (abs(smcta)+abs(smctb))
    if (factor .ne. r0) resnor = resnor/factor
!
    if (resnor .le. tol) then
! Check validity of 2-vector return to edge
        aux1 = consta
        aux2 = constb
        aux3 = constb
        s1 = pstrs1-aux1*dgama-aux2*dgamb
        s2 = pstrs2-aux2*dgama-aux1*dgamb
        s3 = pstrs3-aux3*(dgama+dgamb)
        edge = r1
        p = (s1+s2+s3)*r1d3
!
!         if (epflag) then
!             write(6, '(A)') '!'
!             write(6, '(A)') '!!! RETURN TO EDGE !!!'
!             write(6, '(A)') '!'
!             write(6, '(4(A,E15.8))') &
!             '! * S1=',s1,' S2=',s2,' S3=',s3
!         endif
!
        goto 70
    end if
! failure of stress update procedure
    sufail = r2
    codret = 1
!        if (epflag) write(6,'(A)')'LCMOHR :: Projection Failure on edge'
!
! (-_-) goto 999 if a problem occured
    goto 999
50  continue
! ==================================================================
!
! 3. Apply multi-vector return mapping to APEX
!
! ==================================================================
    ifplas = r3
!
    smcta = pstrs1
    smctb = pstrs2
    smctc = pstrs3
    phia = smcta-cohe
    phib = smctb-cohe
    phic = smctc-cohe
    consta = r4g*r1d3+bulk
    constb = -r2g*r1d3+bulk
    drvaa = consta*r1d3/r2g/bulk
    drvbb = constb*r1d3/r2g/bulk
    drvab = drvaa+drvbb
!
! Compute Plastic Multipliers
    dgama = drvab*phia-drvbb*(phib+phic)
    dgamb = drvab*phib-drvbb*(phia+phic)
    dgamc = drvab*phic-drvbb*(phib+phia)
!
!     if (epflag) then
!         write(6, '(A)') '!'
!         write(6, '(A)') '!!! RETURN TO APEX !!!'
!         write(6, '(A)') '!'
!         write(6, '(4(A,E15.8))')&
!         '! * CONSTA =',consta,' CONSTB =',constb
!         write(6, '(4(A,E15.8))')&
!         '! * DGAMA =',dgama,' DGAMB =',dgamb,' DGAMC =',dgamc
!         write(6, '(4(A,E15.8))')&
!         '! * PHIA =',phia,' PHIB =',phib,' PHIC =',phic
!     endif
! Compute new residual
    phia = smcta-consta*dgama-constb*(dgamb+dgamc)-cohe
    phib = smctb-consta*dgamb-constb*(dgama+dgamc)-cohe
    phic = smctc-consta*dgamc-constb*(dgamb+dgama)-cohe
!
! Check convergence
    resnor = (abs(phia)+abs(phib)+abs(phic))
    factor = (abs(smcta)+abs(smctb)+abs(smctc))
    if (factor .ne. r0) resnor = resnor/factor
!
    if (resnor .le. tol) then
! Set initial guess for volumetric plastic strain increment DEPV
        res = cohe-pt
! volumetric plastic strain
        ddepv = -res/bulk
! final volumetric stress
        p = pt-bulk*ddepv
!
! update principal stresses
        s1 = p
        s2 = p
        s3 = p
        apex = r1
!
        goto 70
    end if
!
! Projection to apex failure
    sufail = r3
    codret = 1
!
! (-_-) goto 999 if a problem occured
    goto 999
! ==================================================================
!
! 4. END OF INTEGRATION: CONTINUE IF NO PROBLEM OCCURED
!
! ==================================================================
70  continue
! ------------------------------------------------------
!
! 4.1. update internal variables and output stress components
!
! ------------------------------------------------------
    tr(:) = r0
    tr(ii) = s1
    tr(jj) = s3
    tr(mm) = s2
    call bptobg(tr, stres, eigprj)
!
! INCREMENT de deformation plastique dans la base cartesienne
    tr(:) = r0
    tr(ii) = dgama
    tr(jj) = dgamc
    tr(mm) = dgamb
    call bptobg(tr, strain, eigprj)
!
! Attention: On re-projete les contraintes dans la base de Voigt
    do i = mmax+1, nmax
        stres(i) = sqr*stres(i)
        strain(i) = sqr*strain(i)
    end do
!
! Update internal variables
    if (apex .eq. r1) then
        vip(1) = vim(1)+p/bulk
        vip(2) = vim(2)+r2*r1d3* &
                 sqrt(dgama*dgama+dgamb*dgamb+dgamc*dgamc- &
                      dgama*dgamb-dgama*dgamc-dgamb*dgamc)
    else
        vip(1) = vim(1)+(dgama+dgamb)
        vip(2) = vim(2)+r2*r1d3* &
                 sqrt(dgama*dgama+dgamb*dgamb-dgama*dgamb)
    end if
!
! ------------------------------------------------------
!
! 4.2. Elastic step: update stress using linear elastic law
!
! ------------------------------------------------------
75  continue
    if (ifplas .eq. r0) then
!
!         if (epflag) then
!             write(6,'(A)')'!'
!             write(6,'(A)')'!!!'
!             write(6,'(A)')'!!!    ELASTIC STATE    !!!'
!             write(6,'(A)')'!!!'
!         endif
!
        do i = 1, mmax
            stres(i) = strest(i)
            stres(mmax+i) = sqr*strest(mmax+i)
        end do
!
! Update internal variables
        vip(1) = vim(1)
        vip(2) = vim(2)
!
        do i = 1, 6
            vip(3+i) = vim(3+i)
        end do
    else
!
        do i = 1, 6
            vip(3+i) = vim(3+i)+strain(i)
        end do
!
    end if
! ------------------------------------------------------
!
! 4.3. Update algorithmic variables
!      and Check criterium before exit
!
! ------------------------------------------------------
    vip(3) = ifplas
!
! recalcul du seuil en sortie
! ---------------------------
    if (ifplas .gt. r0) then
        res = s1-cohe
    else
        res = pstrs1-cohe
    end if
!
!
    if (cohe .ne. r0) res = res/abs(cohe)
!
    if (res .gt. tol) then
!
        sufail = r4
        codret = 1
        goto 999
!
    end if
!
!
    if (option(1:9) .eq. 'RIGI_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
! --------------------------------------
!
!   Compute Elastic Stiffness Matrix
!
! --------------------------------------
        tr(2) = rprops(1)
        tr(3) = rprops(2)
!
        call mctgel(dsidep, tr)
!
        if (ifplas .eq. r0) goto 888
!
! -----------------------------------------
!
! Plastic Stiffness Matrix in the principal base
! Compute DPSTRS needed to compute the general stiffness matrix
!
! =========================================
!
!           Tridimensionnal case
!
! =========================================
!
        if (sufail .eq. r0 .and. tridim) then
!
            call ratg3d(stres, strait, rprops, dsidep, ii, &
                        jj, mm, edge, apex, codret)
!
            if (codret .eq. 1) goto 999
!
! =========================================
!
!            Bidimensionnal case
!
! =========================================
!
        else if (sufail .eq. r0) then
!
!
            call ratg2d(stres, strait, rprops, dsidep, ii, &
                        jj, mm, edge, apex, outofp)
!
! Endif SUFAIL
        end if
! Endif OPTION
    end if
!
888 continue
    codret = 0
999 continue
!
end subroutine
