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
subroutine nm1dco(fami, kpg, ksp, option, imate, &
                  materi, e, sigm, epsm, deps, &
                  vim, sigp, vip, dsde, crildc, &
                  codret)
!
!
    implicit none
! ----------------------------------------------------------------------
!          PLASTICITE VON MISES ISOTROPE BILINEAIRE MONODIM
!          ON PEUT AVOIR T0 DIFF TREF
!
!
! IN  T        : TEMPERATURE PLUS
! IN  TM       : TEMPERATURE MOINS
! IN  E        : MODULE D EG
! IN  ET       : PENTE D ECROUISSAGE
! IN  ALPH     : COEF DILAT THERMIQUE
! IN  SY       : LIMITE D ELASTICITE INITIALE
!
! IN  SIGM    : CONTRAINTE AU TEMPS MOINS
!               UTILISE UNIQUEMENT POUR EVALUER DSDEM
! IN  DEPS    : DEFORMATION  TOTALE PLUS - DEFORMATION TOTALE MOINS
! IN  EPSPM   : DEFORMATION  PLASTIQUE MOINS
! IN  PM      : DEFORMATION  PLASTIQUE CUMULEE MOINS
!
! OUT SIGP     : CONTRAINTES PLUS
! OUT EPSP    : DEFORMATION  PLASTIQUE PLUS
! OUT P       : DEFORMATION  PLASTIQUE CUMULEE PLUS
! OUT DSDE    : DSIG/DEPS
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    real(kind=8) :: sigm, deps, pm, vim(*), vip(*), epspm, corrm
    real(kind=8) :: sigp, dsde, resi, crildc(*)
    character(len=16) :: option
    character(len=*) :: fami, materi
    integer(kind=8) :: imate, codret, kpg, ksp
!     ------------------------------------------------------------------
!     VARIABLES LOCALES
!     ------------------------------------------------------------------
    real(kind=8) :: epsm
    integer(kind=8) :: codres(1)
    real(kind=8) :: e, sy, dc, v, k, m
    real(kind=8) :: epsilf, epsd, epsc, d, p, epsp, ecr, fplas, indi
    real(kind=8) :: dfds, dfpds, dfdecr, difecr, lambp, fd, var1
    real(kind=8) :: var2, var3, rv, fini, fplas2, b, val(1)
    aster_logical :: dconv, pconv, melas
    integer(kind=8) :: iter, itemax, i, j, ibid
    pm = vim(1)
    epspm = vim(1)
    d = vim(2)
    codret = 0
    indi = 0.d0
!
!
! --- CARACTERISTIQUES ECROUISSAGE LINEAIRE
    call rcvalb(fami, kpg, ksp, '+', imate, &
                materi, 'CORR_ACIER', 0, ' ', [0.d0], &
                1, 'D_CORR', val, codres, 1)
    dc = val(1)
    call rcvalb(fami, kpg, ksp, '+', imate, &
                materi, 'CORR_ACIER', 0, ' ', [0.d0], &
                1, 'ECRO_K', val, codres, 1)
    k = val(1)
    call rcvalb(fami, kpg, ksp, '+', imate, &
                materi, 'CORR_ACIER', 0, ' ', [0.d0], &
                1, 'ECRO_M', val, codres, 1)
    m = val(1)
    call rcvalb(fami, kpg, ksp, '+', imate, &
                materi, 'CORR_ACIER', 0, ' ', [0.d0], &
                1, 'SY', val, codres, 1)
    sy = val(1)
    call rcvalb(fami, kpg, ksp, '+', imate, &
                materi, 'ELAS', 0, ' ', [0.d0], &
                1, 'NU', val, codres, 1)
    v = val(1)
!
! --- PARAMETRES DE CONVERGENCE
    resi = crildc(3)
    itemax = nint(crildc(1))
!
    call rcvarc('F', 'CORR', '-', fami, kpg, &
                ksp, corrm, ibid)
    if (corrm .le. 15.d0) then
        epsc = 2.345d-1-(1.11d-2*corrm)
    else
        epsc = 5.1d-2-(6.d-4*corrm)
    end if
!       ENDIF
!
!    DEFORMATION PLASTIQUE DE DEBUT D'ENDOMMAGMENT
    epsd = 0.8d0*epsc
!    RV LE PARAMETRE QUI DEPEND DU TAUX DE TRIAXIALITE
    var1 = 1.d0+v
    var2 = 1.d0-(2.d0*v)
    var3 = ((1.d0/3.d0)**2.d0)
    rv = (((2.d0/3.d0)*var1)+(3.d0*var2*var3))
    epsilf = epsm+deps
    epsp = epspm
    p = pm
    sigp = sigm
    dconv = .false.
    melas = (option .eq. 'RIGI_MECA_ELAS') .or.&
     &      (option .eq. 'FULL_MECA_ELAS')
    if ((option .eq. 'FULL_MECA') .or. (option .eq. 'RAPH_MECA')) then
!
        iter = 0
        do i = 1, itemax
            if (.not. dconv) then
                iter = iter+1
!
!    *******ELASTICITE********************
                sigp = e*(epsilf-epsp)
!JMP         SIGP =(1.D0-D)* E*(EPSILF-EPSP)
!
                ecr = k*(p**(1.d0/m))
                fini = ((abs(sigp)/(1.d0-d))-ecr-sy)
                fplas = fini
                if (fini .le. 0.d0) then
                    vip(3) = 0.d0
                    dconv = .true.
                else
                    vip(3) = 1.d0
                    pconv = .false.
!
!    ******PLASTICITE**********************
                    do j = 1, itemax
                        if (.not. pconv) then
                            dfds = (1.d0/(1.d0-d))
                            dfpds = (1.d0/(1.d0-d))
                            dfdecr = -1.d0
                            difecr = ((k/m)*((sigp/((1.d0-d)*k))-(sy/k))**(1.d0-m))
                            lambp = (fplas/((dfds*e*dfpds)-(dfdecr*difecr)))
                            epsp = epsp+lambp*dfpds
                            p = p+(lambp/(1.d0-d))
                            sigp = sigp-((e*lambp)/(1.d0-d))
                            ecr = k*(p**(1.d0/m))
                            fplas = ((abs(sigp)/(1.d0-d))-ecr-sy)
                            pconv = (abs(fplas/sy) .le. resi)
                        else
                            goto 141
                        end if
                    end do
141                 continue
                    if (j .ge. itemax) then
                        call utmess('I', 'MODELISA5_40')
                        codret = 1
                        goto 999
                    end if
                end if
!
            end if
!
!    *****ENDOMMAGEMENT*********************
            fd = epsp-epsd
            if (fd .le. 0.d0) then
                dconv = .true.
                goto 142
            else
                d = (dc*((rv*epsp)-epsd))/(epsc-epsd)
                fplas2 = ((abs(sigp)/(1.d0-d))-ecr-sy)
                if (iter .ne. 1) then
                    fd = fplas2
!  CHANGEMENT DE CRITERE
!              DCONV = (ABS(FD/FDINI) .LE. RESI)
                    dconv = (abs(fd/sy) .le. resi)
                    if (dconv) goto 142
                end if
            end if
            if (d .gt. 0.99d0) then
                dconv = .true.
                sigp = 0.d0
                goto 142
            end if
        end do
142     continue
!
        if (i .ge. itemax) then
            call utmess('I', 'MODELISA5_41')
            codret = 1
            goto 999
        end if
        vip(1) = p
        vip(2) = d
    end if
!
    if (option .eq. 'RIGI_MECA_TANG') then
        p = vim(1)
        d = vim(2)
        indi = vim(3)
    else if (option .eq. 'FULL_MECA') then
        indi = vip(3)
    end if
!
    if (indi .lt. 0.5d0) then
        dsde = e
    else
        if (d .gt. 0.d0) then
            if (option .eq. 'RIGI_MECA_TANG') then
                b = e+k/m*(1.d0-d)*p**((1.d0/m)-1.d0)-sigm/(1.d0-d)*(dc*rv/(epsc-epsd))
                dsde = e/b*(k/m*(1.d0-d)*p**((1.d0/m)-1.d0)-sigm/(1.d0-d)*(dc*rv/(epsc-epsd)))
            else
                b = e+k/m*(1.d0-d)*p**((1.d0/m)-1.d0)-sigp/(1.d0-d)*(dc*rv/(epsc-epsd))
                dsde = e/b*(k/m*(1.d0-d)*p**((1.d0/m)-1.d0)-sigp/(1.d0-d)*(dc*rv/(epsc-epsd)))
            end if
        else
            if (option .eq. 'RIGI_MECA_TANG') then
                dsde = ( &
                       ( &
                       k*(1.d0/m)*(p**((1.d0/m)-1.d0)))/(1.d0+(((k*(1.d0/m))/e)*(p**((1.d0/m)-1&
                       &.d0))) &
                       ) &
                       )
            else
                dsde = ( &
                       ( &
                       k*(1.d0/m)*(p**((1.d0/m)-1.d0)))/(1.d0+(((k*(1.d0/m))/e)*(p**((1.d0/m)-1&
                       &.d0))) &
                       ) &
                       )
            end if
        end if
!
!     CAS RIGI_MECA_ELAS ET FULL_MECA_ELAS AVEC ENDOMMAGEMENT
    end if
    if (melas) dsde = (1.d0-d)*dsde
!
999 continue
end subroutine
