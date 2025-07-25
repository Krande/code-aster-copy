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

subroutine vpvers(eigsol, modes, checksd)
! VERIFICATION DE LA COHERENCE DES PARAMETRES SOLVEUR MODAL ET DES DONNEES/OBJETS SOUS-JACENTS
!
! RQ . SI OPTION='PLUS_GRANDE' ON MODIFIE LES VALEURS RAIDE/MASSE DE LA SD EIGENSOLVER
!      VIA VPECRI.
! -------------------------------------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
#include "asterf_types.h"
#include "asterc/cheksd.h"
#include "asterc/r8prem.h"
#include "asterfort/ajlagr.h"
#include "asterfort/assert.h"
#include "asterfort/exisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/utmess.h"
#include "asterfort/vpcrea.h"
#include "asterfort/vpecri.h"
#include "asterfort/vplecs.h"
#include "asterfort/vrrefe.h"
!
! --- INPUT
!
    character(len=8), intent(in) :: modes
    character(len=19), intent(in) :: eigsol
    aster_logical, intent(in) :: checksd
!
! --- OUTPUT
! None
!
! --- INPUT/OUTPUT
! None
!
! --- VARIABLES LOCALES
!
    integer(kind=8)           :: ibid, iret
    real(kind=8)      :: alpha, eps, freq1, rbid
    character(len=1)  :: appr
    character(len=8)  :: method
    character(len=9)  :: typevp
    character(len=14) :: matra, matrb, matrc
    character(len=16) :: modrig, optiof, typeqz, typres
    character(len=19) :: amor, masse, masse0, raide, raide0, numedd
    character(len=24) :: valk(2), k24buff
    aster_logical :: lc, lkr, lns, lpg, lqz
!
! -----------------------
! --- CORPS DE LA ROUTINE
! -----------------------
!
!
! --  INITS.
    eps = 1.d+4*r8prem()
!
! --  LECTURE DES PARAMETRES MODAUX
    call vplecs(eigsol, alpha_=alpha, freq1_=freq1, &
                appr_=appr, method_=method, typevp_=typevp, matra_=matra, &
                matrb_=matrb, matrc_=matrc, modrig_=modrig, optiof_=optiof, &
                typeqz_=typeqz, typres_=typres, amor_=amor, masse_=masse, raide_=raide, &
                lc_=lc, lkr_=lkr, lns_=lns, lpg_=lpg, lqz_=lqz)
!
! --  REGLES D'EXCLUSION
!
! --  OPTION DE CALCUL
    if ((optiof(1:4).eq.'TOUT').and.(.not.lqz)) call utmess('F', 'ALGELINE5_65', sk='CALC_'//typevp)

! --  MODES RIGIDES
    if ((modrig(1:11) .eq. 'MODE_RIGIDE') .and. (method(1:8) .ne. 'TRI_DIAG')) &
        call utmess('F', 'ALGELINE2_65')

! --  QEP
    if (lc) then
        valk(1) = matra
        valk(2) = matrc
        if (lpg) call utmess('F', 'ALGELINE5_82')
        if (optiof(1:5) .eq. 'BANDE') call utmess('F', 'ALGELINE2_66', nk=2, valk=valk)
        if (((appr .eq. 'I') .or. (appr .eq. 'C')) .and. &
            ((abs(freq1) .lt. eps) .or. (optiof(1:11) .eq. 'PLUS_PETITE'))) then
            call utmess('F', 'ALGELINE2_67')
        end if
        if (modrig(1:11) .eq. 'MODE_RIGIDE') call utmess('F', 'ALGELINE2_68', nk=2, valk=valk)
        if ((method(1:8) .eq. 'SORENSEN') .and. &
            ((abs(freq1) .lt. eps) .or. (optiof(1:11) .eq. 'PLUS_PETITE'))) then
            call utmess('F', 'ALGELINE2_71')
        end if
        if (method(1:6) .eq. 'JACOBI') call utmess('F', 'ALGELINE5_64', sk=matrc)
    end if
!
! --  MATRICE K COMPLEXE
    if (.not. lkr) then
        valk(1) = matra
        valk(2) = matrc
        if (lpg) call utmess('F', 'ALGELINE5_82')
        if ((method(1:8) .ne. 'SORENSEN') .and. (.not. lqz)) call utmess('F', 'ALGELINE2_69')
        if (optiof(1:5) .eq. 'BANDE') call utmess('F', 'ALGELINE2_66', nk=2, valk=valk)
        if (abs(freq1) .lt. eps) call utmess('F', 'ALGELINE2_70')
        if (modrig(1:11) .eq. 'MODE_RIGIDE') call utmess('F', 'ALGELINE2_68', nk=2, valk=valk)
        if (typres(1:10) .eq. 'MODE_FLAMB') call utmess('F', 'ALGELINE2_46', sk=matra)
    end if
!
! --  METHODE='SORENSEN'
    if (method(1:8) .eq. 'SORENSEN') then
        if ((abs(alpha) .lt. (1.2d0*eps)) .or. (abs(alpha) .gt. (0.83d0-eps))) then
            call utmess('E', 'ALGELINE2_64')
        end if
    end if
!
! --  METHODE='QZ'
    if (lqz) then
        if ((typeqz(1:5) .eq. 'QZ_QR') .and. &
            ((typres(1:10) .eq. 'FLAMBEMENT') .or. lc .or. (.not. lkr))) then
            valk(1) = matra
            valk(2) = matrc
            call utmess('F', 'ALGELINE5_60', nk=2, valk=valk)
        end if
    end if
!
! -- MATRICE K ET/OU M ET OU C NON SYMETRIQUE(S)
    if (lns) then
        if (lpg) call utmess('F', 'ALGELINE5_82')
        if ((.not. lqz) .and. (method(1:8) .ne. 'SORENSEN')) call utmess('F', 'ALGELINE5_69')
        if (.not. lkr) call utmess('F', 'ALGELINE5_70', sk=matra)
        if (optiof(1:5) .eq. 'BANDE') call utmess('F', 'ALGELINE4_39')
        if (modrig(1:11) .eq. 'MODE_RIGIDE') call utmess('F', 'ALGELINE4_40')
        if (typres(1:10) .eq. 'MODE_FLAMB') call utmess('F', 'ALGELINE4_41')
    end if
!
! --  TRAITEMENT PARTICULIER LIE A L'OPTION PLUS GRANDE
    if (lpg) then
        raide0 = masse
        masse0 = raide
        raide = 'MATR'
        masse = 'MATM'
        call ajlagr(masse0, raide0, raide)
        call mtdefs(masse, masse0, 'V', ' ')
        call mtcmbl(1, 'R', [1.d0], masse0, masse, &
                    'LAGR', ' ', 'ELIM=')
! --  ON MODIFIE QUELQUES VALEURS DE LA SD EIGENSOLVER
        k24buff = raide
        call vpecri(eigsol, 'K', 2, k24buff, rbid, ibid)
        k24buff = masse
        call vpecri(eigsol, 'K', 3, k24buff, rbid, ibid)
    end if
!
! --  COMPATIBILITE DES MODES (DONNEES ALTEREES)
    call exisd('MATR_ASSE', raide, iret)
    if (iret .ne. 0) then
        call dismoi('NOM_NUME_DDL', raide, 'MATR_ASSE', repk=numedd)
    else
        numedd = ''
    end if
    if (lpg) then
        call vpcrea(0, modes, raide0, amor, masse0, numedd, iret)
    else
        call vpcrea(0, modes, masse, amor, raide, numedd, iret)
    end if
!
! --  VERIFICATION DES "REFE"
    call vrrefe(masse, raide, iret)
    if (iret .gt. 0) then
        valk(1) = raide
        valk(2) = masse
        call utmess('F', 'ALGELINE2_58', nk=2, valk=valk)
    end if
    if (lc) then
        call vrrefe(raide, amor, iret)
        if (iret .gt. 0) then
            valk(1) = raide
            valk(2) = amor
            call utmess('F', 'ALGELINE2_58', nk=2, valk=valk)
        end if
    end if
!
    if (checksd) then
        call cheksd(eigsol, 'SD_EIGENSOLVER', iret)
        ASSERT(iret .eq. 0)
    end if
!
!     FIN DE VPVERS
!
end subroutine
