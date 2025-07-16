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
subroutine pmvtgt(option, carcri, deps2, sigp, vip, &
                  nbvari, epsilo, varia, matper, dsidep, &
                  smatr, sdeps, ssigp, svip, iret)
! person_in_charge: jean-michel.proix at edf.fr
!-----------------------------------------------------------------------
!     OPERATEUR CALC_POINT_MAT : MATRICE TANGENTE PAR PERTURBATION
!     RESSEMBLE A TGVERI MAIS SANS ELEMENT FINI
!-----------------------------------------------------------------------
! ----------------------------------------------------------------------
! VAR OPTION NOM DE L'OPTION DE CALCUL
!             IN  : CELLE UTILISEE PAR CALC_POINT_MAT
!             OUT : 'RAPH_MECA' SI BOUCLE, 'FULL_MECA' SI FIN DE BOUCLE
! IN  CARCRI  : CARCRI(1) = type de matrice tangente
!               0 : ANALYTIQUE, on ne passe pas ici
!               1 : PERTURBATION, on calcule Ktgte (FULL_MECA)
!               2 : VERIFICATION, on calcule Ktgte (FULL_MECA) + Kpertu
!               CARCRI(7) = valeur de la perturbation
! IN  DEPS2   : DEFORMATIONS
! IN  SIGP    : CONTRAINTES
! IN  VIP     : VARIABLES INTERNES
! IN  NBVARI  : Nombre de variables internes
! VAR EPSILO  : VALEUR DE LA PERTURBATION, A GARDER
! VAR VARIA   : TABLEAU DES VARIATIONS
! VAR MATPER  : MATRICE TANGENTE PAR PERTURBATION
! VAR SMATR   : SAUVEGARDE MATRICE TANGENTE
! VAR SDEPS   : SAUVEGARDE DEFORMATIONS
! VAR SSIGP   : SAUVEGARDE CONTRAINTES
! VAR SVIP    : VARIABLES INTERNES
! OUT IRET  SI IRET = 0 -> FIN, SINON -> BOUCLE
! ----------------------------------------------------------------------
    implicit none
#include "jeveux.h"
#include "asterc/r8miem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/r8inir.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    character(len=16) :: option
    integer(kind=8) :: iret, nbvari
    real(kind=8) :: carcri(*), deps2(6), sigp(6), matper(36), dsidep(6, 6)
    real(kind=8) :: sdeps(6), ssigp(6), vip(nbvari), svip(nbvari), smatr(36)
    real(kind=8) :: varia(2*36)
!
!
    character(len=24) :: matra, matrc
    integer(kind=8) :: ematra, ematrc, exi, i, j, indi, nvar, init, pos
    real(kind=8) :: v, epsilo, fp, fm, pertu, maxeps
    blas_int :: b_incx, b_incy, b_n
    save init, pos
    data matra/'PYTHON.TANGENT.MATA'/
    data matrc/'PYTHON.TANGENT.MATC'/
    data init, pos/1, 0/
! ----------------------------------------------------------------------
!
    call jemarq()
!
!     Calcul de la matrice TGTE par PERTURBATION
!
    iret = 0
    if (abs(carcri(2)) .lt. 0.1d0) then
        goto 999
    end if
    if (option(1:9) .eq. 'RIGI_MECA') then
        goto 999
    end if
!
! --  INITIALISATION (PREMIER APPEL)
!
    if (init .eq. 1) then
!
!       PERTURBATION OU VERIFICATION => FULL_MECA
        if (option .ne. 'FULL_MECA') then
            goto 999
        end if
!
!       CALCUL de la valeur de la perturbation
        maxeps = 0.d0
        do i = 1, 6
            maxeps = max(maxeps, abs(deps2(i)))
        end do
        pertu = carcri(7)
        epsilo = pertu*maxeps
        if (epsilo .lt. r8miem()) then
            call utmess('A', 'ALGORITH11_86')
            goto 999
        end if
!
!      ARCHIVAGE DES VALEURS DE REFERENCE
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deps2, b_incx, sdeps, b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp, b_incx, ssigp, b_incy)
        b_n = to_blas_int(nbvari)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vip, b_incx, svip, b_incy)
!       ARCHIVAGE DE LA MATRICE TANGENTE COHERENTE
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dsidep, b_incx, smatr, b_incy)
!      PREPARATION DES ITERATIONS
        option = 'RAPH_MECA'
        iret = 1
        init = 0
        pos = 0
!
    end if
!
! -- TRAITEMENT DES VARIATIONS
!
!
!    SAUVEGARDE DE LA FORCE INTERIEURE PERTURBEE
!
    nvar = int((pos+1)/2)
!
    if (nvar .gt. 0) then
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigp, b_incx, varia(1+(pos-1)*6), b_incy)
    end if
!
    pos = pos+1
    nvar = int((pos+1)/2)
    indi = 1-2*mod(pos, 2)
!
    if (nvar .le. 6) then
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sdeps, b_incx, deps2, b_incy)
        deps2(nvar) = sdeps(nvar)+indi*epsilo
!
!      INITIALISATION DES CHAMPS 'E'
        call r8inir(6, 0.d0, sigp, 1)
        iret = 1
        goto 999
    end if
!
!    CALCUL DE LA MATRICE TANGENTE
!
    do i = 1, 6
        do j = 1, 6
            fm = varia((2*j-2)*6+i)
            fp = varia((2*j-1)*6+i)
            v = (fp-fm)/(2*epsilo)
            matper((i-1)*6+j) = v
        end do
    end do
!
!    MENAGE POUR ARRET DE LA ROUTINE
!
    iret = 0
    init = 1
    option = 'FULL_MECA'
!
!    RETABLISSEMENT DE LA SOLUTION
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sdeps, b_incx, deps2, b_incy)
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, ssigp, b_incx, sigp, b_incy)
    b_n = to_blas_int(nbvari)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, svip, b_incx, vip, b_incy)
!
!     PERTURBATION => SAUVEGARDE DE LA MATRICE CALCULEE PAR
!     DIFFERENCES FINIES COMME MATRICE TANGENTE
!
    if (abs(carcri(2)-1.d0) .lt. 0.1d0) then
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, matper, b_incx, dsidep, b_incy)
!
!     VERIFICATION
!
    else if (abs(carcri(2)-2.d0) .lt. 0.1d0) then
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, smatr, b_incx, dsidep, b_incy)
!
!      CREATION DES OBJETS
!      CE N'EST PAS LA PREMIERE FOIS QU'ON CALCULE LA MATRICE TANGENTE
!      -> ON NE CONSERVE QUE LE DERNIER CALCUL (EN COURS)
        call jeexin(matra, exi)
        if (exi .ne. 0) then
            call jedetr(matra)
            call jedetr(matrc)
        end if
        call wkvect(matra, 'G V R', 36, ematra)
        call wkvect(matrc, 'G V R', 36, ematrc)
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, smatr, b_incx, zr(ematra), b_incy)
        b_n = to_blas_int(36)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, matper, b_incx, zr(ematrc), b_incy)
!         CALL JELIBE(MATRA)
!         CALL JELIBE(MATRC)
    end if
!
999 continue
!
    call jedema()
end subroutine
