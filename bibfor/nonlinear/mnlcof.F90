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

subroutine mnlcof(imat, numdrv, matdrv, xcdl, parcho, &
                  adime, xvecu0, xtang, ninc, nd, &
                  nchoc, h, hf, ordman, xups, &
                  xfpnla, lbif, nextr, epsbif)
    implicit none
!
!
!     MODE_NON_LINE CALCUL DES COEFFICIENTS DE LA SERIE ENTIERE
!     -    -   -               -- -
! ----------------------------------------------------------------------
! IN      IMAT   : I(2) : DESCRIPTEUR DES MATRICES :
!                           - IMAT(1) => MATRICE DE RAIDEUR
!                           - IMAT(2) => MATRICE DE MASSE
! IN      NUMDRV : K14  : NUME_DDL_GENE DE LA MATRICE JACOBIENNE
! IN      MATDRV : K19  : NOM DE  LA MATRICE JACOBIENNE
! IN      XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN      PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN      ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN      XVECU0 : K14  : NOM DU VECTEUR INITIAL (U0)
! IN      XTANG  : K14  : NOM DU VECTEUR TANGENT (U1)
! IN      NINC   : I    : NOMBRE D'INCONNUES DU SYSTEME
! IN      ND     : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN      NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN      H      : I    : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN      HF     : I    : NOMBRE D'HARMONIQUES POUR LA FORCE
! IN      ORDMAN : I    : ORDRE DE LA MAN
! OUT     XUPS   : K14  : NOM DU VECTEUR QUI CONTIENT LES COEFFICIENTS
!                         DE LA SERIE ENTIERE DE LA VARIABLE
! OUT     XFPNLA : K14  : NOM DU VECTEUR QUI CONTIENT LE SECOND MEMBRE
!                         SUPPLEMENTAIRE
! OUT     LBIF   : L    : SI TRUE ALORS IL EXISTE UNE BIFURCATION
! IN      NEXTR  : I    : NOMBRE DE TERMES A PRENDRE EN COMPTE POUR LA GESTION DES BIF.
! IN      EPSBIF : R8   : ERREUR DE PRECISION POUR LA GESTION DES BIF.
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnldrv.h"
#include "asterfort/mnlqnl.h"
#include "asterfort/resoud.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer :: imat(2), ninc, nd, nchoc, h, hf, ordman, nextr
    real(kind=8) :: epsbif
    character(len=14) :: numdrv, xcdl, parcho, adime, xvecu0, xtang, xups, xfpnla
    character(len=19) :: matdrv
    aster_logical :: lbif
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    integer :: ivecu0, iups, itang, p, r, ivecu1, ivecu2, iqnl, iret
    integer :: ifpnla
    character(len=14) :: xvecu1, xvecu2, xqnl
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES POUR EXTRACTION GEOMETRIQUE
! ----------------------------------------------------------------------
    integer :: k
!
    character(len=19) :: solveu
    real(kind=8) :: nvec, nratio, necar, ac, nudom
    complex(kind=8) :: cbid
    real(kind=8), pointer :: alpha(:) => null()
    real(kind=8), pointer :: ecar(:) => null()
    real(kind=8), pointer :: fpnl(:) => null()
    real(kind=8), pointer :: ratio(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
    call jemarq()
    solveu = '&&OP0061.SOLVEUR'
! ----------------------------------------------------------------------
! --- ON RECUPERE LE VECTEUR INITIAL, TANGENT ET LA SERIE
! ----------------------------------------------------------------------
    call jeveuo(xvecu0, 'L', ivecu0)
    call jeveuo(xups, 'E', iups)
    call jeveuo(xtang, 'L', itang)
    lbif = .false.
! ----------------------------------------------------------------------
! --- ON INSERE LE VECTEUR INITIAL
! ----------------------------------------------------------------------
    call dcopy(ninc, zr(ivecu0), 1, zr(iups), 1)
! ----------------------------------------------------------------------
! --- ON INSERE LE VECTEUR CORRESPONDANT AU PREMIER ORDRE DE LA SERIE
! ----------------------------------------------------------------------
    call dcopy(ninc, zr(itang), 1, zr(iups+ninc), 1)
! ----------------------------------------------------------------------
! --- ON CALCUL LA MATRICE JACOBIENNE
! ----------------------------------------------------------------------
    call mnldrv(.false._1, imat, numdrv, matdrv, xcdl, &
                parcho, adime, xvecu0, zr(itang), ninc, &
                nd, nchoc, h, hf)
! ----------------------------------------------------------------------
! --- CALCUL DES ORDRES P=2,...,ORDMAN
! ----------------------------------------------------------------------
    AS_ALLOCATE(vr=fpnl, size=ninc)
    xvecu1 = '&&MNLCOF.VECU1'
    call wkvect(xvecu1, 'V V R', ninc, ivecu1)
    xvecu2 = '&&MNLCOF.VECU2'
    call wkvect(xvecu2, 'V V R', ninc, ivecu2)
    xqnl = '&&MNLCOF.VECQ'
    call wkvect(xqnl, 'V V R', ninc-1, iqnl)
    do p = 2, ordman
!       REMISE A ZERO DU SECOND MEMBRE
        call dscal(ninc, 0.d0, fpnl, 1)
!       CALCUL DU SECOND MEMBRE
        do r = 1, p-1
!         VECU1 = UPS(:,R)
            call dcopy(ninc, zr(iups+r*ninc), 1, zr(ivecu1), 1)
!         VECUI = UPS(:,P-R)
            call dcopy(ninc, zr(iups+(p-r)*ninc), 1, zr(ivecu2), 1)
!         CALCULE DE Q(UPS(:,R),UPS(:,P-R))
            call mnlqnl(imat, xcdl, parcho, adime, xvecu1, &
                        xvecu2, ninc, nd, nchoc, h, &
                        hf, xqnl)
!         CALCUL DE FPNL(1:NEQ)=FPNL(1:NEQ)-Q(SYS,UPS(:,R),UPS(:,P-R))
            call daxpy(ninc-1, -1.d0, zr(iqnl), 1, fpnl, &
                       1)
        end do
        fpnl(ninc) = 0.d0
! ---   RESOLUTION DU SYSTEME LINEAIRE UPS(:,P) = K\FPNL
        call resoud(matdrv, ' ', solveu, ' ', 1, &
                    ' ', ' ', 'v', fpnl, [cbid], &
                    ' ', .false._1, 0, iret)
        call dcopy(ninc, fpnl, 1, zr(iups+p*ninc), 1)
    end do
! ----------------------------------------------------------------------
! --- REMISE A ZERO DES VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call dscal(ninc-1, 0.d0, zr(iqnl), 1)
    call dscal(ninc, 0.d0, zr(ivecu1), 1)
    call dscal(ninc, 0.d0, zr(ivecu2), 1)
! ----------------------------------------------------------------------
! --- EXTRACTION SERIE GEOMETRIQUE
! ----------------------------------------------------------------------
    AS_ALLOCATE(vr=alpha, size=nextr)
    AS_ALLOCATE(vr=ratio, size=nextr)
    AS_ALLOCATE(vr=ecar, size=nextr-1)
    do k = 1, nextr
        call dscal(ninc, 0.d0, zr(ivecu1), 1)
        alpha(k) = ddot(ninc, zr(iups+(ordman-k+2-1)*ninc), 1, zr(iups+(ordman-k+1-1)*ninc), 1)
        alpha(k) = alpha(k)/ddot(ninc, zr(iups+(ordman-k+2-1)*ninc), &
                                 1, zr(iups+(ordman-k+2-1)*ninc), 1)
        call dcopy(ninc, zr(iups+(ordman-k+1-1)*ninc), 1, zr(ivecu1), 1)
        call daxpy(ninc, -alpha(k), zr(iups+(ordman-k+2-1)*ninc), 1, zr(ivecu1), &
                   1)
        nvec = dnrm2(ninc, zr(ivecu1), 1)
        ratio(k) = nvec/dnrm2(ninc, zr(iups+(ordman-k+1-1)*ninc), 1)
        if (k .gt. 1) then
            ecar(k-1) = (alpha(k-1)-alpha(k))/alpha(k-1)
        end if
    end do
    nratio = dnrm2(nextr, ratio, 1)
    necar = dnrm2(nextr-1, ecar, 1)
    if (nratio .lt. epsbif .and. necar .lt. epsbif) then
        lbif = .true.
        call dscal(ninc, 0.d0, zr(ivecu1), 1)
        ac = ddot(ninc, zr(iups+ordman*ninc), 1, zr(iups+(ordman-1)*ninc), 1)
        ac = ac/ddot(ninc, zr(iups+ordman*ninc), 1, zr(iups+ordman*ninc), 1)
        call dcopy(ninc, zr(iups+ordman*ninc), 1, zr(ivecu1), 1)
        call dscal(ninc, ac**ordman, zr(ivecu1), 1)
        nudom = dnrm2(ninc, zr(ivecu1), 1)
        do k = 1, ordman
            call daxpy(ninc, -((1.d0/ac)**k), zr(ivecu1), 1, zr(iups+k*ninc), &
                       1)
        end do
    end if
!
    call dscal(ninc, 0.d0, zr(ivecu1), 1)
! ----------------------------------------------------------------------
! --- CALCUL DU SECOND MEMBRE SUPPLEMENTAIRE POUR LE CALCUL DE LA
! ---                                                     LONGUEUR D'ARC
! ----------------------------------------------------------------------
    call jeveuo(xfpnla, 'E', ifpnla)
    call dscal(ninc-1, 0.d0, zr(ifpnla), 1)
! ---   CALCULE DE FPNLA = FPNLA - Q(SYS,UPS(:,R),UPS(:,P+1-R))
    do r = 1, ordman
        call dscal(ninc-1, 0.d0, zr(iqnl), 1)
        call dscal(ninc, 0.d0, zr(ivecu1), 1)
        call dscal(ninc, 0.d0, zr(ivecu2), 1)
! ---   VECU1 = UPS(:,R+1)
        call dcopy(ninc, zr(iups+r*ninc), 1, zr(ivecu1), 1)
! ---   VECU2 = UPS(:,ORDMAN+2-R)
        call dcopy(ninc, zr(iups+(ordman+1-r)*ninc), 1, zr(ivecu2), 1)
! ---   Q(VECU1,VECU2)
        call mnlqnl(imat, xcdl, parcho, adime, xvecu1, &
                    xvecu2, ninc, nd, nchoc, h, &
                    hf, xqnl)
!       AJOUT DES DEUX VECTEURS DANS XFPNLA
        call daxpy(ninc-1, -1.d0, zr(iqnl), 1, zr(ifpnla), &
                   1)
    end do
! ----------------------------------------------------------------------
! --- DESTRUCTION DES VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call jedetr(xvecu1)
    call jedetr(xvecu2)
    AS_DEALLOCATE(vr=fpnl)
    call jedetr(xqnl)
    AS_DEALLOCATE(vr=alpha)
    AS_DEALLOCATE(vr=ratio)
    AS_DEALLOCATE(vr=ecar)
!
    call jedema()
!
end subroutine
