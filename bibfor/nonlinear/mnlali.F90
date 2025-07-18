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
subroutine mnlali(reprise, modini, imat, xcdl, parcho, &
                  adime, ninc, nd, nchoc, h, &
                  hf, ampl, xvect, lnm, num_ordr)
    implicit none
!
!
!       MODE NON LINEAIRE - ALGORITHME - INITIALISATION
!       -         -        -             --           -
! ----------------------------------------------------------------------------
!
! INITIALISATION DU POINT DE DEPART DE LA MAN
! ----------------------------------------------------------------------------
! IN   REPRISE  : L    : INDIQUE SI ON EFFECTUE UNE REPRISE OU NON
! IN   MODINI   : K8   : RESULTAT POUR LA REPRISE (SI REPRISE=TRUE)
! IN   IMAT     : I(2) : DESCRIPTEUR DES MATRICES :
!                       - IMAT(1) => MATRICE DE RAIDEUR
!                       - IMAT(2) => MATRICE DE MASSE
! IN   XCDL     : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN   PARCHO   : K14  : SD PARAMETRE DES CONTACTEURS
! IN   ADIME    : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN   NINC     : I    : NOMBRE D INCONNUES DU SYSTEME
! IN   ND       : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN   NCHOC    : I    : NOMBRE DE CONTACTEURS
! IN   H        : I    : NOMBRE D'HARMONIQUES POUR U
! IN   HF       : I    : NOMBRE D'HARMONIQUES POUR F
! IN   AMPL     : R8   : AMPLITUDE DE DEPART
! OUT  XVECT    : K14  : NOM DU VECTEUR A INITIALISER
! IN   LNM      : K8   : NOM DU CONCEPT MODE LINE POUR INITIALISER
! IN   NUM_ORDR : K8   : NUMERO D'ORDRE DU CONCEPT MODE LINE POUR INITIALISER
! ----------------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8depi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlbil.h"
#include "asterfort/mnlcir.h"
#include "asterfort/mnluil.h"
#include "asterfort/rsadpa.h"
#include "asterfort/vprecu.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "blas/idamax.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    aster_logical :: reprise
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf, num_ordr
    character(len=8) :: modini
    character(len=14) :: parcho, adime, xcdl, xvect
    real(kind=8) :: ampl
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    character(len=19) :: matrma, matrig, clnm
    character(len=8) :: lnm
    character(len=14) :: xdep1, xdep2, xtemp
    character(len=24) :: typmod
    character(len=16) :: k16bid
    real(kind=8) :: omega, lambda, ampref, alpha, eta, jeu
    integer(kind=8) :: ivect, neq, neqv, nbmode, nbpari, nbparr, nbpark, ilnm, ifreq
    integer(kind=8) :: iadim, i, j, k, icdl
    integer(kind=8) :: neqs, ijmax, nddlx, nddly
    integer(kind=8) :: nt, nddl, ht, idep1, idep2, itemp
    real(kind=8), pointer :: raid(:) => null()
    real(kind=8), pointer :: reg(:) => null()
    character(len=8), pointer :: type(:) => null()
    integer(kind=8), pointer :: vneqs(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: orig(:) => null()
    integer(kind=8), pointer :: vnddl(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
! ----------------------------------------------------------------------
! --- RECUPERATION DES NOMS ET DE LA TAILLE DES
! ---                                   MATRICES DE RIGIDITE ET DE MASSE
! ----------------------------------------------------------------------
    matrig = zk24(zi(imat(1)+1)) (1:19)
    matrma = zk24(zi(imat(2)+1)) (1:19)
    neq = zi(imat(1)+2)
! ----------------------------------------------------------------------
! --- QUELQUES VALEURS UTILES
! ----------------------------------------------------------------------
    lambda = 0.d0
    call jeveuo(adime, 'L', iadim)
    call jeveuo(xcdl, 'L', icdl)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.REG', 'L', vr=reg)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', ijmax)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.ORIG', 'L', vr=orig)
! ----------------------------------------------------------------------
! --- RECUPERATION ET MISE A ZERO DU VECTEUR D'INITIALISATION
! ----------------------------------------------------------------------
    call jeveuo(xvect, 'E', ivect)
! ----------------------------------------------------------------------
! --- INITIALISATION
! ----------------------------------------------------------------------
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(ivect), b_incx)
! ----------------------------------------------------------------------
! --- RECUPERATION DE LA FREQUENCE PROPRE DU MODE LINEAIRE
! ----------------------------------------------------------------------
    if (reprise) then
        call rsadpa(modini, 'L', 1, 'FREQ', 2, &
                    0, sjv=ifreq, styp=k16bid)
    else
        call rsadpa(lnm, 'L', 1, 'FREQ', num_ordr, &
                    0, sjv=ifreq, styp=k16bid)
    end if
    omega = r8depi()*zr(ifreq)/zr(iadim-1+3)
! ----------------------------------------------------------------------
! --- RECUPERATION DU MODE D'INITIALISATION
! ----------------------------------------------------------------------
    clnm = '&&MNLALI.RECUP     '
    if (reprise) then
        call vprecu(modini, 'DEPL', -1, [0], clnm, &
                    0, ' ', ' ', ' ', ' ', &
                    neqv, nbmode, typmod, nbpari, nbparr, &
                    nbpark)
    else
        call vprecu(lnm, 'DEPL', -1, [0], clnm, &
                    0, ' ', ' ', ' ', ' ', &
                    neqv, nbmode, typmod, nbpari, nbparr, &
                    nbpark)
    end if
    call jeveuo(clnm, 'L', ilnm)
!
! ----------------------------------------------------------------------
! --- COPIE DU MODE PROPRE DANS LE VECTEUR D'INITIALISATION
! ----------------------------------------------------------------------
    if (reprise) then
        ht = (nbmode-1)/2
        do j = 1, nbmode
            i = 0
            do k = 1, neq
                if (zi(icdl-1+k) .eq. 0) then
                    i = i+1
                    if (h .gt. ht .and. j .gt. ht+1) then
                        zr(ivect-1+(h-ht+j-1)*nd+i) = zr(ilnm-1+(j-1)*neq+k)
                    else if (h .lt. ht .and. j .gt. (h+1) .and. j .le. (2*h+1)) &
                        then
                        zr(ivect-1+(j-1)*nd+i) = zr(ilnm-1+(ht-h+j-1)*neq+k)
                    else
                        zr(ivect-1+(j-1)*nd+i) = zr(ilnm-1+(j-1)*neq+k)
                    end if
                end if
            end do
        end do
    else
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(ivect-1+nd+i) = zr(ilnm-1+k+(num_ordr-1)*neq)
            end if
        end do
    end if
! ----------------------------------------------------------------------
! --- ADIMENSIONNEMENT
! ----------------------------------------------------------------------
    if (reprise) then
        b_n = to_blas_int(nd*(2*h+1))
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/zr(ijmax), zr(ivect), b_incx)
    else
! --- MISE A L'ECHELLE PAR L'AMPLITUDE DE DEPART
!        iamax=idamax(nd,zr(ivect+nd),1)
!        ampref=zr(ivect-1+nd+iamax)
        ampref = 1.d0
        b_n = to_blas_int(nd)
        b_incx = to_blas_int(1)
        call dscal(b_n, ampl/ampref, zr(ivect+nd), b_incx)
    end if
! ----------------------------------------------------------------------
! --- REMPLISSAGE DES FORCES DE CHOCS ET DES VECTEURS AUXILIAIRES
! ----------------------------------------------------------------------
    xdep1 = '&&MNLALI.DEP1'
    xdep2 = '&&MNLALI.DEP2'
    xtemp = '&&MNLALI.TEMP'
    call wkvect(xdep1, 'V V R', 2*h+1, idep1)
    call wkvect(xdep2, 'V V R', 2*h+1, idep2)
    call wkvect(xtemp, 'V V R', ninc, itemp)
    nt = int(2**int(dlog(2.d0*dble(hf)+1.d0)/dlog(2.d0)+1.d0))
    neqs = 0
    do i = 1, nchoc
! ---   ON RECUPERE LES PARAMETRES DE CHOCS
        alpha = raid(i)/zr(iadim)
        eta = reg(i)
        jeu = vjeu(i)/zr(ijmax)
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0/jeu, zr(ivect-1+nddl), b_incx, zr(idep1), &
                       b_incy)
            call mnlbil(zr(idep1), omega, alpha, eta, h, &
                        hf, nt, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)))
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep2), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect-1+nddlx), b_incx, zr(idep1), b_incy)
            zr(idep1) = zr(idep1)-orig(3*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/jeu, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(ivect-1+nddly), b_incx, zr(idep2), b_incy)
            zr(idep2) = zr(idep2)-orig(3*(i-1)+2)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 1.d0/jeu, zr(idep2), b_incx)
            call mnlcir(xdep1, xdep2, omega, alpha, eta, &
                        h, hf, nt, xtemp)
            b_n = to_blas_int(4*(2*hf+1))
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(itemp), b_incx, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)), b_incy)
        else if (type(i) (1:4) .eq. 'PLAN') then
            nddl = vnddl(6*(i-1)+1)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idep1), b_incx)
            b_n = to_blas_int(2*h+1)
            b_incx = to_blas_int(nd)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0/jeu, zr(ivect-1+nddl), b_incx, zr(idep1), &
                       b_incy)
            call mnluil(zr(idep1), omega, alpha, eta, h, &
                        hf, nt, zr(ivect+nd*(2*h+1)+neqs*(2*hf+1)))
        end if
        neqs = neqs+vneqs(i)
    end do
!
! ----------------------------------------------------------------------
! --- REMPLISSAGE DES PARAMETRES
! ----------------------------------------------------------------------
! --- ON REMPLI LA PARTIE DU VECTEUR CORRESPONDANT A GAMMA1=LAMBDA*OMEGA
    zr(ivect+ninc-4) = lambda*omega
! --- ON REMPLI LA PARTIE DU VECTEUR CORRESPONDANT A GAMMA2=OMEGA*OMEGA
    zr(ivect+ninc-3) = omega*omega
! --- ON REMPLI LA PARTIE DU VECTEUR CORRESPONDANT A LAMBDA
    zr(ivect+ninc-2) = lambda
! --- ON REMPLI LA PARTIE DU VECTEUR CORRESPONDANT A OMEGA
    zr(ivect+ninc-1) = omega
! ----------------------------------------------------------------------
! --- DESTRUCTION DE LA SD_RESULTAT CONTENANT LE(S) MODE(S) LINEAIRE(S)
! ----------------------------------------------------------------------
!
!
!    call jedetr(clnm)
    call jedetr(xdep1)
    call jedetr(xdep2)
    call jedetr(xtemp)
    call jedema()
!
end subroutine
