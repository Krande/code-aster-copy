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
subroutine mnlqnl(imat, xcdl, parcho, adime, xvec1, &
                  xvec2, ninc, nd, nchoc, h, &
                  hf, xqnl)
    implicit none
!
!
!     MODE_NON_LINE PARTIE "QUADRATIQUE" NON-LINEAIRE
!     -    -                -            -   -
! ----------------------------------------------------------------------
!
! REGROUPE LES TERMES NON-LINEAIRES DU PROBLEME A RESOUDRE (MBH)
! ----------------------------------------------------------------------
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  IMAT   : I(2) : DESCRIPTEUR DES MATRICES :
!                       - IMAT(1) => MATRICE DE RAIDEUR
!                       - IMAT(2) => MATRICE DE MASSE
! IN  XCDL   : K14  : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  XVEC1  : K14  : NOM DU PREMIER VECTEUR SOLUTION
! IN  XVEC2  : K14  : NOM DU SECOND VECTEUR SOLUTION
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I    : NOMBRE DE DEGRES DE LIBERTE ACTIFS
! IN  NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN  H      : I    : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN  HF     : I    : NOMBRE D'HARMONIQUES POUR LA FORCE
! OUT XQNL   : I    : NOM DU VECTEUR DES TERMES NON-LINEAIRES
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mnlaft.h"
#include "asterfort/mrmult.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer :: imat(2), ninc, nd, nchoc, h, hf
    character(len=14) :: xcdl, parcho, adime, xvec1, xvec2, xqnl
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: kk, alpha, jeu
    integer :: puismax, nt, neq, ivec1, ivec2, icdl, iqnl, k, ivtp1, ivtp2
    integer :: ivtp3, ivtp4, ivtp5, nddl, j, i
    integer :: iadim, neqs, nddlx, nddly
    real(kind=8), pointer :: vtep6(:) => null()
    character(len=8), pointer :: type(:) => null()
    integer, pointer :: vneqs(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    real(kind=8), pointer :: jeumax(:) => null()
    integer, pointer :: vnddl(:) => null()
!
    call jemarq()
! ----------------------------------------------------------------------
! --- INITIALISATION DES VARIABLES POUR UTILISER MNLAFT
! ----------------------------------------------------------------------
! --- NT EST LA TAILLE DU VECTEUR AUQUEL ON APPLIQUE LA FFT
    puismax = int(dlog(4.d0*dble(hf)+1.d0)/dlog(2.d0)+1.d0)
    nt = 2**puismax
! ----------------------------------------------------------------------
! --- RECUPERATION DU NOM DE LA MATRICE ET TAILLE DE LA MATRICE
! ----------------------------------------------------------------------
    neq = zi(imat(1)+2)
! ----------------------------------------------------------------------
! --- RECUPERATION POINTEUR DE XVEC1, XVEC2, Q(XVEC1,XVEC2) ET XCDL
! ----------------------------------------------------------------------
    call jeveuo(xvec1, 'L', ivec1)
    call jeveuo(xvec2, 'L', ivec2)
    call jeveuo(xcdl, 'L', icdl)
    call jeveuo(adime, 'L', iadim)
    call jeveuo(xqnl, 'E', iqnl)
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.NDDL', 'L', vi=vnddl)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', vr=jeumax)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call dscal(ninc-1, 0.d0, zr(iqnl), 1)
! ----------------------------------------------------------------------
! --- EQUATION DE LA DYNAMIQUE
! ----------------------------------------------------------------------
! --- ON MET XSK PUIS XCK DANS LE VECTEUR QNL
    call dcopy(nd*h, zr(ivec2+(h+1)*nd), 1, zr(iqnl+nd), 1)
    call dcopy(nd*h, zr(ivec2+nd), 1, zr(iqnl+nd*(h+1)), 1)
!     QNL(COS) = -K*GAMA1*XSK
    call dscal(nd*h, -zr(ivec1-1+ninc-3), zr(iqnl+nd), 1)
    do k = 1, h
        call dscal(nd, dble(k), zr(iqnl-1+k*nd+1), 1)
    end do
! --- QNL(SIN) =  K*GAMA1*XCK
    call dscal(nd*h, zr(ivec1-1+ninc-3), zr(iqnl+nd*(h+1)), 1)
    do k = 1, h
        call dscal(nd, dble(k), zr(iqnl-1+(k+h)*nd+1), 1)
    end do
! --- CREATION DE 2 VECTEURS TEMPORAIRES
    call wkvect('&&MNLQNL.VTEP1', 'V V R', neq*(2*h), ivtp1)
    call wkvect('&&MNLQNL.VTEP2', 'V V R', neq*(2*h), ivtp2)
! --- VECTEMP1 = X_K DE MEME TAILLE QUE LE NBRE D'EQUATION
    do j = 1, 2*h
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(ivtp1-1+(j-1)*neq+k) = zr(ivec2-1+j*nd+i)
            end if
        end do
    end do
! --- VECTEMP2 = M*VECTEMP1
    call mrmult('ZERO', imat(2), zr(ivtp1), zr(ivtp2), 2*h, &
                .false._1)
! --- VECTEMP3 = VECTEMP2 (ON ELIMINE LES DDLS NON ACTIFS)
    call wkvect('&&MNLQNL.VTEP3', 'V V R', nd*(2*h), ivtp3)
    do j = 1, 2*h
        i = 0
        do k = 1, neq
            if (zi(icdl-1+k) .eq. 0) then
                i = i+1
                zr(ivtp3-1+(j-1)*nd+i) = zr(ivtp2-1+(j-1)*neq+k)/zr( &
                                         iadim-1+2)
            end if
        end do
    end do
! --- QNL = QNL - (K^2)*GAMMA2*VECTEMP3
    do k = 1, h
        kk = dble(k)**2
        call daxpy(nd, -kk*zr(ivec1-1+ninc-2), zr(ivtp3-1+(k-1)*nd+1), 1, zr(iqnl-1+k*nd+1), &
                   1)
        call daxpy(nd, -kk*zr(ivec1-1+ninc-2), zr(ivtp3-1+(h+k-1)*nd+1), 1, &
                   zr(iqnl-1+(h+k)*nd+1), 1)
    end do
! ----------------------------------------------------------------------
! --- EQUATIONS SUPPLEMENTAIRES
! ----------------------------------------------------------------------
!      CALL JEVEUO(PARCHO//'.ORIG','L',IORIG)
    call wkvect('&&MNLQNL.VTEP4', 'V V R', 2*hf+1, ivtp4)
    call wkvect('&&MNLQNL.VTEP5', 'V V R', 2*hf+1, ivtp5)
    AS_ALLOCATE(vr=vtep6, size=2*hf+1)
    neqs = 0
    do i = 1, nchoc
        alpha = raid(i)/zr(iadim-1+1)
        jeu = vjeu(i)/jeumax(1)
!        WRITE(6,*) 'JEUV',JEU
        if (type(i) (1:7) .eq. 'BI_PLAN') then
            nddl = vnddl(6*(i-1)+1)
!          WRITE(6,*) 'NDDLV',NDDL
! ---     -F*Z
            call mnlaft(zr(ivec1-1+nd*(2*h+1)+neqs*(2*hf+1)+1), &
                        zr(ivec2-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1), hf, nt, &
                        zr(iqnl+nd*(2*h+1)+neqs*(2*hf+1)))
            call dscal(2*hf+1, -1.d0, zr(iqnl+nd*(2*h+1)+neqs*(2*hf+1)), 1)
! ---     -(F/ALPHA-XG)*(F/ALPHA-XG))
!           VECTEMP4=F1/ALPHA - XG
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
            call daxpy(2*hf+1, 1.d0/alpha, zr(ivec1-1+nd*(2*h+1)+neqs*(2*hf+1)+1), 1, zr(ivtp4), &
                       1)
!           CSTE & COS
            call daxpy(h+1, -1.d0/jeu, zr(ivec1-1+nddl), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, -1.d0/jeu, zr(ivec1-1+nd*(h+1)+nddl), nd, zr(ivtp4+hf+1), &
                       1)
!           VECTEMP5=F2/ALPHA - XG
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
            call daxpy(2*hf+1, 1.d0/alpha, zr(ivec2-1+nd*(2*h+1)+neqs*(2*hf+1)+1), 1, zr(ivtp5), &
                       1)
!           CSTE & COS
            call daxpy(h+1, -1.d0/jeu, zr(ivec2-1+nddl), nd, zr(ivtp5), &
                       1)
!           SIN
            call daxpy(h, -1.d0/jeu, zr(ivec2-1+nd*(h+1)+nddl), nd, zr(ivtp5+hf+1), &
                       1)
!          WRITE(6,*) 'VECT5',ZR(IVTP5:IVTP5+2*HF)
            call mnlaft(zr(ivtp4), zr(ivtp5), hf, nt, zr(iqnl-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1))
            call dscal(2*hf+1, -1.d0, zr(iqnl-1+nd*(2*h+1)+(neqs+1)*(2*hf+1)+1), 1)
        else if (type(i) (1:6) .eq. 'CERCLE') then
            nddlx = vnddl(6*(i-1)+1)
            nddly = vnddl(6*(i-1)+2)
! ---     FX*R - FN*(UX/JEU)
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec2-1+nddlx), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec2-1+nd*(h+1)+nddlx), nd, zr(ivtp4+hf+1), &
                       1)
!           FN*(UX/JEU)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
            call mnlaft(zr(ivec1+nd*(2*h+1)+(neqs+3)*(2*hf+1)), zr(ivtp4), hf, nt, zr(ivtp5))
!           FX*R
            call mnlaft(zr(ivec1+nd*(2*h+1)+neqs*(2*hf+1)), &
                        zr(ivec2+nd*(2*h+1)+(neqs+2)*(2*hf+1)), hf, nt, &
                        zr(iqnl+nd*(2*h+1)+neqs*(2*hf+1)))
            call daxpy(2*hf+1, -1.d0, zr(ivtp5), 1, zr(iqnl+nd*(2*h+1)+neqs*(2*hf+1)), &
                       1)
! ---     FY*R - FN*(UY/JEU)
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec2-1+nddly), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec2-1+nd*(h+1)+nddly), nd, zr(ivtp4+hf+1), &
                       1)
!           FN*(UY/JEU)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
            call mnlaft(zr(ivec1+nd*(2*h+1)+(neqs+3)*(2*hf+1)), zr(ivtp4), hf, nt, zr(ivtp5))
!           FY*R
            call mnlaft(zr(ivec1+nd*(2*h+1)+(neqs+1)*(2*hf+1)), &
                        zr(ivec2+nd*(2*h+1)+(neqs+2)*(2*hf+1)), hf, nt, &
                        zr(iqnl+nd*(2*h+1)+(neqs+1)*(2*hf+1)))
            call daxpy(2*hf+1, -1.d0, zr(ivtp5), 1, zr(iqnl+nd*(2*h+1)+(neqs+1)*(2*hf+1)), &
                       1)
! ---     R*R - (UX/JEU)^2 - (UY/JEU)^2
!          - (UY/JEU)^2
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec1-1+nddly), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec1-1+nd*(h+1)+nddly), nd, zr(ivtp4+hf+1), &
                       1)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec2-1+nddly), nd, zr(ivtp5), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec2-1+nd*(h+1)+nddly), nd, zr(ivtp5+hf+1), &
                       1)
            call dscal(2*hf+1, 0.d0, vtep6, 1)
            call mnlaft(zr(ivtp4), zr(ivtp5), hf, nt, vtep6)
            call daxpy(2*hf+1, -1.d0, vtep6, 1, zr(iqnl+nd*(2*h+1)+(neqs+2)*(2*hf+1)), &
                       1)
!         - (UX/JEU)^2
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec1-1+nddlx), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec1-1+nd*(h+1)+nddlx), nd, zr(ivtp4+hf+1), &
                       1)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
!           CSTE & COS
            call daxpy(h+1, 1.d0/jeu, zr(ivec2-1+nddlx), nd, zr(ivtp5), &
                       1)
!           SIN
            call daxpy(h, 1.d0/jeu, zr(ivec2-1+nd*(h+1)+nddlx), nd, zr(ivtp5+hf+1), &
                       1)
            call dscal(2*hf+1, 0.d0, vtep6, 1)
            call mnlaft(zr(ivtp4), zr(ivtp5), hf, nt, vtep6)
            call daxpy(2*hf+1, -1.d0, vtep6, 1, zr(iqnl+nd*(2*h+1)+(neqs+2)*(2*hf+1)), &
                       1)
!          + R^2
            call dscal(2*hf+1, 0.d0, vtep6, 1)
            call mnlaft(zr(ivec1+nd*(2*h+1)+(neqs+2)*(2*hf+1)), &
                        zr(ivec2+nd*(2*h+1)+(neqs+2)*(2*hf+1)), hf, nt, vtep6)
            call daxpy(2*hf+1, 1.d0, vtep6, 1, zr(iqnl+nd*(2*h+1)+(neqs+2)*(2*hf+1)), &
                       1)
! ---     (FN/ALPHA - R)*FN
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
            call daxpy(2*hf+1, -1.d0, zr(ivec1+nd*(2*h+1)+(neqs+2)*(2*hf+1)), 1, zr(ivtp4), &
                       1)
            call daxpy(2*hf+1, 1.d0/alpha, zr(ivec1+nd*(2*h+1)+(neqs+3)*(2*hf+1)), 1, zr(ivtp4), &
                       1)
            call dcopy(2*hf+1, zr(ivec2+nd*(2*h+1)+(neqs+3)*(2*hf+1)), 1, zr(ivtp5), 1)
            call mnlaft(zr(ivtp4), zr(ivtp5), hf, nt, zr(iqnl+nd*(2*h+1)+(neqs+3)*(2*hf+1)))
        else if (type(i) (1:4) .eq. 'PLAN') then
            nddl = vnddl(6*(i-1)+1)
! ---     (F/ALPHA - XG)*F
            call dscal(2*hf+1, 0.d0, zr(ivtp4), 1)
            call dscal(2*hf+1, 0.d0, zr(ivtp5), 1)
!            call jxveri(' ', ' ')
!           (F/ALPHA - XG)
            call daxpy(2*hf+1, 1.d0/alpha, zr(ivec1+nd*(2*h+1)+neqs*(2*hf+1)), 1, zr(ivtp4), &
                       1)
!           CSTE & COS
            call daxpy(h+1, -1.d0/jeu, zr(ivec1-1+nddl), nd, zr(ivtp4), &
                       1)
!           SIN
            call daxpy(h, -1.d0/jeu, zr(ivec1-1+nd*(h+1)+nddl), nd, zr(ivtp4+hf+1), &
                       1)
!           F
            call daxpy(2*hf+1, 1.d0, zr(ivec2+nd*(2*h+1)+neqs*(2*hf+1)), 1, zr(ivtp5), &
                       1)
!
            call mnlaft(zr(ivtp4), zr(ivtp5), hf, nt, zr(iqnl+nd*(2*h+1)+neqs*(2*hf+1)))
        end if
        neqs = neqs+vneqs(i)
!        WRITE(6,*) 'NEQS',NEQS
    end do
! ----------------------------------------------------------------------
! --- AUTRES EQUATIONS
! ----------------------------------------------------------------------
! --- QNL(NINC-3) = -LAMBDA*OMEGA
    zr(iqnl+ninc-4) = -1.d0*zr(ivec1+ninc-2)*zr(ivec2+ninc-1)
! --- QNL(NINC-2) =-OMEGA*OMEGA
    zr(iqnl+ninc-3) = -1.d0*zr(ivec1+ninc-1)*zr(ivec2+ninc-1)
! --- EQUATION DE PHASE
    zr(iqnl+ninc-2) = 0.d0
    do k = 1, h
        zr(iqnl+ninc-2) = zr(iqnl+ninc-2)+k*zr(ivec1-1+ninc)*zr(ivec2-1+(h+k)*nd+1)
    end do
! ----------------------------------------------------------------------
! --- DESTRUCTION DES VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
    call jedetr('&&MNLQNL.VTEP1')
    call jedetr('&&MNLQNL.VTEP2')
    call jedetr('&&MNLQNL.VTEP3')
    call jedetr('&&MNLQNL.VTEP4')
    call jedetr('&&MNLQNL.VTEP5')
    AS_DEALLOCATE(vr=vtep6)
!
    call jedema()
!
end subroutine
