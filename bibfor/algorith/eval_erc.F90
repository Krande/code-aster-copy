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
subroutine eval_erc(baseno, dynam1, vecterc, nommes, matobs, &
                    obsdim, ifreq, omega, alpha, cout_fon, &
                    terme_uv)
!
!
    implicit none
! ----------------------------------------------------------------------
!
!  ROUTINE LIEE A L'OPERATEUR CALC_ERC_DYN
!
!  ROUTINE DE REMPLISSAGE DES VALEURS DU SECOND MEMBRE ERC EN MODALE
!  ON CONSTRUIT LE PRODUIT coef_alpha*H^T*G*tilde(u)
! ----------------------------------------------------------------------
! IN  : BASENO        : NOM COMMUN POUR LES OBJETS JEVEUX A CREER
! IN  : DYNAM1        : NOM DE L'OBJET JEVEUX CONTENANT LA PREMIERE MATRICE
!                       DYNAMIQUE DU BLOC DIAGONAL DE LA MATRICE D'ERC
! IN  : VECT_ERC      : VECTEUR DU 2ND MEMBRE ASSOCIE AU PB MATRICIEL D'ERC
! IN  : NOMMES        : NOM DU CONCEPT JEVEUX CONTENANT LA MESURE
! IN  : MATOBS        : LISTE DES NOMS DES OBJETS JEVEUX DEFINISSANT LA MATRICE
!                       D'OBSERVATION. LA MATRICE EST STOCKEE EN SPARSE SOUS LE
!                       FORMAT COO (FILE,COLONNE,VALEUR)
! IN  : OBSDIM        : TABLEAU DONNANT LES INFORMATIONS DIMENSIONNELLES DE LA
!                       MATRICE D'OBSERVATION (DIM_FILE,DIM_COLONNE,NOMBRE_DE_
!                       VALEURS_NONNULLES)
! IN  : I_FREQ        : NUMERO D'ORDRE DE LA MESURE ASSOCIE A LA FREQ EN COURS
! IN  : OMEGA         : PULSATION ASSOCIEE A LA FREQUENCE EN COURS
! IN  : ALPHA         : PARAMETRE ALPHA DE LA FONCTIONNELLE D'ERC
! OUT : COUT_FON      : VALEUR DE LA FONCTION COUT
! OUT : TERME_UV      : PART DE LA VALEUR DE LA FONCTION COUT ASSOCIE AUX CHAMPS D'ERR
! ----------------------------------------------------------------------!
! ----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jedetr.h"
#include "asterfort/jexnum.h"
#include "asterfort/getvtx.h"
#include "asterfort/getvid.h"
#include "asterfort/mrmult.h"
#include "asterfort/r8inir.h"
#include "asterfort/wkvect.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "blas/dspmv.h"
#include "blas/ddot.h"
!
    character(len=8), intent(in) :: baseno, nommes
    character(len=19), intent(in) :: dynam1
    integer(kind=8), intent(in) :: ifreq, obsdim(3)
    real(kind=8), intent(in) :: vecterc(*), omega, alpha
    real(kind=8), intent(out) :: cout_fon, terme_uv
    character(len=24), intent(in) :: matobs(3)
    character(len=4) :: type_mes
    character(len=8) :: mnorme
    character(len=11) :: bl11
    character(len=19) :: nom_objev_mes
    integer(kind=8) :: occ, i_tach, i_mes, iaux1, iaux2, iaux3, iaux4, iobsfil, iobscol
    integer(kind=8) :: iobsval, ii, n_fil, n_col, idesc, nvect_mes, nvale_norme, ivale_norm, lmat
    real(kind=8) :: coef_mes, coef_alpha, terme_obs
    aster_logical :: isdiag
    blas_int :: b_incx, b_incy, b_n
!
    bl11 = '           '
!
! --- CALCUL DU PRODUIT (H*u-\tilde{u})
! --- --- recuperation de la mesure
    call getvtx(' ', 'CHAMP_MESURE', scal=type_mes)
    if (type_mes .eq. 'DEPL') then
        occ = 1
        coef_mes = 1.0d0
    else if (type_mes .eq. 'VITE') then
        occ = 2
        coef_mes = 1.0d0/omega
    else
        occ = 3
        coef_mes = -1.0d0/(omega*omega)
    end if
!
    call jeveuo(jexnum(nommes//'           .TACH', occ), 'L', i_tach)
    nom_objev_mes = zk24(i_tach+ifreq-1) (1:19)
    call jeveuo(nom_objev_mes//'.VALE', 'L', i_mes)
!  on recopie la mesure dans un vecteur de travail avec signe oppose
    call wkvect(baseno//'.EVALF_AUX1.VAL', 'V V R', obsdim(1), iaux1)
    b_n = to_blas_int(obsdim(1))
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(i_mes), b_incx, zr(iaux1), b_incy)
    b_n = to_blas_int(obsdim(1))
    b_incx = to_blas_int(1)
    call dscal(b_n, -1.d0*coef_mes, zr(iaux1), b_incx)
! --- recuperation de la matrice d'observation
    call jeveuo(matobs(1), 'L', iobsfil)
    call jeveuo(matobs(2), 'L', iobscol)
    call jeveuo(matobs(3), 'L', iobsval)
!   calcul du produit sachant que iaux a ete initialise avec -\tilde{u}
    do ii = 1, obsdim(3)
        n_fil = zi(iobsfil-1+ii)
        n_col = zi(iobscol-1+ii)
        zr(iaux1-1+n_fil) = zr(iaux1-1+n_fil)+zr(iobsval-1+ii)*vecterc(obsdim(2)+n_col)
    end do
! --- FIN CALCUL DU PRODUIT (H*u-\tilde{u}) qui se trouve dans iaux1
!
    isdiag = .true._1
!
! --- RECUPERATION DE LA MATRICE NORME
    call getvid(' ', 'MATR_NORME', scal=mnorme)
    call jeveuo(mnorme//bl11//'.DESC', 'L', idesc)
    nvect_mes = zi(idesc+1)
    nvale_norme = nvect_mes
    if (zi(idesc+2) .eq. 2) isdiag = .false._1
    if (.not. isdiag) nvale_norme = (nvect_mes*(nvect_mes+1))/2
!
    if (nvect_mes .ne. obsdim(1)) then
        call utmess('F', 'ALGORITH9_71')
    end if
!
    call jeveuo(jexnum(mnorme//bl11//'.VALM', 1), 'L', ivale_norm)
!
! --- PRODUIT  (alpha/(1-alpha))*(H*u-\tilde{u})^T*G*(H*u-\tilde{u})
!
! --- on cree un deuxieme vecteur de travail  iaux2
    call wkvect(baseno//'.EVALF_AUX2.VAL', 'V V R', nvect_mes, iaux2)
    call r8inir(nvect_mes, 0.d0, zr(iaux2), 1)
!
! --- premier produit matrice vecteur   coef_alpha*G*(H*u-\tilde{u})
    coef_alpha = alpha/(1.0d0-alpha)
    if (isdiag) then
!
        do ii = 1, nvect_mes
            zr(iaux2-1+ii) = zr(ivale_norm-1+ii)*zr(iaux1-1+ii)*coef_alpha
        end do
    else
        b_n = to_blas_int(nvect_mes)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dspmv('u', b_n, coef_alpha, zr(ivale_norm), zr(iaux1), &
                   b_incx, 0.d0, zr(iaux2), b_incy)
    end if
!   on finalise par le prduit de iaux1 et iaux2 pour avoir le produit final du terme d'obs
    b_n = to_blas_int(nvect_mes)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    terme_obs = ddot(b_n, zr(iaux1), b_incx, zr(iaux2), b_incy)
! --- PRODUIT DU BLOC DYNAM1 AVEC LE CHAMP (u-v)
    call wkvect(baseno//'.EVALF_AUX3.VAL', 'V V R', obsdim(2), iaux3)
    b_n = to_blas_int(obsdim(2))
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vecterc(1), b_incx, zr(iaux3), b_incy)
    call wkvect(baseno//'.EVALF_AUX4.VAL', 'V V R', obsdim(2), iaux4)
    call r8inir(obsdim(2), 0.d0, zr(iaux4), 1)
    call jeveuo(dynam1//'.&INT', 'L', lmat)
    call mrmult('ZERO', lmat, zr(iaux3), zr(iaux4), 1, &
                .false._1)
    b_n = to_blas_int(obsdim(2))
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    terme_uv = 0.5d0*ddot(b_n, zr(iaux3), b_incx, zr(iaux4), b_incy)
    cout_fon = terme_uv+terme_obs
!     NETOYAGE DES OBJETS JEVEUX TEMPORAIRES
    call jedetr(baseno//'.EVALF_AUX1.VAL')
    call jedetr(baseno//'.EVALF_AUX2.VAL')
    call jedetr(baseno//'.EVALF_AUX3.VAL')
    call jedetr(baseno//'.EVALF_AUX4.VAL')
!
end subroutine
