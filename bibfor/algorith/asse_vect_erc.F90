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
subroutine asse_vect_erc(baseno, nom_vect_erc, nommes, matobs, obsdim, &
                         alpha, n_ordre_mes, omega)
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
! IN  : NOM_VECT_ERC  : NOM DE L'OBJET JEVEUX DU 2ND MEMBRE ASSOCIE
!                       AU PB MATRICIEL D'ERC
! IN  : NOM_NUME_ERC  : NOM DE L'OBJET JEVEUX DU NUME_DDL CREE ASSOCIE
!                       AU PB MATRICIEL D'ERC
! IN  : NOMMES        : NOM DU CONCEPT JEVEUX CONTENANT LA MESURE
! IN  : MATOBS        : LISTE DES NOMS DES OBJETS JEVEUX DEFINISSANT LA MATRICE
!                       D'OBSERVATION. LA MATRICE EST STOCKEE EN SPARSE SOUS LE
!                       FORMAT COO (FILE,COLONNE,VALEUR)
! IN  : OBSDIM        : TABLEAU DONNANT LES INFORMATIONS DIMENSIONNELLES DE LA
!                       MATRICE D'OBSERVATION (DIM_FILE,DIM_COLONNE,NOMBRE_DE_
!                       VALEURS_NONNULLES)
! IN  : ALPHA         : PARAMETRE ALPHA DE LA FONCTIONNELLE D'ERC
! IN  : N_ORDRE_MES   : NUMERO D'ORDRE DE LA MESURE ASSOCIE A LA FREQ EN COURS
! IN  : OMEGA         : PULSATION ASSOCIEE A LA FREQUENCE EN COURS
! ----------------------------------------------------------------------!
! ----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/getvtx.h"
#include "asterfort/getvid.h"
#include "asterfort/r8inir.h"
#include "asterfort/wkvect.h"
#include "asterfort/utmess.h"
#include "asterfort/jedetr.h"
#include "blas/dcopy.h"
#include "blas/dspmv.h"
!
    character(len=4) :: type_mes
    character(len=8), intent(in) :: nommes, baseno
    character(len=8) :: mnorme
    character(len=11) :: bl11
    character(len=19), intent(in) :: nom_vect_erc
    character(len=19) :: nom_objev_mes
    integer(kind=8), intent(in) :: obsdim(3), n_ordre_mes
    integer(kind=8) :: occ, i_tach, i_mes, ii, idesc, nvect_mes, nvale_norme, ivale_norm
    integer(kind=8) :: iaux1, iaux2, iobsfil, iobscol, iobsval, ivecterc, n_fil, n_col
    real(kind=8), intent(in) :: alpha, omega
    real(kind=8) :: coeff_alpha, coef_mes
    character(len=24), intent(in) :: matobs(3)
    logical :: isdiag
    blas_int :: b_incx, b_incy, b_n
!
    bl11 = '           '
!
! --- RECUPERATION DE LA MESURE
    call getvtx(' ', 'CHAMP_MESURE', scal=type_mes)
    if (type_mes .eq. 'DEPL') then
        occ = 1
        coef_mes = 1.0d0
    else if (type_mes .eq. 'VITE') then
        occ = 2
        coef_mes = 1.0d0/omega
    else
        occ = 3
        coef_mes = 1.0d0/(-omega*omega)
    end if
!
    call jeveuo(jexnum(nommes//'           .TACH', occ), 'L', i_tach)
    nom_objev_mes = zk24(i_tach+n_ordre_mes-1) (1:19)
    call jeveuo(nom_objev_mes//'.VALE', 'L', i_mes)
!
!
    isdiag = .true.
! --- RECUPERATION DE LA MATRICE NORME
    call getvid(' ', 'MATR_NORME', scal=mnorme)
    call jeveuo(mnorme//bl11//'.DESC', 'L', idesc)
    nvect_mes = zi(idesc+1)
    nvale_norme = nvect_mes
    if (zi(idesc+2) .eq. 2) isdiag = .false.
    if (.not. isdiag) nvale_norme = (nvect_mes*(nvect_mes+1))/2
!
    if (nvect_mes .ne. obsdim(1)) then
        call utmess('F', 'ALGORITH9_71')
    end if
!
    call jeveuo(jexnum(mnorme//bl11//'.VALM', 1), 'L', ivale_norm)
!
! --- PRODUIT   coef_alpha*G*tilde(u)
    coeff_alpha = -2.0d0*alpha/(1.0d0-alpha)*coef_mes
! --- ON RECOPIE LA MESURE DANS UN VECTEUR DE TRAVAIL iaux1
    call wkvect(baseno//'.VECAUX1ERC.VAL', 'V V R', nvect_mes, iaux1)
    b_n = to_blas_int(nvect_mes)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(i_mes), b_incx, zr(iaux1), b_incy)
! --- ON CREE UN (PETIT) VECTEUR DE TRAVAIL  iaux2
    call wkvect(baseno//'.VECAUX2ERC.VAL', 'V V R', nvect_mes, iaux2)
    call r8inir(nvect_mes, 0.d0, zr(iaux2), 1)
!
! --- PREMIER PRODUIT MATRICE VECTEUR   coef_alpha*G*tilde(u)
!
    if (isdiag) then
!
        do ii = 1, nvect_mes
            zr(iaux2-1+ii) = coeff_alpha*zr(ivale_norm-1+ii)*zr(iaux1-1+ii)
        end do
!
    else
        b_n = to_blas_int(nvect_mes)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dspmv('u', b_n, coeff_alpha, zr(ivale_norm), zr(iaux1), &
                   b_incx, 0.d0, zr(iaux2), b_incy)
    end if
! --- RECUPERATION DE LA MATRICE D'OBSERVATION
    call jeveuo(matobs(1), 'L', iobsfil)
    call jeveuo(matobs(2), 'L', iobscol)
    call jeveuo(matobs(3), 'L', iobsval)
!   --- RECUPERATION ET PRECONDITIONNEMENT DU VECTEUR ERC
    call jeveuo(nom_vect_erc//'.VALE', 'L', ivecterc)
    call r8inir(2*obsdim(2), 0.d0, zr(ivecterc), 1)
!
! --- FINALISATION DU PRODUIT   H^T* VECT_AUX_2 (coef_alpha*G*tilde(u))
!
    do ii = 1, obsdim(3)
        n_fil = zi(iobsfil-1+ii)
        n_col = zi(iobscol-1+ii)
        zr(ivecterc-1+obsdim(2)+n_col) = zr( &
                                         ivecterc-1+obsdim(2)+n_col)+zr(iobsval-1+ii)*zr(iaux2-1+&
                                         &n_fil &
                                         )
    end do
!     NETOYAGE DES OBJETS JEVEUX TEMPORAIRES
    call jedetr(baseno//'.VECAUX1ERC.VAL')
    call jedetr(baseno//'.VECAUX2ERC.VAL')
!
end subroutine
