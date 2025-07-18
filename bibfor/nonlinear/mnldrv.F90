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
subroutine mnldrv(lcal, imat, numdrv, matdrv, xcdl, &
                  parcho, adime, xvect, vecplu, ninc, &
                  nd, nchoc, h, hf)
    implicit none
!
!
!     MODE_NON_LINE -- MATRICE JACOBIENNE
!     -    -                -            -   -
! ----------------------------------------------------------------------
!
! CALCUL LA MATRICE JACOBIENNE POUR UN CERTAIN VECTEUR SOLUTION
! ----------------------------------------------------------------------
! IN  LCAL   : L         : SI .TRUE. ALORS ON RECALCUL LA MATRICE
!                                    SINON ON REMPLACE SEULEMENT LA
!                                           DERNIERE LIGNE DE LA MATRICE
! IN  IMAT   : I(2)      : DESCRIPTEUR DES MATRICES :
!                            - IMAT(1) => MATRICE DE RAIDEUR
!                            - IMAT(2) => MATRICE DE MASSE
! IN  NUMDRV : K14       : NUME_DDL_GENE DE LA MATRICE JACOBIENNE
! IN  MATDRV : K19       : NOM DE  LA MATRICE JACOBIENNE
! IN  XCDL   : K14       : INDICE DES CONDITIONS AUX LIMITES
! IN  PARCHO : K14       : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14       : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  XVECT  : K14       : NOM DU VECTEUR SOLUTION
! IN  VECPLU : R8(NINC)  : NOM DU VECTEUR SOLUTION
! IN  NINC   : I         : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I         : NOMBRE DE DEGRES DE LIBERTE
! IN  NCHOC  : I         : NOMBRE DE CONTACTEURS
! IN  H      : I         : NOMBRE D'HARMONIQUES POUR LE DEPLACEMENT
! IN  HF     : I         : NOMBRE D'HARMONIQUES POUR LA FORCE
! ----------------------------------------------------------------------
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/juveca.h"
#include "asterfort/mnlind.h"
#include "asterfort/mnlldr.h"
#include "asterfort/mnlqd1.h"
#include "asterfort/mnlqd2.h"
#include "asterfort/preres.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    aster_logical :: lcal
    integer(kind=8) :: imat(2), ninc, nd, nchoc, h, hf
    character(len=14) :: numdrv, xcdl, parcho, adime, xvect
    character(len=19) :: matdrv, solveu
    real(kind=8) :: vecplu(ninc)
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    integer(kind=8) :: ivect, ismct, ivat1, ivat2, ivinf, iiinf, ininf
    integer(kind=8) :: ivei, itemp, idrvj, ninf, ndrdv, ind, i, j, ival1, ival2
    integer(kind=8) :: ismhc, ismdi, iret, ibid, irefa, ndrva, nzmk, nind
    character(len=14) :: xei, xtemp, xdrvj, xiinf, xsmct, xvat1, xvat2, xvinf
    character(len=14) :: xninf, numedd, xtemp2
    character(len=8) :: kbid
    character(len=19) :: matk
    real(kind=8) :: eps, cle, res
    integer(kind=8) :: neq, itep2
    integer(kind=8), pointer :: smde(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
    call jemarq()
! ----------------------------------------------------------------------
! --- ON RECUPERE LE VECTEUR SOLUTION
! ----------------------------------------------------------------------
    if (lcal) then
        call jeveuo(xvect, 'L', ivect)
! ----------------------------------------------------------------------
! --- HEURISTIQUE POUR OBTENIR LA TAILLE APPROXIMATIVE DES VALEURES
! ---                                NON NULLES DE LA MATRICE JACOBIENNE
! ----------------------------------------------------------------------
        matk = zk24(zi(imat(1)+1)) (1:19)
        call jelira(matk//'.VALM', 'NMAXOC', nzmk, kbid)
        call jelira(matk//'.VALM', 'LONMAX', nzmk, kbid)
        call jeexin(matdrv//'.VALM', iret)
        if (iret .eq. 0) then
            ndrva = (2*h+1)*nzmk+nchoc*(2*hf+1)+3*((nchoc*(2*hf+1))**2)/4+2*nd*(2*h+1)+nchoc*(2*h&
                    &f+1)*nchoc*(2*h+1)
        else
            call jelira(jexnum(matdrv//'.VALM', 1), 'LONMAX', ndrva, kbid)
            ndrva = 101*ndrva/100
        end if
! ----------------------------------------------------------------------
! --- CREATION DE VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
! --- POUR LE STOCKAGE MORSE
        xsmct = '&&MNLDRV.SMHC'
        call wkvect(xsmct, 'V V I', ndrva, ismct)
        xvat1 = '&&MNLDRV.VAL1'
        call wkvect(xvat1, 'V V R', ndrva, ivat1)
        xvat2 = '&&MNLDRV.VAL2'
        call wkvect(xvat2, 'V V R', ndrva, ivat2)
        xvinf = '&&MNLDRV.VINF'
        call wkvect(xvinf, 'V V R', ndrva, ivinf)
        xiinf = '&&MNLDRV.IINF'
        call wkvect(xiinf, 'V V R', ndrva, iiinf)
        xninf = '&&MNLDRV.NINF'
        call wkvect(xninf, 'V V I', ninc, ininf)
! --- POUR LE CALCUL DE LA MATRICE JACOBIENNE
        xei = '&&MNLDRV.EI'
        call wkvect(xei, 'V V R', ninc, ivei)
        xtemp = '&&MNLDRV.TEMP'
        call wkvect(xtemp, 'V V R', ninc-1, itemp)
        xtemp2 = '&&MNLDRV.TEMP2'
        call wkvect(xtemp2, 'V V R', ninc-1, itep2)
        xdrvj = '&&MNLDRV.DRVJ'
        call wkvect(xdrvj, 'V V R', ninc, idrvj)
!
! ----------------------------------------------------------------------
! --- PARAMETRES
! ----------------------------------------------------------------------
        eps = 0.d0
        ninf = 0
        ndrdv = 0
        res = 0
! ----------------------------------------------------------------------
! --- CALCUL DE LA MATRICE JACOBIENNE
! ----------------------------------------------------------------------
        call jeveuo(numdrv//'.SMOS.SMDI', 'E', ismdi)
        call dismoi('NOM_NUME_DDL', matk, 'MATR_ASSE', repk=numedd)
        neq = zi(imat(1)+2)
        do j = 1, ninc
! ---     MISE A ZERO DES VECTEURS TEMPORAIRES
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(idrvj), b_incx)
            b_n = to_blas_int(ninc-1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itemp), b_incx)
            b_n = to_blas_int(ninc-1)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(itep2), b_incx)
            b_n = to_blas_int(ninc)
            b_incx = to_blas_int(1)
            call dscal(b_n, 0.d0, zr(ivei), b_incx)
! ---     CREATION DU VECTEUR DE BASE CANONIQUE (E_J)
            zr(ivei-1+j) = 1.d0
! ---     CALCUL DE L(E_J)
            call mnlldr(j, imat, neq, ninc, nd, &
                        nchoc, h, hf, parcho, xcdl, &
                        adime, xtemp)
            b_n = to_blas_int(ninc-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, zr(itemp), b_incx, zr(idrvj), &
                       b_incy)
! ---     CALCUL DE Q(V,E_J)
            call mnlqd2(j, imat, neq, ninc, nd, &
                        nchoc, h, hf, parcho, xcdl, &
                        adime, xvect, xtemp)
            b_n = to_blas_int(ninc-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, zr(itemp), b_incx, zr(idrvj), &
                       b_incy)
! ---     CALCUL DE Q(E_J,V)
            call mnlqd1(j, imat, neq, ninc, nd, &
                        nchoc, h, hf, parcho, xcdl, &
                        adime, xvect, xtemp)
            b_n = to_blas_int(ninc-1)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, zr(itemp), b_incx, zr(idrvj), &
                       b_incy)
            zr(idrvj-1+ninc) = 1.d0
! ---     STOCKAGE MORSE
! ---     ON STOCKE LES VALEURS NON NULLES DE LA PARTIE TRIANGULAIRE
! ---     INFERIEURE DE LA MATRICE
            do i = j+1, ninc
                if (abs(zr(idrvj-1+i)) .gt. eps) then
                    ninf = ninf+1
! ---         ON VERIFIE QUE LA TAILLE APPROXIMATIVE DES VECTEURS
! ---         TEMPORAIRES EST SUFFISANTE
                    if (ninf .gt. ndrva) then
                        ndrva = 110*ndrva/100
                        call juveca(xsmct, ndrva)
                        call jeveuo(xsmct, 'E', ismct)
                        call juveca(xvat1, ndrva)
                        call jeveuo(xvat1, 'E', ivat1)
                        call juveca(xvat2, ndrva)
                        call jeveuo(xvat2, 'E', ivat2)
                        call juveca(xvinf, ndrva)
                        call jeveuo(xvinf, 'E', ivinf)
                        call juveca(xiinf, ndrva)
                        call jeveuo(xiinf, 'E', iiinf)
                    end if
! ---         1/ ON CREE UNE CLE UNIQUE A L'AIDE DES INDICES DE LIGNE ET
! ---         DE COLONNES DE LA MATRICE
                    zr(iiinf-1+ninf) = ninc*i+j
! ---         2/ ON RECUPERE LA VALEUR ASSOCIE A CETTE CLE
                    zr(ivinf-1+ninf) = zr(idrvj-1+i)
                end if
            end do
! ---     3/ ON INDIQUE LE NOMBRE DE TERMES NON NULS POUR CETTE COLONNE
            zi(ininf-1+j) = ninf
! ---     ON STOCKE LES VALEURS (NON NULLES) DE LA MATRICE DANS LES
! ---     VECTEURS TEMPORAIRES
            do i = 1, j
! ---       LA DIAGONALE EST TOUJOURS SOTCKEE (POUR LA PARTIE
! ---       SUPERIEURE ET INFERIEURE)
                if (i .eq. j) then
                    ndrdv = ndrdv+1
! ---         ON VERIFIE QUE LA TAILLE APPROXIMATIVE DES VECTEURS
! ---         TEMPORAIRES EST SUFFISANTE
                    if (ndrdv .gt. ndrva) then
                        ndrva = 110*ndrva/100
                        call juveca(xsmct, ndrva)
                        call jeveuo(xsmct, 'E', ismct)
                        call juveca(xvat1, ndrva)
                        call jeveuo(xvat1, 'E', ivat1)
                        call juveca(xvat2, ndrva)
                        call jeveuo(xvat2, 'E', ivat2)
                        call juveca(xvinf, ndrva)
                        call jeveuo(xvinf, 'E', ivinf)
                        call juveca(xiinf, ndrva)
                        call jeveuo(xiinf, 'E', iiinf)
                    end if
                    zi(ismct-1+ndrdv) = i
                    zr(ivat1-1+ndrdv) = zr(idrvj-1+i)
                    zr(ivat2-1+ndrdv) = zr(idrvj-1+i)
! ---       ON TRAITE LE CAS OU UNE VALEUR DANS LA PARTIE SUP. EST !=0
                else if (abs(zr(idrvj-1+i)) .gt. eps) then
                    ndrdv = ndrdv+1
! ---         ON VERIFIE QUE LA TAILLE APPROXIMATIVE DES VECTEURS
! ---         TEMPORAIRES EST SUFFISANTE
                    if (ndrdv .gt. ndrva) then
                        ndrva = 110*ndrva/100
                        call juveca(xsmct, ndrva)
                        call jeveuo(xsmct, 'E', ismct)
                        call juveca(xvat1, ndrva)
                        call jeveuo(xvat1, 'E', ivat1)
                        call juveca(xvat2, ndrva)
                        call jeveuo(xvat2, 'E', ivat2)
                        call juveca(xvinf, ndrva)
                        call jeveuo(xvinf, 'E', ivinf)
                        call juveca(xiinf, ndrva)
                        call jeveuo(xiinf, 'E', iiinf)
                    end if
! ---         ON STOCKE DANS LES VECTEURS TEMPORAIRES
                    zi(ismct-1+ndrdv) = i
                    zr(ivat1-1+ndrdv) = zr(idrvj-1+i)
! ---         ON S'OCCUPE DE LA PARTIE INF.
                    cle = ninc*j+i
                    if (i .eq. 1) then
                        nind = zi(ininf-1+i)
                        call mnlind(nind, 0, cle, zr(iiinf), ind)
                    else
                        nind = zi(ininf-1+i)-zi(ininf-1+i-1)
                        call mnlind(nind, zi(ininf-1+i-1), cle, zr(iiinf-1+zi(ininf-1+i-1)+1), &
                                    ind)
                    end if
                    if (ind .lt. 0) then
                        zr(ivat2-1+ndrdv) = 0.d0
                    else
                        zr(ivat2-1+ndrdv) = zr(ivinf-1+ind)
                    end if
! ---       ON TRAITE LE CAS OU UNE VALEUR DANS LA PARTIE SUP. EST =0
                else
! ---         ON VERIFIE SI LA VALEUR DANS LA PARTIE INF. EST =0 AUSSI
                    cle = ninc*j+i
                    if (i .eq. 1) then
                        nind = zi(ininf-1+i)
                        call mnlind(nind, 0, cle, zr(iiinf), ind)
                    else
                        nind = zi(ininf-1+i)-zi(ininf-1+i-1)
                        call mnlind(nind, zi(ininf-1+i-1), cle, zr(iiinf-1+zi(ininf-1+i-1)+1), &
                                    ind)
                    end if
                    if (ind .gt. 0) then
                        ndrdv = ndrdv+1
! ---         ON VERIFIE QUE LA TAILLE APPROXIMATIVE DES VECTEURS
! ---         TEMPORAIRES EST SUFFISANTE
                        if (ndrdv .gt. ndrva) then
                            ndrva = 110*ndrva/100
                            call juveca(xsmct, ndrva)
                            call jeveuo(xsmct, 'E', ismct)
                            call juveca(xvat1, ndrva)
                            call jeveuo(xvat1, 'E', ivat1)
                            call juveca(xvat2, ndrva)
                            call jeveuo(xvat2, 'E', ivat2)
                            call juveca(xvinf, ndrva)
                            call jeveuo(xvinf, 'E', ivinf)
                            call juveca(xiinf, ndrva)
                            call jeveuo(xiinf, 'E', iiinf)
                        end if
                        zi(ismct-1+ndrdv) = i
                        zr(ivat1-1+ndrdv) = zr(idrvj-1+i)
                        zr(ivat2-1+ndrdv) = zr(ivinf-1+ind)
                    end if
                end if
            end do
            zi(ismdi-1+j) = ndrdv
        end do
        call jeveuo(numdrv//'.SMOS.SMDE', 'E', vi=smde)
        smde(2) = ndrdv
! ---   DESTRUCTION DES CHAMPS A REMPLIR
        call jedetr(matdrv//'.VALM')
        call jedetr(numdrv//'.SMOS.SMHC')
! ---   AU CAS OU LA MATRICE A ETE FACTORISEE
        call jeveuo(matdrv//'.REFA', 'E', irefa)
        zk24(irefa-1+8) = ' '
! ---   REMPLISSAGE CHAMP SMDI, SMHC ET VALE
        call jecrec(matdrv//'.VALM', 'V V R', 'NU', 'DISPERSE', 'VARIABLE', &
                    2)
        call jeecra(jexnum(matdrv//'.VALM', 1), 'LONMAX', ndrdv, ' ')
        call jeecra(jexnum(matdrv//'.VALM', 2), 'LONMAX', ndrdv, ' ')
        call jeveuo(jexnum(matdrv//'.VALM', 1), 'E', ival1)
        call jeveuo(jexnum(matdrv//'.VALM', 2), 'E', ival2)
        call wkvect(numdrv//'.SMOS.SMHC', 'V V S', ndrdv, ismhc)
        b_n = to_blas_int(ndrdv)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivat1), b_incx, zr(ival1), b_incy)
        b_n = to_blas_int(ndrdv)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ivat2), b_incx, zr(ival2), b_incy)
        do i = 1, ndrdv
            zi4(ismhc-1+i) = zi(ismct-1+i)
        end do
! ----------------------------------------------------------------------
! --- DESTRUCTION DES VECTEURS TEMPORAIRES
! ----------------------------------------------------------------------
        call jedetr(xsmct)
        call jedetr(xvat1)
        call jedetr(xvat2)
        call jedetr(xvinf)
        call jedetr(xninf)
        call jedetr(xiinf)
        call jedetr(xei)
        call jedetr(xtemp)
        call jedetr(xtemp2)
        call jedetr(xdrvj)
! ----------------------------------------------------------------------
! --- ON REMPLACE LA DERNIERE LIGNE PAR LE VECTEUR VECPLU
! ----------------------------------------------------------------------
    else
        call jeveuo(matdrv//'.REFA', 'E', irefa)
        zk24(irefa-1+8) = ' '
        call jeveuo(numdrv//'.SMOS.SMDI', 'L', ismdi)
        call jeveuo(jexnum(matdrv//'.VALM', 1), 'E', ival1)
        call jeveuo(jexnum(matdrv//'.VALM', 2), 'E', ival2)
    end if
    zr(ival1-1+zi(ismdi-1+ninc)) = vecplu(ninc)
    b_n = to_blas_int(ninc)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vecplu, b_incx, zr(ival2-1+zi(ismdi-1+ninc-1)+1), b_incy)
! ----------------------------------------------------------------------
! --- FACTORISATION DE LA MATRICE
! ----------------------------------------------------------------------
    solveu = '&&OP0061.SOLVEUR'
    call preres(solveu, 'V', iret, ' ', matdrv, &
                ibid, -9999)
!
    call jedema()
!
end subroutine
