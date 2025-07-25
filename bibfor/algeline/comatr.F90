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
subroutine comatr(option, typev, nbproc, rang, vnconv, &
                  dim1i, dim2i, vecti, dim1r, dim2r, &
                  vectr, dim1c, dim2c, vectc)
!     COMMUNICATION VIA LE COMMUNICATEUR MPI COURANT D'UNE MATRICE SOIT
!     REELLE, SOIT ENTIERE, SOIT DE CHAR*, SOIT COMPLEXE.
!     ROUTINE CREE POUR LES BESOINS DU PARALLELISME MPI DANS LES MACROS.
!     EN INPUT: SEULES LES VNCONV(RANG) COLONNES (SI OPTION='S') ET
!              LIGNES (SI OPTION='T') SONT SIGNIFIANTES POUR LE PROCES
!              SUS COURANT.
!     EN OUTPUT: TOUS LES PROCESSUS RECUPERENT LA MEME MATRICE. ELLE EST
!              COMPOSEE DES NBPROC PAQUETS DE VNCONV(I) COLONNES
!              (RESP. LIGNES) DE TOUS LES PROCESSUS I. CHAQUE PAQUET EST
!              RANGE DS LA MATRICE PAR ORDRE DE RANG CROISSANT:
!              EN PREMIER LES VNCONV(1) COLONNES OU LIGNES DU PROCESSUS
!              DE RANG 0, PUIS LES VNCONV(2) DE CELUI DE RANG 1...
! ======================================================================
! IN  OPTION  : K1  : 'S' POUR STANDARD , 'T' POUR TRANSPOSE.
! IN  TYPEV   : K1  : 'R' POUR REEL (AVEC VECTR), 'I' POUR ENTIER (AVEC
!                    VECTI) ET 'C' POUR COMPLEXE (VECTC).
! IN  NBPROC  : IS  : ENTIER CORRESPONDANT AU NBRE DE PROCESSUS MPI.
! IN  RANG    : IS  : ENTIER CORRESPONDANT AU RANG DU PROCESSUS MPI.
! IN  DIM1   : IS  : VECTEURS DU NBRE DE LIGNES DE LA MATRICE CONSIDEREE
! IN  DIM2   : IS  : IDEM NBRE DE COLONNES.
! IN  VNCONV  : IS  : VECTEUR DE NBPROC ENTIERS CORRESPONDANT AUX
!                     DECALAGES PAR PROC.
! IN/OUT VECTI: IS  : MATRICE D'ENTIERS A COMMUNIQUER (DIM1I X DIM2I)
! IN/OUT VECTR: R8  : MATRICE REELLE A COMMUNIQUER    (DIM1R X DIM2R)
! IN/OUT VECTC: C8  : MATRICE COMPLEXE A COMMUNIQUER  (DIM1C X DIM2C)
! ======================================================================
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
! PARAMETRES D'APPEL
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/vecinc.h"
#include "asterfort/vecini.h"
#include "asterfort/vecint.h"
#include "blas/dnrm2.h"
    integer(kind=8) :: nbproc, rang, dim1i, dim2i, dim1r, dim2r, dim1c, dim2c
    integer(kind=8) :: vnconv(nbproc), vecti(dim1i, *)
    real(kind=8) :: vectr(dim1r, *)
    complex(kind=8) :: vectc(dim1c, *)
    character(len=1) :: option, typev
!
! VARIABLES LOCALES
    integer(kind=8) :: nconv, nconvg, i, j, idecal, iaux1, izero, idim1, idim2
    integer(kind=8) :: ifm, niv
    real(kind=8) :: rzero, res
    complex(kind=8) :: czero, dcmplx
    aster_logical :: ldebug
    blas_int :: b_incx, b_n
!
! --- INIT.
    call jemarq()
    call infniv(ifm, niv)
    izero = 0
    rzero = 0.d0
    czero = dcmplx(0.d0, 0.d0)
    ldebug = .false.
!      LDEBUG=.TRUE.
!
! ----------------------------------------------------------------------
! --- VERIF PARAMETRES INPUT
!-----------------------------------------------------------------------
    ASSERT((option .eq. 'S') .or. (option .eq. 'T'))
    ASSERT((typev .eq. 'R') .or. (typev .eq. 'I') .or. (typev .eq. 'C'))
    ASSERT((nbproc .ge. 1) .and. (rang .ge. 0) .and. (rang+1 .le. nbproc))
!
    if (typev .eq. 'I') then
        idim1 = dim1i
        idim2 = dim2i
    else if (typev .eq. 'R') then
        idim1 = dim1r
        idim2 = dim2r
    else if (typev .eq. 'C') then
        idim1 = dim1c
        idim2 = dim2c
    else
        ASSERT(.false.)
    end if
!
! ----------------------------------------------------------------------
! --- CALCULS PRELIMINAIRES
!-----------------------------------------------------------------------
! --- NCONV:  NBRE DE PREMIERES COLONNES A DECALER
! --- NCONVG: SOMME DE DECALAGES
! --- IDECAL: DECALAGE POUR LE PROC COURANT
    nconv = vnconv(rang+1)
    nconvg = 0
    idecal = 0
    do i = 1, nbproc
        ASSERT(vnconv(i) .ge. 0)
        if ((i-1) .lt. rang) idecal = idecal+vnconv(i)
        nconvg = nconvg+vnconv(i)
    end do
    if (option .eq. 'S') then
        ASSERT(idim2 .eq. nconvg)
    else if (option .eq. 'T') then
        if (idim1 .ne. nconvg) then
            ASSERT(.false.)
        end if
    end if
!
! --- VERIF INIT.
    if (ldebug) then
        write (ifm, *) 'INITIALISATION***************************'
        if ((typev .eq. 'R') .and. (option .eq. 'S')) then
            do j = 1, idim2
                b_n = to_blas_int(idim1)
                b_incx = to_blas_int(1)
                res = dnrm2(b_n, vectr(1, j), b_incx)
                write (ifm, *) j, res
            end do
        else if ((typev .eq. 'R') .and. (option .eq. 'T')) then
! --- ON NE FAIT QU'IMPRIMER LES TERMES CAR CERTAINS SONT EN 1.E+308
            do i = 1, idim1
                write (ifm, *) i, (vectr(i, j), j=1, idim2)
            end do
        else
            write (ifm, *) '! ATTENTION: DEBUG OPTION NON PRISE EN COMPTE !'
        end if
    end if
!
! ----------------------------------------------------------------------
! --- COMMUNICATIONS PROPREMENTS DITES
!-----------------------------------------------------------------------
! --- STEP 1:
! --- POUR LE PROCESSUS COURANT, ON INITIALISE LA FIN DE LA MATRICE
! --- A ZERO: COMME SEULS LES NCONV PREMIERES COLONNES (RESP. LIGNES)
! --- SONT SIGNIFIANTES.
    if (option .eq. 'S') then
        iaux1 = idim1*(idim2-nconv)
    else
        iaux1 = idim1-nconv
    end if
    if ((option .eq. 'S') .and. (iaux1 .gt. 0)) then
!
        if (typev .eq. 'R') then
            call vecini(iaux1, rzero, vectr(1, nconv+1))
        else if (typev .eq. 'I') then
            call vecint(iaux1, izero, vecti(1, nconv+1))
        else if (typev .eq. 'C') then
            call vecinc(iaux1, czero, vectc(1, nconv+1))
        end if
!
    else if ((option .eq. 'T') .and. (iaux1 .gt. 0)) then
!
        if (typev .eq. 'R') then
            do j = 1, idim2
                call vecini(iaux1, rzero, vectr(nconv+1, j))
            end do
        else if (typev .eq. 'I') then
            do j = 1, idim2
                call vecint(iaux1, izero, vecti(nconv+1, j))
            end do
        else if (typev .eq. 'C') then
            do j = 1, idim2
                call vecinc(iaux1, czero, vectc(nconv+1, j))
            end do
        end if
!
    end if
!
! --- VERIF STEP 1.
    if (ldebug) then
        write (ifm, *) 'STEP 1***************************'
        if ((typev .eq. 'R') .and. (option .eq. 'S')) then
            do j = 1, idim2
                b_n = to_blas_int(idim1)
                b_incx = to_blas_int(1)
                res = dnrm2(b_n, vectr(1, j), b_incx)
                write (ifm, *) j, res
            end do
        else if ((typev .eq. 'R') .and. (option .eq. 'T')) then
            do i = 1, idim1
                write (ifm, *) i, (vectr(i, j), j=1, idim2)
            end do
        else
            write (ifm, *) '! ATTENTION: DEBUG OPTION NON PRISE EN COMPTE !'
        end if
    end if
!
! --- STEP 2:
! --- ON DECALE LES NCONV PREMIERES LIGNES OU COLONNES POUR LES METTRE
! --- BIEN EN PLACE DS LE BUFFER DE COMMUNICATION. ON DECALE EN COMMEN
! --- CANT PAR LES COLONNES OU LES LIGNES LES PLUS ELOIGNEES DE MANIERE
! --- A NE PAS ECRASER DE DONNEES.
    if ((option .eq. 'S') .and. (idecal .gt. 0)) then
!
        if (typev .eq. 'R') then
            do j = nconv, 1, -1
                do i = 1, idim1
                    vectr(i, j+idecal) = vectr(i, j)
                end do
            end do
        else if (typev .eq. 'I') then
            do j = nconv, 1, -1
                do i = 1, idim1
                    vecti(i, j+idecal) = vecti(i, j)
                end do
            end do
        else if (typev .eq. 'C') then
            do j = nconv, 1, -1
                do i = 1, idim1
                    vectc(i, j+idecal) = vectc(i, j)
                end do
            end do
        end if
!
    else if ((option .eq. 'T') .and. (idecal .gt. 0)) then
!
        if (typev .eq. 'R') then
            do j = 1, idim2
                do i = nconv, 1, -1
                    vectr(i+idecal, j) = vectr(i, j)
                end do
            end do
        else if (typev .eq. 'I') then
            do j = 1, idim2
                do i = nconv, 1, -1
                    vecti(i+idecal, j) = vecti(i, j)
                end do
            end do
        else if (typev .eq. 'C') then
            do j = 1, idim2
                do i = nconv, 1, -1
                    vectc(i+idecal, j) = vectc(i, j)
                end do
            end do
        end if
!
    end if
!
! --- VERIF STEP2.
    if (ldebug) then
        write (ifm, *) 'STEP 2***************************'
        if ((typev .eq. 'R') .and. (option .eq. 'S')) then
            do j = 1, idim2
                b_n = to_blas_int(idim1)
                b_incx = to_blas_int(1)
                res = dnrm2(b_n, vectr(1, j), b_incx)
                write (ifm, *) j, res
            end do
        else if ((typev .eq. 'R') .and. (option .eq. 'T')) then
            do i = 1, idim1
                write (ifm, *) i, (vectr(i, j), j=1, idim2)
            end do
        else
            write (ifm, *) '! ATTENTION: DEBUG OPTION NON PRISE EN COMPTE !'
        end if
    end if
!
! --- STEP 3:
! --- ON ANNULE LES IDECAL PREMIERES COLONNES OU LIGNES
    if ((option .eq. 'S') .and. (idecal .gt. 0)) then
!
        iaux1 = idim1*idecal
        if (typev .eq. 'R') then
            call vecini(iaux1, rzero, vectr(1, 1))
        else if (typev .eq. 'I') then
            call vecint(iaux1, izero, vecti(1, 1))
        else if (typev .eq. 'C') then
            call vecinc(iaux1, czero, vectc(1, 1))
        end if
!
    else if ((option .eq. 'T') .and. (idecal .gt. 0)) then
!
        if (typev .eq. 'R') then
            do j = 1, idim2
                call vecini(idecal, rzero, vectr(1, j))
            end do
        else if (typev .eq. 'I') then
            do j = 1, idim2
                call vecint(idecal, izero, vecti(1, j))
            end do
        else if (typev .eq. 'C') then
            do j = 1, idim2
                call vecinc(idecal, czero, vectc(1, j))
            end do
        end if
!
    end if
!
! --- VERIF STEP3.
    if (ldebug) then
        write (ifm, *) 'STEP 3***************************'
        if ((typev .eq. 'R') .and. (option .eq. 'S')) then
            do j = 1, idim2
                b_n = to_blas_int(idim1)
                b_incx = to_blas_int(1)
                res = dnrm2(b_n, vectr(1, j), b_incx)
                write (ifm, *) j, res
            end do
        else if ((typev .eq. 'R') .and. (option .eq. 'T')) then
            do i = 1, idim1
                write (ifm, *) i, (vectr(i, j), j=1, idim2)
            end do
        else
            write (ifm, *) '! ATTENTION: DEBUG OPTION NON PRISE EN COMPTE !'
        end if
    end if
!
! --- STEP 4 FINAL:
! --- ON COMMUNIQUE TOUTE LA MATRICE
    iaux1 = idim1*idim2
    if (typev .eq. 'R') then
        call asmpi_comm_vect('MPI_SUM', 'R', nbval=iaux1, vr=vectr(1, 1))
    else if (typev .eq. 'I') then
        call asmpi_comm_vect('MPI_SUM', 'I', nbval=iaux1, vi=vecti(1, 1))
    else if (typev .eq. 'C') then
        call asmpi_comm_vect('MPI_SUM', 'C', nbval=iaux1, vc=vectc(1, 1))
    end if
!
! --- VERIF FINALIZATION.
    if (ldebug) then
        write (ifm, *) 'FINALISATION***************************'
        if ((typev .eq. 'R') .and. (option .eq. 'S')) then
            do j = 1, idim2
                b_n = to_blas_int(idim1)
                b_incx = to_blas_int(1)
                res = dnrm2(b_n, vectr(1, j), b_incx)
                write (ifm, *) j, res
            end do
        else if ((typev .eq. 'R') .and. (option .eq. 'T')) then
            do i = 1, idim1
                write (ifm, *) i, (vectr(i, j), j=1, idim2)
            end do
        else
            write (ifm, *) '! ATTENTION: DEBUG OPTION NON PRISE EN COMPTE !'
        end if
    end if
    call jedema()
end subroutine
