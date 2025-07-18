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
subroutine cfpcdi(resoco, neq, nbliai, tole, epsipc, &
                  mu, apcoef, apddl, appoin, inliac, &
                  matass, solveu, premax, ssgrad, ssgrpr)
!
    implicit none
#include "jeveux.h"
#include "asterfort/caladu.h"
#include "asterfort/calatm.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/r8inir.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dscal.h"
    character(len=24) :: resoco
    integer(kind=8) :: neq, nbliai, apddl(*), appoin(*), inliac(*), premax
    real(kind=8) :: apcoef(*), ssgrad(*), ssgrpr(*), mu(*)
    real(kind=8) :: tole, epsipc
    character(len=19) :: matass, solveu
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (RESOLUTION - GCP)
!
! PRECONDITIONNEMENT DE L'ALGORITHME DU GRADIENT CONJUGUE PROJETE
!
! ----------------------------------------------------------------------
!
! RESOLUTION D'UN PROBLEME ANNEXE A DEPLACEMENT IMPOSE SUR LES NOEUDS
! EFFECTIVEMENT EN CONTACT.
!
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
! IN  NEQ    : NOMBRE D'EQUATIONS DU SYSTEME
! IN  NBLIAI : NOMBRE DE LIAISONS DE CONTACT
! IN  TOLE   : TOLERANCE DE DETECTION DE MU NUL
! IN  EPSIPC : TOLERANCE SOLVEUR ITERATIF
! IN  MATASS : NOM DE LA MATRICE DU PREMIER MEMBRE ASSEMBLEE
! IN  PREMAX : NOMBRE MAXI D'ITERATIONS DU SOLVEUR ITERATIF
! IN  SOLVEU : SD SOLVEUR
! IN  SSGRAD : SOUS-GRADIENT NON-PRECONDITIONNE
! OUT SSGRPR : SOUS-GRADIENT PRECONDITIONNE
!
!
!
!
    integer(kind=8) :: ifm, niv
    real(kind=8) :: numer, denom, conver, alpha
    real(kind=8) :: numerp, numerm, beta
    real(kind=8) :: convm, coef
    integer(kind=8) :: iliac, iliai, jdecal, nbddl, iterat, nbliac
    character(len=24) :: cncin0, secmbr, ddelt, pcresi, pcdire, pcdepl
    integer(kind=8) :: jsecmb, jddelt, jpcres, jpcdir, jpcdep
    character(len=19) :: k19bla
    complex(kind=8) :: c16bid
    parameter(coef=1.d-2)
    integer(kind=8) :: iret
    blas_int :: b_incx, b_incy, b_n
    c16bid = dcmplx(0.d0, 0.d0)
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    k19bla = ' '
!
! --- ACCES AUX CHAMPS DE TRAVAIL
!
    cncin0 = resoco(1:14)//'.CIN0'
    secmbr = resoco(1:14)//'.SECM'
    ddelt = resoco(1:14)//'.DDEL'
    call jeveuo(secmbr(1:19)//'.VALE', 'E', jsecmb)
    call jeveuo(ddelt(1:19)//'.VALE', 'E', jddelt)
    pcresi = resoco(1:14)//'.PCRS'
    pcdire = resoco(1:14)//'.PCDR'
    pcdepl = resoco(1:14)//'.PCUU'
    call jeveuo(pcresi, 'E', jpcres)
    call jeveuo(pcdire, 'E', jpcdir)
    call jeveuo(pcdepl, 'E', jpcdep)
!
! --- INITIALISATION DE DELTA
!
    call r8inir(neq, 0.d0, zr(jddelt), 1)
!
! --- COMPTAGE DU NOMBRE DE LIAISONS REELLEMENT ACTIVES
!
    nbliac = 0
    do iliai = 1, nbliai
        if ((mu(iliai) .gt. tole) .or. (ssgrad(iliai) .gt. epsipc)) then
            nbliac = nbliac+1
            inliac(nbliac) = iliai
        end if
    end do
!
! --- SI AUCUNE LIAISON ACTIVE ON SORT CAR
! --- LE PRECONDITIONNEUR EST INUTILE
!
    if (nbliac .eq. 0) then
        b_n = to_blas_int(nbliai)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, ssgrad, b_incx, ssgrpr, b_incy)
        if (niv .ge. 2) then
            write (ifm, *) '<CONTACT><CALC> PAS DE '//&
     &                  'PRECONDITIONNEMENT (PAS DE LIAISONS ACTIVES)'
        end if
        goto 120
    end if
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT><CALC> PRECONDITIONNEUR DIRICHLET'
        write (ifm, 9010) nbliac, epsipc
    end if
!
! --- NOMBRE D'ITERATIONS MAX
!
    if (premax .eq. 0) premax = 2*nbliac
!
! --- MISE A ZERO DES VECTEURS DE TRAVAIL
!
    call r8inir(neq, 0.d0, zr(jpcdep), 1)
    call r8inir(nbliai, 0.d0, zr(jpcres), 1)
    iterat = 1
!
! ======================================================================
! =========================== BOUCLE PRINCIPALE ========================
! ======================================================================
!
20  continue
!
! --- NOUVELLE VALEUR DU GRADIENT
!
    do iliac = 1, nbliac
        iliai = inliac(iliac)
        jdecal = appoin(iliai)
        nbddl = appoin(iliai+1)-appoin(iliai)
!       RESIDU=A.UU-SSGRAD(ACT)
        call caladu(neq, nbddl, apcoef(1+jdecal), apddl(1+jdecal), zr(jpcdep), &
                    zr(jpcres-1+iliac))
        zr(jpcres-1+iliac) = zr(jpcres-1+iliac)-ssgrad(iliai)
    end do
!
! --- TEST DE CONVERGENCE
!
    conver = -1.d0
    do iliac = 1, nbliac
        conver = max(conver, abs(zr(jpcres-1+iliac)))
    end do
    if (niv .ge. 2) then
        if (iterat .eq. 1) convm = 10*conver/coef
        if (conver .lt. (coef*convm)) then
            write (ifm, 9000) iterat, conver
            convm = conver
        end if
    end if
!
! --- ON A CONVERGE
!
    if (conver .lt. epsipc) then
        if (niv .ge. 2) then
            write (ifm, 9020) iterat, conver
        end if
        goto 90
    end if
!
! --- NOUVELLE DIRECTION DE RECHERCHE
! --- DIRECH=RESIDU+BETA*DIRECH
!
    if (iterat .eq. 1) then
        b_n = to_blas_int(nbliac)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        numerp = ddot(b_n, zr(jpcres), b_incx, zr(jpcres), b_incy)
        b_n = to_blas_int(nbliac)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(jpcres), b_incx, zr(jpcdir), b_incy)
    else
        numerm = numerp
        b_n = to_blas_int(nbliac)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        numerp = ddot(b_n, zr(jpcres), b_incx, zr(jpcres), b_incy)
        beta = numerp/numerm
        b_n = to_blas_int(nbliac)
        b_incx = to_blas_int(1)
        call dscal(b_n, beta, zr(jpcdir), b_incx)
        b_n = to_blas_int(nbliac)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, zr(jpcres), b_incx, zr(jpcdir), &
                   b_incy)
    end if
!
! --- CALCUL DU SECOND MEMBRE
! --- AT.DIRECH
!
    call r8inir(neq, 0.d0, zr(jsecmb), 1)
    do iliac = 1, nbliac
        iliai = inliac(iliac)
        jdecal = appoin(iliai)
        nbddl = appoin(iliai+1)-appoin(iliai)
        call calatm(neq, nbddl, zr(jpcdir-1+iliac), apcoef(1+jdecal), apddl(1+jdecal), &
                    zr(jsecmb))
    end do
!
! --- RESOLUTION
! --- DU=K-1*(AT.DIRECH)
!
    call resoud(matass, k19bla, solveu, cncin0, 0, &
                secmbr, ddelt, 'V', [0.d0], [c16bid], &
                k19bla, .true._1, 0, iret)
    call jeveuo(ddelt(1:19)//'.VALE', 'E', jddelt)
!
! --- PAS D'AVANCEMENT
!
    b_n = to_blas_int(nbliac)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    numer = ddot(b_n, zr(jpcres), b_incx, zr(jpcres), b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    denom = ddot(b_n, zr(jddelt), b_incx, zr(jsecmb), b_incy)
    alpha = numer/denom
!
    if (alpha .lt. 0.d0) then
        call utmess('F', 'CONTACT_7')
    end if
!
! --- ACTUALISATION DU SOUS GRADIENT ET DU DEPLACEMENT
!
    do iliac = 1, nbliac
        iliai = inliac(iliac)
        ssgrpr(iliai) = ssgrpr(iliai)+alpha*zr(jpcdir-1+iliac)
    end do
!
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, -alpha, zr(jddelt), b_incx, zr(jpcdep), &
               b_incy)
!
!
! --- ON A ATTEINT LE NOMBRE D'ITERATION MAXIMAL
    if (iterat .ge. premax) goto 80
!
!
! --- ON N A PAS CONVERGE MAIS IL RESTE DES ITERATIONS A FAIRE
    iterat = iterat+1
    goto 20
!
80  continue
!
!     ON A DEPASSE LE NOMBRE D'ITERATIONS MAX
    if (niv .ge. 2) then
        write (ifm, 9000) iterat, conver
        call utmess('I', 'CONTACT_3', si=premax)
    end if
!
!
90  continue
!
! ======================================================================
! ============================= ON A CONVERGE ==========================
! ======================================================================
!
!     LES CRITERES DE CONVERGENCE SONT DECALES ENTRE L'APPELANT
!     ET CETTE ROUTINE. DU COUP, ON PEUT ENTRER ICI ET S'APERCEVOIR
!     QUE L'ON A RIEN A FAIRE. DANS CE CAS, ON RECOPIE.
    if (iterat .eq. 1) then
        do iliac = 1, nbliac
            iliai = inliac(iliac)
            ssgrpr(iliai) = zr(jpcres-1+iliac)
        end do
    end if
!
!     ON REPROJETE LE SOUS-GRADIENT PRECONDITIONNE POUR
!     ASSURER LA POSITIVITE DES MULTIPLICATEURS
    b_n = to_blas_int(nbliai)
    b_incx = to_blas_int(1)
    call dscal(b_n, -1.d0, ssgrpr, b_incx)
    do iliai = 1, nbliai
        if (mu(iliai) .le. tole) then
            ssgrpr(iliai) = max(ssgrpr(iliai), 0.d0)
        end if
    end do
!
!
120 continue
!
    call jedema()
!
9000 format(' <CONTACT><CALC> PRECONDITIONNEUR : ITERATION =', i6,&
&        ' RESIDU =', 1pe12.5)
9010 format(' <CONTACT><CALC> PRECONDITIONNEUR : ', i6,&
&        ' LIAISON ACTIVES, CRITERE DE CONVERGENCE =', 1pe12.5)
9020 format(' <CONTACT><CALC> PRECONDITIONNEUR : ITERATION =', i6,&
&        ' RESIDU =', 1pe12.5, ' => CONVERGENCE')
end subroutine
