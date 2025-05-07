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

subroutine dintels(cequi, ht, bw, enrobi, enrobs, &
                   scmaxi, scmaxs, ssmax, uc, &
                   ntot, dnsinf, dnssup, nrd, mrd, ndemi) bind(C)
!______________________________________________________________________
!
!      DINTELS

!      CONSTRUCTION DU DIAGRAMME D'INTERACTION D'UNE SECTION
!      FERRAILLEE - VERIFICATION D'UN FERRAILLAGE EXISTANT
!      CRITERE = LIMITATION DES CONTRAINTES (ELS)
!
!      I CEQUI     COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I HT        HAUTEUR DE LA SECTION
!      I BW        LARGEUR DE LA SECTION
!      I ENROBI    ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS    ENROBAGE DES ARMATURES SUPERIEURES
!      I SCMAXI    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE INF
!      I SCMAXS    CONTRAINTE DE COMPRESSION MAXI DU BETON EN FIBRE SUP
!      I SSMAX     CONTRAINTE MAXI DE L'ACIER DE FLEXION
!      I UC        UNITE DES CONTRAINTES :
!                     UC = 0 CONTRAINTES EN Pa
!                     UC = 1 CONTRAINTES EN MPa
!      I DNSINF    DENSITE DE L'ACIER INFERIEUR
!      I DNSSUP    DENSITE DE L'ACIER SUPERIEUR
!
!      I NTOT      DIMENSIONS DES VECTEURS
!      O NRD       VECTEUR DES EFFORTS NORMAUX RESISTANTS (DIAG INTER)
!      O MRD       VECTEUR DES MOMENTS RESISTANTS (DIAG INTER)
!
!______________________________________________________________________
!
!
    use, intrinsic :: iso_c_binding
    implicit none
!
    real(c_double), intent(in) :: cequi
    real(c_double), intent(in) :: ht
    real(c_double), intent(in) :: bw
    real(c_double), intent(in) :: enrobi
    real(c_double), intent(in) :: enrobs
    real(c_double), intent(in) :: scmaxi
    real(c_double), intent(in) :: scmaxs
    real(c_double), intent(in) :: ssmax
    integer(c_int64_t), intent(in) :: uc
    integer(c_int64_t), intent(inout) :: ntot
    real(c_double), intent(in), optional :: dnsinf
    real(c_double), intent(in), optional :: dnssup
    real(c_double), intent(out), optional :: nrd(1:ntot)
    real(c_double), intent(out), optional :: mrd(1:ntot)
    integer(c_int64_t), intent(out), optional :: ndemi

!-----------------------------------------------------------------------
!!!!VARIABLES DE CALCUL
!-----------------------------------------------------------------------

    real(kind=8) :: d, d0, dneg, d0neg
    real(kind=8) :: unite_pa
    real(kind=8) :: X, scmax, scmaxneg
    real(kind=8) :: alpha_12, alpha_12neg, alpha
    real(kind=8) :: ScSUP, ScINF, Delta
    real(kind=8) :: SsSUP, SsINF, Ncc, Mcc
    integer(kind=8) :: N_ET, N_PC, N_EC, N_PCAC, N_PCACN, k
    integer(kind=8) :: N_ECN

    real(kind=8), allocatable :: N_P1(:), M_P1(:)
    real(kind=8), allocatable :: N_P2(:), M_P2(:)
    real(kind=8), allocatable :: N_P3(:), M_P3(:)
    real(kind=8), allocatable :: N_P4(:), M_P4(:)
    real(kind=8), allocatable :: N_P5(:), M_P5(:)
    real(kind=8), allocatable :: N_P6(:), M_P6(:)

!   Paramètres de calcul
    if (uc .eq. 0) then
        unite_pa = 1.d-6
    else if (uc .eq. 1) then
        unite_pa = 1.d0
    end if

    d = ht-enrobi
    d0 = enrobs
    scmax = scmaxs
    dneg = ht-enrobs
    d0neg = enrobi
    scmaxneg = scmaxi

    alpha_12 = 1.0d0/(1.0d0+(ssmax/cequi)/scmax)
    alpha_12neg = 1.0d0/(1.0d0+(ssmax/cequi)/scmaxneg)

!   Initialisation des entiers

    N_ET = 11
    N_PC = 101
    N_EC = ceiling(10.d0*(scmaxneg*unite_pa))+1
    N_PCAC = ceiling((N_PC-1)*(ht/d))+1
    N_PCACN = ceiling((N_PC-1)*(ht/dneg))+1
    N_ECN = ceiling(10.d0*(scmax*unite_pa))+1
    k = 1
    if (ntot < 0) then
        ntot = N_ET+N_PCAC+N_EC+N_ECN+N_PCACN+N_ET
        if (present(ndemi)) then
            ndemi = N_ET+N_PCAC+N_EC
        end if
        return
    end if
    if (.not. present(dnsinf) .or. .not. present(dnssup)) then
        write (6, *) "SyntaxError: dnsinf and dnssup are required"
        return
    end if
    if (.not. present(nrd) .or. .not. present(mrd)) then
        write (6, *) "SyntaxError: nrd and mrd are required"
        return
    end if
    nrd = 0.0d0
    mrd = 0.0d0

!-----------------------------------------------------------------------
!Traitement des différents cas (Pivots A / B / C)
!-----------------------------------------------------------------------

!ET = Entièrement Tendue
!PC = Partiellement Comprimée
!EC = Entièrement Comprimée

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Positif
!---------------------------------------------------------------

    allocate (N_P1(N_ET))
    allocate (M_P1(N_ET))

    do k = 1, N_ET

        SsINF = -ssmax
        ScSUP = -(ssmax/cequi)*(1.d0-0.1d0*(k-1.d0))
        SsSUP = ((SsINF/cequi-ScSUP)*(d0/d)+ScSUP)*cequi
        ScINF = (SsINF/cequi-ScSUP)*(ht/d)+ScSUP
        N_P1(k) = dnsinf*SsINF+dnssup*SsSUP
        M_P1(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)

    end do

!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Positif
!--------------------------------------------------------------------------

    allocate (N_P2(N_PCAC))
    allocate (M_P2(N_PCAC))

    do k = 1, N_PCAC

        if (k .lt. N_PCAC) then
            alpha = (k-1.d0)/(N_PC-1.d0)
        else
            alpha = ht/d
        end if
        X = alpha*d
        if (alpha .lt. alpha_12) then
            SsINF = -ssmax
            ScSUP = (X/(d-X))*(ssmax/cequi)
        else
            ScSUP = scmax
            SsINF = scmax*cequi*(1.d0-d/X)
        end if
        ScINF = ScSUP+(SsINF/cequi-ScSUP)*(ht/d)
        SsSUP = (ScSUP+(SsINF/cequi-ScSUP)*(d0/d))*cequi
        Ncc = ScSUP*0.5d0*X*bw
        Mcc = ScSUP*((1.d0/4.d0)*ht*X-(1.d0/6.d0)*X*X)*bw
        N_P2(k) = dnsinf*SsINF+dnssup*SsSUP+Ncc
        M_P2(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)+Mcc

    end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Positif
!-------------------------------------------------------------------

    allocate (N_P3(N_EC))
    allocate (M_P3(N_EC))

    do k = 1, N_EC

        ScINF = scmaxneg*((k-1.d0)/(N_EC-1.d0))
        ScSUP = scmax
        Delta = ScSUP-ScINF
        if (abs(Delta) .gt. epsilon(Delta)) then
            X = (ScSUP/(ScSUP-ScINF))*ht
            alpha = X/d
        else
            alpha = -1000.d0
        end if
        Ncc = 0.5d0*(ScSUP+ScINF)*ht*bw
        Mcc = (1.d0/12.d0)*(ScSUP-ScINF)*ht*ht*bw
        SsSUP = (ScSUP+(ScINF-ScSUP)*(d0/ht))*cequi
        SsINF = (ScSUP+(ScINF-ScSUP)*(d/ht))*cequi
        N_P3(k) = dnsinf*SsINF+dnssup*SsSUP+Ncc
        M_P3(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)+Mcc

    end do

!Traitement en PIVOT C - Entièrement Comprimée (EC) + Moment Negatif
!-------------------------------------------------------------------

    allocate (N_P4(N_ECN))
    allocate (M_P4(N_ECN))

    do k = 1, N_ECN

        ScSUP = scmax*(1.d0-(k-1.d0)/(N_ECN-1.d0))
        ScINF = scmaxneg
        Delta = ScSUP-ScINF
        if (abs(Delta) .gt. epsilon(Delta)) then
            X = (ScSUP/(ScSUP-ScINF))*ht
            alpha = X/d
        else
            alpha = -1000.d0
        end if
        Ncc = 0.5d0*(ScSUP+ScINF)*ht*bw
        Mcc = (1.d0/12.d0)*(ScINF-ScSUP)*ht*ht*bw
        SsSUP = (ScINF+(ScSUP-ScINF)*(dneg/ht))*cequi
        SsINF = (ScINF+(ScSUP-ScINF)*(d0neg/ht))*cequi
        N_P4(k) = dnsinf*SsINF+dnssup*SsSUP+Ncc
        M_P4(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)-Mcc

    end do

!Traitement en PIVOT A et B - Partiellement Comprimée (PC) + Moment Negatif
!--------------------------------------------------------------------------

    allocate (N_P5(N_PCACN))
    allocate (M_P5(N_PCACN))

    do k = 1, N_PCACN

        if (k .lt. N_PCACN) then
            alpha = ht/dneg-(k-1.d0)/(N_PC-1.d0)
        else
            alpha = 0.d0
        end if
        X = alpha*dneg
        if (alpha .lt. alpha_12neg) then
            SsSUP = -ssmax
            ScINF = (X/(dneg-X))*(ssmax/cequi)
        else
            ScINF = scmaxneg
            SsSUP = scmaxneg*cequi*(1.d0-dneg/X)
        end if
        ScSUP = ScINF+(SsSUP/cequi-ScINF)*(ht/dneg)
        SsINF = (ScINF+(SsSUP/cequi-ScINF)*(d0neg/dneg))*cequi
        Ncc = ScINF*0.5d0*X*bw
        Mcc = ScINF*((1.d0/4.d0)*ht*X-(1.d0/6.d0)*X*X)*bw
        N_P5(k) = dnsinf*SsINF+dnssup*SsSUP+Ncc
        M_P5(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)-Mcc

    end do

!Traitement en PIVOT A - Entièrement Tendu (ET) + Moment Negatif
!---------------------------------------------------------------
    allocate (N_P6(N_ET))
    allocate (M_P6(N_ET))

    do k = 1, N_ET

        SsSUP = -ssmax
        ScINF = -(ssmax/cequi)*0.1d0*(k-1.d0)
        SsINF = ((SsSUP/cequi-ScINF)*(d0neg/dneg)+ScINF)*cequi
        ScSUP = (SsSUP/cequi-ScINF)*(ht/dneg)+ScINF
        N_P6(k) = dnsinf*SsINF+dnssup*SsSUP
        M_P6(k) = -dnsinf*SsINF*(d-0.5d0*ht)+dnssup*SsSUP*(0.5d0*ht-d0)

    end do

!-----------------------------------------------------------------------
!Fin de Traitement des différents cas
!-----------------------------------------------------------------------
    if (N_ET+N_PCAC+N_EC+N_ECN+N_PCACN+N_ET > ntot) then
        write (6, *) "IndexError: ntot argument must be greater than", &
            N_ET+N_PCAC+N_EC+N_ECN+N_PCACN+N_ET
    else
        do k = 1, N_ET
            nrd(k) = N_P1(k)
            mrd(k) = M_P1(k)
        end do
        do k = 1, N_PCAC
            nrd(k+N_ET) = N_P2(k)
            mrd(k+N_ET) = M_P2(k)
        end do
        do k = 1, N_EC
            nrd(k+N_ET+N_PCAC) = N_P3(k)
            mrd(k+N_ET+N_PCAC) = M_P3(k)
        end do
        do k = 1, N_ECN
            nrd(k+N_ET+N_PCAC+N_EC) = N_P4(k)
            mrd(k+N_ET+N_PCAC+N_EC) = M_P4(k)
        end do
        do k = 1, N_PCACN
            nrd(k+N_ET+N_PCAC+N_EC+N_ECN) = N_P5(k)
            mrd(k+N_ET+N_PCAC+N_EC+N_ECN) = M_P5(k)
        end do
        do k = 1, N_ET
            nrd(k+N_ET+N_PCAC+N_EC+N_ECN+N_PCACN) = N_P6(k)
            mrd(k+N_ET+N_PCAC+N_EC+N_ECN+N_PCACN) = M_P6(k)
        end do
    end if

    deallocate (N_P1)
    deallocate (N_P2)
    deallocate (N_P3)
    deallocate (N_P4)
    deallocate (N_P5)
    deallocate (N_P6)
    deallocate (M_P1)
    deallocate (M_P2)
    deallocate (M_P3)
    deallocate (M_P4)
    deallocate (M_P5)
    deallocate (M_P6)

end subroutine
