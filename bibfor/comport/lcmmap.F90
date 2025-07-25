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
subroutine lcmmap(fami, kpg, ksp, mult_comp, mod, &
                  imat, nmat, angmas, pgl, materd, &
                  materf, matcst, nbcomm, cpmono, ndt, &
                  ndi, nr, nvi, nfs, nsg, &
                  nhsr, numhsr, hsr)
! aslint: disable=W1504
    implicit none
!       POLYCRISTAL : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                    NB DE CMP DIRECTES/CISAILLEMENT , NB VAR. INTERNES
!
!       OBJETS DE STOCKAGE DES COMPORTEMENTS:
!           MATER(*,1) = E , NU , ALPHA
!           MATER(*,2) = Fractions Volumiques et angles de chaque phase
!                      + COEFFICIENT DE CHAQUE COMPORTEMENT MONOCRSITAL
!                        pour chaque famille de systèmes de glissement
!                        à la température actuelle (MATERF)
!                       et à la température précédente (MATERD)
!           NBCOMM = indices des coefficents de chaque comportement
!                    dans MATER(*,2)
!           CPMONO = noms des différentes "briques" de comportement
!
!      STRUCTURE DES OBJETS CREES
!
!           MATER(*,2) : NBMONO = Nombre de monocristaux
!                        indice debut premier monocristal
!                        indice debut deuxième monocristal
!..............
!                        indice debut denier monocristal
!                        indice des paramètes localisation
!                        Fv et 3 angles par phase
!           pour chaque monocristal différent
!                 par famille de système de glissement
!                    nb coef écoulement + coef,
!                    nb coef écrou isot + coef,
!                    nb coef ecou cine + coef
!                        puis 2 (ou plus) paramètres localisation
!
!           CPMONO(*) : nom de la methode de localisation
!                 puis, pour chaque matériau différent
!                 nom du monocristal, nombre de familles SG, et,
!                    par famille de système de glissement
!                       Nom de la famille
!                       Nom du materiau
!                       Nom de la loi d'écoulement
!                       Nom de la loi d'écrouissage isotrope
!                       Nom de la loi d'écrouissage cinématique
!                       Nom de la loi d'élasticité
!           NBCOMM(*,3) :
!                        Colonne 1      Colonne 2      Colonne3
!                    _____________________________________________
!
!            Ligne 1     Nb phases      Nb var.int.   Nb monocristaux
!                                                     différents
!   pour chaque phase g  Num ligne g    Ind CPMONO    ind frac vol MATER
!   ..................
!   pour chaque phase
!   pour la localisation  indice coef   nb param      0
!   phase g              nb fam g       0            NVIg
!                ... et pour chaque famille de système de glissement :
!             famille 1  ind coef       ind coef      ind coef
!                        ecoulement     ecr iso       ecr cin
!    .....
!         (ind signifie l'indice dans MATER(*,2)
!                    _____________________________________________
!                    VARIABLES INTERNES :
!                    Evp(6)+Norme(Evp)+
!                       Nphase * (Betag (6) ou Evpg(6))
!                       Nphase*(Nsyst*(ALPHAsg,GAMMAsg,Psg)
!                    dernière : indicateur
!                        ou s désigne le SYSTEME DE GLISSEMENT
!                        ou g désigne le "grain" ou la phase
!
!       ----------------------------------------------------------------
!       IN  FAMI   : FAMILLE DES POINTS DE GAUSS
!           KPG    : POINT DE GAUSS
!           KSP    : SOUS-POINT DE GAUSS
!           COMP   :  GRANDEUR COMPOR
!           IMAT   :  ADRESSE DU MATERIAU CODE
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION  MAXIMUM DE MATER
!      ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
!       OUT MATERD :  COEFFICIENTS MATERIAU A T
!           PGL    : MATRICE DE PASSAGE GLOBAL LOCAL
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!                     MATER(*,1) = CARACTERISTIQUES   ELASTIQUES
!                     MATER(*,2) = CARACTERISTIQUES   PLASTIQUES
!           MATCST :  'OUI' SI  MATERIAU A T = MATERIAU A T+DT
!                     'NON' SINON
!           NBCOMM : POSITION DES COEF POUR CHAQUE LOI DE CHAQUE SYSTEME
!           CPMONO : NOMS DES LOIS POUR CHAQUE FAMILLE DE SYSTEME
!
!           NDT    :  NB TOTAL DE COMPOSANTES TENSEURS
!           NDI    :  NB DE COMPOSANTES DIRECTES  TENSEURS
!           NR     :  NB DE COMPOSANTES SYSTEME NL
!           NVI    :  NB DE VARIABLES INTERNES
!     ----------------------------------------------------------------
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/d1ma3d.h"
#include "asterfort/dmat3d.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lcmaec.h"
#include "asterfort/lcmaei.h"
#include "asterfort/lcmafl.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
    integer(kind=8) :: nfs, nsg, kpg, ksp, nmat, ndt, ndi, nr, nvi, nbcomm(nmat, 3)
    integer(kind=8) :: nhsr
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), hook(6, 6)
    real(kind=8) :: hsr(nsg, nsg, nhsr)
    real(kind=8) :: kooh(6, 6), tbsysg
    real(kind=8) :: epsi, angmas(3), pgl(3, 3), hookf(6, 6)
    real(kind=8) :: valres(nmat), ms(6), ng(3), q(3, 3), lg(3)
    character(len=8) :: mod, nomc(14)
    integer(kind=8) :: cerr(14), itbint, nbtbsy, nbsysi, imonoi, imonor, numhsr(nhsr)
    character(len=3) :: matcst
    character(len=*) :: fami
    character(len=16) :: mult_comp
    character(len=16) :: nmater, necoul, necris, necrci, nomfam
    character(len=16) :: phenom, compk, compi, compr, monoi, monor
    character(len=24) :: cpmono(5*nmat+1)
    integer(kind=8) :: i, imat, nbfsys, ifa, j, nbmono, nbsys, nbsyst, idecal
    integer(kind=8) :: nbphas, icompk, icompi, icompr, dimk, tabicp(nmat), nvloc
    integer(kind=8) :: indmat, indcp, imono, nbval, indloc, indcom, iphas, nbfam
    integer(kind=8) :: numono, nvintg, idmono, nbval1, nbval2, nbval3, nbcoef
    blas_int :: b_incx, b_incy, b_n
    common/tbsysg/tbsysg(900)
!     ----------------------------------------------------------------
    call jemarq()
!
! -   NB DE COMPOSANTES / VARIABLES INTERNES -------------------------
!
    if (mod(1:2) .eq. '3D') then
        ndt = 6
        ndi = 3
    else if (mod(1:6) .eq. 'D_PLAN' .or. mod(1:4) .eq. 'AXIS') then
        ndt = 6
        ndi = 3
    else if (mod(1:6) .eq. 'C_PLAN') then
        ndt = 6
        ndi = 3
    end if
    call r8inir(2*nmat, 0.d0, materd, 1)
    call r8inir(2*nmat, 0.d0, materf, 1)
!
!
!
!     LA DERNIERE VARIABLE INTERNE EST L'INDICATEUR PLASTIQUE
    nr = nvi+ndt-1
    compk = mult_comp(1:8)//'.CPRK'
    compi = mult_comp(1:8)//'.CPRI'
    compr = mult_comp(1:8)//'.CPRR'
    call jeveuo(compk, 'L', icompk)
    call jeveuo(compi, 'L', icompi)
    call jeveuo(compr, 'L', icompr)
    nbphas = zi(icompi-1+2)
    dimk = zi(icompi-1+5+3*nbphas)
    nvloc = zi(icompi-1+5+3*nbphas+1)
!
    itbint = 0
!
    do i = 1, nmat
        do j = 1, 3
            nbcomm(i, j) = 0
        end do
    end do
!
    nbmono = zi(icompi-1+4)
    do i = 1, dimk
        cpmono(i) = zk24(icompk-1+i)
    end do
!           MATER(*,2) : Nombre de monocristaux
!                        indice debut premier monocristal
!                        indice debut deuxième monocristal
!..............
!                        indice debut denier monocristal
!                        indice des paramètes localisation
!                        Fv et 3 angles par phase
!           pour chaque monocristal différent
!                 par famille de système de glissement
!                    coef écoulement, coef écrou isot, coef ecou cine
!                        puis 2 (ou plus) paramètres localisation
    materd(1, 2) = nbmono
    materf(1, 2) = nbmono
    indmat = 1+nbmono+1
    do i = 1, 4*nbphas
        materd(indmat+i, 2) = zr(icompr-1+i)
        materf(indmat+i, 2) = zr(icompr-1+i)
    end do
    indmat = 1+nbmono+1+4*nbphas
    materd(1+1, 2) = indmat+1
    materf(1+1, 2) = indmat+1
    indcp = 3
!     Boucle sur le nombre de monocristaux
    do imono = 1, nbmono
        read (cpmono(indcp), '(I24)') nbfsys
!        LESCTURE EVENTUELLE DE LA TABLE DE SYST. GLISS
        monoi = cpmono(indcp-1) (1:8)//'.CPRI'
        monor = cpmono(indcp-1) (1:8)//'.CPRR'
        call jeveuo(monoi, 'L', imonoi)
        ASSERT(nbfsys .eq. zi(imonoi-1+5))
        itbint = zi(imonoi-1+4)
        nbsyst = zi(imonoi-1+8)
        nbtbsy = 0
        do ifa = 1, nbfsys
            nbsysi = zi(imonoi-1+8+ifa)
            nbtbsy = nbtbsy+nbsysi
        end do
        if (nbtbsy .ne. 0) then
            call r8inir(900, 0.d0, tbsysg, 1)
            call jeveuo(monor, 'L', imonor)
!           TABLE CONTENANT LES SYSTEMES
            b_n = to_blas_int(6*nbtbsy+12)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(imonor), b_incx, tbsysg, b_incy)
        else
            tbsysg(1) = 0.d0
        end if
!        TABLE CONTENANT LA MATRICE D'INTERACTION
        if (itbint .eq. 1) then
            idecal = 0
            if (nbtbsy .eq. 0) then
                call jeveuo(monor, 'L', imonor)
            end if
            idecal = nint(zr(imonor+1))
            do i = 1, nbsyst
                do j = 1, nbsyst
                    hsr(i, j, imono) = zr(imonor-2+idecal+nbsyst*(i-1)+j)
                end do
            end do
        end if
!
        do ifa = 1, nbfsys
            nomfam = cpmono(indcp+5*(ifa-1)+1) (1:16)
            call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                        ng, lg, 0, q)
            nmater = cpmono(indcp+5*(ifa-1)+2) (1:16)
            necoul = cpmono(indcp+5*(ifa-1)+3) (1:16)
            necris = cpmono(indcp+5*(ifa-1)+4) (1:16)
            necrci = cpmono(indcp+5*(ifa-1)+5) (1:16)
!           NOMBRE DE MATRICE D'INTERACTION DIFFERENTES
!           COEFFICIENTS MATERIAUX LIES A L'ECOULEMENT
            call lcmafl(fami, kpg, ksp, '-', nmater, &
                        imat, necoul, nbval, valres, nmat, &
                        itbint, nfs, nsg, hsr(1, 1, imono), nbsys)
            materd(indmat+1, 2) = nbval
            materf(indmat+1, 2) = nbval
            indmat = indmat+1
            do i = 1, nbval
                materd(indmat+i, 2) = valres(i)
            end do
            call lcmafl(fami, kpg, ksp, '+', nmater, &
                        imat, necoul, nbval, valres, nmat, &
                        itbint, nfs, nsg, hsr(1, 1, imono), nbsys)
            do i = 1, nbval
                materf(indmat+i, 2) = valres(i)
            end do
            indmat = indmat+nbval
!           COEFFICIENTS MATERIAUX LIES A L'ECROUISSAGE CINEMATIQUE
            call lcmaec(fami, kpg, ksp, '-', nmater, &
                        imat, necrci, nbval, valres, nmat)
            materd(indmat+1, 2) = nbval
            materf(indmat+1, 2) = nbval
            indmat = indmat+1
            do i = 1, nbval
                materd(indmat+i, 2) = valres(i)
            end do
            call lcmaec(fami, kpg, ksp, '+', nmater, &
                        imat, necrci, nbval, valres, nmat)
            do i = 1, nbval
                materf(indmat+i, 2) = valres(i)
            end do
            indmat = indmat+nbval
!           COEFFICIENTS MATERIAUX LIES A L'ECROUISSAGE ISOTROPE
            call lcmaei(fami, kpg, ksp, '-', nmater, &
                        imat, necris, necoul, nbval, valres, &
                        nmat, itbint, nfs, nsg, hsr(1, 1, imono), &
                        ifa, nomfam, nbsys)
            materd(indmat+1, 2) = nbval
            materf(indmat+1, 2) = nbval
            indmat = indmat+1
            do i = 1, nbval
                materd(indmat+i, 2) = valres(i)
            end do
            call lcmaei(fami, kpg, ksp, '+', nmater, &
                        imat, necris, necoul, nbval, valres, &
                        nmat, itbint, nfs, nsg, hsr(1, 1, imono), &
                        ifa, nomfam, nbsys)
            do i = 1, nbval
                materf(indmat+i, 2) = valres(i)
            end do
            indmat = indmat+nbval
!
        end do
        tabicp(imono) = indcp
        indcp = indcp+5*nbfsys+1+2
!        INDICE DU DEBUT DU MONO SUIVANT DANS MATER
        materd(1+imono+1, 2) = indmat+1
        materf(1+imono+1, 2) = indmat+1
    end do
!     Paramètres de la loi de localisation
    indloc = indmat+1
    do i = 1, nvloc
        materd(indmat+i, 2) = zr(icompr-1+4*nbphas+i)
        materf(indmat+i, 2) = zr(icompr-1+4*nbphas+i)
    end do
    nbcoef = indmat+nvloc
!
!  FIN remplissage de MATER(*,2)
!
    call rccoma(imat, 'ELAS', 1, phenom, cerr(1))
    if (phenom .eq. 'ELAS') then
!
! -    ELASTICITE ISOTROPE
!
        nomc(1) = 'E       '
        nomc(2) = 'NU      '
        nomc(3) = 'ALPHA   '
!
! -     RECUPERATION MATERIAU A TEMPD (T)
!
        call rcvalb(fami, kpg, ksp, '-', imat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomc(1), materd(1, 1), cerr(1), 1)
        call rcvalb(fami, kpg, ksp, '-', imat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, nomc(3), materd(3, 1), cerr(3), 0)
        if (cerr(3) .ne. 0) materd(3, 1) = 0.d0
        materd(nmat, 1) = 0
!
! -     RECUPERATION MATERIAU A TEMPF (T+DT)
!
        call rcvalb(fami, kpg, ksp, '+', imat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomc(1), materf(1, 1), cerr(1), 1)
        call rcvalb(fami, kpg, ksp, '+', imat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, nomc(3), materf(3, 1), cerr(3), 0)
        if (cerr(3) .ne. 0) materf(3, 1) = 0.d0
        materf(nmat, 1) = 0
!
    else if (phenom .eq. 'ELAS_ORTH') then
! -    ELASTICITE ORTHOTROPE
! -     MATRICE D'ELASTICITE ET SON INVERSE A TEMPD(T)
!
        call dmat3d(fami, imat, r8vide(), '-', kpg, &
                    ksp, angmas, hook)
        call d1ma3d(fami, imat, r8vide(), '-', kpg, &
                    ksp, angmas, kooh)
!
!         termes  SQRT(2) qui ne sont pas mis dans DMAT3D
        do j = 4, 6
            do i = 1, 6
                hook(i, j) = hook(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                hook(i, j) = hook(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 4, 6
            do i = 1, 6
                kooh(i, j) = kooh(i, j)/sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                kooh(i, j) = kooh(i, j)/sqrt(2.d0)
            end do
        end do
        do i = 1, 6
            do j = 1, 6
                materd(6*(j-1)+i, 1) = hook(i, j)
                materd(36+6*(j-1)+i, 1) = kooh(i, j)
            end do
        end do
        materd(nmat, 1) = 1
        nomc(1) = 'ALPHA_L'
        nomc(2) = 'ALPHA_T'
        nomc(3) = 'ALPHA_N'
        call rcvalb(fami, kpg, ksp, '-', imat, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    3, nomc, materd(73, 1), cerr, 0)
        if (cerr(1) .ne. 0) materd(73, 1) = 0.d0
        if (cerr(2) .ne. 0) materd(74, 1) = 0.d0
        if (cerr(3) .ne. 0) materd(75, 1) = 0.d0
!
! -     MATRICE D'ELASTICITE ET SON INVERSE A A TEMPF (T+DT)
        call dmat3d(fami, imat, r8vide(), '+', kpg, &
                    ksp, angmas, hookf)
        call d1ma3d(fami, imat, r8vide(), '+', kpg, &
                    ksp, angmas, kooh)
!       termes  SQRT(2) qui ne sont pas mis dans DMAT3D
        do j = 4, 6
            do i = 1, 6
                hookf(i, j) = hookf(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                hookf(i, j) = hookf(i, j)*sqrt(2.d0)
            end do
        end do
        do j = 4, 6
            do i = 1, 6
                kooh(i, j) = kooh(i, j)/sqrt(2.d0)
            end do
        end do
        do j = 1, 6
            do i = 4, 6
                kooh(i, j) = kooh(i, j)/sqrt(2.d0)
            end do
        end do
        do i = 1, 6
            do j = 1, 6
                materf(6*(j-1)+i, 1) = hookf(i, j)
                materf(36+6*(j-1)+i, 1) = kooh(i, j)
            end do
        end do
        materf(nmat, 1) = 1
        call rcvalb(fami, kpg, ksp, '+', imat, &
                    ' ', phenom, 0, ' ', [0.d0], &
                    3, nomc, materf(73, 1), cerr, 0)
        if (cerr(1) .ne. 0) materf(73, 1) = 0.d0
        if (cerr(2) .ne. 0) materf(74, 1) = 0.d0
        if (cerr(3) .ne. 0) materf(75, 1) = 0.d0
    else
        call utmess('F', 'ALGORITH4_65', sk=phenom)
    end if
!
!     Remplissage de NBCOMM : Boucle sur le nombre de phases
!     3 : Nombre de familles de systèmes de glissement phase g
!     Nb phases
    nbcomm(1, 1) = zi(icompi-1+2)
!     Nb var. int. total
    nbcomm(1, 2) = zi(icompi-1+3)
!     Nb materiaux (monocristaux) différents
    nbcomm(1, 3) = zi(icompi-1+4)
!     2 : Numéro du matériau phase g
!     3 : indice fraction volumique  dans MATER
!     Indices de la première phase    dans NBCOMM
    indcom = nbphas+1
!     Ligne précisant l'indice des coef localisation
    indcom = indcom+1
    nbcomm(indcom, 1) = indloc
    nbcomm(indcom, 2) = nvloc
    nbcomm(2, 1) = indcom+1
    indcom = indcom+1
    indcp = 1
    do iphas = 1, nbphas
!         Indice de la fraction volumique dans MATER
        nbcomm(1+iphas, 3) = 1+nbcomm(1, 3)+1+4*(iphas-1)+1
        nbfam = zi(icompi-1+4+3*(iphas-1)+1)
        numono = zi(icompi-1+4+3*(iphas-1)+2)
        numhsr(iphas) = numono
        nvintg = zi(icompi-1+4+3*(iphas-1)+3)
!        Indice du monocristal dans CPMONO
        nbcomm(1+iphas, 2) = tabicp(numono)
        nbcomm(indcom, 1) = nbfam
!        Indice début monocristal dans CPMONO
        nbcomm(indcom, 3) = nvintg
!        Debut du monocristal
        idmono = nint(materd((1+numono), 2))
!        Indice coef ecoulement famille ifa dans MATER
        do ifa = 1, nbfam
            nbval1 = nint(materd(idmono, 2))
            nbcomm(indcom+ifa, 1) = idmono+1
            nbval2 = nint(materd((idmono+1+nbval1), 2))
            nbcomm(indcom+ifa, 2) = idmono+1+nbval1+1
            nbval3 = nint(materd((idmono+1+nbval1+1+nbval2), 2))
            nbcomm(indcom+ifa, 3) = idmono+1+nbval1+1+nbval2+1
            idmono = idmono+3+nbval1+nbval2+nbval3
        end do
        indcom = indcom+nbfam+1
        if (iphas .lt. nbphas) then
!           Indice du monocristal (phase) dans NBCOMM
            nbcomm(1+iphas+1, 1) = indcom
        end if
    end do
!     Nombre total de COEF
    if (nbcoef .gt. nmat) then
        call utmess('F', 'COMPOR2_6')
    end if
    if (nbphas .gt. nmat) then
        call utmess('F', 'COMPOR2_6')
    end if
    nbcomm(nmat, 3) = nbcoef
!
! -   MATERIAU CONSTANT ?
    matcst = 'OUI'
    epsi = r8prem()
    do i = 1, nmat
        if (abs(materd(i, 1)-materf(i, 1)) .gt. epsi*materd(i, 1)) then
            matcst = 'NON'
            goto 999
        end if
    end do
    do i = 1, nmat
        if (abs(materd(i, 2)-materf(i, 2)) .gt. epsi*materd(i, 2)) then
            matcst = 'NON'
            call utmess('F', 'COMPOR1_28')
            goto 999
        end if
    end do
!
999 continue
    ASSERT(nmat .gt. indcom)
    ASSERT((5*nmat+1) .gt. dimk)
    ASSERT(nmat .gt. (indmat+nvloc))
!
!     ON STOCKE A LA FIN LE NOMBRE TOTAL DE COEF MATERIAU
    materd(nmat, 2) = nbcomm(nmat, 3)
    materf(nmat, 2) = nbcomm(nmat, 3)
!
    call jedema()
end subroutine
