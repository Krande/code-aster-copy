! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

subroutine dis_choc_frot_syme(DD, icodma, ulp, xg, klv, kgv, &
                              dvl, dpe, dve, Predic, force, varmo, varpl)
!
    use te0047_type
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/diraidklv.h"
#include "asterfort/diklvraid.h"
#include "asterfort/in_liste_entier.h"
#include "asterfort/rcvala.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpsgl.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: DD
    integer         :: icodma
    real(kind=8)    :: ulp(*), dvl(*)
    real(kind=8)    :: dpe(*), dve(*)
    real(kind=8)    :: klv(*), xg(*), kgv(*)
    real(kind=8)    :: varmo(*), varpl(*), force(*)
    aster_logical   :: Predic
!
! --------------------------------------------------------------------------------------------------
!
!     RELATION DE COMPORTEMENT "DIS_CHOC"
!
! --------------------------------------------------------------------------------------------------
! in :
!       icodma  : adresse du materiau code
!       ulp     : deplacement
!       xg      : coordonnees des noeuds repere global
!       dvl     : vitesse
!       dpe     : déplacement d'entrainement
!       dve     : vitesse d'entrainement
!       varmo   : variables internes (temps moins)
!       kgv     : matrice de raideur repère global          (triangulaire supérieure)
! in/out :
!       klv     : matrice de raideur symétrique initiale    (triangulaire supérieure)
!       klv     : matrice de raideur symétrique             (triangulaire supérieure)
! out :
!       force   : efforts
!       varpl   : variables internes (temps plus)
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    integer, parameter  :: nbre1 = 9
    real(kind=8)        :: valre1(nbre1)
    integer             :: codre1(nbre1)
    character(len=12)   :: nomre1(nbre1)
!   Index des variables internes
    integer, parameter :: idepx = 1, idepy = 2, idepz = 3, iidic = 4, idepyp = 5, idepzp = 6
    integer, parameter :: ifx = 7, ify = 8, ifz = 9, icalc = 10
!   État du discret : adhérent, glissant, décollé
    integer, parameter :: EtatAdher = 0, EtatGliss = 1, EtatDecol = 2
    integer, parameter :: EnVitesse = 1, EnPlasticite = 2
!
    integer      :: ii
    real(kind=8) :: xl(6), xd(3), auxil(6), raide(6), rignor, rigtan, depxyz(3), vitxyz(3)
    real(kind=8) :: coulom, dist12, psca, vit123(3)
    real(kind=8) :: vitt, fort
!
    integer         :: axes(3), ContactInGlobal, TestOK, TestNOK
    real(kind=8)    :: ldp(3), ldm(3), SigneAxe(3)
    real(kind=8)    :: ldpglob(3), forceglob(3), raideglob(6)
    aster_logical   :: IsEnfonce
!
!   ContactInGlobal :
!       0 : discret normal
!       1 : discret repère global
!       2 : discret repère global en diagonale
    integer, parameter :: ReperLocal = 0, ReperGlobal = 1, ReperBizarre = 2
!
    character(len=32) :: messak(3)
!
    data nomre1/'RIGI_NOR', 'RIGI_TAN', 'AMOR_NOR', 'AMOR_TAN', 'COULOMB', &
        'DIST_1', 'DIST_2', 'JEU', 'CONTACT'/
! ----------------------------------------------------------------------
!
!   Définition des parametres
    xl(:) = 0.d0; auxil(:) = 0.d0; xd(:) = 0.d0
!   Coordonnees dans le repere local
    if (DD%ndim .eq. 3) then
        call utpvgl(DD%nno, 3, DD%pgl, xg, xl)
    else
        call ut2vgl(DD%nno, 2, DD%pgl, xg, xl)
    end if
!   Raideurs du discret
!       ==> Elles sont surchargées par celles du matériau
    call diraidklv(DD%nomte, raide, klv)
    !
    valre1(:) = 0.0
    valre1(1) = raide(1)
!   Caractéristiques du matériau
    call rcvala(icodma, ' ', 'DIS_CONTACT', 0, ' ', [0.0d0], &
                nbre1, nomre1, valre1, codre1, 0, nan='NON')
    rignor = valre1(1)
    rigtan = valre1(2)
    coulom = valre1(5)
    ContactInGlobal = nint(valre1(9))
    if (.not. in_liste_entier(ContactInGlobal, [ReperLocal, ReperGlobal])) then
        messak(1) = 'DIS_CONTACT'
        messak(2) = 'DIS_CHOC'
        messak(3) = '"1D"|"COIN_2D"'
        call utmess('F', 'DISCRETS_35', nk=3, valk=messak)
    end if
    !
    ! Si ContactInGlobal [ReperGlobal, ReperBizarre]
    !   Prise en compte de          : RIGI_NOR, REPERE
    !   Pas de prise en compte de   : DIST_1, DIST_2, JEU, RIGI_TAN, AMOR_NOR, AMOR_TAN, COULOMB
    !       OK    =0  [ 1,                      9 ]
    !       NOOK  =1  [    2, 3, 4, 5, 6, 7, 8,   ]
    if (in_liste_entier(ContactInGlobal, [ReperGlobal, ReperBizarre])) then
        TestOK = codre1(1)+codre1(9)
        TestNOK = codre1(2)+codre1(3)+codre1(4)+codre1(5)+codre1(6)+codre1(7)+codre1(8)
        if (TestOK .ne. 0 .or. TestNOK .ne. 7) then
            messak(1) = 'DIS_CONTACT'
            messak(2) = 'DIS_CHOC'
            messak(3) = '"COIN_2D"'
            call utmess('F', 'DISCRETS_36', nk=3, valk=messak)
        end if
    end if
!
    varpl(icalc) = EnVitesse
!   Élément avec 2 noeuds
    if (DD%nno .eq. 2) then
        dist12 = valre1(6)+valre1(7)
        ! Vitesse tangente
        vit123(:) = 0.0
        vit123(2) = dvl(2+DD%nc)-dvl(2)
        if (DD%ndim .eq. 3) then
            vit123(3) = dvl(3+DD%nc)-dvl(3)
        end if
        ! ------------------------------------------------------------------------------------------
        ! Détermination du plan du discret : géométrie initiale
        ldm(1:3) = xg(4:6)-xg(1:3); axes(:) = 0; SigneAxe(:) = 0.0
        if (in_liste_entier(ContactInGlobal, [ReperGlobal, ReperBizarre])) then
            ! Plan du discret     : [ axes(1), axes(2) ]
            ! Axe perpendiculaire : axes(3)
            if (abs(ldm(1)) <= r8prem()) then
                ! Plan YZ, vect ↑ X
                axes(:) = [2, 3, 1]
            else if (abs(ldm(2)) <= r8prem()) then
                ! Plan XZ, vect ↑ Y
                axes(:) = [1, 3, 2]
            else if (abs(ldm(3)) <= r8prem()) then
                ! Plan XY, vect ↑ Z
                axes(:) = [1, 2, 3]
            else
                ! <F> Le discret n'est pas plan
                messak(1) = 'DIS_CONTACT'
                messak(2) = 'DIS_CHOC (cas symétrique)'
                call utmess('F', 'DISCRETS_33', nk=2, valk=messak)
            end if
            if (abs(ldm(axes(1))) <= r8prem() .or. &
                abs(ldm(axes(2))) <= r8prem()) then
                ! <F> Le discret est suivant un axe
                messak(1) = 'DIS_CONTACT'
                messak(2) = 'DIS_CHOC (cas symétrique)'
                call utmess('F', 'DISCRETS_34', nk=2, valk=messak)
            end if
            do ii = 1, 3
                SigneAxe(axes(ii)) = LeSigne(ldm(axes(ii)))
            end do
        end if
        ! Instant '-' : Longueur du discret, repère GLOBAL
        ldm(1:3) = ldm(1:3)+DD%ugm(4:6)-DD%ugm(1:3)
        ! Instant '+' : Longueur du discret, repère GLOBAL
        ldp(1:3) = ldm(1:3)+DD%dug(4:6)-DD%dug(1:3)
        ! ------------------------------------------------------------------------------------------
        ! Longueur du discret, repère LOCAL
        do ii = 1, DD%ndim
            xd(ii) = xl(DD%ndim+ii)-xl(ii)
        end do
        ! Déplacement d'entrainement
        call dcopy(DD%ndim, dpe(1), 1, auxil, 1)
        call dcopy(DD%ndim, dpe(1+DD%nc), 1, auxil(4), 1)
        depxyz(:) = 0.0
        depxyz(1) = (xd(1)+ulp(1+DD%nc)-ulp(1)+auxil(4)-auxil(1)-dist12)-r8prem()
        depxyz(2) = (xd(2)+ulp(2+DD%nc)-ulp(2)+auxil(5)-auxil(2))
        if (DD%ndim .eq. 3) then
            depxyz(3) = xd(3)+ulp(3+DD%nc)-ulp(3)+auxil(6)-auxil(3)
        end if
        ! Vitesse d'entrainement
        call dcopy(DD%ndim, dve(1), 1, auxil, 1)
        call dcopy(DD%ndim, dve(1+DD%nc), 1, auxil(4), 1)
        ! Vitesse tangente
        vitxyz(:) = 0.0
        vitxyz(2) = vit123(2)+auxil(5)-auxil(2)
        if (DD%ndim .eq. 3) then
            vitxyz(3) = vit123(3)+auxil(6)-auxil(3)
        end if
        ! ------------------------------------------------------------------------------------------
        force(1:3) = 0.0; forceglob(:) = 0.0; raideglob(:) = 0.0
        IsEnfonce = ASTER_FALSE
        if (in_liste_entier(ContactInGlobal, [ReperGlobal, ReperBizarre])) then
            if ((ldp(axes(1))*SigneAxe(axes(1)) <= 0.0) .and. &
                (ldp(axes(2))*SigneAxe(axes(2)) <= 0.0)) then
                IsEnfonce = ASTER_TRUE
                ! Les 2 ldp sont <= 0.0
                !   On garde le plus petit pour le calcul de l'effort
                !   Pas de déplacement, ni d'effort dans le plan perpendiculaire
                ldpglob(:) = 0.0
                if (abs(ldp(axes(1))) <= abs(ldp(axes(2)))) then
                    ldpglob(axes(1)) = ldp(axes(1))
                    ldpglob(axes(2)) = abs(ldpglob(axes(1)))*LeSigne(ldp(axes(2)))
                    forceglob(axes(1)) = rignor*ldpglob(axes(1))
                    raideglob(axes(1)) = rignor
                else
                    ldpglob(axes(1)) = abs(ldp(axes(2)))*LeSigne(ldp(axes(1)))
                    ldpglob(axes(2)) = ldp(axes(2))
                    forceglob(axes(2)) = rignor*ldpglob(axes(2))
                    raideglob(axes(2)) = rignor
                end if
                call utpvgl(1, 3, DD%pgl, ldpglob, depxyz)
                if (ContactInGlobal == ReperGlobal) then
                    call diklvraid(DD%nomte, kgv, raideglob)
                    ! Matrice de raideur dans le repère local
                    ! Force dans le repère local
                    if (DD%ndim .eq. 3) then
                        call utpsgl(DD%nno, DD%nc, DD%pgl, kgv, klv)
                        call utpvgl(DD%nno, DD%nc, DD%pgl, &
                                    forceglob, force)
                    else
                        call ut2mgl(DD%nno, DD%nc, DD%pgl, kgv, klv)
                        call ut2vgl(DD%nno, DD%nc, DD%pgl, &
                                    forceglob, force)
                    end if
                end if
            end if
        else
            if (depxyz(1) <= 0.0) IsEnfonce = ASTER_TRUE
        end if
        if (IsEnfonce) then
            if (ContactInGlobal == ReperGlobal) then
                ! Plus de calcul à faire si fait dans le repère global
                varpl(ifx) = force(1)
                varpl(ify) = force(2)
                varpl(ifz) = force(3)
                varpl(iidic) = EtatAdher
            else
                force(1) = rignor*depxyz(1)
                if (force(1) .gt. 0.0) force(1) = 0.0
                psca = varmo(ify)*vitxyz(2)+varmo(ifz)*vitxyz(3)
                if ((psca .ge. 0.0) .and. (nint(varmo(iidic)) .eq. EtatGliss)) then
                    vitt = (vitxyz(2)**2+vitxyz(3)**2)**0.5d0
                    force(2) = 0.0
                    force(3) = 0.0
                    if (vitt .gt. 0.00) then
                        force(2) = -coulom*force(1)*vitxyz(2)/vitt
                        force(3) = -coulom*force(1)*vitxyz(3)/vitt
                    end if
                    varpl(iidic) = EtatGliss
                else
                    force(2) = rigtan*(depxyz(2)-varmo(idepy))+varmo(ify)
                    force(3) = rigtan*(depxyz(3)-varmo(idepz))+varmo(ifz)
                    varpl(iidic) = EtatAdher
                    fort = (force(2)**2+force(3)**2)**0.5d0
                    if (fort .gt. abs(coulom*force(1))) then
                        vitt = (vitxyz(2)**2+vitxyz(3)**2)**0.5d0
                        force(2) = 0.0
                        force(3) = 0.0
                        if (vitt .gt. 0.0) then
                            force(2) = -coulom*force(1)*vitxyz(2)/vitt
                            force(3) = -coulom*force(1)*vitxyz(3)/vitt
                            varpl(iidic) = EtatGliss
                        end if
                    end if
                end if
                varpl(ifx) = force(1)
                varpl(ify) = force(2)
                varpl(ifz) = force(3)
                !
                if (abs(rigtan) .gt. r8prem()) then
                    varpl(idepyp) = depxyz(2)-varmo(ify)/rigtan
                    varpl(idepzp) = depxyz(3)-varmo(ifz)/rigtan
                end if
                !
                force(2) = force(2)+raide(2)*(ulp(2+DD%nc)-ulp(2))
                if (DD%ndim .eq. 3) then
                    force(3) = force(3)+raide(3)*(ulp(3+DD%nc)-ulp(3))
                end if
                ! Actualisation de la matrice de raideur
                raide(1) = rignor
                call diklvraid(DD%nomte, klv, raide)
            end if
        else
            varpl(ifx) = 0.0
            varpl(ify) = 0.0
            varpl(ifz) = 0.0
            varpl(iidic) = EtatDecol
            varpl(idepyp) = 0.0
            varpl(idepzp) = 0.0
            klv(1:78) = 0.0
        end if
!
!   Élément avec 1 noeud
    else
        if (ContactInGlobal .ne. ReperLocal) then
            messak(1) = 'DIS_CONTACT'
            messak(2) = 'DIS_CHOC (cas 1 noeud)'
            messak(3) = '"LOCAL"'
            call utmess('F', 'DISCRETS_35', nk=3, valk=messak)
        end if
        dist12 = valre1(8)-valre1(6)
        ! Vitesse tangente
        vit123(:) = 0.0
        vit123(2) = dvl(2)
        if (DD%ndim .eq. 3) then
            vit123(3) = dvl(3)
        end if
        call dcopy(DD%ndim, dpe(1), 1, auxil, 1)
        depxyz(:) = 0.0
        depxyz(1) = ulp(1)+dist12+auxil(1)
        depxyz(2) = ulp(2)+auxil(2)
        if (DD%ndim .eq. 3) then
            depxyz(3) = ulp(3)+auxil(3)
        end if
        call dcopy(DD%ndim, dve(1), 1, auxil, 1)
        ! Vitesse tangente
        vitxyz(:) = 0.0
        vitxyz(2) = vit123(2)+auxil(2)
        if (DD%ndim .eq. 3) then
            vitxyz(3) = vit123(3)+auxil(3)
        end if
        force(1:3) = 0.0
        if (depxyz(1) .le. 0.0d0) then
            force(1) = rignor*depxyz(1)
            if (force(1) .gt. 0.0) force(1) = 0.0
            psca = varmo(ify)*vitxyz(2)+varmo(ifz)*vitxyz(3)
            if ((psca .ge. 0.0) .and. (nint(varmo(iidic)) .eq. EtatGliss)) then
                vitt = (vitxyz(2)**2+vitxyz(3)**2)**0.5d0
                force(2) = 0.0
                force(3) = 0.0
                if (vitt .gt. 0.0) then
                    force(2) = -coulom*force(1)*vitxyz(2)/vitt
                    force(3) = -coulom*force(1)*vitxyz(3)/vitt
                end if
                varpl(iidic) = EtatGliss
            else
                force(2) = rigtan*(depxyz(2)-varmo(idepy))+varmo(ify)
                force(3) = rigtan*(depxyz(3)-varmo(idepz))+varmo(ifz)
                varpl(iidic) = EtatAdher
                fort = (force(2)**2+force(3)**2)**0.5d0
                if (fort .gt. abs(coulom*force(1))) then
                    vitt = (vitxyz(2)**2+vitxyz(3)**2)**0.5d0
                    force(2) = 0.0
                    force(3) = 0.0
                    if (vitt .gt. 0.0) then
                        force(2) = -coulom*force(1)*vitxyz(2)/vitt
                        force(3) = -coulom*force(1)*vitxyz(3)/vitt
                        varpl(iidic) = EtatGliss
                    end if
                end if
            end if
            varpl(ifx) = force(1)
            varpl(ify) = force(2)
            varpl(ifz) = force(3)
            !
            if (abs(rigtan) .gt. r8prem()) then
                varpl(idepyp) = depxyz(2)-varmo(ify)/rigtan
                varpl(idepzp) = depxyz(3)-varmo(ifz)/rigtan
            end if
            !
            force(2) = force(2)+raide(2)*ulp(2)
            if (DD%ndim .eq. 3) then
                force(3) = force(3)+raide(3)*ulp(3)
            end if
            ! Actualisation de la matrice de raideur
            raide(1) = rignor
            call diklvraid(DD%nomte, klv, raide)
        else
            varpl(ifx) = 0.0
            varpl(ify) = 0.0
            varpl(ifz) = 0.0
            varpl(iidic) = EtatDecol
            varpl(idepyp) = 0.0
            varpl(idepzp) = 0.0
            klv(1:78) = 0.0
        end if
    end if
    varpl(idepx) = depxyz(1)
    varpl(idepy) = depxyz(2)
    varpl(idepz) = depxyz(3)
!
!   Prédiction en dynamique, on retourne les efforts précédents
!       Si on passe d'une formulation 'plastique' à une en 'vitesse'
!       On le fait à la fin, la raideur doit être mise comme il faut
    if (Predic .and. (nint(varmo(icalc)) .eq. EnPlasticite)) then
        ! Les efforts précédents
        force(1) = varmo(ifx)
        force(2) = varmo(ify)
        force(3) = varmo(ifz)
        ! On remet les varpl comme il faut. Elles ont peut-être été modifiées
        varpl(ifx) = varmo(ifx)
        varpl(ify) = varmo(ify)
        varpl(ifz) = varmo(ifz)
    end if
!
contains

    function LeSigne(xx)
        implicit none
#include "asterc/r8prem.h"
        real(kind=8) :: LeSigne
        real(kind=8), intent(in) :: xx
!
        if (xx >= r8prem()) then
            LeSigne = 1.0
        else if (xx <= -r8prem()) then
            LeSigne = -1.0
        else
            LeSigne = 0.0
        end if
    end function LeSigne
end subroutine
