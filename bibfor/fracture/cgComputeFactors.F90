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
! person_in_charge: nicolas.pignet at edf.fr
!
subroutine cgComputeFactors(cgField, cgTheta, cgStat)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/imprsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cgComputeLayers.h"
#include "jeveux.h"
!
    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
    type(CalcG_stat), intent(inout) :: cgStat

    integer(kind=8) :: i
    integer(kind=8) :: nbel, iret, jcnsl

    real(kind=8) :: d, xm, ym, zm, xn, yn, zn, eps, alpha
    real(kind=8) :: valpar(1), valres_i, valres_s

    character(len=8), parameter :: licmp(6) = ['MODULE  ', 'DIR_X   ', 'DIR_Y   ', 'DIR_Z   ', &
                                               'ABSC_CUR', 'LONG    ']
    character(len=8) :: nompar(1)
    character(len=24) :: cnstet
    real(kind=8) :: theta0, start, finish
    real(kind=8), pointer :: v_theta(:) => null()
    real(kind=8), pointer :: v_coor(:) => null()
    real(kind=8), pointer :: v_base(:) => null()
    real(kind=8), pointer :: v_absc(:) => null()
    real(kind=8), pointer :: v_numc(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute Theta Field in 2D and 3D
!
!    Contenu de theta: structure de données stockant les informations
!            permettant de calculer le champ theta aux points de Gauss
!            d'un élémentdans les te.
!            Le champ stocké dans jeveux sous le nom
!            cgTheta%theta_factors//"_CHAM_THETA_FACT".
!
!            Les composantes de ce champ sont, pour chaque noeud :
!               - 1: la valeur de la fonction theta0(r) pour ce noeud,
!                    où r est la distance entre le noeud et son projeté
!                    sur le fond de fissure.
!               - 2 à 4: le vecteur t donnant la direction et le sens de
!                        theta pour ce noeud, associé au projeté du noeud
!                        sur le fond de fissure.
!               - 5: l'abscisse curviligne du projeté du noeud sur le
!                    fond de fissure. (3D uniquement)
!               - 6: longueur de la fissure. (3D uniquement)
!
!    Remplissage de la variable cgTheta%nb_theta_field
!
! --------------------------------------------------------------------------------------------------
!
    call cpu_time(start)
!
    call jemarq()
!
    eps = 1.d-06
!
    call cgTheta%getCoorNodes(v_coor)
    call cgTheta%getBaseLoc(v_base)
    call cgTheta%getAbscurv(v_absc)
!
    call dismoi('NB_NO_MAILLA', cgTheta%mesh, 'MAILLAGE', repi=nbel)
!
    if (cgField%level_info > 1) then
        call utmess('I', 'RUPTURE3_2')
    end if

    ! Vérifications spécifiques en 3D
    if (cgField%ndim .eq. 3) then
        ! VERIFICATION PRESENCE NB_POINT_FOND
        if (cgTheta%nb_point_fond .ne. 0) then
            ! INTERDICTION D AVOIR NB_POINT_FOND AVEC
            ! DISCTRETISATION =  LEGENDRE
            if (cgTheta%discretization .eq. 'LEGENDRE') then
                call utmess('F', 'RUPTURE1_73')
            end if
        end if
    end if

!
! allocation de la structure de données temporaire de theta

!   N.B.: le champ theta contient seulement les ndimte champs utiles
!     * stockage des champs simples utilisés pour le remplissages
!
    cnstet = '&&CNSTET_CHAM'
    call cnscre(cgTheta%mesh, 'THET_R', 6, licmp, 'V', cnstet)
    call jeveuo(cnstet(1:19)//'.CNSL', 'E', jcnsl)
    call jeveuo(cnstet(1:19)//'.CNSV', 'E', vr=v_theta)

!   CALCUL DES NUMEROS DE COUCHE EN FONCTION DE LA CONNECTIVITE
    if (cgTheta%radius_type .eq. 'NB_COUCHE') then
        AS_ALLOCATE(vr=v_numc, size=nbel)
        call cgComputeLayers(cgField, cgTheta, cnstet, v_numc)
    end if

!   BOUCLE SUR LES NOEUDS M COURANTS DU MAILLAGE
    do i = 1, nbel

        if (cgTheta%radius_type .ne. 'NB_COUCHE') then
            !       CALCUL DE LA FONCTION DETERMINANT LA NORME DE THETA EN
            !       FONTION DE R_INF, R_SUP ET LA DISTANCE DU NOEUD AU FRONT
            !
            !       COORDONNEES DU NOEUD COURANT M
            xm = v_coor((i-1)*3+1)
            ym = v_coor((i-1)*3+2)
            zm = 0.d0
            !
            if (cgField%ndim .eq. 3) then
                zm = v_coor((i-1)*3+3)
            end if
            !
            !       COORDONNEES DU PROJETE N DE CE NOEUD SUR LE FRONT DE FISSURE
            xn = v_base(3*cgField%ndim*(i-1)+1)
            yn = v_base(3*cgField%ndim*(i-1)+2)
            zn = 0.d0
            if (cgField%ndim .eq. 3) then
                zn = v_base(3*cgField%ndim*(i-1)+3)
            end if
            d = sqrt((xn-xm)*(xn-xm)+(yn-ym)*(yn-ym)+(zn-zm)*(zn-zm))
            if ((cgTheta%radius_type .eq. 'R') .or. (cgTheta%radius_type .eq. 'DEFAUT')) then
                alpha = (d-cgTheta%r_inf)/(cgTheta%r_sup-cgTheta%r_inf)
            else if (cgTheta%radius_type .eq. 'R_FO') then
                nompar(1) = 'ABSC'
                valpar(1) = v_absc(i)
                call fointe('FM', cgTheta%r_inf_fo, 1, nompar, valpar, valres_i, iret)
                call fointe('FM', cgTheta%r_sup_fo, 1, nompar, valpar, valres_s, iret)
                alpha = (d-valres_i)/(valres_s-valres_i)
            else
                ASSERT(.FALSE.)
            end if
        else
            !       CALCUL DE LA FONCTION DETERMINANT LA NORME DE THETA EN
            !       FONTION DE NB_COUCHE_INF, NB_COUCHE_SUP ET LE NUMERO DE COUCHE DU NOEUD
            d = v_numc(i)
            alpha = (d-cgTheta%nb_couche_inf)/(cgTheta%nb_couche_sup-cgTheta%nb_couche_inf)
        end if

!       calcul de theta0 du noeud i
        if ((abs(alpha) .le. eps) .or. (alpha .lt. 0)) then
            theta0 = 1.d0
        else if ((abs(alpha-1) .le. eps) .or. ((alpha-1) .gt. 0)) then
            theta0 = 0.d0
        else
            theta0 = 1.d0-alpha
        end if

!       stockage de la norme de theta, évaluée au noeud i
        zl(jcnsl-1+6*(i-1)+1) = ASTER_TRUE
        v_theta((i-1)*6+1) = theta0

!       stockage de t, la direction de theta, pour le noeud i issue de BASLOC
        zl(jcnsl-1+6*(i-1)+2) = ASTER_TRUE
        zl(jcnsl-1+6*(i-1)+3) = ASTER_TRUE
        zl(jcnsl-1+6*(i-1)+4) = ASTER_TRUE
        v_theta((i-1)*6+2) = v_base(3*cgField%ndim*(i-1)+cgField%ndim+1)
        v_theta((i-1)*6+3) = v_base(3*cgField%ndim*(i-1)+cgField%ndim+2)
        if (cgField%ndim .eq. 2) then
            v_theta((i-1)*6+4) = 0.d0
        else
            v_theta((i-1)*6+4) = v_base(3*cgField%ndim*(i-1)+cgField%ndim+3)
        end if
!
!       stockage de l'abscisse curviligne s, pour le noeud i en 3D
        zl(jcnsl-1+6*(i-1)+5) = ASTER_TRUE
        if (cgField%ndim .eq. 2) then
            v_theta((i-1)*6+5) = 0.d0
        else
            v_theta((i-1)*6+5) = v_absc(i)
        end if
!
!       stockage de la longueur de la fissure, utile en 3D seulement
        zl(jcnsl-1+6*(i-1)+6) = ASTER_TRUE
        if (cgField%ndim .eq. 2) then
            v_theta((i-1)*6+6) = 0.d0
        else
            v_theta((i-1)*6+6) = cgTheta%lonfis
        end if
    end do

    if (cgTheta%radius_type .eq. 'NB_COUCHE') then
        AS_DEALLOCATE(vr=v_numc)
    end if
!
!   ALLOCATION DES OBJETS POUR STOCKER LE VRAI CHAMP_NO THETA_FACTORS
!
    cgTheta%theta_factors = cgTheta%theta_factors(1:8)//'_CHAM_THETA_FACT'
    call cnscno(cnstet, ' ', 'OUI', 'G', cgTheta%theta_factors, 'F', iret)
    call detrsd('CHAM_NO_S', cnstet)
!    call imprsd('CHAMP', cgTheta%theta_factors, 6, ' VECTEUR THETA FACTORS')
!
    call jedema()
!
    call cpu_time(finish)
    cgStat%cgThetaFact = finish-start
end subroutine
