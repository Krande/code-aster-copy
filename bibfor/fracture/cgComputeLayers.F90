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
subroutine cgComputeLayers(cgField, cgTheta, cnstet, v_numc)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"

#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cncinv.h"
#include "asterfort/jexatr.h"
#include "asterfort/jenuno.h"
#include "asterfort/conare.h"
#include "asterfort/char8_to_int.h"

#include "jeveux.h"
!
    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
    character(len=24), intent(in) :: cnstet
    real(kind=8), intent(inout), pointer :: v_numc(:)

    integer(kind=8) :: nbel, i

    integer(kind=8), pointer :: fondNoeudNume(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: v_l_no_cm(:) => null()
    integer(kind=8), pointer :: v_l_no_cp(:) => null()

    real(kind=8), pointer :: v_theta(:) => null()
    character(len=8), pointer   :: fondNoeud(:) => null()

    integer(kind=8) :: nume, num_fr, num_c, num_no, iar, iatyma, ima, na, nbar, ndime
    integer(kind=8) :: ino1, ino2, ino3, i_node
    integer(kind=8) :: nno1, nno2, nno3, numac
    integer(kind=8) :: nno_cm, nno_cp
    integer(kind=8) :: adra, nbmaca, ar(12, 3)
    integer(kind=8) :: ibid, ityp
    integer(kind=8) :: jdrvlc, jcncin, jconx2
    character(len=24) :: cnxinv
    character(len=8) :: type

!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute Layers Field in 2D
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
!    Dans cette routine, calcul du numero de couche
!
! --------------------------------------------------------------------------------------------------
!

    call jeveuo(cnstet(1:19)//'.CNSV', 'E', vr=v_theta)

    !Recuperation de la connectivite

    call jeveuo(cgTheta%mesh//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(cgTheta%mesh//'.CONNEX', 'LONCUM'), 'L', jconx2)

    call jeveuo(cgTheta%mesh//'.TYPMAIL', 'L', iatyma)

    !Recuperation de la connectivite inverse
    cnxinv = '&&'//cgTheta%mesh//'.CNXINV'
    ibid = 0
    call cncinv(cgTheta%mesh, [ibid], 0, 'V', cnxinv)
    call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', jdrvlc)
    call jeveuo(jexnum(cnxinv, 1), 'L', jcncin)

    !Nombre elements
    call dismoi('NB_NO_MAILLA', cgTheta%mesh, 'MAILLAGE', repi=nbel)

    AS_ALLOCATE(vi=v_l_no_cp, size=nbel)
    AS_ALLOCATE(vi=v_l_no_cm, size=nbel)

    !Initialisation de tous les numeros de couche a nb_couche_sup
    do i = 1, nbel
        v_numc(i) = cgTheta%nb_couche_sup
        v_l_no_cp(i) = 0
        v_l_no_cm(i) = 0
    end do

    !Recuperation des noeuds du maillage du front
    call jeveuo(cgTheta%crack//'.FOND.NOEU', 'L', vk8=fondNoeud)
    AS_ALLOCATE(vi=fondNoeudNume, size=cgTheta%nb_fondNoeud)
    do i_node = 1, cgTheta%nb_fondNoeud
        nume = char8_to_int(fondNoeud(i_node))
        fondNoeudNume(i_node) = nume
    end do

    !Boucle sur les noeuds du front
    do num_fr = 1, cgTheta%nb_fondNoeud

        !Initialisation du numero de couche pour le noeud du front
        nno_cm = 1
        v_l_no_cm(1) = fondNoeudNume(num_fr)
        v_numc(fondNoeudNume(num_fr)) = 0

        !Boucle sur les couches
        do num_c = 1, cgTheta%nb_couche_sup

            !Nombre de noeuds de la couche suivante : initialisation
            nno_cp = 0

            !Boucle sur les noeuds de la couche precedente
            do num_no = 1, nno_cm

                !Noeud courant
                na = v_l_no_cm(num_no)

                ! Mailles connectees au noeud courant
                adra = zi(jdrvlc-1+na)
                ! Nombre de mailles connectees
                nbmaca = zi(jdrvlc-1+na+1)-zi(jdrvlc-1+na)

                do ima = 1, nbmaca
                    !NUMERO DE LA MAILLE
                    numac = zi(jcncin-1+adra+ima-1)

                    ityp = iatyma-1+numac
                    call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
                    call dismoi('DIM_TOPO', type, 'TYPE_MAILLE', repi=ndime)

                    !       On ne traite pas les mailles de bord
                    if (ndime .eq. cgField%ndim) then
                        !
                        ! Connectivite des arretes de la maille
                        call conare(type, ar, nbar)

                        !
                        !       BOUCLE SUR LE NOMBRE D'ARETES DE LA MAILLE NUMAC
                        do iar = 1, nbar

                            !Noeuds extremites de l'arrete
                            ino1 = ar(iar, 1)
                            nno1 = connex(zi(jconx2+numac-1)+ino1-1)
                            ino2 = ar(iar, 2)
                            nno2 = connex(zi(jconx2+numac-1)+ino2-1)

                            if (na .ne. nno2) then
                                if (num_c .lt. v_numc(nno2)) then
                                    v_numc(nno2) = num_c
                                    nno_cp = nno_cp+1
                                    v_l_no_cp(nno_cp) = nno2
                                end if
                            end if
                            if (na .ne. nno1) then
                                if (num_c .lt. v_numc(nno1)) then
                                    v_numc(nno1) = num_c
                                    nno_cp = nno_cp+1
                                    v_l_no_cp(nno_cp) = nno1
                                end if
                            end if

                            !Eventuels noeuds milieux
                            if (ar(iar, 3) .ne. 0) then
                                ino3 = ar(iar, 3)
                                nno3 = connex(zi(jconx2+numac-1)+ino3-1)
                                if (abs(v_numc(nno1)-v_numc(nno2)) .le. 1e-8) then
                                    v_numc(nno3) = v_numc(nno1)
                                else
                                    v_numc(nno3) = min(num_c-0.5, v_numc(nno3))
                                end if
                            end if

                        end do

                    end if

                end do

            end do

            if (nno_cp .gt. 0) then
                nno_cm = nno_cp
                v_l_no_cm = v_l_no_cp
            else
                nno_cm = 0
            end if

        end do

    end do

    AS_DEALLOCATE(vi=v_l_no_cp)
    AS_DEALLOCATE(vi=v_l_no_cm)
    AS_DEALLOCATE(vi=fondNoeudNume)

end subroutine
