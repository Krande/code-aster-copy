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
subroutine rco3d_clcrela(ligrel, noma, nb_pairs, nbnocot, &
                         list_total_no_co, map_noco_pair, map_noco_nbelem, &
                         map_noco_nbnoco, resuelem, fonrez, lisrel)
    !
    use raco3d_module
    !
    implicit none
    !
#include "jeveux.h"
#include "asterfort/afrela.h"
#include "asterfort/assert.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/digdel.h"
#include "asterfort/imprel.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/int_to_char8.h"

    character(len=19), intent(in) :: ligrel, resuelem, lisrel
    character(len=8), intent(in) :: noma
    integer(kind=8), intent(in) :: nb_pairs, nbnocot
    integer(kind=8), intent(in) :: map_noco_pair(:, :, :)
    integer(kind=8), intent(in) :: map_noco_nbnoco(:, :, :)
    integer(kind=8), intent(in) :: map_noco_nbelem(:, :)
    integer(kind=8), pointer, intent(in) :: list_total_no_co(:)
    character(len=*), intent(in) :: fonrez

!  SUBROUTINE: rco3d_addrela
!
!  DESCRIPTION:
!  This subroutine adds the linear relations imposed by the 3D-SHELL CONNECTION.
!
!  INPUT PARAMETERS:
!  -----------------
!  ligrel           - IN    - K19    - Name of the LIGREL (3D-SHELL CONNECTION)
!
!  noma             - IN    - K8     - Name of the mesh
!
!  nb_pairs         - IN    - I      - Total number of element pairs involved in the connection.
!
!  nbnocot          - IN    - I      - Total number of shell nodes.
!
!  list_total_no_co - IN    - PTR    - Array containing the indices of shell nodes.
!                                     Dimensions: (nbnocot)
!
!  map_noco_pair    - IN   - I(*)    - Mapping of each shell node to the pairs
!                                      that contain the node.
!                                      Dimensions: (9, nbnocot, nb_pairs)
!
!  map_noco_nbnoco  - IN   - I(*)    - Mapping of each shell node to the number
!                                       of shell nodes in the pair.
!                                     Dimensions: (9, nbnocot, nb_pairs)
!
!  map_noco_nbelem  - IN   - I(*)    - Mapping of each shell node to the total number of pairs that
!                                     contain the node.
!                                     Dimensions: (9, nbnocot)
!
!  fonrez           - IN    - R      - 'REEL'
!
!  resuelem         - IN    - K19    - Name of the elementary matrices (`resuelem`) used to compute
!                                    the linear realtion .
!
!  lisrel           - IN    - K19    - Name of the SD LISTE_RELA structure, which holds
!                                     the list of linear relations.
!
! ------------------------------------------------------------------------------

    complex(kind=8) :: betac
    real(kind=8) :: beta
    character(len=24) :: noeuma
    character(len=16) :: motfac
    character(len=8) :: betaf, dofs(6), nomnoe
    integer(kind=8) :: nbterm, i, j, k, l
    character(len=4) :: typval, typcoe
    character(len=8), pointer :: lisddl(:) => null()
    character(len=8), pointer :: lisno(:) => null()
    integer(kind=8), pointer :: repe_type(:) => null()
    real(kind=8), pointer :: repe_defi(:) => null()
    complex(kind=8), pointer :: coec(:) => null()
    real(kind=8), pointer :: coer(:) => null()
    integer(kind=8), pointer :: v_desc(:) => null()
    integer(kind=8), pointer :: v_list_no_pair(:) => null()
    integer(kind=8) :: iret, nb_gr, ncomp, mode
    integer(kind=8) :: jv_liel, num_pair, nno
    type(PointerContainer), allocatable :: grel_ptr(:), resu_ptr(:)
    integer(kind=8) :: nnco, nn3d, idx, nddl, row_index, ico, i3d
    integer(kind=8) :: dofco, dof3d, idof, jdof
    aster_logical :: found, check, lcolle

    call jemarq()
    ! Fill the linear relations

    ! --- VECTEUR DU NOM DES NOEUDS
    AS_ALLOCATE(vk8=lisno, size=NB_NDDL_MAX*nb_pairs)
    ! --- VECTEUR DU NOM DES DDLS
    AS_ALLOCATE(vk8=lisddl, size=NB_NDDL_MAX*nb_pairs)
    ! --- VECTEUR DES COEFFICIENTS REELS
    AS_ALLOCATE(vr=coer, size=NB_NDDL_MAX*nb_pairs)
    ! --- VECTEUR DES COEFFICIENTS COMPLEXES
    AS_ALLOCATE(vc=coec, size=NB_NDDL_MAX*nb_pairs)
    ! --- VECTEUR DES DIRECTIONS DES DDLS A CONTRAINDRE
    AS_ALLOCATE(vr=repe_defi, size=3*NB_NDDL_MAX*nb_pairs)
    ! --- VECTEUR DES DIMENSIONS DE CES DIRECTIONS
    AS_ALLOCATE(vi=repe_type, size=NB_NDDL_MAX*nb_pairs)

    ! --- TYPE DES VALEURS AU SECOND MEMBRE DES RELATIONS
    typval = fonrez
    ! --- TYPE DES VALEURS DES COEFFICIENTS DES RELATIONS
    typcoe = 'REEL'
    ! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UNE FONCTION
    betaf = '&FOZERO'
    ! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UN REEL
    beta = 0.0d0
    ! --- VALEUR DU SECOND MEMBRE DES RELATIONS QUAND C'EST UN COMPLEXE
    betac = (0.0d0, 0.0d0)
    !
    coer = 0.0d0
    ! Dofs names
    dofs(1) = 'DX'
    dofs(2) = 'DY'
    dofs(3) = 'DZ'
    dofs(4) = 'DRX'
    dofs(5) = 'DRY'
    dofs(6) = 'DRZ'
    !
    noeuma = noma//'.NOMNOE'
    !
    motfac = 'LIAISON_ELEM'
    !
    !
    call jeexin(resuelem//'.DESC', iret)
    !
    call jeveuo(resuelem//'.DESC', 'L', vi=v_desc)
    nb_gr = v_desc(2)
    !
    allocate (grel_ptr(nb_gr))
    allocate (resu_ptr(nb_gr))
    !
    nbterm = 0
    dofco = 6
    dof3d = 3
    coer = 0.0d0
    lisno = ''
    lisddl = ''
    lcolle = .false.
    !
    do i = 1, nb_gr
        call jeveuo(jexnum(ligrel(1:19)//'.LIEL', i), 'L', vi=grel_ptr(i)%iptr)
        call jeveuo(jexnum(resuelem//'.RESL', i), 'L', vr=resu_ptr(i)%rptr)
    end do

    ! Assembly
    do j = 1, nbnocot
        do idof = 1, 6
            do i = 1, nb_gr
                mode = v_desc(2+i)
                ncomp = digdel(mode)
                do k = 1, map_noco_nbelem(i, j)
                    nnco = map_noco_nbnoco(i, j, k)
                    ! jv_liel indice de l element dans le grel
                    jv_liel = map_noco_pair(i, j, k)
                    num_pair = -grel_ptr(i)%iptr(jv_liel)
                    call jelira(jexnum(ligrel//'.NEMA', num_pair), 'LONMAX', nno)
                    nn3d = nno-nnco-1
                    nddl = 6*nnco+3*nn3d
                    call jeveuo(jexnum(ligrel//'.NEMA', num_pair), 'L', vi=v_list_no_pair)
                    check = .false.
                    do l = 1, nnco
                        if (v_list_no_pair(l) .eq. list_total_no_co(j)) then
                            idx = l
                            check = .true.
                            exit
                        end if
                    end do
                    if (.not. check) then
                        ASSERT(.false.)
                    end if

                    !
                    row_index = (jv_liel-1)*ncomp+6*(idx-1)*nddl+(idof-1)*nddl
                    do ico = 1, nnco
                        !call jenuno(jexnum(noeuma, v_list_no_pair(ico)), nomnoe)
                        nomnoe = int_to_char8(v_list_no_pair(ico), lcolle, noma, 'NOEUD')
                        do jdof = 1, dofco
                            found = .false.
                            do l = 1, nbterm
                                if ((lisno(l) .eq. nomnoe) .and. (lisddl(l) .eq. dofs(jdof))) then
                                    coer(l) = coer(l)+ &
                                              resu_ptr(i)%rptr(row_index+6*(ico-1)+jdof)
                                    found = .true.
                                    exit
                                end if
                            end do
                            if (.not. found) then
                                nbterm = nbterm+1
                                lisno(nbterm) = nomnoe
                                lisddl(nbterm) = dofs(jdof)
                                coer(nbterm) = resu_ptr(i)%rptr(row_index+6*(ico-1)+jdof)
                            end if
                        end do
                    end do
                    do i3d = nnco+1, nno-1
                        !call jenuno(jexnum(noeuma, v_list_no_pair(i3d)), nomnoe)
                        nomnoe = int_to_char8(v_list_no_pair(i3d), lcolle, noma, 'NOEUD')
                        do jdof = 1, dof3d
                            found = .false.
                            do l = 1, nbterm
                                if ((lisno(l) .eq. nomnoe) .and. (lisddl(l) .eq. dofs(jdof))) then
                                    coer(l) = coer(l) &
                                              +resu_ptr(i)%rptr( &
                                              row_index+6*nnco+3*(i3d-nnco-1)+jdof)
                                    found = .true.
                                    exit
                                end if
                            end do
                            if (.not. found) then
                                nbterm = nbterm+1
                                lisno(nbterm) = nomnoe
                                lisddl(nbterm) = dofs(jdof)
                                coer(nbterm) = &
                                    resu_ptr(i)%rptr(row_index+6*nnco+3*(i3d-nnco-1)+jdof)
                            end if
                        end do
                    end do
                end do
            end do

            call afrela(coer, coec, lisddl, lisno, repe_type, &
                        repe_defi, nbterm, beta, betac, betaf, &
                        typcoe, typval, 0.d0, lisrel)
            call imprel(motfac, nbterm, coer, lisddl, lisno, &
                        beta, 1.0e-16)
            nbterm = 0
            coer = 0.0d0
            lisno = ''
            lisddl = ''
        end do
    end do

    deallocate (grel_ptr)
    deallocate (resu_ptr)
    AS_DEALLOCATE(vk8=lisno)
    AS_DEALLOCATE(vk8=lisddl)
    AS_DEALLOCATE(vr=coer)
    AS_DEALLOCATE(vc=coec)
    AS_DEALLOCATE(vr=repe_defi)
    AS_DEALLOCATE(vi=repe_type)

    call jedema()

end subroutine
