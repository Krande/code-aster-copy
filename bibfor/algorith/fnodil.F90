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
subroutine fnodil(option, typmod, ds_dil, ndim, &
                  nnos, nnom, npg, nddl, dimdef, iw, vff, &
                  vffb, idff, idffb, geomi, compor, &
                  sief, fint)
!
    use dil_type
    use bloc_fe_module, only: add_fint, prod_sb
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmbeps.h"

!
    character(len=8), intent(in)    :: typmod(*)
    character(len=16), intent(in)   :: option, compor(*)
    type(dil_modelisation)         :: ds_dil
    integer(kind=8), intent(in)             :: ndim, nnos, nnom, npg, nddl, dimdef
    integer(kind=8), intent(in)             :: iw, idff, idffb
    real(kind=8), intent(in)        :: geomi(ndim, nnos+nnom)
    real(kind=8), intent(in)        :: vff(nnos+nnom, npg), vffb(nnos, npg)
    real(kind=8), intent(in)        :: sief(dimdef*npg)
    real(kind=8), intent(out)       :: fint(nddl)
!
! ----------------------------------------------------------------------
!     BUT:  CALCUL  DES OPTIONS FORC_NODA
!           EN PETITES DEFORMATIONS
!           POUR SECOND GRADIENT DE DILATATION : XXXX_DIL
! ----------------------------------------------------------------------
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNOS    : NOMBRE DE NOEUDS SOMMETS
! IN  NNOM    : NOMBRE DE NOEUDS MILIEUX
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  NDDL    : NOMBRE DE DEGRES DE LIBERTE D'UN ELEMENT ENRICHI
! IN  IW      : PTR. POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME DE DEPLACEMENT
! IN  VFFB    : VALEUR  DES FONCTIONS DE FORME DE A ET LAMBDA
! IN  IDFF    : PTR. DERIVEE DES FONCTIONS DE FORME DE DEPLACEMENT ELEMENT DE REF.
! IN  IDFFB   : PTR. DERIVEE DES FONCTIONS DE FORME DE A ET LAMBDA ELEMENT DE REF.
! IN  GEOMI   : COORDONNEES DES NOEUDS (CONFIGURATION INITIALE)
! IN  COMPOR  : COMPORTEMENT
! IN  SIGMG    : CONTRAINTES GENERALISEES
! IN  VIM     : VARIABLES INTERNES
! IN SIGPG    : CONTRAINTES GENERALIEES
! OUT FINT    : FORCES INTERIEURES

! ----------------------------------------------------------------------
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8), parameter :: vrac2(6) = [1.d0, 1.d0, 1.d0, rac2, rac2, rac2]
    real(kind=8), dimension(6), parameter  :: kron = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
! ----------------------------------------------------------------------
    aster_logical :: axi
    integer(kind=8)       :: g, n, i
    integer(kind=8)       :: xu(ndim, nnos+nnom), xg(1, nnos), xp(1, nnos)
    integer(kind=8)       :: cod(npg)
    integer(kind=8)       :: nnu, nng, nnp, ndu, ndg, ndp, neu, neg, nep
    real(kind=8)  :: r, dff(nnos+nnom, ndim), dffb(nnos, ndim), poids
    real(kind=8)  :: butmp(2*ndim, ndim, nnos+nnom)
    real(kind=8)  :: bu(2*ndim+1, ndim, nnos+nnom)
    real(kind=8)  :: bg(1+ndim, 1, nnos), bp(1, 1, nnos)
    real(kind=8)  :: siefu(2*ndim+1), siefg(1+ndim), siefp(1)
    real(kind=8)  :: kr(2*ndim)

! --- INITIALISATION ---
    axi = ASTER_FALSE

    nnu = nnos+nnom
    nng = nnos
    nnp = nnos
    ndu = ndim
    ndg = 1
    ndp = 1
    neu = 1+2*ndim
    neg = 1+ndim
    nep = 1

    kr = kron(1:2*ndim)

    fint = 0
    cod = 0

    ! tableaux de reference bloc (depl,gonf,pres) -> numero du ddl
    forall (i=1:ndu, n=1:nng) xu(i, n) = (n-1)*(ndu+ndg+ndp)+i
    forall (i=1:ndp, n=1:nnp) xp(i, n) = (n-1)*(ndu+ndg+ndp)+ndu+i
    forall (i=1:ndg, n=1:nng) xg(i, n) = (n-1)*(ndu+ndg+ndp)+ndu+ndp+i
    forall (i=1:ndu, n=nng+1:nnu) xu(i, n) = (ndu+ndg+ndp)*nng+(n-1-nng)*ndu+i

    gauss: do g = 1, npg

        ! -----------------------!
        !  ELEMENTS CINEMATIQUES !
        ! -----------------------!

        ! Calcul des derivees des fonctions de forme P1
        call dfdmip(ndim, nnos, axi, geomi, g, iw, vffb(1, g), idffb, r, poids, dffb)

        ! Calcul des derivees des fonctions de forme P2
        call dfdmip(ndim, nnu, axi, geomi, g, iw, vff(1, g), idff, r, poids, dff)

        ! Assemblage de la matrice BU
        call nmbeps(axi, r, vff(:, g), dff, butmp)
        bu = 0.0
        bu(1:2*ndim, :, :) = butmp
        bu(2*ndim+1, :, :) = butmp(1, :, :)+butmp(2, :, :)+butmp(3, :, :)

        ! Assemblage des matrices BG et BP
        bg = 0.0
        bg(1, 1, :) = vffb(:, g)
        bg(2:neg, 1, :) = transpose(dffb)
        bp = 0.0
        bp(1, 1, :) = vffb(:, g)

        ! Correction matrices B
        if (option .eq. 'REFE_FORC_NODA') then
            bu = abs(bu)
            bg = abs(bg)
            bp = abs(bp)
        end if

        ! Contraintes generalisees EF par bloc
        if (ds_dil%inco) then
            siefu(1:2*ndim) = sief(dimdef*(g-1)+1:dimdef*(g-1)+2*ndim)*vrac2(1:2*ndim) &
                              -sief(dimdef*(g-1)+2*ndim+1)*kr
        else
            siefu(1:2*ndim) = sief(dimdef*(g-1)+1:dimdef*(g-1)+2*ndim)*vrac2(1:2*ndim)
        end if
        siefu(2*ndim+1) = sief(dimdef*(g-1)+2*ndim+1)
        siefg = sief(dimdef*(g-1)+2*ndim+2:dimdef*(g-1)+3*ndim+2)
        siefp = sief(dimdef*g)

        ! Forces interieures au point de Gauss g
        call add_fint(fint, xu, poids*prod_sb(siefu, bu))
        call add_fint(fint, xg, poids*prod_sb(siefg, bg))
        call add_fint(fint, xp, poids*prod_sb(siefp, bp))

    end do gauss

end subroutine
