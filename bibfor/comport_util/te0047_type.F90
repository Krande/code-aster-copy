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
! person_in_charge: jean-luc.flejou at edf.fr
!
module te0047_type
!
    use Behaviour_module, only: behaviourOption
    implicit none
!
#include "asterf_types.h"
!
    type :: te0047_dscr
        ! Ce qu'il faut sauvegarder en fonction de l'option :
        !       matrice Tangente
        !       vecteur Force
        !       vecteur Sigma
        !       vecteur Variables internes
        !   lMatr       :   FULL_MECA*             RIGI_MECA*
        !   lVect       :   FULL_MECA*  RAPH_MECA  RIGI_MECA_TANG
        !   lSigm       :   FULL_MECA*  RAPH_MECA  RIGI_MECA_TANG
        !   lVari       :   FULL_MECA*  RAPH_MECA
        !   lPred       :                          RIGI_MECA_TANG
        ! Si on fait uniquement une prédiction
        !   lMatrPred   :                          RIGI*
        !
        !   rigi        :   (1:4)=='RIGI'   |   (1:4)=='FULL'
        !   resi        :   (1:4)=='RAPH'   |   (1:4)=='FULL'
        !
        !   nomte   : nom terme élémentaire
        !   ndim    : dimension de l'espace
        !   nbt     : nombre de terme dans la matrice
        !   nno     : nombre de noeuds de l'élément
        !   nc      : nombre de composante par noeud
        !   neq     : nombre de terme de la matrice complète
        !   syme    : symetrique=1, non-symetrique=2
        !
        !   ugm     : déplacement moins dans le repère global
        !   dug     : incrément de déplacement dans le repère global
        !
        !   ulm     : déplacement moins dans le repère local
        !   dul     : incrément de déplacement dans le repère local
        !
        !   pgl     : matrice de passage global vers local
        !
        aster_logical       :: lVect = ASTER_FALSE
        aster_logical       :: lMatr = ASTER_FALSE
        aster_logical       :: lVari = ASTER_FALSE
        aster_logical       :: lSigm = ASTER_FALSE
        aster_logical       :: lPred = ASTER_FALSE
        aster_logical       :: lMatrPred = ASTER_FALSE
        aster_logical       :: lTraceDbg = ASTER_FALSE
        !
        aster_logical       :: rigi = ASTER_FALSE
        aster_logical       :: resi = ASTER_FALSE
        !
        integer             :: syme = 0
        !
        character(len=16)   :: option = ''
        character(len=16)   :: nomte = ''
        character(len=16)   :: rela_comp = ''
        character(len=16)   :: defo_comp = ''
        character(len=16)   :: type_comp = ''
        !
        integer             :: ndim
        integer             :: nbt
        integer             :: nno
        integer             :: nc
        integer             :: neq
        !
        real(kind=8)        :: ulm(12), ugm(12)
        real(kind=8)        :: dul(12), dug(12)
        real(kind=8)        :: pgl(3, 3)
        real(kind=8)        :: TempsPlus, TempsMoins
        real(kind=8)        :: Dilatation
        real(kind=8)        :: DeltaDilatation
        !
    end type te0047_dscr
!
! -------------------------------------------------------------------------------------------------
!
    public :: te0047_dscr
    public :: getDiscretInformations, te0047_dscr_write
!
    private
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
!
! -------------------------------------------------------------------------------------------------
!&<
contains
    ! ---------------------------------------------------------------------------------------------
    !
    ! Récupération d'informations sur les discrets
    !
    ! ---------------------------------------------------------------------------------------------
    ! IN
    !   nomte  : nom terme élémentaire
    !   option : nom de l'option
    !   syme   : 1 ou 2
    ! OUT
    !   tout le reste
    ! ---------------------------------------------------------------------------------------------
    subroutine getDiscretInformations(D)
        type(te0047_dscr), intent(inout) :: D
        !
        integer       :: icompo, itype, ideplm, ideplp, jtempsm, jtempsp, ii
        integer       :: imate, codret, cod_res(1), igeom
        real(kind=8)  :: temp_moins, temp_plus, temp_refe, val_res(1), xl2
        aster_logical :: IsOk
        !
        IsOk = (D%nomte .ne. '') .and. (D%option .ne. '') .and. &
               ((D%syme .eq. 1) .or. (D%syme .eq. 2))
        ASSERT(IsOk)
        !
        call infted(D%nomte, D%syme, D%nbt, D%nno, D%nc, D%ndim, itype)
        !
        D%neq = D%nno*D%nc
        ! Récupération des infos concernant les comportements :
        !     rela_comp   zk16(icompo-1+RELA_NAME)   NOM_DU_COMPORTEMENT
        !     nbvar       zk16(icompo-1+NVAR)        nbvar = read (zk16(icompo-1+NVAR),'(i16)')
        !     defo_comp   zk16(icompo-1+DEFO)        PETIT   PETIT_REAC  GROT_GDEP
        !     type_comp   zk16(icompo-1+INCRELAS)    COMP_ELAS   COMP_INCR
        !
        ! Properties of behaviour
        call jevech('PCOMPOR', 'L', icompo)
        D%rela_comp = zk16(icompo-1+RELA_NAME)
        D%defo_comp = zk16(icompo-1+DEFO)
        D%type_comp = zk16(icompo-1+INCRELAS)
        !
        ! Si c'est ELAS on calcule la dilatation thermique :
        !   Dilatation       = alpha*( T+ - Tref )*xl
        !   DeltaDilatation  = alpha*( T+ - T-   )*xl
        D%DeltaDilatation = 0.0
        D%Dilatation = 0.0
        if ( D%nno .eq. 1 ) goto 100
        if ( D%rela_comp.eq.'ELAS' ) then
            ! Température
            !   Si une température n'existe pas ==> pas de dilatation
            call rcvarc(' ', 'TEMP', '+',   'RIGI', 1, 1, temp_plus,  codret)
            if ( codret.ne.0 ) goto 100
            call rcvarc(' ', 'TEMP', '-',   'RIGI', 1, 1, temp_moins, codret)
            if ( codret.ne.0 ) goto 100
            call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, 1, temp_refe,  codret)
            if ( codret.ne.0 ) goto 100
            ! Matériau
            !   Il peut ne pas exister sur les discrets
            call tecach('ONN', 'PMATERC', 'L', codret, iad=imate)
            if ( codret .ne.0 ) goto 100
            call rcvalb('RIGI', 1, 1, '+', zi(imate), ' ', 'ELAS', &
                        0, ' ', [0.0d0], 1, 'ALPHA', val_res, cod_res, 0)
            !
            if ( cod_res(1).ne.0 ) goto 100
            call jevech('PGEOMER', 'L', igeom)
            igeom = igeom-1
            if ( D%ndim.eq.3 ) then
                xl2 = ( zr(igeom+4) - zr(igeom+1) )**2
                xl2 = ( zr(igeom+5) - zr(igeom+2) )**2 + xl2
                xl2 = ( zr(igeom+6) - zr(igeom+3) )**2 + xl2
            else if ( D%ndim.eq.2 ) then
                xl2 = ( zr(igeom+3) - zr(igeom+1) )**2
                xl2 = ( zr(igeom+4) - zr(igeom+2) )**2 + xl2
            else
                xl2 = 0.0
            endif
            D%DeltaDilatation = val_res(1)*(temp_plus-temp_moins)*sqrt(xl2)
            D%Dilatation      = val_res(1)*(temp_plus-temp_refe)*sqrt(xl2)
        endif
100     continue
        !
        ! Select objects to construct from option name
        !   lVari       : (1:9)'RAPH_MECA'  (1:9)'FULL_MECA'
        !   lSigm       : (1:9)'RAPH_MECA'  (1:9)'FULL_MECA'       'RIGI_MECA_TANG'
        !   lVect       : (1:9)'RAPH_MECA'  (1:9)'FULL_MECA'       'RIGI_MECA_TANG'
        !   lMatr       :                   (1:9)'FULL_MECA'  (1:9)'RIGI_MECA'
        !   lPred       :                                          'RIGI_MECA_TANG'
        !   lMatrPred   :                                     (1:4)'RIGI'
        call behaviourOption(D%option, zk16(icompo), D%lMatr, D%lVect, &
                             D%lVari, D%lSigm, codret)
        D%lMatrPred = D%option(1:4) .eq. 'RIGI'
        D%lPred = D%option .eq. 'RIGI_MECA_TANG'
        !   RIGI_MECA_TANG ->        DSIDEP        -->  RIGI
        !   FULL_MECA      ->  SIGP  DSIDEP  VARP  -->  RIGI  RESI
        !   RAPH_MECA      ->  SIGP          VARP  -->        RESI
        D%rigi = (D%option(1:4) .eq. 'RIGI') .or. (D%option(1:4) .eq. 'FULL')
        D%resi = (D%option(1:4) .eq. 'RAPH') .or. (D%option(1:4) .eq. 'FULL')
        !
        ! Déplacements dans le repère global :
        !     ugm = déplacement précédent
        !     dug = incrément de déplacement
        D%ugm(:) = 0.0
        D%dug(:) = 0.0
        D%ulm(:) = 0.0
        D%dul(:) = 0.0
        ! Instants
        D%TempsPlus  = 0.0
        D%TempsMoins = 0.0
        !
        ! Récupération : Déplacements, Instants à + et -
        !   Si déplacement - n'existe pas les autres non plus
        call tecach('NNN', 'PDEPLMR', 'L', codret, iad=ideplm)
        if ( codret .eq. 0 ) then
            ! call jevech('PDEPLMR', 'L', ideplm)
            call jevech('PDEPLPR', 'L', ideplp)
            do ii = 1, D%neq
                D%ugm(ii) = zr(ideplm+ii-1)
                D%dug(ii) = zr(ideplp+ii-1)
            end do
            !
            ! Temps + et - , calcul de dt
            call jevech('PINSTMR', 'L', jtempsm)
            call jevech('PINSTPR', 'L', jtempsp)
            D%TempsPlus  = zr(jtempsp)
            D%TempsMoins = zr(jtempsm)
        endif
        !
    end subroutine getDiscretInformations
    !
    !
    subroutine te0047_dscr_write(D)
        type(te0047_dscr), intent(in) :: D
        !
        write (*, '(5(A16))') [character(16) ::'option', 'nomte', 'rela_comp', &
                               'type_comp', 'defo_comp']
        write (*, '(5(A16))') D%option, D%nomte, D%rela_comp, D%type_comp, D%defo_comp
        !
        write (*, '(A)') 'Vect Matr Vari Sigm Pred MPred '
        write (*, '(6(L1,4X))') D%lVect, D%lMatr, D%lVari, D%lSigm, D%lPred, D%lMatrPred
        write (*, '(A)') 'resi rigi'
        write (*, '(2(L1,4X))') D%resi, D%rigi
    end subroutine te0047_dscr_write
!&>
end module te0047_type
