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

subroutine cfinal(ds_contact, nbliac)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cftabl.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8) :: nbliac
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODE DISCRETE - ALGORITHME)
!
! ACTIVATION DES LIAISONS INITIALES
!
! ----------------------------------------------------------------------
!
! In  ds_contact       : datastructure for contact management
! I/O NBLIAC : NOMBRE DE LIAISONS ACTIVES
!
!
    aster_logical :: liaact, liaexi
    real(kind=8) :: jeuini, jeumin
    integer(kind=8) :: posit, ajliai, spliai, indic, nbliac_init
    integer(kind=8) :: nbliai
    integer(kind=8) :: iliai, iliac
    aster_logical :: lgcp, lgliss
    character(len=1) :: typeaj
    character(len=19) :: liac
    integer(kind=8) :: jliac
    character(len=24) :: jeuite, jeux
    integer(kind=8) :: jjeuit, jjeux
    character(len=24) :: numlia
    integer(kind=8) :: jnumli
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- LECTURE DES STRUCTURES DE DONNEES DE CONTACT
!
    liac = ds_contact%sdcont_solv(1:14)//'.LIAC'
    call jeveuo(liac, 'L', jliac)
    jeuite = ds_contact%sdcont_solv(1:14)//'.JEUITE'
    jeux = ds_contact%sdcont_solv(1:14)//'.JEUX'
    call jeveuo(jeuite, 'L', jjeuit)
    call jeveuo(jeux, 'L', jjeux)
    numlia = ds_contact%sdcont_solv(1:14)//'.NUMLIA'
    call jeveuo(numlia, 'L', jnumli)
!
! --- INITIALISATIONS
!
    jeumin = r8prem()
    posit = 0
    typeaj = 'A'
    spliai = 0
    ajliai = 0
    nbliac_init = nbliac
!
! --- PARAMETRES
!
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    lgcp = cfdisl(ds_contact%sdcont_defi, 'CONT_GCP')
    lgliss = cfdisl(ds_contact%sdcont_defi, 'CONT_DISC_GLIS')
!
! --- DETECTION DES COUPLES DE NOEUDS INTERPENETRES
!
    do iliai = 1, nbliai
!
! ----- JEU SANS CORRECTION DU CONTACT
!
        jeuini = zr(jjeux+3*(iliai-1)+1-1)
!
! ----- LIAISON ACTIVEE ?
!
        liaact = .false.
        if (lgcp) then
            liaact = .true.
        else
            if (jeuini .lt. jeumin) then
                liaact = .true.
            else
                liaact = .false.
            end if
        end if
!
! ----- LIAISON GLISSIERE -> TOUTES LES LIAISONS SONT ACTIVEES
!
        if (lgliss) then
            liaact = .true.
        end if
!
! ----- LA LIAISON EXISTE-T-ELLE DEJA ?
!
        liaexi = .false.
        do iliac = 1, nbliac_init
            if (zi(jliac-1+iliac) .eq. iliai) then
                liaexi = .true.
            end if
        end do
!
! ----- INDICE DE LA NOUVELLE LIAISON ACTIVE
!
        if (liaact) then
            if (lgcp) then
                posit = iliai
            else
                posit = nbliac+1
            end if
        end if
!
! ----- ACTIVATION DE LA LIAISON DE CONTACT
!
        if (liaact) then
            call cftabl(indic, nbliac, ajliai, spliai, &
                        ds_contact%sdcont_solv, typeaj, posit, &
                        iliai)
        end if
!
    end do
!
    call jedema()
!
end subroutine
