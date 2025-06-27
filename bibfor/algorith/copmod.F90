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
subroutine copmod(base, bmodr, bmodz, champ, numer, &
                  nbmodes, nequa)
    implicit none
!
!                              Function
!     _____________________________________________________________
!    | Extract, inside a temporary work vector, the fields from a  |
!    | modal basis concept while changing or not their numbering   |
!    |_____________________________________________________________|
!
! ---------
! Examples: call copmod ( base, bmodr = zr(jbase) )
! --------- call copmod ( base, bmodr = zr(jbase), numer = nuddl )
!           call copmod ( base, bmodz = zc(jbase) )
!           call copmod ( base, bmodr = zr(jbase), numer = nuddl, champ = 'VITE')
!
!
!                     Description of the input/output arguments
!   _________________________________________________________________________________
!  | in < obl > base      : Entry modal basis (mode_meca)                        [k8]|
!  |            ------                                                               |
!  |out < fac >|basout|   : Output modal basis (mode_meca)                       [r8]|
!  |            ======    : after changing numbering                                 |
!  |out < fac >|bmodr |   : Output work vector containing the copied fields      [r8]|
!  |            ======    : (real fields case)                                       |
!  |out < fac >|bmodz |   : Output work vector containing the copied fields     [c16]|
!  |            ------|   : (complex fields case)                                    |
!  |                  |                                                              |
!  |                   => Validation rule : One of these two output vectors must     |
!  |                                        be given                                 |
!  |---------------------------------------------------------------------------------|
!  | in < fac > champ     : Field type to copy (default = 'DEPL')                [k*]|
!  | in < fac > numer     : - If given, the name of the nume_ddl concept         [k*]|
!  |                      :   or that of the nume_equa_user giving the new numbering      |
!  |                      :   of the copied fields                                   |
!  |                      : - If absent, the numbering is unchanged                  |
!  | in < fac > nbmodes   : The number of modes to be copied                      [i]|
!  |                      : (default = total number of modes in base)              |
!  | in < fac > nequa     : The number of equations in each vector                [i]|
!  |                      : (default = determine automatically from the nume_ddl)    |
!  |_________________________________________________________________________________|
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copy_field_with_numbering.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/idensd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
#include "blas/zcopy.h"
!   ___________________________________________________________________
!
!  - 0 - INITIALISATIONS DIVERSES
!   ___________________________________________________________________
!
!     0.1 - DECLARATION DES VARIABLES D'ENTREE/SORTIE
!
    character(len=8), intent(in) :: base
    real(kind=8), optional, intent(out) :: bmodr(*)
    complex(kind=8), optional, intent(out) :: bmodz(*)
    character(len=*), optional, intent(in) :: champ
    character(len=*), optional, intent(in) :: numer
    integer(kind=8), optional, intent(in) :: nbmodes
    integer(kind=8), optional, intent(in) :: nequa
!
!     0.2 - DECLARATION DES VARIABLES LOCALES
!
    character(len=1) :: typc, typbase
    character(len=4) :: docu
    aster_logical :: modnum, chnoeud, is_nume_equa_user, r2zbase
    integer(kind=8) :: i, iret, neq, nbmode
    integer(kind=8) :: jdeeq, jval
    character(len=14) :: nume_ddl_user, nume_ddl_field
    character(len=16) :: champ2
    character(len=19) :: nume_equa_user, nume_equa_field
    character(len=19) :: nomcha, tmpcha
    character(len=24) :: maill1, maill2, valk(4), valcha
    character(len=24), pointer :: refe(:) => null()
    blas_int :: b_incx, b_incy, b_n
!
!     0.3 - ACTUALISATION DE LA VALEUR DE LA MARQUE COURANTE
!
    call jemarq()
!
!     0.4 - TEST DES ARGUMENTS D'ENTREES, ET ATTRIBUTION DES VALEURS PAR DEFAUT
!
    typc = 'R'
    ASSERT(UN_PARMI2(bmodr, bmodz))
    if (present(bmodz)) typc = 'C'
!
    neq = 0
    !
    champ2 = 'DEPL'
    if (present(champ)) champ2 = champ
    !
! Determine if input "numer" is a nume_ddl (14 caracters) or a nume_equa_user (19 caracters)
! Retrieve the underlying nume_equa_user: "nume_equa_user"
    nume_equa_user = ' '
    is_nume_equa_user = .false.
!
    if (present(numer)) then
        call jeexin(numer(1:14)//'.NUME.NEQU', iret)
        if (iret .ne. 0) then
            nume_ddl_user = numer(1:14)
            nume_equa_user = nume_ddl_user//'.NUME'
        else
            nume_equa_user = numer
        end if
        is_nume_equa_user = .true.
    end if
    if (present(nequa)) neq = nequa
!
    chnoeud = .true.
    if (champ2(6:7) .eq. 'EL') chnoeud = .false.
!
!   --- RECUPERATION/VERIFICATION DU NOMBRE D'EQUATIONS RENSEIGNE
!   --- 1. CHAMP AUX NOEUDS : PAR RAPPORT A L'INFORMATION DANS LE NUME_DDL
    if (chnoeud) then
        !
! If optional argument nequa is given, dimension should be checked
! Check with reference nume_ddl if provided
        if (present(numer)) then
            call dismoi('NB_EQUA', nume_equa_user, 'NUME_EQUA', repi=neq)
            if (present(nequa)) then
                ASSERT(nequa .eq. neq)
            end if
        else
! Check with nume_ddl in resu_dyna (meant to be reference) otherwise
            call dismoi('NUME_DDL', base, 'RESU_DYNA', repk=nume_ddl_field, arret='C', &
                        ier=iret)
            if (iret .eq. 0) then
                call jeexin(nume_ddl_field//'.NUME.NEQU', iret)
                if (iret .ne. 0) then
                    call dismoi('NB_EQUA', nume_ddl_field, 'NUME_DDL', repi=neq)
                    if (present(nequa)) then
                        ASSERT(nequa .eq. neq)
                    end if
                end if
            end if
        end if
    else
!   --- 2. CHAMP AUX ELEMENTS : PAR RAPPORT A UN CHAMP DU MEME TYPE DE LA BASE
        call rsexch('F', base, champ2, 1, nomcha, &
                    iret)
        call jelira(nomcha(1:19)//'.CELV', 'LONMAX', neq)
        if (present(nequa)) then
            ASSERT(nequa .eq. neq)
        end if
    end if
    ASSERT(neq .ne. 0)
!   --- FIN DE LA RECUPERATION/VERIFICATION DU NOMBRE D'EQUATIONS
!
    call dismoi('NB_MODES_TOT', base, 'RESULTAT', repi=nbmode)
    if (present(nbmodes)) then
        ASSERT(nbmodes .le. nbmode)
        nbmode = nbmodes
    end if
!
!  ____________________________________________________________________
!
!  - 1 - RECHERCHE DES INFORMATIONS SUR LES CHAMPS DANS LA BASE MODALE
!  ____________________________________________________________________
!
!     1.1 - CHERCHER UN OBJET .REFE DANS UN CHAMP DE LA BASE MODALE
!
!     1.1.1 - RECUPERER LE NOM DE CHAMP DU 1ER NUMERO ORDRE
!
    call rsexch('F', base, champ2, 1, nomcha, &
                iret)
    valcha = nomcha//'.VALE'
    call jeexin(nomcha//'.VALE', iret)
    if (iret .le. 0) valcha = nomcha//'.CELV'
    call jelira(valcha, 'TYPE', cval=typbase)
    r2zbase = (typc == 'C') .and. (typbase == 'R')
!
!     1.1.2 - POUR TRAITER LES CAS AVEC SS-STRUCTURATION, TESTER SI
!             L'OBJET .REFE EXISTE DANS CE CHAMP, SI NON (IRET.EQ.0)
!             RECUPERER LE .REFE DU CHAMP DE DEPLACEMENT
!
    call jeexin(nomcha//'.REFE', iret)
    if (iret .eq. 0) call rsexch('F', base, 'DEPL', 1, nomcha, &
                                 iret)
!
    call jeveuo(nomcha//'.REFE', 'L', vk24=refe)
!
!     1.2 - EXTRAIRE LE NOM DE MAILLAGE .REFE[1] ET DU NUME_DDL .REFE[2]
!
    call dismoi('DOCU', nomcha, 'CHAMP', repk=docu)
!
    if (docu == "CHNO") then
        call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=maill1)
        call dismoi('NUME_EQUA', nomcha, 'CHAM_NO', repk=nume_equa_field)
    else
        maill1 = refe(1) (1:8)
        nume_equa_field = refe(2) (1:19)
    end if
!
!     1.3 - TRAITEMENT DES CAS AVEC UN NUME_EQUA ET NON PAS UN NUME_DDL
!           COMPLET.
!
    if (is_nume_equa_user) then
        call dismoi('NOM_MAILLA', nume_equa_user, 'NUME_EQUA', repk=maill2)
!       --- ON NE FAIT PAS DE TEST DE COMPATIBILITE SUR LES MAILLAGES
!         - SI ON NE DISPOSE PAS DE NUMEDDL COMPLET
!         - IMPORTANT : LE TEST DOIT SE FAIRE QUAND MEME EN DEHORS DE
!                       L'APPEL A COPMOD (VOIR OP0072 PAR EXEMPLE)
    else
        maill2 = maill1
    end if
!
!     1.4 - LIBERER L'OBJET .REFE PARCE QU'ON N'EN A PLUS BESOIN
!
    call jelibe(nomcha(1:19)//'.REFE')
!  ____________________________________________________________________
!
!  - 2 - RECHERCHE DES INFORMATIONS SUR LA NUMEROTATION FINALE
!  ____________________________________________________________________
!
!     2.1 - NOUVELLE NUMEROTATION ? (SUR UN UNIQUE MAILLAGE)
!     2.2 - SI OUI, VERIFIER LA COMPATIB. DES 2 MAILLAGES DES NUME_DDL
!           RESTITUTION SUR SQUELLETE : CAS SPECIAL
!
    call jeexin(maill1(1:8)//'.INV.SKELETON', iret)
    modnum = .false.
    if (is_nume_equa_user .and. (iret .eq. 0)) then
        if (.not. idensd('NUME_EQUA', nume_equa_user, nume_equa_field)) then
            modnum = .true.
            call dismoi('NOM_MAILLA', nume_equa_user, 'NUME_EQUA', repk=maill2)
            if (maill1 .ne. maill2) then
                valk(1) = nume_equa_user
                valk(2) = maill2
                valk(3) = nume_equa_field
                valk(4) = maill1
                call utmess('F', 'ALGORITH12_62', nk=4, valk=valk)
            end if
        end if
    end if
!
!
!     2.3 - RECUPERER L'OBJET .DEEQ
!
    if (modnum) then
        call jeveuo(nume_equa_user//'.DEEQ', 'L', jdeeq)
    else
        call jeveuo(nume_equa_field//'.DEEQ', 'L', jdeeq)
    end if
!  ____________________________________________________________________
!
!  - 3 - RECOPIE DES CHAMPS ET MODIFICATION DE LA NUMER. SI NECESSAIRE
!  ____________________________________________________________________
!
!     3.1 - BOUCLE SUR LES MODES DE LA BASE
    do i = 1, nbmode
!       3.1.1 - EXTRAIRE LE NOM DU CHAMP D'INTERET (NOMCHA)
        call rsexch('F', base, champ2, i, nomcha, &
                    iret)
!
!       3.1.2 - NOUVELLE NUMER.? BASE IN REELLE et BASE OUT COMPLEXE?
!               ALORS CREER UN NOUVEAU CHAMP TEMPORAIRE
!               AVEC LA BONNE NUMEROTATION
        if (docu(1:4) == 'CHNO' .and. (modnum .or. r2zbase)) then
            tmpcha = '&&COPMOD.CHAMP'
            call copy_field_with_numbering(nomcha, tmpcha, maill2, nume_equa_user, 'V', &
                                           typc=typc, nequa=neq)
            nomcha = tmpcha
        end if
!
!       3.1.3 - OBTENIR L'OBJET DES VALEURS DU CHAMP (.VALE OU .CELV)
!               POUR LES CHAM_NO ET CHAM_ELEM RESPECTIVEMENT
        valcha = nomcha(1:19)//'.VALE'
        call jeexin(nomcha(1:19)//'.VALE', iret)
        if (iret .le. 0) valcha = nomcha(1:19)//'.CELV'
        call jeveuo(valcha, 'L', jval)
!
!       3.1.4 - COPIER LES VALEURS DU CHAMP DANS LE VECTEUR DE SORTIE
        if (typc .ne. 'C') then
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jval), b_incx, bmodr((i-1)*neq+1), b_incy)
        else
            b_n = to_blas_int(neq)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call zcopy(b_n, zc(jval), b_incx, bmodz((i-1)*neq+1), b_incy)
        end if
!
!       3.1.5 - MENAGE ET LIBERATION DE LA MEMOIRE SELON LE BESOIN
        call jelibe(valcha)
        if (docu(1:4) == 'CHNO' .and. (modnum .or. r2zbase)) then
            call detrsd('CHAMP', tmpcha)
        end if
!
!       3.1.6 - ANNULER LES DDL DE LAGRANGE S'IL S'AGIT DES CHAMPS DE
!               DEPLACEMENTS
        if (champ2 .eq. 'DEPL') then
            if (typc .ne. 'C') then
                call zerlag(neq, zi(jdeeq), vectr=bmodr((i-1)*neq+1))
            else
                call zerlag(neq, zi(jdeeq), vectz=bmodz((i-1)*neq+1))
            end if
        end if
!
    end do
!     FIN DE LA BOUCLE (3.1) SUR LES MODES
!  ____________________________________________________________________
!
    call jedema()
end subroutine
