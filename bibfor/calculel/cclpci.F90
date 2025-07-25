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

subroutine cclpci(option, modele, resuin, resuou, mater, mateco, &
                  carael, ligrel, numord, nbpain, lipain, &
                  lichin, codret)
!
    use HHO_precalc_module, only: hhoAddInputField
!
    implicit none
!     --- ARGUMENTS ---
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbpain, numord, codret, nb_in_maxi
    character(len=8) :: modele, resuin, resuou, mater, mateco, carael
    character(len=8) :: lipain(*)
    character(len=16) :: option
    character(len=24) :: lichin(*), ligrel
!  CALC_CHAMP - DETERMINATION LISTE DE PARAMETRES ET LISTE DE CHAMPS IN
!  -    -                     -        -                      -      -
! ----------------------------------------------------------------------
!
! IN  :
!   OPTION  K16  NOM DE L'OPTION A CALCULER
!   MODELE  K8   NOM DU MODELE
!   RESUIN  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT IN
!   RESUOU  K8   NOM DE LA STRUCUTRE DE DONNEES RESULTAT OUT
!   MATER   K8   NOM DU MATERIAU
!   CARAEL  K8   NOM DU CARAELE
!   LIGREL  K24  NOM DU LIGREL
!   NUMORD  I    NUMERO D'ORDRE COURANT
!
! OUT :
!   NBPAIN  I    NOMBRE DE PARAMETRES IN
!   LIPAIN  K8*  LISTE DES PARAMETRES IN
!   LICHIN  K8*  LISTE DES CHAMPS IN
!   CODRET  I    CODE RETOUR (0 SI OK, 1 SINON)
! ----------------------------------------------------------------------
! person_in_charge: nicolas.sellenet at edf.fr
!
    integer(kind=8) :: opt, iaopds, iaoplo, iapara, nparin, ipara, opt2, ierd
    integer(kind=8) :: decal
    character(len=8) :: noma
    character(len=16) :: optio2
    character(len=19) :: nochin
!
    call jemarq()
!
    codret = 0
!
    if (option(6:9) .eq. 'NOEU') then
        nparin = 0
    else
        call jenonu(jexnom('&CATA.OP.NOMOPT', option), opt)
        call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopds)
        call jeveuo(jexnum('&CATA.OP.LOCALIS', opt), 'L', iaoplo)
        call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iapara)
!
        nparin = zi(iaopds-1+2)
        nbpain = 0
    end if
!
!     BOUCLE SUR LES PARAMETRES DE L'OPTION
    do ipara = 1, nparin
        nochin = ' '
!
        nbpain = nbpain+1
        lipain(nbpain) = zk8(iapara+ipara-1)
!
        optio2 = zk24(iaoplo+3*ipara-2) (1:16)
!
!       CAS OU CE PARAM EST UNE OPTION OU UN CHAMP DANS LA
!       SD RESULTAT
        call jenonu(jexnom('&CATA.OP.NOMOPT', optio2), opt2)
        if ((opt2 .ne. 0) .or. (zk24(iaoplo+3*ipara-3) .eq. 'RESU')) then
            if (zk24(iaoplo+3*ipara-1) .eq. 'NP1') then
                decal = 1
            else if (zk24(iaoplo+3*ipara-1) (1:3) .eq. 'NM1') then
                decal = -1
            else
                decal = 0
            end if
            call rsexch(' ', resuin, optio2, numord+decal, nochin, &
                        ierd)
            if (ierd .ne. 0) then
                call rsexch(' ', resuou, optio2, numord+decal, nochin, &
                            ierd)
            end if
!
            if (ierd .ne. 0) then
                if ((option .eq. optio2)) then
!             CAS OU UN CHAMP DEPEND DE LUI MEME A L'INSTANT N-1
!             EXEMPLE : ENDO_ELGA
                    if (zk24(iaoplo+3*ipara-1) .eq. 'NM1T') then
                        nochin = '&&CALCOP.INT_0'
                        call alchml(ligrel, optio2, lipain(nbpain), 'V', nochin, &
                                    ierd, ' ')
                        if (ierd .gt. 0) then
                            call utmess('A', 'CALCCHAMP_19', sk=option)
                            goto 10
                        end if
                    else
                        call rsexch(' ', resuou, optio2, numord+decal, nochin, &
                                    ierd)
                        call alchml(ligrel, optio2, lipain(nbpain), 'G', nochin, &
                                    ierd, ' ')
                        if (ierd .gt. 0) then
                            call utmess('A', 'CALCCHAMP_19', sk=option)
                            goto 10
                        end if
                        call rsnoch(resuou, optio2, numord+decal)
                    end if
                else
                    nochin = ' '
                end if
            end if
!       CAS OU CE PARAM EST UN OBJET DU MAILLAGE
        else if (zk24(iaoplo+3*ipara-3) .eq. 'MAIL') then
            call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
            nochin = noma//zk24(iaoplo+3*ipara-2)
!       CAS OU CE PARAM EST UN OBJET DU MODELE
        else if (zk24(iaoplo+3*ipara-3) .eq. 'MODL') then
            nochin = modele//zk24(iaoplo+3*ipara-2)
!       CAS OU CE PARAM EST UN OBJET DU CARA_ELEM
        else if (zk24(iaoplo+3*ipara-3) .eq. 'CARA') then
            nochin = carael//zk24(iaoplo+3*ipara-2)
!       CAS OU CE PARAM EST UN OBJET PARTICULIER SUR LA VOLATILE
        else if (zk24(iaoplo+3*ipara-3) .eq. 'VOLA') then
            nochin = zk24(iaoplo+3*ipara-2)
!       CAS OU CE PARAM EST UN OBJET DU CHAMMAT
        else if (zk24(iaoplo+3*ipara-3) .eq. 'CHMA') then
            nochin = mater//zk24(iaoplo+3*ipara-2)
!       CAS OU CE PARAM EST UN OBJET DU MATECODE
        else if (zk24(iaoplo+3*ipara-3) .eq. 'MACO') then
            nochin = mateco//zk24(iaoplo+3*ipara-2)
        end if
        lichin(nbpain) = nochin
10      continue
    end do
!
    nb_in_maxi = nbpain+3
    call hhoAddInputField(modele, nb_in_maxi, lichin, lipain, nbpain)
!
    call jedema()
!
end subroutine
