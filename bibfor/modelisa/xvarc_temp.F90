! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine xvarc_temp(novarc, evouch, evol, prolga, proldr, finst,&
                      nboccv, carte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/xtmafi.h"
    character(len=8), intent(in) :: novarc, evouch, evol, finst
    integer, intent(in) :: nboccv
    character(len=16), intent(in) :: prolga, proldr
    character(len=19), intent(in) :: carte
! person_in_charge: sam.cuvilliez at edf.fr
! ----------------------------------------------------------------------
!
!     AFFE_MATERIAU / AFFE_VARC
!
!     -> cas particulier du chainage thermo-mecanique avec X-FEM
!
!     but : modifier dans la carte chmat//'.TEMP    .2' le nom
!           symbolique associe a la variable de commande TEMP.
!           'TEMP' doit etre remplace par 'TEMP_ELGA'.
!           Cette modification affecte les mailles qui portent
!           des elements enrichis dans le modele.
!
!     in novarc : nom de la varc de l'occurence courante de AFFE_VARC
!     in evouch : 'EVOL' / 'CHAMP' / 'VIDE'
!     in evol   : nom de la sd resultat eventuellement renseignee
!                 sous le MCS EVOL
!     in prolga : valeur du MCS PROL_GAUCHE
!     in proldr : valeur du MCS PROL_DROITE
!     in finst  : nom de la sd fonction eventuellement renseignee
!                 sous le MCS FONC_INST
!     in nboccv : nombre d'occurence du MCF AFFE_VARC
!     in carte  : carte chmat//'.TEMP    .2'
!
!     "out"     : ecrire dans carte si les conditions sont reunies
!
! ----------------------------------------------------------------------
    integer :: iret, nfiss, nbmx, jmax, ndim
    integer :: ibid, nbord, iord, icode_ini
    character(len=8) :: modein, modevo, noma
    character(len=19) :: chamel, resu19
    character(len=24) :: mesmai, lismai, ligrch
    integer, pointer :: vordr(:) => null()
    integer, pointer :: vcode(:) => null()
    character(len=8), pointer :: fiss(:) => null()
    character(len=8), pointer :: p_mod_ther(:) => null()
    character(len=16), pointer :: vale(:) => null()
    character(len=24), pointer :: vcelk(:) => null()
    character(len=24), pointer :: vligr(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
! ----------------------------------------------------------------------
! --- verifications prealables
! ----------------------------------------------------------------------
!
!   on sort s'il ne s'agit pas la variable de commande TEMP
    if (novarc .ne. 'TEMP') goto 999
!
!   on sort si le MCS EVOL n'est pas renseigne
    if (evouch .ne. 'EVOL') goto 999
!
!   on sort si evol ne contient pas le champ de nom symbolique TEMP_ELGA
    resu19 = evol
    call jeveuo(resu19//'.ORDR', 'L', vi=vordr)
    call jelira(resu19//'.ORDR', 'LONUTI', nbord)
    ASSERT(nbord .ge. 1)
    AS_ALLOCATE(vi=vcode, size=nbord)
    do iord=1,nbord
        call rsexch(' ', evol, 'TEMP_ELGA', vordr(iord), chamel, vcode(iord))
    enddo
!   verif de coherence (pb par ex. si on a fait une mauvaise utilisation de CREA_RESU)
    icode_ini = vcode(1)
    do iord=1,nbord
        ASSERT(vcode(iord) .eq. icode_ini)
    enddo
    AS_DEALLOCATE(vi=vcode)
    if (icode_ini .ne. 0) goto 999
!
! ----------------------------------------------------------------------
! --- TEMP_ELGA est present -> on s'est assure que l'utilisateur veut
! --- faire du chainage thermo-mecanique avec un resultat thermique xfem.
! --- Il faut faire des verifications supplementaires
! ----------------------------------------------------------------------
!
!   evol doit faire reference a une et une seule sd_modele "modevo"
!   (pb par ex. si on a fait une mauvaise utilisation de CREA_RESU)
!   -> on recupere modevo via le ligrel de definition des cham_elem
!      TEMP_ELGA, ligrel qui doit lui aussi etre unique
    AS_ALLOCATE(vk24=vligr, size=nbord)
    do iord=1,nbord
        call rsexch(' ', evol, 'TEMP_ELGA', vordr(iord), chamel, ibid)
        call jeveuo(chamel//'.CELK', 'L', vk24=vcelk)
        ligrch = vcelk(1)
        call exisd('LIGREL', ligrch, iret)
        ASSERT(iret .eq. 1)
        vligr(iord) = ligrch
    enddo
    ligrch = vligr(1)
    do iord=1,nbord
        ASSERT(vligr(iord) .eq. ligrch)
    enddo
    AS_DEALLOCATE(vk24=vligr)
    modevo = ligrch(1:8)
    call exisd('MODELE', modevo, iret)
    ASSERT(iret .eq. 1)
!
!   modevo est necessairement un modele xfem
    call exixfe(modevo, iret)
    ASSERT(iret .eq. 1)
!
!   dans ce cas on ne peut avoir qu'une seule occurence de AFFE_VARC
    if (nboccv .ne. 1) call utmess('F', 'XFEM_96')
!
!   dans ce cas le MCS MODELE devient obligatoire
    call getvid(' ', 'MODELE', scal=modein, nbret=iret)
    if (iret .ne. 1) call utmess('F', 'XFEM_97')
!
!   enfin on s'assure que ce modele "modein" a ete cree par MODI_MODELE_XFEM
!   avec le MCS MODELE_THER == modevo
    call jeexin(modein//'.MODELE_THER', iret)
    ASSERT(iret .ne. 0)
    call jeveuo(modein//'.MODELE_THER', 'L', vk8=p_mod_ther)
    ASSERT(p_mod_ther(1) .eq. modevo)
!
! ----------------------------------------------------------------------
! --- recuperation des mailles portant des EF (principaux) enrichis
! ----------------------------------------------------------------------
!
!   recup de la dimension du probleme
    call dismoi('DIM_GEOM', modein, 'MODELE', repi=ndim)
    ASSERT( (ndim .eq. 2) .or. (ndim .eq. 3) )
!
    lismai = '&&XVARCT.NUM_MAILLES'
    mesmai = '&&XVARCT.MES_MAILLES'
!
    call dismoi('NOM_MAILLA', modein, 'MODELE', repk=noma)
    call dismoi('NB_FISS_XFEM', modein, 'MODELE', repi=nfiss)
    call jeveuo(modein//'.FISS', 'L', vk8=fiss)
!
    call xtmafi(ndim, fiss, nfiss, lismai, mesmai, nbmx, model=modein)
    call jeveuo(lismai, 'L', jadr=jmax)
!
! ----------------------------------------------------------------------
! --- modification dans la carte du nom symbolique du champ de
! --- temperature pour ces mailles : 'TEMP' -> 'TEMP_ELGA'
! ----------------------------------------------------------------------
!
    call jeveuo(carte//'.VALV', 'E', vk16=vale)
!
    vale(1) = 'TEMP'
    vale(2) = 'EVOL'
    vale(3) = evol
    vale(4) = 'TEMP_ELGA'
    vale(5) = prolga
    vale(6) = proldr
    vale(7) = finst
!
    call nocart(carte, 3, 7, mode='NUM', nma=nbmx, limanu=zi(jmax))
!
    call jedetr(mesmai)
    call jedetr(lismai)
!
999 continue
!
    call jedema()
!
end subroutine
