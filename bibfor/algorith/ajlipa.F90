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

subroutine ajlipa(modelz, base, kdis)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvis.h"
#include "asterfort/gnoms3.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"

    character(len=*), intent(in) :: modelz
    character(len=1), intent(in) :: base
    character(len=24), intent(in) :: kdis
! ----------------------------------------------------------------------
!  But :
!     Creation (ou modification) de la sd_partition d'un modele
!     (Commandes AFFE_MODELE et MODI_MODELE)
!  Remarques :
!     * la sd n'est creee que dans le cas du parallelisme mpi distribue
!     * il faut appeler cette routine apres adalig si cette derniere
!       est appelee (cas de op0018)
! ----------------------------------------------------------------------
    character(len=8) :: modele, mopart, valk(3), nomres, partsd
    character(len=16) :: typres, nomcom
    character(len=19) :: ligrmo
    character(len=24) :: k24b
    integer(kind=8) :: i, rang, nbproc, ifm, niv, ibid, nbsd, nbma
    integer(kind=8) :: nmpp, nmp0, nmp0af, ico, nbpro1, krang, nmp1
    integer(kind=8) :: iexi, iad
    integer(kind=8) :: icobis, dist0, jnumsd, vali(3), nbmamo, ima
    integer(kind=8) ::  jprti, jprtk
    aster_logical :: plein0, exi_partsd
    integer(kind=8), pointer :: fdim(:) => null()
    character(len=8), pointer :: fref(:) => null()
    integer(kind=8), pointer :: maille(:) => null()
    mpi_int :: mrank, msize
    data k24b/' '/

!-----------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)

! ----------------------------------------------------------------------
!   --  Verifications et initialisations
! ----------------------------------------------------------------------
    modele = modelz
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
    call dismoi('PARTITION', ligrmo, 'LIGREL', repk=partsd)

    call getres(nomres, typres, nomcom)
    ASSERT(nomres .eq. modele)
    ASSERT(nomcom .eq. 'AFFE_MODELE' .or. nomcom .eq. 'MODI_MODELE')

!   -- s'il existe deja une partition, on la detruit :
!   --------------------------------------------------
    call exisd('PARTITION', partsd, iexi)
    if (iexi .gt. 0) then
        ASSERT(nomcom .eq. 'MODI_MODELE')
        call detrsd('PARTITION', partsd)
    end if
!
!   -- s'il n'y a pas d'elements finis dans le modele :
!   ---------------------------------------------------
    call jeexin(ligrmo//'.LIEL', iexi)
    if (iexi .eq. 0) goto 999

!   -- s'il n'a qu'un seul proc, il n'y a rien a faire :
!   ----------------------------------------------------
    nbproc = 1
    rang = 0
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
    if (nbproc .le. 1) goto 999

!   -- si le modele n'a pas de mailles, il n'y a rien a faire :
!   -----------------------------------------------------------
    call jeexin(ligrmo//'.TYFE', iexi)
    if (iexi .eq. 0) goto 999

!   -- si l'utilisateur ne veut pas de distribution des calculs,
!      il n'y a rien a faire :
!   ------------------------------------------------------------
    if (kdis .eq. 'CENTRALISE') goto 999

!   -- la partsd est a fournir si 'SOUS_DOMAINE'
!   -------------------------------------------------------------------
    exi_partsd = kdis .eq. 'SOUS_DOMAINE'

! ----------------------------------------------------------------------
!   Lecture des mot-cles et verifications supplementaires
!   Creation de la partsd
! ----------------------------------------------------------------------
    dist0 = 0

!   -- Creation de la sd_partition :
!   ----------------------------------------------------
    if (partsd == ' ') then
        partsd = "PART"
        call gnoms3(partsd, 5, 8, ".PRTK")
        call jeveuo(ligrmo//'.LGRF', 'E', iad)
        zk8(iad-1+2) = partsd
    end if
    call jeveuo(ligrmo//'.TYFE', 'L', vi=maille)
    call jelira(ligrmo//'.TYFE', 'LONMAX', nbma)
    call wkvect(partsd//'.PRTI', base//' V I', 1, jprti)
    zi(jprti-1+1) = nbproc
    call wkvect(partsd//'.PRTK', base//' V K24', 1, jprtk)
    zk24(jprtk-1+1) = kdis
    if (kdis(1:5) .eq. 'MAIL_') then
        call wkvect(partsd//'.NUPR', base//' V I', nbma+1, jnumsd)
        zi(jnumsd-1+nbma+1) = nbproc

!       nbmamo : nbre de mailles du modele
        nbmamo = 0
        do ima = 1, nbma
            zi(jnumsd-1+ima) = -999
            if (maille(ima) .ne. 0) nbmamo = nbmamo+1
        end do
    end if

!   -- Recuperations des mot-cles :
!   -------------------------------
!   uniquement pour DISTRIBUTION/METHODE='MAIL_*'
    call getvis('DISTRIBUTION', 'CHARGE_PROC0_MA', iocc=1, scal=dist0, nbret=ibid)

!   -- Verification pour le cas du partitionnement avec partsd :
!   -----------------------------------------------------------------
    if (exi_partsd) then
        call jeveuo(partsd//'.FREF', 'L', vk8=fref)
        mopart = fref(1)
        if (modele .ne. mopart) then
            valk(1) = partsd(1:8)
            valk(2) = modele
            valk(3) = mopart
            call utmess('F', 'PARTITION1_17', nk=3, valk=valk)
        end if
    end if

!   -- Verifications sur le nombre de mailles ou de sous-domaines :
!      par rapport au nombre de processeurs
!   ---------------------------------------------------------------
    if (exi_partsd) then
        call jeveuo(partsd//'.FDIM', 'L', vi=fdim)
        nbsd = fdim(1)
!       il faut au moins un sd par proc hors proc0
        if (((nbsd-dist0) .lt. (nbproc-1)) .and. (dist0 .gt. 0)) then
            call utmess('F', 'PARTITION1_99')
        end if
        if ((nbsd .lt. nbproc) .and. (dist0 .eq. 0)) then
            vali(1) = nbsd
            vali(2) = nbproc
            call utmess('F', 'PARTITION1_1', ni=2, vali=vali)
        end if
    else if (kdis(1:5) .eq. 'MAIL_') then
!       il faut au moins une maille par proc
        if (nbmamo .lt. nbproc) then
            vali(1) = nbmamo
            vali(2) = nbproc
            call utmess('F', 'PARTITION1_93', ni=2, vali=vali)
        end if
    end if

! ----------------------------------------------------------------------
!   Remplissage de la sd
! ----------------------------------------------------------------------

    if (kdis .eq. 'MAIL_DISPERSE') then
!   ---------------------------------------
!       -- le proc 0 a une charge differente des autres (dist0) :
!       nmpp nbre de mailles par proc (a la louche)
        nmpp = max(1, nbmamo/nbproc)
!       nmp0 nbre de mailles affectees au proc0 (a la louche)
        nmp0 = (dist0*nmpp)/100

!       -- affectation des mailles aux differents procs :
        nmp0af = 0
        ico = 0
        nbpro1 = nbproc
        plein0 = .false.
        do ima = 1, nbma
            if (maille(ima) .eq. 0) goto 40
            ico = ico+1
            krang = mod(ico, nbpro1)
            if (plein0) krang = krang+1
            if (krang .eq. 0) nmp0af = nmp0af+1
            zi(jnumsd-1+ima) = krang
            if (nmp0af .eq. nmp0) then
                plein0 = .true.
                nbpro1 = nbproc-1
            end if
40          continue
        end do

    else if (kdis .eq. 'MAIL_CONTIGU') then
!   --------------------------------------
!       nmp0 nbre de mailles affectees au proc0 :
        nmpp = max(1, nbmamo/nbproc)
        nmp0 = (dist0*nmpp)/100
        nmp1 = ((nbmamo-nmp0)/(nbproc-1))+1

!       -- affectation des mailles aux differents procs :
!          on affecte les 1eres mailles au proc0 puis les autres
!          aux autres procs.
        nmpp = nmp0
        krang = 0
        ico = 0
        do ima = 1, nbma
            if (maille(ima) .eq. 0) goto 50
            ico = ico+1
!         -- on change de proc :
            if (ico .gt. nmpp) then
                ico = 1
                nmpp = nmp1
                krang = krang+1
            end if
            zi(jnumsd-1+ima) = krang
50          continue
        end do

!       -- on verifie que toutes les mailles sont distribuees :
        ico = 0
        icobis = 0
        do i = 1, nbma
            if (zi(jnumsd-1+i) .ge. 0) ico = ico+1
            if (zi(jnumsd-1+i) .eq. rang) icobis = icobis+1
        end do
        ASSERT(ico .eq. nbmamo)

    else if (kdis .eq. 'GROUP_ELEM' .or. kdis .eq. 'SOUS_DOMAINE') then
!   ----------------------------------------------------------------
!       -- il n'y a rien a faire !
!       La regle pour les calculs elementaires et les assemblages est :
!       quelque soit le ligrel (modele, charge, ....) :
!       le grel igrel est traite par le processeur
!       de rang=mod(igrel,nbproc)

    else
        ASSERT(.false.)
    end if

999 continue

    call jedema()
end subroutine
