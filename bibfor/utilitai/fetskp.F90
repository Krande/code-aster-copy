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

subroutine fetskp(mod, meth, nbpart)
!
!    - FONCTION REALISEE:
!       - CREATION DU GRAPHE D'ENTREE DU PARTITIONNEUR
!       - APPEL A METIS OU EXECUTION DE SCOTCH
!       - CREATION DE NOUVEAUX GROUPES DE MAILLES
!----------------------------------------------------------------------
! person_in_charge: jacques.pellet at edf.fr
!
    implicit none
    character(len=8), intent(in) :: mod, meth
    integer(kind=8), intent(in) :: nbpart
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/fetsco.h"
#include "asterc/gpmetis_aster.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/creaco.h"
#include "asterfort/creagm.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbmama, idco, nbmato, renum2, nbma, nomsdm, masd
    integer(kind=8) :: nbmasd, id, co, renum
    integer(kind=8) :: numsdm, nmap, i, ima
    integer(kind=8) :: iocc, nocc, ifm, niv, nblien, renum3, idma
    integer(kind=8) :: rang, nbproc, versco, n1, n2, n3, ier, iaux, iaux2
    integer(kind=8) :: iexi, numsd0
    real(kind=8) :: tmps(7)
    character(len=8) :: ma
    character(len=8) :: kersco
    character(len=24) :: k24b, liel
    character(len=19) :: ligrel
    integer(kind=4), pointer :: vedlo(:) => null()
    integer(kind=8), pointer :: vrenum1(:) => null()
    integer(kind=4), pointer :: vvelo(:) => null()
    mpi_int :: mrank, msize
!
!---------------------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
!
    nbproc = 1
    rang = 0
    call asmpi_info(rank=mrank, size=msize)
    rang = to_aster_int(mrank)
    nbproc = to_aster_int(msize)
!
! ********************************************************************
!                       CREATION DU GRAPHE
!
! ------- ON RECUPERE LES DONNEES DU MAILLAGE OU DU MODELE
!
    call dismoi('NOM_MAILLA', mod, 'MODELE', repk=ma)
    ASSERT(.not. isParallelMesh(ma))

    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbmato)
    AS_ALLOCATE(vi=vrenum1, size=nbmato)

    nbmato = 0

    call dismoi('NOM_LIGREL', mod, 'MODELE', repk=ligrel)
    liel = ligrel//".LIEL"
    call jeexin(liel, iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'PARTITION_2')
    end if
    call jelira(liel, 'NMAXOC', nocc)
    do iocc = 1, nocc
        call jelira(jexnum(liel, iocc), 'LONMAX', nbma)
        nbmato = nbmato+nbma-1
    end do
    call wkvect('&&FETSKP.RENUM', 'V V I', nbmato, renum)
    id = 1
    do iocc = 1, nocc
        call jelira(jexnum(liel, iocc), 'LONMAX', nbma)
        call jeveuo(jexnum(liel, iocc), 'L', idma)
        do ima = 1, nbma-1
            zi(renum-1+id) = zi(idma-1+ima)
! ----- ON VERIFIE QUE LE MODELE NE CONTIENT PAS DE MAILLES TARDIVES
! ----- QUI SONT ESSENTIELLEMENT DES NOEUDS A CE STADE
            if (zi(idma-1+ima) .lt. 0) then
                call utmess('F', 'PARTITION_3')
            end if
            vrenum1(zi(idma-1+ima)) = id
            id = id+1
        end do
    end do
!
! ------- CREATION DE LA CONNECTIVITE DES MAILLES
!
    call creaco(nbmato, ma, nblien)
    if (nblien .eq. 0) call utmess('F', 'PARTITION_4')
!
! ------ ON RECUPERE LES TABLEAUX CONSTRUITS DANS CREACO
!
    call jeveuo('&&FETSKP.RENUM2', 'L', renum2)
    call jeveuo('&&FETSKP.RENUM3', 'L', renum3)
    call jeveuo('&&FETSKP.CO', 'L', co)
    call jeveuo('&&FETSKP.IDCO', 'L', idco)
    call jeveuo('&&FETSKP.NBMAMA', 'L', nbmama)
!
! -------  UTILISATION DE CONTRAINTES
!
    AS_ALLOCATE(vi4=vvelo, size=nbmato)
    AS_ALLOCATE(vi4=vedlo, size=nblien)
    do ima = 1, nbmato
        vvelo(ima) = 1
    end do
    do i = 1, nblien
        vedlo(i) = 1
    end do
!
!!
    AS_DEALLOCATE(vi=vrenum1)
    call jedetr('&&FETSKP.NBMAMA')
!
!
! ********************************************************************
!                       LANCEMENT DU LOGICIEL
!
!
!     ************** LANCEMENT DE SCOTCH
!
    if (meth(1:6) .eq. 'SCOTCH') call wkvect('&&FETSKP.NMAP', 'V V S', nbmato, nmap)
    if ((meth(1:6) .eq. 'SCOTCH') .and. (rang .eq. 0)) then
        if (niv .ge. 2) then
            call uttcpu('CPU.FETSKP', 'INIT', ' ')
            call uttcpu('CPU.FETSKP', 'DEBUT', ' ')
        end if
        write (ifm, *) ' '
        write (ifm, *) '***************** SCOTCH *****************'
        write (ifm, *) ' '
        write (ifm, *) ' '
        write (ifm, *) ' * LE NOMBRE DE MAILLES    :', nbmato
        write (ifm, *) ' * LE NOMBRE DE CONNEXIONS :', nblien
        write (ifm, *) ' '
        call fetsco(nbmato, nblien, zi4(co), zi4(idco), nbpart, &
                    zi4(nmap), vedlo(1), vvelo(1), versco, ier)
        n1 = versco/10000
        n2 = (versco-n1*10000)/100
        n3 = versco-n1*10000-n2*100
        kersco(1:8) = '........'
        write (kersco(1:2), '(I2)') n1
        write (kersco(4:5), '(I2)') n2
        write (kersco(7:8), '(I2)') n3
        if (niv .ge. 2) then
            call uttcpu('CPU.FETSKP', 'FIN', ' ')
            call uttcpr('CPU.FETSKP', 7, tmps)
            write (ifm, *) ' * TEMPS DE PARTITIONNEMENT  :', tmps(7)
            write (ifm, *) ' '
        end if
        write (ifm, *) '********** FIN SCOTCH ', kersco, ' *********'
        if (ier .ne. 0) then
            call utmess('F', 'UTILITAI_56', si=ier)
        end if
        write (ifm, *) ' '
!
!     ************** LANCEMENT DE METIS
!
    else if (rang .eq. 0) then
        call wkvect('&&FETSKP.TMPNUMSDM', 'V V I', nbmato, numsd0)
        if (meth .eq. 'METIS  ') then
            call gpmetis_aster(nbmato, nblien, zi4(idco), zi4(co), nbpart, zi(numsd0))
        end if
    end if
!
    AS_DEALLOCATE(vi4=vedlo)
    AS_DEALLOCATE(vi4=vvelo)
!
!
! ********************************************************************
!                    CREATION DES GROUPES DE MAILLES
!
!
    call wkvect('&&FETSKP.NUMSDM', 'V V I', nbmato, numsdm)
    call wkvect('&&FETSKP.NBMASD', 'V V I', nbpart, nbmasd)
!
! ------- LECTURE DU RESULTAT DU PARTITONNEUR
!
    if (meth(1:6) .ne. 'SCOTCH') then
        if (rang .eq. 0) then
            do ima = 1, nbmato
                zi(numsdm-1+zi(renum2-1+ima)) = zi(numsd0-1+ima)
            end do
            call jedetr('&&FETSKP.TMPNUMSDM')
        end if
        k24b = '&&FETSKP.NUMSDM'
        call asmpi_comm_jev('BCAST', k24b)
        do ima = 1, nbmato
            iaux = zi(renum2-1+ima)
            iaux2 = zi(numsdm-1+iaux)
            zi(nbmasd+iaux2) = zi(nbmasd+iaux2)+1
        end do
    else
        k24b = '&&FETSKP.NMAP'
        call asmpi_comm_jev('BCAST', k24b)
        do ima = 1, nbmato
            zi(numsdm-1+zi(renum2-1+ima)) = zi4(nmap-1+ima)
            zi(nbmasd+zi(numsdm-1+zi(renum2-1+ima))) = zi(nbmasd+zi( &
                                                          numsdm-1+zi(renum2-1+ima)))+1
        end do
        call jedetr(k24b)
    end if
!
! ------- CREATION DES GROUP_MA
!
    call creagm(nbmato, nbpart, ma, masd)
!
    call jeveuo('&&FETSKP.NOMSDM', 'L', nomsdm)
    if (niv .gt. 1) then
        do i = 1, nbpart
            write (ifm, *) 'LE SOUS DOMAINE ', zk24(nomsdm-1+i), ' CONTIENT ' &
                , zi(nbmasd-1+i), ' MAILLES '
        end do
    end if

    call jedetr('&&FETSKP.TEMP')
    call jedetr('&&FETSKP.RENUM2')
    call jedetr('&&FETSKP.RENUM3')
    call jedetr('&&FETSKP.IDMASD')
    call jedetr('&&FETSKP.NOMSDM')
    call jedetr('&&FETSKP.MASD')
    call jedetr('&&FETSKP.NBMASD')
    call jedetr('&&FETSKP.NUMSDM')
    call jedetr('&&FETSKP.RENUM')
    call jedetr('&&FETSKP.CO')
    call jedetr('&&FETSKP.IDCO')

    call jedema()
end subroutine
