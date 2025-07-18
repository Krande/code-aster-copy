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
subroutine patrma(llist1, llist2, t, nbtymx, nomma, &
                  llistt, ntypm)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/pacoor.h"
#include "asterfort/padtma.h"
#include "asterfort/panbno.h"
#include "asterfort/parotr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    character(len=*) :: llist1, llist2, llistt
    character(len=8) :: nomma
    real(kind=8) :: t(3)
    integer(kind=8) :: nbtymx, ntypm
!     BUT: TRIER DE 2 LISTES DE MAILLES GENEREES PAR LA ROUTINE PALIMA
!     A FIN DE GENERER UN OJB METTANT EN VIS A VIS LES 2 LISTES VIA
!     LA TRANSLATION  T
!     LLISTT : OJB XD V I NUM()
!     LLISTT(I) : VECTEUR D'ENTIER DE DIM=3+NBMA*(2+2*NBNTOT)
!          V(1) = NUMTMC NUMERO DU TYPE_MAIL DES MAILLES DE COUPLAGE
!          V(2) = NBMA : NBRE DE COUPLAGE
!          V(3) = NBNTOT : NOMBRE DE NOEUDS DES MAILLES EN VIS A VIS
!          POUR IC = 1,NBMA
!    V(3+(IC-1)*(2+2*NBNTOT)+1)= NUMA1 NUM. DE LA 1ERE MAILLE DU COUPLE
!    V(3+(IC-1)*(2+2*NBNTOT)+2)= NUMA2 NUM. DE LA 2EME MAILLE DU COUPLE
!                        EN VIS A VIS AVEC NUMA1 VIA LA TRANSFORMATION
!           POUR INO = 1,NBNTOT
! V(3+(IC-1)*(2+2*NBNTOT)+2+2*(INO-1)+1)=N1(INO)NUM.DU NO. INO DE NUMA1
! V(3+(IC-1)*(2+2*NBNTOT)+2+2*(INO-1)+1)= N2(N1(INO)) NUM.DU NOEUD DE
!                   NUMA2 EN VIS A VIS AVEC LE NOEUD N1(INO)
!
! ARGUMENTS D'ENTREE:
! IN   LLIST1 K24 : NOM DE LA 1ERE LISTE
! IN   LLIST2 K24 : NOM DE LA 2EME LISTE
! IN   NOMMA  K8  : NOM DU MAILLAGE
! IN   T      R(3): COORDONNEES DE LA TRANSLATION
! IN   NBTYMX I   : NBRE MAX DE TYPE_MAILLE DIFFERENT (< OU = NBGREL)
! OUT  LLISTT K24 : NOM DE LA LISTE TRIEE RESULTAT
! OUT  NTYPM  I   : NBRE DE TYPE_MAIL DIFFERENT
!
! ROUTINES APPELEES:
    aster_logical :: fintyp
    integer(kind=8) :: nbnott(3)
    integer(kind=8) :: vali
    character(len=8) :: nomma1, nomma2, nomma3
    character(len=24) :: l1, l2, lt, coor1, coor2, couple, wcpl, connex, biject
    character(len=24) :: valk(3)
    real(kind=8) :: centre(3)
    real(kind=8) :: mrot(3, 3)
!-----------------------------------------------------------------------
    integer(kind=8) :: i1, iageom, idbij, idcoo1, idcoo2, idcopl, idl1
    integer(kind=8) :: idl2, idlt, idno1, idno2, idtyp, idwcpl, iftyp
    integer(kind=8) :: ima, ima2, ino, iret, ityp, ityp1, j
    integer(kind=8) :: j2, jdeb, jfin, jtyp, k, lonmx, nbma, err, tot_error
    integer(kind=8) :: nbma1, nbma2, nbmaty, nbntot, numa1, numa2, numa3
!
    real(kind=8) :: d, dmin
!-----------------------------------------------------------------------
    data centre/0.d0, 0.d0, 0.d0/
    data mrot/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
! --- DEBUT
    call jemarq()
    l1 = llist1
    l2 = llist2
    lt = llistt
    call jeexin(lt, iret)
    if (iret .ne. 0) call jedetr(lt)
    connex = nomma//'.CONNEX'
    coor1 = '&&PATRMA.COOR1'
    coor2 = '&&PATRMA.COOR2'
    wcpl = '&&PATRMA.WCOUPLE'
    couple = '&&PATRMA.COUPLE'
    biject = '&&PATRMA.BIJECT'
    call jeveuo(nomma//'.COORDO    .VALE', 'L', iageom)
!
!
!
    call jeveuo(l1, 'L', idl1)
    nbma1 = zi(idl1)
    call jeveuo(l2, 'L', idl2)
    nbma2 = zi(idl2)
    if (nbma1 .ne. nbma2) then
        valk(1) = l1
        valk(2) = l2
        call utmess('F', 'MODELISA6_22', nk=2, valk=valk)
    end if
    nbma = nbma1
    call wkvect(biject, 'V V I', nbma, idbij)
!
!     CREATION DE LLISTT
!
    call jecrec(lt, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbtymx)
    fintyp = .false.
    ntypm = 0
    jdeb = 1
    jfin = nbma
    idtyp = 1
    iftyp = nbma
    ityp1 = -1
    do i1 = 1, nbma
        numa1 = zi(idl1+2*i1-1)
        nomma1 = int_to_char8(numa1)
        ityp = zi(idl1+2*i1)
        if (ityp .ne. ityp1) then
            if (ityp1 .ne. -1) then
                call jedetr(coor1)
                call jedetr(coor2)
                call jedetr(wcpl)
                call jedetr(couple)
            end if
            ima = 0
            ntypm = ntypm+1
            ityp1 = ityp
            call panbno(ityp1, nbnott)
            nbntot = nbnott(1)+nbnott(2)+nbnott(3)
! --- NUM DU TYPE_MAILLE DE COUPLAGE : NUMTMC
!
            jdeb = 1
            jfin = nbma
            fintyp = .false.
            idtyp = 1
            iftyp = nbma
            call wkvect(coor1, 'V V R', 3*nbnott(1), idcoo1)
            call wkvect(coor2, 'V V R', 3*nbnott(1), idcoo2)
            call wkvect(wcpl, 'V V I', nbntot, idwcpl)
            call wkvect(couple, 'V V I', nbntot, idcopl)
        end if
        call parotr(nomma, iageom, numa1, nbnott(1), centre, &
                    mrot, t, zr(idcoo1))
        dmin = 999999999999.d0
        tot_error = 0
        do j = jdeb, jfin
            err = 0
            jtyp = zi(idl2+2*j)
            if (jtyp .ne. ityp1) then
                if (fintyp) then
                    iftyp = j-1
                    goto 21
                else
                    idtyp = j+1
                    goto 2
                end if
            else
                fintyp = .true.
                ima2 = zi(idl2+2*j-1)
                nomma2 = int_to_char8(ima2)
                call pacoor(nomma, ima2, nbnott(1), zr(idcoo2))
                call padtma(zr(idcoo1), zr(idcoo2), nbnott, zi(idwcpl), d, &
                            err)
                tot_error = tot_error+err
                if (d .lt. dmin .and. err .eq. 0) then
                    dmin = d
                    numa2 = ima2
                    j2 = j
                    do ino = 1, nbntot
                        zi(idcopl-1+ino) = zi(idwcpl-1+ino)
                    end do
                end if
            end if
2           continue
        end do
        if (tot_error .eq. jfin) then
            call utmess('F', 'MODELISA6_8')
        end if
        if (idtyp .eq. nbma+1) then
            vali = ityp
            call utmess('F', 'MODELISA8_89', si=vali)
        end if
21      continue
        ima = ima+1
        jdeb = idtyp
        jfin = iftyp
        nbmaty = jfin-jdeb+1
        fintyp = .false.
        if (ima .eq. 1) then
            if (ntypm .gt. 1) call jelibe(jexnum(lt, ntypm-1))
            call jecroc(jexnum(lt, ntypm))
            lonmx = 3+nbmaty*(2+2*nbntot)
            call jeecra(jexnum(lt, ntypm), 'LONMAX', lonmx)
            call jeveuo(jexnum(lt, ntypm), 'E', idlt)
            zi(idlt) = ityp1
            zi(idlt+1) = nbmaty
            zi(idlt+2) = nbntot
            idlt = idlt+2
        end if
        zi(idlt+1) = numa1
        zi(idlt+2) = numa2
        if (zi(idbij-1+j2) .eq. 0) then
            zi(idbij-1+j2) = numa1
        else
            nomma1 = int_to_char8(numa1)
            nomma2 = int_to_char8(numa2)
            numa3 = zi(idbij+j2-1)
            nomma3 = int_to_char8(numa3)
            valk(1) = nomma1
            valk(2) = nomma3
            valk(3) = nomma2
            call utmess('F', 'MODELISA8_90', nk=3, valk=valk)
        end if
        call jeveuo(jexnum(connex, numa1), 'L', idno1)
        call jeveuo(jexnum(connex, numa2), 'L', idno2)
        do k = 1, nbntot
            zi(idlt+2*k+1) = zi(idno1-1+k)
            zi(idlt+2*k+2) = zi(idno2-1+zi(idcopl-1+k))
        end do
        idlt = idlt+2+2*nbntot
    end do
    call jelibe(jexnum(lt, ntypm))
    call jedetr(biject)
    call jedetr(coor1)
    call jedetr(coor2)
    call jedetr(wcpl)
    call jedetr(couple)
    call jedema()
end subroutine
