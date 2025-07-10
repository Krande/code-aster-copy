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
subroutine gmeelt(imod, nbtyma, nomail, nbnoma, nuconn, &
                  nbmail, nbgrou)
    implicit none
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/codnop.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: imod, nbtyma, nbmail, nbnoma(19), nuconn(19, 32)
    character(len=8) :: nomail(*)
    integer(kind=8), intent(in) :: nbgrou
!
!      GMEELT --   ECRITURE DES MAILLES ET DES GROUP_MA VENANT
!                  D'UN FICHIER .GMSH DANS LE FICHIER .MAIL
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NBTYMA         IN    I         NOMBRE  DE TYPES DE MAILLES
!    NOMAIL(*)      IN    K8        TABLEAU DES NOMS DES TYPES DE MAILLE
!    NBMAIL         IN    I         NOMBRE TOTAL DE MAILLES
!    NUCONN         IN    I         PASSAGE DE LA NUMEROTATION DES NDS
!                                     D'UNE MAILLE : ASTER -> GMSH
!    NBGROU         IN    I         NOMBRE DE GROUPES
!
!
!
!
    integer(kind=8) :: neu2(32), ier, i, ij, nte, ima, ityp, nbno, inum, nbnoas
    integer(kind=8) :: idiv, ino, irest, k, l, maxmai, jgrmai, jgr, ima1
    integer(kind=8) :: vali(2)
    character(len=1) :: prfnoe, prfmai
    character(len=8) :: chgrou, chtab(32), chmail, k8bid
    character(len=12) :: chenti
    integer(kind=8), pointer :: nbnma(:) => null()
    integer(kind=8), pointer :: nbtym(:) => null()
    integer(kind=8), pointer :: noma(:) => null()
    integer(kind=8), pointer :: numa(:) => null()
    integer(kind=8), pointer :: typma(:) => null()
    integer(kind=8), pointer :: indma(:) => null()
    integer(kind=8), pointer :: nbmag(:) => null()
    character(len=8), pointer :: nomgr(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    prfmai = 'M'
    prfnoe = 'N'
    chgrou = '        '
    chmail = '        '
    k8bid = '        '
    chenti = 'NBOBJ=      '
!
    do i = 1, 32
        chtab(i) = '        '
    end do
!
! --- RECUPERATION DES OBJETS DE TRAVAIL :
!     ----------------------------------
    call jeveuo('&&PREGMS.NUMERO.MAILLES', 'L', vi=numa)
    call jeveuo('&&PREGMS.TYPE.MAILLES', 'L', vi=typma)
    call jeveuo('&&PREGMS.NBNO.MAILLES', 'L', vi=nbnma)
    call jeveuo('&&PREGMS.CONNEC.MAILLES', 'L', vi=noma)
    call jeveuo('&&PREGMS.NBTYP.MAILLES', 'L', vi=nbtym)
    if (nbgrou > 0) then
        call jeveuo('&&PREGMS.NBMA.GROUP_MA', 'L', vi=nbmag)
        call jeveuo('&&PREGMS.NUMERO.GROUP_MA', 'L', vi=indma)
        call jeveuo('&&PREGMS.NOMS.GROUP_MA', 'L', vk8=nomgr)
    end if
!
! --- ECRITURE DES MAILLES :
!     --------------------
    do nte = 1, nbtyma
!
        if (nbtym(nte) .eq. 0) goto 20
        call codent(nbtym(nte), 'G', chenti(7:12))
!
! ---   ECRITURE DE LA DATE :
!       -------------------
        write (unit=imod, fmt='(A,3X,A,3X,A)') nomail(nte),&
     &    'NOM=INDEFINI', chenti
!
        ij = 0
!
! --- ON VERIFIE QUE LE NOMBRE MAX DE MAILLES N'EST PAS ATTEINT
!     LA LIMITE EST DE 9 999 999 MAILLES
!
        if (nbmail .ge. 10000000) then
            vali(1) = nbmail
            call utmess('E', 'PREPOST6_43', si=vali(1))
        end if
!
! ---   BOUCLE SUR LES MAILLES :
!       ----------------------
        do ima = 1, nbmail
            ityp = typma(ima)
            nbno = nbnma(ima)
            if (ityp .eq. nte) then
                inum = numa(ima)
                call codnop(chmail, prfmai, 1, 1)
                call codent(inum, 'G', chmail(2:8))
!
                nbnoas = nbnoma(nte)
!
                do ino = 1, nbnoas
                    neu2(ino) = noma(1+ij+nuconn(nte, ino)-1)
                    call codnop(chtab(ino), prfnoe, 1, 1)
                    call codent(neu2(ino), 'G', chtab(ino) (2:8))
                end do
!
                idiv = int(nbnoas/8)
                irest = mod(nbnoas, 8)
!
                if (irest .ne. 0) then
                    write (imod, 202) chmail, (chtab(i), i=1, nbnoas)
                else
                    do k = 1, idiv
                        l = 8*(k-1)
                        if (idiv .eq. 1) then
                            write (imod, '(A,8(1X,A))') chmail, (chtab(i) &
                                                                 , i=1+l, 8+l)
                        else
                            write (imod, '(8X,8(1X,A))') (chtab(i), i=1+ &
                                                          l, 8+l)
                        end if
                    end do
                end if
            end if
            ij = ij+nbno
        end do
!
        write (imod, '(A)') 'FINSF'
        write (imod, '(A)') '%'
!
20      continue
    end do
!
202 format(a, 8(1x, a), /, (8x, 8(1x, a)))
!
! --- ECRITURE DES GROUP_MA :
!     ---------------------
    ier = 0
!
    maxmai = 0
    do i = 1, nbgrou
        maxmai = max(maxmai, nbmag(i))
    end do
!
! --- SI IL N Y A AU MOINS UN GROUPE :
!     ------------------------------
    if (maxmai .ne. 0) then
!
        call wkvect('&&PREGMS.GRMA.MAILLES', 'V V K8', maxmai, jgrmai)
!
! --- BOUCLE SUR LES GROUPES DE MAILLES :
!     ---------------------------------
        do i = 1, nbgrou
            chgrou = nomgr(i)
            write (imod, '(A,4X,2A)') 'GROUP_MA', 'NOM=', chgrou
            call jeveuo(jexnum('&&PREGMS.LISTE.GROUP_MA', i), 'E', jgr)
            do k = 1, nbmag(i)
                call codnop(chmail, prfmai, 1, 1)
                call codent(zi(jgr+k-1), 'G', chmail(2:8))
                zk8(jgrmai+k-1) = chmail
            end do
!
! ---   ECRITURE DES MAILLES DU GROUPE DE MAILLES COURANT :
!       -------------------------------------------------
            write (imod, '(8(2X,A))') (zk8(jgrmai+k-1), k=1, nbmag(1+i- &
                                                                   1))
!
            write (imod, '(A)') 'FINSF'
            write (imod, '(A)') '%'
!
! --- DANS LE CAS D'UN POINT ECRITURE D'UN GROUPNO
! ---  LE GROUPE DE MAILLE CONTIENT ALORS UNE SEULE MAILLE POI1
!
            if (nbmag(i) .eq. 1) then
                call jeveuo(jexnum('&&PREGMS.LISTE.GROUP_MA', i), 'E', jgr)
                ima1 = zi(jgr)
                ij = 0
                do ima = 1, nbmail
                    inum = numa(ima)
                    nbno = nbnma(ima)
                    if (inum .eq. ima1) then
                        if (nbno .eq. 1) then
                            write (imod, '(A,4X,2A)') 'GROUP_NO', 'NOM=', &
                                chgrou
                            neu2(ino) = noma(ij+1)
                            call codnop(chtab(ino), prfnoe, 1, 1)
                            call codent(neu2(ino), 'G', chtab(ino) (2:8))
                            write (imod, '((2X,A))') chtab(ino)
                            write (imod, '(A)') 'FINSF'
                            write (imod, '(A)') '%'
                            goto 90
                        end if
                    end if
                    ij = ij+nbno
                end do
            end if
90          continue
        end do
!
    end if
!
    if (ier .ne. 0) then
        call utmess('F', 'PREPOST_60')
    end if
!
    call jedema()
!
!.============================ FIN DE LA ROUTINE ======================
end subroutine
