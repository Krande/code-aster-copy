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

subroutine verecy(intf, numd, numg, nbsec, prec, &
                  distrf)
!    P. RICHARD     DATE 13/12/91
!-----------------------------------------------------------------------
!  BUT:       < VERIFICATION REPETITIVITE CYCLIQUE>
! aslint: disable=
    implicit none
!
!  VERIFICATION DE LA REPETITIVITE CYCLIQUE SUR LE MAILLAGE ET LA
!  DEFINITION DES INTERFACES
!
!-----------------------------------------------------------------------
!
! INTF     /I/: NOM UTILISATEUR DE L'INTERF_DYNA
! NUMD     /I/: NUMERO DE L'INTERFACE DE DROITE
! NUMG     /I/: NUMERO DE L'INTERFACE DE GAUCHE
! NBSEC    /I/: NOMBRE DE SECTEUR
! PREC     /R/: PRECISION DE RECHERCHE DE PROXIMITE
! DISTRF   /R/: DISTANCE DE REFERENCE
!
!
#include "asterf_types.h"
#include "asterc/r8pi.h"
#include "jeveux.h"
#include "asterfort/bmnoin.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: vali(2)
!
!
!
    character(len=6) :: pgc
    character(len=24) :: valk(3)
    character(len=8) :: intf, kbid, mailla, nomnod, nomnog, nomnj
    character(len=50) :: diag
    aster_logical :: ordre
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, j, jnode, llintg
    integer(kind=8) :: llista, llistb, ltnd, ltng, nbd, nbg, nbpbax
    integer(kind=8) :: nbpbr, nbpbse, nbpbto, nbpbvt, nbsec, numd, numg
    integer(kind=8) :: nunod, nunog
    real(kind=8) :: crit, difr, difz, dist, distj, distr, distrf
    real(kind=8) :: distrj, distz, distzj, pi, prec, pvdif, rd
    real(kind=8) :: rg, teta, xd, xg, yd, yg, zd
    real(kind=8) :: zg, zpv, zpvref
    real(kind=8), pointer :: vale(:) => null()
!-----------------------------------------------------------------------
    data pgc/'VERECY'/
!-----------------------------------------------------------------------
!
    distrj = 0.d0
    distzj = 0.d0
!
    call jemarq()
    pi = r8pi()
!
!--------VERIFICATION NOMBRE DE NOEUDS INTERFACES DROITE ET GAUCHE------
!
    kbid = ' '
    call bmnoin(' ', intf, kbid, numd, 0, &
                [0], nbd)
    kbid = ' '
    call bmnoin(' ', intf, kbid, numg, 0, &
                [0], nbg)
!
!
    if (nbg .ne. nbd) then
        vali(1) = nbd
        vali(2) = nbg
        call utmess('E', 'ALGORITH16_50', ni=2, vali=vali)
    end if
!
!
!--------------------VERIFICATION REPETITIVITE GEOMETRIQUE--------------
!
!
    call dismoi('NOM_MAILLA', intf, 'INTERF_DYNA', repk=mailla)
!
!
!
    call wkvect('&&'//pgc//'.NOEUD.DROITE', 'V V I', nbd, ltnd)
    call wkvect('&&'//pgc//'.NOEUD.GAUCHE', 'V V I', nbg, ltng)
!
    kbid = ' '
    call bmnoin(' ', intf, kbid, numd, nbd, &
                zi(ltnd), ibid)
    kbid = ' '
    call bmnoin(' ', intf, kbid, numg, nbg, &
                zi(ltng), ibid)
!
    call jeveuo(mailla//'.COORDO    .VALE', 'L', vr=vale)
!
    teta = 2.d0*pi/nbsec
!
!     --- CONSTITUTION DE LISTA ET LISTB :
!         LE IEME NOEUD DE L'INTERFACE DROITE A POUR VIS-A-VIS
!         LE ZI(LISTA-1+I) EME NOEUD DE L'INTERFACE GAUCHE
!         RECIPROQUEMENT LE NOEUD DE POSITION J DE L'INTERFACE GAUCHE
!         EST LE VIS-A-VIS DU NOEUD DE POSITION ZI(LISTB-1+J) DE
!         L'INTERFACE DROITE.
    call wkvect('&&'//pgc//'.LISTA', 'V V I', nbd, llista)
    call wkvect('&&'//pgc//'.LISTB', 'V V I', nbd, llistb)
    nbpbax = 0
    nbpbr = 0
    nbpbse = 0
    nbpbvt = 0
    ordre = .true.
    do i = 1, nbd
!     --- BOUCLE SUR LES NOEUDS DE L'INTERFACE DROITE ---
        nunod = zi(ltnd+i-1)
        nomnod = int_to_char8(nunod)
!
        xd = vale(1+3*(nunod-1))
        yd = vale(1+3*(nunod-1)+1)
        zd = vale(1+3*(nunod-1)+2)
        rd = sqrt(xd*xd+yd*yd)
!
!       RECHERCHE DU NOEUD J (GAUCHE) LE PLUS PROCHE DE I (DROITE)
        do j = 1, nbd
!       --- BOUCLE SUR LES NOEUDS DE L'INTERFACE GAUCHE ---
            nunog = zi(ltng+j-1)
            nomnog = int_to_char8(nunog)
            xg = vale(1+3*(nunog-1))
            yg = vale(1+3*(nunog-1)+1)
            zg = vale(1+3*(nunog-1)+2)
            rg = sqrt(xg*xg+yg*yg)
            distr = abs(rd-rg)
            distz = abs(zd-zg)
            if (j .eq. 1 .or. (distr .le. distrj .and. distz .le. distzj)) then
!          --- CRITERE : RAYON ET HAUTEUR Z LES PLUS PROCHES ---
                distrj = distr
                distzj = distz
                distj = sqrt(distr*distr+distz*distz)
                jnode = j
                nomnj = nomnog
            else if (distr .le. distrj .or. distz .le. distzj) then
!          --- SI UN SEUL CRITERE EST BON, ON COMPARE LES DISTANCES ---
                dist = sqrt(distr*distr+distz*distz)
                if (dist .lt. distj) then
!
                    distrj = distr
                    distzj = distz
                    distj = dist
                    jnode = j
                    nomnj = nomnog
                end if
            end if
        end do
        zi(llista-1+i) = jnode
        if (zi(llistb-1+jnode) .ne. 0) then
!       --- CAS OU JNODE EST DEJA UN VIS-A-VIS ---
            nunog = zi(ltng+zi(llistb-1+jnode)-1)
            nomnog = int_to_char8(nunog)
            valk(1) = nomnj
            valk(2) = nomnod
            valk(3) = nomnog
            call utmess('F', 'ALGORITH16_51', nk=3, valk=valk)
        end if
        zi(llistb-1+jnode) = i
!       SI JNODE EST DIFFERENT DE I, C'EST QUE LES NOEUDS D'INTERFACE
!       ONT ETE DONNES DANS UN ORDRE DE NON CORRESPONDANCE
        if (jnode .ne. i) ordre = .false.
        nunog = zi(ltng+jnode-1)
        nomnog = int_to_char8(nunog)
        xg = vale(1+3*(nunog-1))
        yg = vale(1+3*(nunog-1)+1)
        zg = vale(1+3*(nunog-1)+2)
!
! VERIFICATION OZ AXE REPETITIVITE
!
        difz = abs(zd-zg)
        if (distrf .lt. 0.d0) then
!       --- DISTANCE DE REFERENCE NON CONNUE
            crit = prec*1.d-2*max(abs(zd), abs(zg))
        else
            crit = prec*distrf
        end if
        if (difz .gt. crit) then
            nbpbax = nbpbax+1
            vali(1) = i
            valk(1) = nomnod
            valk(2) = nomnog
            call utmess('E', 'ALGORITH16_52', nk=2, valk=valk, si=vali(1))
        end if
!
!      VERIFICATION RAYON
!
        rd = ((xd**2)+(yd**2))**0.5d0
        rg = ((xg**2)+(yg**2))**0.5d0
!
        difr = abs(rd-rg)
        crit = prec*distrf
        if (distrf .lt. 0.d0) then
!       --- DISTANCE DE REFERENCE NON CONNUE
            crit = prec*1.d-2*max(rd, rg)
        else
            crit = prec*distrf
        end if
        if (difr .gt. crit) then
            nbpbr = nbpbr+1
            vali(1) = i
            valk(1) = nomnod
            valk(2) = nomnog
            call utmess('E', 'ALGORITH16_53', nk=2, valk=valk, si=vali(1))
        end if
!
!  VERIFICATION SENS ANGLE
!
        zpv = (xd*yg)-(yd*xg)
        if (zpv .lt. 0.d0) then
            nbpbse = nbpbse+1
            vali(1) = i
            valk(1) = nomnod
            valk(2) = nomnog
            call utmess('E', 'ALGORITH16_54', nk=2, valk=valk, si=vali(1))
        end if
!
! VERIFICATION VALEUR ANGLE
!
        zpvref = (sin(teta)*rd*rg)
        pvdif = abs(zpvref-abs(zpv))
        crit = zpvref*prec
        if (pvdif .gt. crit) then
            nbpbvt = nbpbvt+1
            vali(1) = i
            valk(1) = nomnod
            valk(2) = nomnog
            call utmess('E', 'ALGORITH16_55', nk=2, valk=valk, si=vali(1))
        end if
!
    end do
!
!
    nbpbto = nbpbax+nbpbr+nbpbse+nbpbvt
!
    if (nbpbto .eq. 0) then
        call utmess('I', 'ALGORITH16_56')
        diag = ' '
    else if (nbpbax .eq. nbd) then
        diag = ' AXE DE REPETITIVITE DIFFERENT DE 0Z      '
    else if (nbpbvt .eq. nbd) then
        diag = 'NOMBRE DE SECTEURS DONNE ERRONE          '
    else if (nbpbse .eq. nbd) then
        diag = 'INVERSION INTERFACE DROITE ET GAUCHE'
    else if (nbpbr .eq. nbpbto) then
        diag = 'INTERFACES DROITE ET GAUCHE NON COMPATIBLES'
    else
        diag = ' PAS DE DIAGNOSTIC SIMPLE TROUVE'
    end if
!
    if (.not. ordre .and. diag .eq. ' ') then
!     --- LES NOEUDS NE SONT PAS EN VIS-A-VIS ---
!         ON REGARDE D'ABORD SI LE TRI EST PLAUSIBLE
        do i = 1, nbd
            if (zi(llistb-1+zi(llista-1+i)) .ne. i) then
                diag = 'TRI DES NOEUDS IMPOSSIBLE'
                goto 40
            end if
        end do
40      continue
!
        call utmess('A', 'ALGORITH16_57')
        call jeveuo(jexnum(intf//'.IDC_LINO', numg), 'E', llintg)
!    --- ON ORDONNE LES NOEUDS DE LLINTG SUIVANT LLISTA
        do i = 1, nbd
!        --- RECOPIE DE LLINT2 DANS LLISTB
            zi(llistb-1+i) = zi(llintg-1+i)
        end do
        do i = 1, nbd
            zi(llintg-1+i) = zi(llistb-1+zi(llista-1+i))
        end do
!
    end if
!
!     --- DESTRUCTION OBJETS SUR VOLATILE
    call jedetr('&&VERECY.LISTA')
    call jedetr('&&VERECY.LISTB')
    call jedetr('&&VERECY.NOEUD.DROITE')
    call jedetr('&&VERECY.NOEUD.GAUCHE')
!
    if (diag .ne. ' ') then
        valk(1) = diag
        call utmess('F', 'ALGORITH16_58', sk=valk(1))
    end if
!
    call jedema()
end subroutine
