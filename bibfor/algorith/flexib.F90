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
subroutine flexib(basmod, nbmod, flex, nl, nc, &
                  numl, numc)
    implicit none
!  P. RICHARD     DATE 09/04/91
!-----------------------------------------------------------------------
!  BUT : CALCULER LA MATRICE DE FLEXIBILITE RESIDUELLE ASSOCIEE
!        A UN PROBLEME CYCLIQUE AVEC INTERFACE MAC NEAL OU AUCUN
!        (FLEXIBILITE NULLE DANS LE CAS AUCUN)
!
!        SEULE LA SOUS MATRICE RELATIVE AUX DEFORMEES (COLONNES) D'UNE
!        INTERFACE ET AUX DDL D'UNE AUTRE (LIGNES) EST CALCULEE
!
!        POUR LES LIGNES IL EST POSSIBLE DE NE PAS DONNER UNE INTERFACE
!        MAIS DE PRENDRE TOUTES LES LIGNES ( = TOUS LES DDL PHYSIQUES)
!        IL SUFFIT POUR CELA DE DONNER UN NUMERO D'INTERFACE NUML = 0
!-----------------------------------------------------------------------
!
! BASMOD /I/ : NOM UTILISATEUR DE LA BASE MODALE
! NBMOD  /I/ : NOMBRE DE MODES PROPRES UTILISES
! FLEX   /O/ : MATRICE DE FLEXIBILITE RESIDUELLE
! NL     /I/ : NOMBRE DE LIGNES DE LA MATRICE DE FLEXIBILITE
! NC     /I/ : NOMBRE DE COLONNES DE LA MATRICE DE FLEXIBILITE
! NUML   /I/ : NUMERO DE L'INTERFACE DE DDL RELATIFS AUX LIGNES
! NUMC   /I/ : NUMERO DE L'INTERFACE DE DDL RELATIFS AUX COLONNES
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/bmnodi.h"
#include "asterfort/bmradi.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/zerlag.h"
#include "blas/dcopy.h"
!
!
    integer(kind=8) :: nl, nc
    real(kind=8) :: flex(nl, nc)
    character(len=6) :: pgc
    character(len=8) :: basmod, typint, intf, kbid, k8bid
    character(len=19) :: numddl
    character(len=24) :: chamva, noeint
    character(len=24) :: valk
!
!----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, iord, iran, j
    integer(kind=8) :: jj, k, kk, ldkge, ldmge, llcham, lldes
    integer(kind=8) :: llnoc, llnol, ltextc, ltextl, ltorc
    integer(kind=8) :: ltvec, nbmod, nbnoc, nbnol, nbnot, neq
    integer(kind=8) :: numc, numl
    real(kind=8) :: toto, xkgen, xx
    integer(kind=8), pointer :: deeq(:) => null()
    character(len=8), pointer :: idc_type(:) => null()
    blas_int :: b_incx, b_incy, b_n
!-----------------------------------------------------------------------
    data pgc/'FLEXIB'/
!-----------------------------------------------------------------------
!
    call jemarq()
    do i = 1, nl
        do j = 1, nc
            flex(i, j) = 0.d0
        end do
    end do
!
! --- RECUPERATION CONCEPTS AMONT
!
    call dismoi('REF_INTD_PREM', basmod, 'RESU_DYNA', repk=intf, arret='C')
    call dismoi('NUME_DDL', basmod, 'RESU_DYNA', repk=numddl)
!
    if (intf .eq. '        ') then
        valk = basmod
        call utmess('F', 'ALGORITH13_17', sk=valk)
    end if
!
! --- TEST SUR LE TYPE D'INTERFACE
!
    call jeveuo(intf//'.IDC_TYPE', 'L', vk8=idc_type)
    typint = idc_type(numc)
    call jelibe(intf//'.IDC_TYPE')
    if (typint .eq. 'AUCUN   ') goto 999
!
! --- ALLOCATION DES TABLEAUX DE TRAVAIL
!
    call wkvect('&&'//pgc//'.ORDREC', 'V V I', nc, ltorc)
    call wkvect('&&'//pgc//'.EXTRACC', 'V V I', nc, ltextc)
    call wkvect('&&'//pgc//'.EXTRACL', 'V V I', nl, ltextl)
!
! --- RECUPERATION DU NOMBRE DE DDL PHYSIQUES ASSEMBLES
!
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!----ON AJOUT .NUME POUR OBTENIR LE NUME_EQUA
    numddl(15:19) = '.NUME'
    call jeveuo(numddl//'.DEEQ', 'L', vi=deeq)
!
! --- RECUPERATION DU NOMBRE DE NOEUDS DES INTERFACES
!
    noeint = intf//'.IDC_LINO'
!
    if (numl .gt. 0) then
        call jelira(jexnum(noeint, numl), 'LONMAX', nbnol)
        call jeveuo(jexnum(noeint, numl), 'L', llnol)
    else
        nbnol = 0
    end if
!
    call jelira(jexnum(noeint, numc), 'LONMAX', nbnoc)
    call jeveuo(jexnum(noeint, numc), 'L', llnoc)
!
! --- RECUPERATION DU DESCRIPTEUR DES DEFORMEES
!
    call jeveuo(intf//'.IDC_DEFO', 'L', lldes)
    call jelira(intf//'.IDC_DEFO', 'LONMAX', nbnot)
    nbnot = nbnot/3
!
! --- RECUPERATION DES NUMEROS D'ORDRE DES DEFORMEES (COLONNES)
!     ET RANGS DES DDL D'INTERFACE (LIGNES) DANS VECTEUR ASSEMBLE
!
! --- RECUPERATION NUMERO ORDRE DEFORMEES ET RANG DDL POUR COLONNES
!
    kbid = ' '
    call bmnodi(basmod, kbid, '        ', numc, nc, &
                zi(ltorc), ibid)
    call bmradi(basmod, kbid, '        ', numc, nc, &
                zi(ltextc), ibid)
!
! --- RECUPERATION DDL PHYSIQUES POUR LES LIGNES
!
    if (numl .gt. 0) then
        call bmradi(basmod, kbid, '        ', numl, nl, &
                    zi(ltextl), ibid)
    else
        do i = 1, neq
            zi(ltextl+i-1) = i
        end do
    end if
!
    if (numl .gt. 0) then
        call jelibe(jexnum(noeint, numl))
    end if
    call jelibe(jexnum(noeint, numc))
    call jelibe(intf//'.IDC_DEFO')
!
! --- EXTRACTION PARTIE INTERFACE DE FLEXIBILITE
!
    do i = 1, nc
        call wkvect('&&'//pgc//'.VECT', 'V V R', neq, ltvec)
        iord = zi(ltorc+i-1)
        call dcapno(basmod, 'DEPL    ', iord, chamva)
        call jeveuo(chamva, 'L', llcham)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(llcham), b_incx, zr(ltvec), b_incy)
        call zerlag(neq, deeq, vectr=zr(ltvec))
!
        do j = 1, nl
!
!  EXTRACTION DDL
!
            iran = zi(ltextl+j-1)
            xx = zr(ltvec+iran-1)
            flex(j, i) = xx
        end do
        call jelibe(chamva)
        call jedetr('&&'//pgc//'.VECT')
    end do
!
! --- SUPPRESSION CONTRIBUTION STATIQUE DES MODES CONNUS
!
    do i = 1, nbmod
!
        call rsadpa(basmod, 'L', 1, 'RIGI_GENE', i, &
                    0, sjv=ldkge, styp=k8bid)
        xkgen = zr(ldkge)
        call rsadpa(basmod, 'L', 1, 'MASS_GENE', i, &
                    0, sjv=ldmge, styp=k8bid)
!
        call dcapno(basmod, 'DEPL    ', i, chamva)
        call jeveuo(chamva, 'L', llcham)
        call wkvect('&&'//pgc//'.VECT', 'V V R', neq, ltvec)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(llcham), b_incx, zr(ltvec), b_incy)
        call zerlag(neq, deeq, vectr=zr(ltvec))
!
        do j = 1, nc
            do k = 1, nl
                jj = zi(ltextc+j-1)
                kk = zi(ltextl+k-1)
                toto = zr(ltvec+jj-1)*zr(ltvec+kk-1)/xkgen
                flex(k, j) = flex(k, j)-toto
            end do
        end do
        call jelibe(chamva)
        call jedetr('&&'//pgc//'.VECT')
    end do
!
    call jedetr('&&'//pgc//'.ORDREC')
    call jedetr('&&'//pgc//'.EXTRACC')
    call jedetr('&&'//pgc//'.EXTRACL')
!
999 continue
    call jedema()
end subroutine
