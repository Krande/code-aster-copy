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

subroutine vtcmbl(nbcmb, typcst, const, typech, nomch, &
                  typres, chpres, basez)
!     ------------------------------------------------------------------
!     COMBINAISON LINEAIRE DE CHAM_NO OU DE CHAM_ELEM
!     *  LES CHAM_NOS OU CHAM_ELEMS SONT REELS OU COMPLEXES
!     *  LES SCALAIRES SONT REELS OU COMPLEXES
!     -----------------------------------------------------------------
! IN  : NBCOMB : IS  : NOMBRE DE CHAM_GDS A COMBINER
! IN  : TYPCST : K1  : TYPE DES CONSTANTES (R OU C, OU I )
! IN  : CONST  : R8  : TABLEAU DES COEFFICIENTS
! IN  : TYPECH : K1  : TYPE DES CHAM_GDS   (R OU C)
! IN  : NOMCH  : K19 : NOMS DES CHAM_GDS
! IN  : TYPRES : K1  : TYPE DU CHAMP RESULTAT (R OU C)
! IN  : CHPRES : K19 : NOM DU CHAMP RESULTAT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/dismoi.h"
#include "asterfort/copisd.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdchgd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8), intent(in) :: nbcmb
    real(kind=8), intent(in) :: const(*)
    character(len=*), intent(in) :: typcst(*), typech(*), nomch(*), typres, chpres
    character(len=1), intent(in), optional :: basez
!
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: ifm, niv, i
    integer(kind=8) :: icmb, iret, ival, lvale, iconst
    integer(kind=8) :: jdesc, jrefe, jvale
    integer(kind=8) :: kdesc, krefe, kvale
    integer(kind=8) :: nbdesc, nbrefe, nbvale
    integer(kind=8) :: nbdes1, nbref1, nbval1
    real(kind=8) :: dimag
    complex(kind=8) :: c8cst
    character(len=1) :: base
    character(len=4) :: docu, type, type2
    character(len=5) :: refe, desc, vale
    character(len=19) :: ch19, ch19r, nume_equa
    character(len=24) :: noojb
!     ------------------------------------------------------------------
!
    call jemarq()
! RECUPERATION ET MAJ DU NIVEAU D'IMPRESSION
    call infniv(ifm, niv)
!
!-----------------------------------------------------------------------
! --- PHASE D'INITIALISATION
!-----------------------------------------------------------------------
    type = typres
    ch19 = nomch(1)
!
    if (present(basez)) then
        base = basez
    else
        base = 'V'
    end if
!
! CHAM_NO OU CHAM_ELEM ?
    call dismoi("DOCU", ch19, "CHAMP", repk=docu)
!
!
! INIT. DE BASE
    if (docu .eq. 'CHNO') then
        refe = '.REFE'
        vale = '.VALE'
    else if (docu .eq. 'CHML') then
        refe = '.CELK'
        desc = '.CELD'
        vale = '.CELV'
    else
        call utmess('F', 'UTILITAI_21')
    end if
!
!
!
!   PREMIER CHAM_NO A CONCATENER
    ch19 = nomch(1)
!
!   OBTENTION DES ADRESSES ET DES TAILLES DES .DESC, .REFE ET .VALE
!   DU PREMIER CHAM_NO A CONCATENER. ON SUPPOSE QUE
!   TOUS LES CHAM_NOS DE LA LISTE NOMCH SONT HOMOGENES SUR CE POINT.
    call jelira(ch19//vale, 'LONMAX', nbvale)
    call jelira(ch19//refe, 'LONMAX', nbrefe)
    call jeveuo(ch19//refe, 'L', jrefe)

    if (docu == 'CHML') then
        call jelira(ch19//desc, 'LONMAX', nbdesc)
        call jeveuo(ch19//desc, 'L', jdesc)
    else
        nbdesc = 0
    end if

!
!   CONSTRUCTION D'UN CHAM_GD RESULTAT SUR LE MODELE DE NOMCH(1)
    ch19r = chpres
    call jeexin(ch19r//vale, iret)
!
    if (iret .eq. 0) then
        if (docu == 'CHML') then
            call wkvect(ch19r//desc, base//' V I', nbdesc, kdesc)
        end if
        call wkvect(ch19r//vale, base//' V '//type, nbvale, kvale)
        call wkvect(ch19r//refe, base//' V K24', nbrefe, krefe)
    else
        if (docu == 'CHML') then
            call jeveuo(ch19r//desc, 'E', kdesc)
            call jelira(ch19r//desc, 'LONMAX', nbdes1)
        else
            nbdes1 = 0
        end if
        call jeveuo(ch19r//refe, 'E', krefe)
        call jelira(ch19r//refe, 'LONMAX', nbref1)
        call jelira(ch19r//vale, 'LONMAX', nbval1)
!       VERIFICATION DE LA COHERENCE DES DIMENSIONS
        ASSERT(nbdes1 .eq. nbdesc)
        ASSERT(nbref1 .eq. nbrefe)
        ASSERT(nbval1 .eq. nbvale)
    end if
!   RECOPIE DU .DESC ET DU .REFE DU PREMIER CHAM_NO DE LA LISTE
!   DANS CEUX DU CHAM_NO SOLUTION
    do i = 0, nbdesc-1
        zi(kdesc+i) = zi(jdesc+i)
    end do
    do i = 0, nbrefe-1
        zk24(krefe+i) = zk24(jrefe+i)
    end do

    if (docu == "CHNO") then
        call dismoi("TYPE_SCA", ch19r, "CHAM_NO", repk=type2)
        if (type2 .ne. type) then
            noojb = '12345678.NUME000000.PRNO'
            call gnomsd(ch19r, noojb, 14, 19)
            noojb(1:8) = ch19r(1:8)
            nume_equa = noojb(1:19)
            call copisd("NUME_EQUA", base, zk24(krefe-1+2), nume_equa)
            zk24(krefe-1+2) = nume_equa
            call jeveuo(nume_equa//".REFN", 'E', krefe)
            zk24(krefe-1+2) = zk24(krefe-1+2) (1:5)//type(1:1)
        end if
    end if
!
!   CHANGER LA GRANDEUR
    call sdchgd(ch19r, type)
    if (docu == 'CHML') then
        call jeecra(ch19r//desc, 'DOCU', cval=docu)
    else
        call jeecra(ch19r//refe, 'DOCU', cval=docu)
    end if
!
!   VECTEUR RECEPTACLE TEMPORAIRE DE LA COMBINAISON LINEAIRE
    call wkvect('&&VTCMBL.VALE', 'V V '//type, nbvale, lvale)
!
!
!-----------------------------------------------------------------------
! --- BOUCLE SUR LES CHAM_GDS A COMBINER
!-----------------------------------------------------------------------
    iconst = 1
    do icmb = 1, nbcmb
!
!       CHAM_NO A CONCATENER
        ch19 = nomch(icmb)
!
        call jeveuo(ch19//vale, 'L', jvale)
        if (typres(1:1) .eq. 'R') then
            if (typech(icmb) (1:1) .eq. 'R') then
                do ival = 0, nbvale-1
                    zr(lvale+ival) = zr(lvale+ival)+const( &
                                     iconst)*zr(jvale+ival)
                end do
            else
                if (typcst(icmb) (1:1) .eq. 'R') then
                    do ival = 0, nbvale-1
                        zr(lvale+ival) = zr(lvale+ival)+ &
                                         const(iconst)*dble(zc(jvale+ival))
                    end do
                else if (typcst(icmb) (1:1) .eq. 'I') then
                    do ival = 0, nbvale-1
                        zr(lvale+ival) = zr(lvale+ival)+ &
                                         const(iconst)*dimag(zc(jvale+ival))
                    end do
                else
                    type = typcst(icmb) (1:1)
                    call utmess('F', 'PREPOST3_6', sk=type)
                end if
            end if
        else
            if (typech(icmb) (1:1) .eq. 'R') then
                if (typcst(icmb) (1:1) .eq. 'R') then
                    do ival = 0, nbvale-1
                        zc(lvale+ival) = zc(lvale+ival)+ &
                                         const(iconst)*zr(jvale+ival)
                    end do
                else if (typcst(icmb) (1:1) .eq. 'C') then
                    c8cst = dcmplx(const(iconst), const(iconst+1) &
                                   )
                    do ival = 0, nbvale-1
                        zc(lvale+ival) = zc(lvale+ival)+c8cst* &
                                         zr(jvale+ival)
                    end do
                end if
            else
                if (typcst(icmb) (1:1) .eq. 'R') then
                    do ival = 0, nbvale-1
                        zc(lvale+ival) = zc(lvale+ival)+ &
                                         const(iconst)*zc(jvale+ival)
                    end do
                else if (typcst(icmb) (:1) .eq. 'C') then
                    c8cst = dcmplx(const(iconst), const(iconst+1) &
                                   )
                    do ival = 0, nbvale-1
                        zc(lvale+ival) = zc(lvale+ival)+c8cst* &
                                         zc(jvale+ival)
                    end do
                end if
            end if
        end if
        call jelibe(ch19//vale)
        iconst = iconst+1
        if (typcst(icmb) (1:1) .eq. 'C') iconst = iconst+1
    end do
!
!
!   IL EST NECESSAIRE D'ACTUALISER KVALE SI LE RESULTAT EST DANS NOMCH()
    call jeveuo(ch19r//vale, 'E', kvale)
    if (type(1:1) .eq. 'R') then
        do ival = 0, nbvale-1
            zr(kvale+ival) = zr(lvale+ival)
        end do
    else if (type(1:1) .eq. 'C') then
        do ival = 0, nbvale-1
            zc(kvale+ival) = zc(lvale+ival)
        end do
    end if
!
    call jedetr('&&VTCMBL.VALE')
    call jelibe(ch19r//vale)
!
    call jedema()
end subroutine
