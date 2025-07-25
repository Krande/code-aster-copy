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

subroutine assvss(base, vec, vecel, nu, vecpro, &
                  motcle, type, fomult, instap)
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/cordd2.h"
#include "asterfort/crelil.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/fointe.h"
#include "asterfort/infniv.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/ssvalv.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: vec, vecpro, base, nu
    character(len=19) :: vecel
    character(len=4) :: motcle
    character(len=24) :: fomult
    integer(kind=8) :: type
    real(kind=8) :: instap
! ----------------------------------------------------------------------
! OUT K19 VEC   : NOM DU CHAM_NO RESULTAT
!                CHAM_NO ::= CHAM_NO_GD + OBJETS PROVISOIRES POUR L'ASS.
! IN  K* BASE   : NOM DE LA BASE SUR LAQUELLE ON VEUT CREER LE CHAM_NO
! IN  K* VECEL  : VECT_ELEM A ASSEMBLER
! IN  K* NU     : NOM D'UN NUMERO_DDL
! IN  K* VECPRO : NOM D'UN CHAM_NO MODELE(NU OU VECPRO EST OBLIGATOIRE)
! IN  K4 MOTCLE : 'ZERO' OU 'CUMU'
! IN  K24 FOMULT: TABLEAU DE FONCTIONS MULTIPLICATRICES DE CHARGES
! IN  R8 INSTAP : INSTANT D'INTERPOLATION
! IN  I  TYPE   : TYPE DU VECTEUR ASSEMBLE : 1 --> REEL
!                                            2 --> COMPLEXE
!----------------------------------------------------------------------
    character(len=8) :: nomacr, exiele
    character(len=14) :: num2
    integer(kind=8) :: gd, nec, nlili
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, iad1, iadlie, iadnem, iadval
    integer(kind=8) :: ialcha, iamail, iancmp, ianueq, ianulo, iaprol
    integer(kind=8) :: iapsdl, ichar, icmp, iconx1, iconx2, idnequ
    integer(kind=8) :: idprn1, idprn2, idresl, idverf, iec
    integer(kind=8) :: ierd, il, ilim, ilimnu, ilivec, ima
    integer(kind=8) :: inold, iret, jec, k1, lgncmp, n1, nbchar
    integer(kind=8) :: nbecmx, nbelm, nbnoss, nbsma, nbssa, ncmp, ncmpel
    integer(kind=8) :: nddl1, nequa, nm, nmxcmp, nnoe, nugd
!-----------------------------------------------------------------------
    parameter(nbecmx=10)
!
    character(len=1) :: bas
    character(len=8) :: ma, mo, mo2, nogdsi, nogdco, nomcas
    character(len=14) :: nudev
    character(len=19) :: vecas, vprof
    character(len=24) :: knueq, kmaila, k24prn
    character(len=24) :: kvelil, kveref, knequa, kvale
    integer(kind=8) :: icodla(nbecmx), icodge(nbecmx)
    integer(kind=8) :: admodl, lcmodl, ifm, niv
    integer(kind=8) :: jfonct
    real(kind=8) :: rcoef
    character(len=24), pointer :: refe(:) => null()
    integer(kind=8), pointer :: sssa(:) => null()
    character(len=8), pointer :: vnomacr(:) => null()
    integer(kind=8), pointer :: conx(:) => null()
!
! --- DEBUT ------------------------------------------------------------
    call jemarq()
!
    call infniv(ifm, niv)
!
    if (motcle(1:4) .eq. 'ZERO') then
    else if (motcle(1:4) .eq. 'CUMU') then
    else
        call utmess('F', 'ASSEMBLA_8', sk=motcle)
    end if
!
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', lcmodl)
    call jeveuo(jexnum('&CATA.TE.MODELOC', 1), 'L', admodl)
!
    vecas = vec
    bas = base
!
! --- SI LE CONCEPT VECAS EXISTE DEJA,ON LE DETRUIT:
    call detrsd('CHAMP_GD', vecas)
    call wkvect(vecas//'.LIVE', bas//' V K24 ', 1, ilivec)
    zk24(ilivec) = vecel
!
! --- NOMS DES PRINCIPAUX OBJETS JEVEUX LIES A VECAS
    kmaila = '&MAILLA                 '
    kvelil = vecas//'.LILI'
    call dismoi('NOM_MODELE', nu, 'NUME_DDL', repk=mo)
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
!
!
! --- CALCUL D UN LILI POUR VECAS
! --- CREATION D'UN VECAS(1:19).ADNE ET VECAS(1:19).ADLI SUR 'V'
    call crelil('F', 1, zk24(ilivec), kvelil, 'V', &
                kmaila, vecas, gd, ma, nec, &
                ncmp, ilim, nlili, nbelm)
    call jeveuo(vecas(1:19)//'.ADLI', 'E', iadlie)
    call jeveuo(vecas(1:19)//'.ADNE', 'E', iadnem)
    call jeexin(ma(1:8)//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(ma(1:8)//'.CONNEX', 'L', iconx1)
        call jeveuo(jexatr(ma(1:8)//'.CONNEX', 'LONCUM'), 'L', iconx2)
    end if
!
! --- ON SUPPOSE QUE LE LE LIGREL DE &MAILLA EST LE PREMIER DE LILINU
    ilimnu = 1
!
! --- NOMS DES PRINCIPAUX OBJETS JEVEUX LIES A NU
! --- IL FAUT ESPERER QUE LE CHAM_NO EST EN INDIRECTION AVEC UN
!     NUME_EQUA APPARTENANT A UNE NUMEROTATION SINON CA VA PLANTER
!     DANS LE JEVEUO SUR KNEQUA
    nudev = nu
    if (nudev(1:1) .eq. ' ') then
        vprof = vecpro
        call jeveuo(vprof//'.REFE', 'L', vk24=refe)
        nudev = refe(2) (1:14)
    end if
!
!
    call dismoi('NOM_MODELE', nudev, 'NUME_DDL', repk=mo)
    call dismoi('NOM_MAILLA', nudev, 'NUME_DDL', repk=ma)
    call dismoi('NB_NO_SS_MAX', ma, 'MAILLAGE', repi=nbnoss)
!
!     100 EST SUPPOSE ETRE LA + GDE DIMENSION D'UNE MAILLE STANDARD:
    nbnoss = max(nbnoss, 100)
!     -- NUMLOC(K,INO) (K=1,3)(INO=1,NBNO(MAILLE))
    call wkvect('&&ASSVEC.NUMLOC', 'V V I', 3*nbnoss, ianulo)
!
    call dismoi('NOM_GD', nudev, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NB_CMP_MAX', nogdsi, 'GRANDEUR', repi=nmxcmp)
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
    nec = nbec(nugd)
    ncmp = nmxcmp
!
    do i = 1, nbecmx
        icodla(i) = 0
        icodge(i) = 0
    end do
!
!   -- POSDDL(ICMP) (ICMP=1,NMXCMP(GD_SI))
    call wkvect('&&ASSVEC.POSDDL', 'V V I', nmxcmp, iapsdl)
!
    call dismoi('NB_NO_MAILLA', mo, 'MODELE', repi=nm)
!
    call jeexin(ma//'.NOMACR', iret)
    if (iret .gt. 0) then
        call jeveuo(ma//'.NOMACR', 'L', vk8=vnomacr)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nogdsi), 'L', iancmp)
        call jelira(jexnom('&CATA.GD.NOMCMP', nogdsi), 'LONMAX', lgncmp)
        icmp = indik8(zk8(iancmp), 'LAGR', 1, lgncmp)
! on ne trouve pas la composante "LAGR" dans la grandeur
        ASSERT(icmp .ne. 0)
! il est imprévu d avoir la composante "LAGR" au delà de 30
        ASSERT(icmp .le. 30)
!       -- icodla est l'entier code correspondant a la cmp "lagr"
        jec = (icmp-1)/30+1
        icodla(jec) = 2**icmp
    end if
!
    k24prn = nudev//'.NUME.PRNO'
    knueq = nudev//'.NUME.NUEQ'
    knequa = nudev//'.NUME.NEQU'
!
    call jeveuo(k24prn, 'L', idprn1)
    call jeveuo(jexatr(k24prn, 'LONCUM'), 'L', idprn2)
    call jeveuo(knueq, 'L', ianueq)
    call jeveuo(knequa, 'L', idnequ)
    nequa = zi(idnequ)
!
!
    kveref = vecas//'.REFE'
    kvale = vecas//'.VALE'
!
    call jecreo(kveref, bas//' V K24')
    call jeecra(kveref, 'LONMAX', 4)
    call jeveuo(kveref, 'E', idverf)
    call jeecra(kveref, 'DOCU', cval='CHNO')
    zk24(idverf+1) = k24prn(1:14)//'.NUME'
!
!
    if (type .eq. 1) then
        call jecreo(kvale, bas//' V R8')
    else if (type .eq. 2) then
        call jecreo(kvale, bas//' V C16')
    else
        call utmess('F', 'ASSEMBLA_11')
    end if
    call jeecra(kvale, 'LONMAX', nequa)
    call jeveuo(kvale, 'E', iadval)
!
!
    call dismoi('NOM_MODELE', vecel, 'VECT_ELEM', repk=mo2)
    if (mo2 .ne. mo) then
        call utmess('F', 'ASSEMBLA_5')
    end if
!
    call dismoi('EXI_ELEM', mo, 'MODELE', repk=exiele)
    call dismoi('NB_SS_ACTI', vecel, 'VECT_ELEM', repi=nbssa)
!
!
!   -- TRAITEMENT DES SOUS-STRUCTURES
!   ----------------------------------------------------------
    if (nbssa .gt. 0) then
        nomcas = ' '
        call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
        call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
        call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
        call ssvalv('DEBUT', nomcas, mo, ma, 0, &
                    idresl, ncmpel, instap)
        call jelira(vecel//'.RELC', 'NUTIOC', nbchar)
        call jeveuo(fomult, 'L', jfonct)
!
        do ichar = 1, nbchar
            call jenuno(jexnum(vecel//'.RELC', ichar), nomcas)
            call jeveuo(jexnum(vecel//'.RELC', ichar), 'L', ialcha)
            if (zk24(jfonct+ichar-1) (1:8) .eq. '&&CONSTA') then
                rcoef = 1.0d0
            else
                call fointe('F ', zk24(jfonct+ichar-1) (1:8), 1, ['INST'], [instap], &
                            rcoef, ierd)
            end if
            do ima = 1, nbsma
!               -- ON N'ASSEMBLE QUE LES SSS VRAIMENT ACTIVES :
                if (sssa(ima) .eq. 0) goto 70
                if (zi(ialcha-1+ima) .eq. 0) goto 70
                call jeveuo(jexnum(ma//'.SUPMAIL', ima), 'L', iamail)
                call jelira(jexnum(ma//'.SUPMAIL', ima), 'LONMAX', nnoe)
                call ssvalv(' ', nomcas, mo, ma, ima, &
                            idresl, ncmpel, instap)
                nomacr = vnomacr(ima)
                call dismoi('NOM_NUME_DDL', nomacr, 'MACR_ELEM_STAT', repk=num2)
                call jeveuo(nomacr//'.CONX', 'L', vi=conx)
                call jeveuo(jexnum(num2//'.NUME.PRNO', 1), 'L', iaprol)
                il = 0
                do k1 = 1, nnoe
                    n1 = zi(iamail-1+k1)
                    if (n1 .gt. nm) then
                        do iec = 1, nbecmx
                            icodge(iec) = icodla(iec)
                        end do
                    else
                        inold = conx(3*(k1-1)+2)
                        do iec = 1, nec
                            icodge(iec) = zi(iaprol-1+(nec+2)*( &
                                             inold-1)+2+iec)
                        end do
                    end if
!
                    iad1 = zi(idprn1-1+zi(idprn2+ilimnu-1)+(n1- &
                                                            1)*(nec+2))
                    call cordd2(idprn1, idprn2, ilimnu, icodge, nec, &
                                ncmp, n1, nddl1, zi(iapsdl))
!
                    if (type .eq. 1) then
                        do i1 = 1, nddl1
                            il = il+1
                            zr(iadval-1+zi(ianueq-1+iad1+zi(iapsdl-1+i1)-1)) = &
                                zr(iadval-1+zi(ianueq-1+iad1+zi(iapsdl-1+i1)-1))+ &
                                zr(idresl+il-1)*rcoef
                        end do
                    else if (type .eq. 2) then
                        do i1 = 1, nddl1
                            il = il+1
                            zc(iadval-1+zi(ianueq-1+iad1+zi(iapsdl-1+i1)-1)) = &
                                zc(iadval-1+zi(ianueq-1+iad1+zi(iapsdl-1+i1)-1))+ &
                                zc(idresl+il-1)*rcoef
                        end do
                    end if
                end do
70              continue
            end do
        end do
        call ssvalv('FIN', nomcas, mo, ma, 0, &
                    idresl, ncmpel, instap)
    end if
!
!
    call jedetr(vecas//'.LILI')
    call jedetr(vecas//'.LIVE')
    call jedetr(vecas//'.ADNE')
    call jedetr(vecas//'.ADLI')
    call jedetr('&&ASSVEC.POSDDL')
    call jedetr('&&ASSVEC.NUMLOC')
    call jedema()
end subroutine
