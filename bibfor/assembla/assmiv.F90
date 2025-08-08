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

subroutine assmiv(base, vec, nbvec, tlivec, licoef, &
                  nu, type)
    implicit none
!
!
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/asmpi_info.h"
#include "asterfort/assert.h"
#include "asterfort/corddl.h"
#include "asterfort/crelil.h"
#include "asterfort/detrsd.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jaexin.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/nbno.h"
#include "asterfort/parti0.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=1) :: base
    character(len=*) :: vec
    integer(kind=8) :: nbvec
    character(len=*) :: tlivec(nbvec)
    real(kind=8) :: licoef(nbvec)
    character(len=*) :: nu
    integer(kind=8) :: type
! ----------------------------------------------------------------------
!    assemblage "particulier" pour convergence en contraintes generalisees
!    realise le min des vect_elem
!
!    Cette routine est inspiree de la routine assvec
!
! OUT K19 VEC   : NOM DU CHAM_NO RESULTAT
!                CHAM_NO ::= CHAM_NO_GD + OBJETS PROVISOIRES POUR L'ASS.
! IN  K* BASE   : NOM DE LA BASE SUR LAQUELLE ON VEUT CREER LE CHAM_NO
! IN  I  NBVEC  : NOMBRE DE VECT_ELEM A ASSEMBLER DANS VEC
! IN  K* TLIVEC : LISTE DES VECT_ELEM A ASSEMBLER
! IN  R  LICOEF : LISTE DES COEF. MULTIPLICATEURS DES VECT_ELEM
! IN  K* NU     : NOM D'UN NUMERO_DDL
! IN  I  TYPE   : TYPE DU VECTEUR ASSEMBLE : 1 --> REEL
!
! ----------------------------------------------------------------------
    real(kind=8) :: rcoef, r
    character(len=24) :: valk(5)
    integer(kind=8) :: gd, nec, nlili
    integer(kind=8) :: rang, nbproc, iret, ifm, niv
    character(len=8) :: ma, mo, mo2, nogdsi, nogdco, partit
    character(len=14) :: nume_ddl
    character(len=19) :: vecas, vprof, vecel, resu
    character(len=24) :: kmaila, k24prn, knulil, kvelil, kveref, nomli
    character(len=24) :: knequa, kvale
    integer(kind=8) :: admodl, lcmodl, iexi
    aster_logical :: ldist, ldgrel
    integer(kind=8) :: vali(4)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, iad, iad1, iadnem
    integer(kind=8) :: ianulo, iconx2
    integer(kind=8) :: idprn1, idprn2, idverf, iel
    integer(kind=8) :: igr, il, ilim, ilimnu, ilinu, ilive
    integer(kind=8) :: ilivec, imat, iresu, jresl, jvale
    integer(kind=8) :: k1, mode, n1, nbelm, nbnoss, nbresu, ncmp, nb_equa, nb_dof
    integer(kind=8) :: ncmpel, nddl1, nel, nm, nmxcmp, nnoe
    integer(kind=8) :: nugd, numa
    integer(kind=8), pointer :: posddl(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    character(len=24), pointer :: relr(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    character(len=24), pointer :: prtk(:) => null()
    integer(kind=8), pointer :: nueq(:) => null()
    integer(kind=8), pointer :: v_nequ(:) => null()
    integer(kind=8), pointer :: maille(:) => null()
    integer(kind=8), pointer :: adli(:) => null()
    character(len=24), pointer :: refe(:) => null()
    mpi_int :: mrank, msize
!-----------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)

!
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', lcmodl)
    call jeveuo(jexnum('&CATA.TE.MODELOC', 1), 'L', admodl)
!
    vecas = vec
!
! ------------------------------------------------------------------
    ldist = .false.
    ldgrel = .false.
    rang = 0
    nbproc = 1
    nume_ddl = nu
    call dismoi('NOM_MODELE', nu, 'NUME_DDL', repk=mo)
    call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
!
    call parti0(nbvec, tlivec, partit)
    if (partit .ne. ' ') then
        ldist = .true.
        call jeveuo(partit//'.PRTK', 'L', vk24=prtk)
        ldgrel = prtk(1) .eq. 'SOUS_DOMAINE' .or. prtk(1) .eq. 'GROUP_ELEM'
        call asmpi_info(rank=mrank, size=msize)
        rang = to_aster_int(mrank)
        nbproc = to_aster_int(msize)
        if (.not. ldgrel) then
            call jeveuo(partit//'.NUPR', 'L', vi=maille)
        end if
    end if
!
!
! --- SI LE CONCEPT VECAS EXISTE DEJA,ON LE DETRUIT:
    call detrsd('CHAMP_GD', vecas)
    call wkvect(vecas//'.LIVE', base//' V K24 ', nbvec, ilivec)
    do i = 1, nbvec
        zk24(ilivec-1+i) = tlivec(i)
    end do
!
! --- NOMS DES PRINCIPAUX OBJETS JEVEUX LIES A VECAS
    kmaila = '&MAILLA                 '
    kvelil = vecas//'.LILI'
    kveref = vecas//'.REFE'
    kvale = vecas//'.VALE'
!
! --- CREATION DE REFE ET DESC
    call jecreo(kveref, base//' V K24')
    call jeecra(kveref, 'LONMAX', 4)
    call jeecra(kveref, 'LONUTI', 4)
    call jeveuo(kveref, 'E', idverf)
    call jeecra(kveref, 'DOCU', cval='CHNO')
!
! --- CALCUL D UN LILI POUR VECAS
! --- CREATION D'UN VECAS(1:19).ADNE ET VECAS(1:19).ADLI SUR 'V'
    call crelil('F', nbvec, tlivec, kvelil, 'V', &
                kmaila, vecas, gd, ma, nec, &
                ncmp, ilim, nlili, nbelm, nume_=nume_ddl)
!
    call jeveuo(vecas(1:19)//'.ADLI', 'E', vi=adli)
    call jeveuo(vecas(1:19)//'.ADNE', 'E', iadnem)
    call jeexin(ma(1:8)//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(ma(1:8)//'.CONNEX', 'L', vi=connex)
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
    if (nume_ddl(1:1) .eq. ' ') then
        vprof = ' '
        call jeveuo(vprof//'.REFE', 'L', vk24=refe)
        nume_ddl = refe(2) (1:14)
    end if
!
    knequa = nume_ddl//'.NUME.NEQU'
    k24prn = nume_ddl//'.NUME.PRNO'
    knulil = nume_ddl//'.NUME.LILI'
    call jeveuo(nume_ddl//'.NUME.NUEQ', 'L', vi=nueq)

!
! - Get number of equations
!
    call jeveuo(knequa, 'L', vi=v_nequ)
    nb_equa = v_nequ(1)
    nb_dof = v_nequ(2)

!
    call dismoi('NOM_MODELE', nume_ddl, 'NUME_DDL', repk=mo)
    call dismoi('NOM_MAILLA', nume_ddl, 'NUME_DDL', repk=ma)
    call dismoi('NB_NO_SS_MAX', ma, 'MAILLAGE', repi=nbnoss)
!
!     100 EST SUPPOSE ETRE LA + GDE DIMENSION D'UNE MAILLE STANDARD:
    nbnoss = max(nbnoss, 100)
!     -- NUMLOC(K,INO) (K=1,3)(INO=1,NBNO(MAILLE))
    call wkvect('&&ASSMIV.NUMLOC', 'V V I', 3*nbnoss, ianulo)
!
    call dismoi('NOM_GD', nume_ddl, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NB_CMP_MAX', nogdsi, 'GRANDEUR', repi=nmxcmp)
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
    nec = nbec(nugd)
    ncmp = nmxcmp
!
!
!     -- POSDDL(ICMP) (ICMP=1,NMXCMP(GD_SI))
    AS_ALLOCATE(vi=posddl, size=nmxcmp)
!
    call dismoi('NB_NO_MAILLA', mo, 'MODELE', repi=nm)
!
!
! ---  RECUPERATION DE PRNO
    call jeveuo(k24prn, 'L', idprn1)
    call jeveuo(jexatr(k24prn, 'LONCUM'), 'L', idprn2)
!
! ---  RECUPERATION DE NEQUA

!
! ---  REMPLISSAGE DE REFE
    zk24(idverf+1) = k24prn(1:14)//'.NUME'
!
!
! --- ALLOCATION VALE
    ASSERT(type .eq. 1)
    call wkvect(kvale, base//' V R8', nb_equa, jvale)
!
    do i = 1, nb_equa
        zr(jvale+i-1) = r8maem()
    end do
!
!
!   -- REMPLISSAGE DE .VALE
!   ------------------------
    do imat = 1, nbvec
        rcoef = licoef(imat)
        vecel = zk24(ilivec+imat-1) (1:19)
!
        call dismoi('NOM_MODELE', vecel, 'VECT_ELEM', repk=mo2)
        if (mo2 .ne. mo) then
            call utmess('F', 'ASSEMBLA_5')
        end if
!
        call jeexin(vecel//'.RELR', iret)
        if (iret .eq. 0) goto 90
!
        call jeveuo(vecel//'.RELR', 'L', vk24=relr)
        call jelira(vecel//'.RELR', 'LONUTI', nbresu)
        do iresu = 1, nbresu
            resu = relr(iresu) (1:19)
            call jeveuo(resu//'.NOLI', 'L', iad)
            nomli = zk24(iad)
!
            call jenonu(jexnom(kvelil, nomli), ilive)
            call jenonu(jexnom(knulil, nomli), ilinu)
!
            do igr = 1, adli(1+3*(ilive-1))
                if (ldgrel .and. mod(igr, nbproc) .ne. rang) goto 70
!
!             -- IL SE PEUT QUE LE GREL IGR SOIT VIDE :
                call jaexin(jexnum(resu//'.RESL', igr), iexi)
                if (iexi .eq. 0) goto 70
!
                call jeveuo(resu//'.DESC', 'L', vi=desc)
                mode = desc(1+igr+1)
                if (mode .gt. 0) then
                    nnoe = nbno(mode)
                    nel = zi(adli(1+3*(ilive-1)+2)+igr)- &
                          zi(adli(1+3*(ilive-1)+2)+igr-1)-1
                    call jeveuo(jexnum(resu//'.RESL', igr), 'L', jresl)
                    ncmpel = digdel(mode)
!
                    do iel = 1, nel
                        numa = zi(adli(1+3*(ilive-1)+1)-1+ &
                                  zi(adli(1+3*(ilive-1)+2)+igr-1)+iel-1)
                        r = rcoef
!
                        if (ldist .and. .not. ldgrel) then
                            if (numa .gt. 0) then
                                if (maille(numa) .ne. rang) goto 60
                            end if
                        end if
!
                        if (numa .gt. 0) then
                            il = 0
                            do k1 = 1, nnoe
                                n1 = connex(zi(iconx2+numa-1)+k1-1)
                                iad1 = zi(idprn1-1+zi(idprn2+ilimnu-1)+ &
                                          (n1-1)*(nec+2)+1-1)
                                call corddl(admodl, lcmodl, idprn1, idprn2, ilimnu, &
                                            mode, nec, ncmp, n1, k1, &
                                            nddl1, posddl)
                                if (nddl1 .eq. 0) goto 50
                                if (iad1 .eq. 0) then
                                    vali(1) = n1
                                    valk(1) = resu
                                    valk(2) = vecel
                                    valk(3) = nume_ddl
                                    call utmess('F', 'ASSEMBLA_41', nk=3, valk=valk, si=vali(1))
                                end if
!
                                if (iad1 .gt. nb_dof) then
                                    vali(1) = n1
                                    vali(2) = iad1
                                    vali(3) = nb_dof
                                    valk(1) = resu
                                    valk(2) = vecel
                                    call utmess('F', 'ASSEMBLA_42', nk=2, valk=valk, ni=3, &
                                                vali=vali)
                                end if
!
                                if (nddl1 .gt. 100) then
                                    vali(1) = nddl1
                                    vali(2) = 100
                                    call utmess('F', 'ASSEMBLA_43', ni=2, vali=vali)
                                end if
!
                                if (type .eq. 1) then
                                    do i1 = 1, nddl1
                                        il = il+1
                                        zr(jvale-1+nueq(iad1+posddl(i1)-1)) = &
                                            min(zr(jvale-1+nueq(iad1+posddl(i1)-1)), &
                                                zr(jresl+(iel-1)*ncmpel+il-1)*r)
                                    end do
                                end if
50                              continue
                            end do
                        end if
60                      continue
                    end do
                    call jelibe(jexnum(resu//'.RESL', igr))
                end if
70              continue
            end do
        end do
!
90      continue
    end do
    call jedetr(vecas//'.LILI')
    call jedetr(vecas//'.LIVE')
    call jedetr(vecas//'.ADNE')
    call jedetr(vecas//'.ADLI')
    AS_DEALLOCATE(vi=posddl)
    call jedetr('&&ASSMIV.NUMLOC')
!
!   -- REDUCTION + DIFFUSION DE VECAS A TOUS LES PROC
    if (ldist) then
        call asmpi_comm_vect('MPI_MIN', 'R', nbval=nb_equa, vr=zr(jvale))
    end if
!
!
    call jedema()
end subroutine
