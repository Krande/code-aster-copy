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

subroutine detrsd(typesd, nomsd)
!
    use elg_data_module
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/amumph.h"
#include "asterfort/apetsc.h"
#include "asterfort/assde1.h"
#include "asterfort/assert.h"
#include "asterfort/detlsp.h"
#include "asterfort/detrs2.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/codent.h"
#include "asterfort/asmpi_info.h"
!
!
    character(len=*) :: typesd, nomsd
!
! ----------------------------------------------------------------------
!  BUT : DETRUIRE UNE STRUCTURE DE DONNEE DONT ON CONNAIT LE TYPE
!  ATTENTION : QUAND ON UTILISE TYPESD=' ', ON APPELLE LA ROUTINE JEDETC
!              QUI EST TRES COUTEUSE EN CPU.
!  IN   TYPESD : TYPE DE LA STRUCTURE DE DONNEE A DETRUIRE
!          'NUME_DDL'     'MLTF'
!          'MATR_ASSE'    'VECT_ASSE'    'MATR_ASSE_GENE'
!          'MATR_ELEM'    'VECT_ELEM'   'PARTITION'
!          'VARI_COM'     'FONCTION' (POUR LES FONCTIONS OU NAPPES)
!          'TABLE_SDASTER' 'TABLE_CONTAINER'
!          'SOLVEUR'      'EIGENSOLVER'  'CORRESP_2_MAILLA'
!          'CHAM_NO_S'    'CHAM_ELEM_S'
!          'CHAM_NO'      'CHAM_ELEM'  'CARTE'
!          'CHAMP' (CHAPEAU AUX CHAM_NO/CHAM_ELEM/CARTE/RESUELEM)
!          'CHAMP_GD' (CHAPEAU DESUET AUX CHAM_NO/CHAM_ELEM/...)
!          'RESULTAT'  'LIGREL'  'NUAGE'  'MAILLAGE' 'CRITERE'
!          'LISTR8'    'LISTIS'  'LISTE_CHARGE'   'NUML_EQUA'
!          (OU ' ' QUAND ON NE CONNAIT PAS LE TYPE).
!       NOMSD   : NOM DE LA STRUCTURE DE DONNEES A DETRUIRE
!          NUME_DDL(K14),MATR_ASSE(K19),VECT_ASSE(K19)
!          CHAMP(K19), MATR_ELEM(K8), VECT_ELEM(K8), VARI_COM(K14)
!          DEFI_CONT(K16), RESO_CONT(K14), TABLE(K19)
!          CHAM_NO(K19), CHAM_NO_S(K19),CHAM_ELEM(K19),CHAM_ELEM_S(K19)
!          CRITERE(K19), LISTE_RELA(K19), CABL_PRECONT(K8), ...
!
!     RESULTAT:
!     ON DETRUIT TOUS LES OBJETS JEVEUX CORRESPONDANT A CES CONCEPTS.
! ----------------------------------------------------------------------
    complex(kind=8) :: cbid
!
    mpi_int :: mrank, msize
    integer(kind=8) :: iret, iad, long, i, nbch, ibid, nbproc, num
    integer(kind=8) :: ityobj, inomsd, nblg, nbpa, nblp, n1, nbjoin
    integer(kind=8) :: iexi, iexi2
    character(len=8) :: metres, k8
    character(len=4) :: chnbjo
    character(len=12) :: vge
    character(len=14) :: nu, com
    character(len=16) :: typ2sd, corres
    character(len=19) :: champ, matas, table, solveu, fnc, resu, matas2, solveu2
    character(len=19) :: ligrel, nuage, ligret, mltf, stock, k19, matel, liste
    character(len=24) :: typobj, knomsd, nojoin
    aster_logical :: lbid
    character(len=24), pointer :: ltns(:) => null()
    character(len=24), pointer :: relr(:) => null()
    character(len=24), pointer :: refa(:) => null()
    cbid = dcmplx(0.d0, 0.d0)
!
! -DEB------------------------------------------------------------------
!
    call jemarq()
    typ2sd = typesd

    call asmpi_info(rank=mrank, size=msize)
    nbproc = to_aster_int(msize)
!
!     ------------------------------------------------------------------
    if (typ2sd .eq. ' ') then
!     -----------------------
!       TYPE_SD INCONNU => CALL JEDETC => COUT CPU IMPORTANT + DANGER
        call jedetc(' ', nomsd, 1)
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CABL_PRECONT') then
!     ------------------------------------
        k8 = nomsd
        call detrs2('CARTE', k8//'.CHME.SIGIN')
        call detrs2('LISTE_RELA', k8//'.LIRELA')
        call detrs2('L_TABLE', k8)
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CHAM_NO_S') then
!     ------------------------------------
        k19 = nomsd
        call jedetr(k19//'.CNSD')
        call jedetr(k19//'.CNSK')
        call jedetr(k19//'.CNSC')
        call jedetr(k19//'.CNSL')
        call jedetr(k19//'.CNSV')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CHAM_ELEM_S') then
!     --------------------------------------
        k19 = nomsd
        call jedetr(k19//'.CESD')
        call jedetr(k19//'.CESK')
        call jedetr(k19//'.CESC')
        call jedetr(k19//'.CESL')
        call jedetr(k19//'.CESV')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'LISTE_RELA') then
!     --------------------------------------
        k19 = nomsd
        call jedetr(k19//'.RLBE')
        call jedetr(k19//'.RLSU')
        call jedetr(k19//'.RLTC')
        call jedetr(k19//'.RLNO')
        call jedetr(k19//'.RLCO')
        call jedetr(k19//'.RLNT')
        call jedetr(k19//'.RLPO')
        call jedetr(k19//'.RLNR')
        call jedetr(k19//'.RLTV')
        call jedetr(k19//'.RLDD')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'PARTITION') then
!     -------------------------------------------
        k19 = nomsd
        call jedetr(k19//'.PRTI')
        call jedetr(k19//'.PRTK')
        call jedetr(k19//'.NUPR')
        call jedetr(k19//'.FDIM')
        call jedetr(k19//'.FETA')
        call jedetr(k19//'.FREF')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CORRESP_2_MAILLA') then
!     -------------------------------------------
        corres = nomsd
        call jedetr(corres//'.PJXX_K1')
        call jedetr(corres//'.PJEF_NB')
        call jedetr(corres//'.PJEF_NU')
        call jedetr(corres//'.PJEF_M1')
        call jedetr(corres//'.PJEF_CF')
        call jedetr(corres//'.PJEF_TR')
        call jedetr(corres//'.PJEF_CO')
        call jedetr(corres//'.PJEF_EL')
        call jedetr(corres//'.PJEF_MP')
        call jedetr(corres//'.PJNG_I1')
        call jedetr(corres//'.PJNG_I2')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CRITERE') then
!     -----------------------------------
        k19 = nomsd
        call jedetr(k19//'.CRTI')
        call jedetr(k19//'.CRTR')
        call jedetr(k19//'.CRDE')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'FONCTION') then
!     -----------------------------------
        fnc = nomsd
        call jedetr(fnc//'.PARA')
        call jedetr(fnc//'.PROL')
        call jedetr(fnc//'.VALE')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'SOLVEUR') then
!     ----------------------------------
        solveu = nomsd
        call jedetr(solveu//'.SLVI')
        call jedetr(solveu//'.SLVK')
        call jedetr(solveu//'.SLVR')
        call jedetr(solveu//'.SLVO')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'EIGENSOLVER') then
!     ----------------------------------
        solveu = nomsd
        call jedetr(solveu//'.ESVK')
        call jedetr(solveu//'.ESVR')
        call jedetr(solveu//'.ESVI')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'VOISINAGE') then
!     ----------------------------------
        vge = nomsd
        call jedetr(vge//'.PTVOIS')
        call jedetr(vge//'.ELVOIS')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'LIGREL') then
!     ----------------------------------
        ligrel = nomsd
        call jedetr(ligrel//'.LGNS')
        call jedetr(ligrel//'.LIEL')
        call jedetr(ligrel//'.NEMA')
        call jedetr(ligrel//'.LGRF')
        call jedetr(ligrel//'.NBNO')
        call jedetr(ligrel//'.NVGE')
        call jedetr(ligrel//'.PRNM')
        call jedetr(ligrel//'.PRNS')
        call jedetr(ligrel//'.REPE')
        call jedetr(ligrel//'.SSSA')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'LIGRET') then
!     ----------------------------------
        ligret = nomsd
        call jedetr(ligret//'.APMA')
        call jedetr(ligret//'.APNO')
        call jedetr(ligret//'.LIMA')
        call jedetr(ligret//'.LINO')
        call jedetr(ligret//'.LITY')
        call jedetr(ligret//'.MATA')
        call jedetr(ligret//'.MODE')
        call jedetr(ligret//'.NBMA')
        call jedetr(ligret//'.LGRF')
        call jedetr(ligret//'.PHEN')
        call jedetr(ligret//'.POMA')
        call jedetr(ligret//'.PONO')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'SQUELETTE') then
!     ----------------------------------
        k8 = nomsd
        call detrs2('MAILLAGE', k8)
        call jedetr(k8//'.CORRES')
        call jedetr(k8//'         .NOMSST')
        call jedetr(k8//'.ANGL_NAUT')
        call jedetr(k8//'.INV.SKELETON')
        call jedetr(k8//'.TRANS')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'L_TABLE') then
!     ----------------------------------
        k19 = nomsd
        call jeexin(k19//'.LTNS', iret)
        if (iret .eq. 0) goto 70
        call jeveuo(k19//'.LTNS', 'L', vk24=ltns)
        call jelira(k19//'.LTNS', 'LONMAX', n1)
        do i = 1, n1
            call detrs2('TABLE', ltns(i))
        end do
        call jedetr(k19//'.LTNS')
        call jedetr(k19//'.LTNT')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'MAILLAGE') then
!     ----------------------------------
        k8 = nomsd
        call detrs2('CHAM_NO', k8//'.COORDO')
        call detrs2('L_TABLE', k8)
        call jedetr(k8//'           .TITR')
        call jedetr(k8//'.CONNEX')
        call jedetr(k8//'.DIME')
        call jedetr(k8//'.GROUPEMA')
        call jedetr(k8//'.GROUPENO')
        call jedetr(k8//'.NOMACR')
        call jedetr(k8//'.PARA_R')
        call jedetr(k8//'.SUPMAIL')
        call jedetr(k8//'.TYPL')
        call jedetr(k8//'.TYPMAIL')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'MODELE') then
!     ----------------------------------
        k8 = nomsd
        call detrs2('LIGREL', k8//'.MODELE')
        call detrs2('L_TABLE', k8)
!
        call jedetr(k8//'           .TITR')
        call jedetr(k8//'.MAILLE')
        call detrs2('PARTITION', k8//'.PARTSD')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'NUAGE') then
!     ----------------------------------
        nuage = nomsd
        call jedetr(nuage//'.NUAI')
        call jedetr(nuage//'.NUAX')
        call jedetr(nuage//'.NUAV')
        call jedetr(nuage//'.NUAL')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'TABLE_CONTAINER') then
!     -----------------------------------
        table = nomsd
        call jeexin(table//'.TBLP', iret)
        if (iret .ne. 0) then
            call jeveuo(table//'.TBNP', 'L', iad)
            nblg = zi(iad+1)
            nbpa = zi(iad)
            ASSERT(nbpa .ge. 3)
            call jeveuo(table//'.TBLP', 'L', iad)
            call jelira(table//'.TBLP', 'LONMAX', long)
            nblp = long/nbpa
            do i = 1, nbpa
                if (zk24(iad+nblp*(i-1)) (1:10) .eq. 'TYPE_OBJET') then
                    typobj = zk24(iad+nblp*(i-1)+3-1)
                else if (zk24(iad+nblp*(i-1)) (1:6) .eq. 'NOM_SD') then
                    knomsd = zk24(iad+nblp*(i-1)+3-1)
                end if
            end do
            call jeveuo(typobj, 'L', ityobj)
            call jeveuo(knomsd, 'L', inomsd)
            do i = 1, nblg
                if (zk16(ityobj+i-1) (1:9) .eq. 'MATR_ELEM') then
                    call detrs2('MATR_ELEM', zk24(inomsd+i-1))
                else if (zk16(ityobj+i-1) (1:9) .eq. 'VECT_ELEM') then
                    call detrs2('VECT_ELEM', zk24(inomsd+i-1))
                else if (zk16(ityobj+i-1) (1:9) .eq. 'CHAM_ELEM') then
                    call detrs2('CHAM_ELEM', zk24(inomsd+i-1))
                end if
            end do
            do i = 1, long
                call jedetr(zk24(iad-1+i))
            end do
            call jedetr(table//'.TBLP')
            call jedetr(table//'.TBNP')
            call jedetr(table//'.TBBA')
        end if
        call jedetr(table//'.TITR')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'TABLE') then
!     --------------------------------
        table = nomsd
        call jeexin(table//'.TBLP', iret)
        if (iret .ne. 0) then
            call jeveuo(table//'.TBLP', 'L', iad)
            call jelira(table//'.TBLP', 'LONMAX', long)
            do i = 1, long
                call jedetr(zk24(iad-1+i))
            end do
            call jedetr(table//'.TBLP')
            call jedetr(table//'.TBNP')
            call jedetr(table//'.TBBA')
        end if
        call jedetr(table//'.TITR')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'MATR_ASSE_GENE' .or. typ2sd .eq. 'MATR_ASSE') then
!     ---------------------------------------
        matas = nomsd
!
!       -- DESTRUCTION DE L'EVENTUELLE MATRICE RéDUITE (ELIM_LAGR) :
        call jeexin(matas//'.REFA', iexi)
        if (iexi .gt. 0) then
            call jeveuo(matas//'.REFA', 'L', vk24=refa)
            if (refa(19) .ne. ' ') then
                matas2 = refa(19) (1:19)
                call jeexin(matas2//".REFA", iexi2)
                if (iexi2 .gt. 0) then
                    call dismoi('SOLVEUR', matas2, 'MATR_ASSE', repk=solveu2, arret='C', &
                                ier=iret)
                    if (iret .eq. 0) then
                        call detlsp(matas2, solveu2)
                    end if
                    call detrs2('MATR_ASSE', matas2)
                    call elg_gest_data('EFFACE', ' ', matas2, ' ')
                end if
!       -- DESTRUCTION DE L'EVENTUELLE MATRICE PRE CONDITIONNEUR XFEM :
!            if (refa(18) .ne. ' ' .or. refa(16) .ne. ' ') then
!                call xfem_pc_detr(matas)
!            endif
            end if
        end if
!
!       -- DESTRUCTION DE L'EVENTUELLE INSTANCE MUMPS OU PETSC :
!          et des matrices PETSC liées à ELIM_LAGR='OUI'
        call jeexin(matas//'.REFA', iexi)
        if (iexi .gt. 0) then
            call dismoi('METH_RESO', matas, 'MATR_ASSE', repk=metres)
            if (metres .eq. 'MUMPS') then
                call amumph('DETR_MAT', ' ', matas, [0.d0], [cbid], &
                            ' ', 0, ibid, lbid)
            else if (metres .eq. 'PETSC') then
                call apetsc('DETR_MAT', ' ', matas, [0.d0], ' ', &
                            0, ibid, iret)
            end if
        end if
!
        call jedetr(matas//'.CCID')
        call jedetr(matas//'.CCII')
        call jedetr(matas//'.CCLL')
        call jedetr(matas//'.CCVA')
        call jedetr(matas//'.CONL')
        call jedetr(matas//'.DESC')
        call jedetr(matas//'.DIGS')
        call jedetr(matas//'.REFA')
        call jedetr(matas//'.PERM')
        call jedetr(matas//'.UALF')
        call jedetr(matas//'.VALF')
        call jedetr(matas//'.VALM')
        call jedetr(matas//'.WALF')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CHAM_NO') then
!     ----------------------------------
        k19 = nomsd
        call jedetr(k19//'.REFE')
        call jedetr(k19//'.VALE')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CARTE') then
!     ----------------------------------
        k19 = nomsd
        call jedetr(k19//'.DESC')
        call jedetr(k19//'.VALE')
        call jedetr(k19//'.NOMA')
        call jedetr(k19//'.NOLI')
        call jedetr(k19//'.LIMA')
        call jedetr(k19//'.PTMA')
        call jedetr(k19//'.PTMS')
        call jedetr(k19//'.NCMP')
        call jedetr(k19//'.VALV')
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'DOMJOINTS') then
        !     ----------------------------------
        k19 = nomsd
        call jedetr(k19//'.DOMJ')
        call jedetr(k19//'.SEND')
        call jedetr(k19//'.RECV')
        call jedetr(k19//'.PGID')
        call jedetr(k19//'.GCOM')
        !
!     ------------------------------------------------------------------
!
    else if (typ2sd .eq. 'NUME_EQUA_GENE') then
!     ------------------------------------
        k19 = nomsd
        call jedetr(k19//'.DEEQ')
        call jedetr(k19//'.LILI')
        call jedetr(k19//'.NUEQ')
        call jedetr(k19//'.PRNO')
        call jedetr(k19//'.ORIG')
        call jedetr(k19//'.NEQU')
        call jedetr(k19//'.DESC')
        call jedetr(k19//'.REFN')
        call jedetr(k19//'.DELG')
!
    else if (typ2sd .eq. 'NUME_EQUA') then
!     ------------------------------------
        k19 = nomsd
        call jedetr(k19//'.DEEQ')
        call jedetr(k19//'.LILI')
        call jedetr(k19//'.NUEQ')
        call jedetr(k19//'.PRNO')
        call jedetr(k19//'.NEQU')
        call jedetr(k19//'.REFN')
        call jedetr(k19//'.DELG')
        call jedetr(k19//'.NULG')
        call jedetr(k19//'.PDDL')
        call jedetr(k19//'.CRCO')
        call jedetr(k19//'.REFP')

    else if (typ2sd .eq. 'NUML_EQUA') then
!  --------------------------------------
        k19 = nomsd
        call jedetr(k19//'.PRNO')
        call jedetr(k19//'.NEQU')
        call jedetr(k19//'.DELG')
        call jedetr(k19//'.NUEQ')
        call jedetr(k19//'.NULG')
        call jedetr(k19//'.NUGL')
        call jedetr(k19//'.NLGP')
        call jedetr(k19//'.PDDL')
        call jeexin(k19//'.JOIN', iexi)
        if (iexi .gt. 0) then
            call jelira(k19//'.JOIN', 'LONMAX', nbjoin)
            call jedetr(k19//'.JOIN')
            do num = 1, nbjoin
                call codent(num, 'G', chnbjo)
                nojoin = k19//'.'//chnbjo
                call jedetr(nojoin)
            end do
        end if
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'CHAM_ELEM') then
!     ------------------------------------
        k19 = nomsd
        call jedetr(k19//'.CELD')
        call jedetr(k19//'.CELK')
        call jedetr(k19//'.CELV')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'RESUELEM') then
!     ------------------------------------
        k19 = nomsd
        call jedetr(k19//'.DESC')
        call jedetr(k19//'.NOLI')
        call jedetr(k19//'.RESL')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'MLTF') then
!     -----------------------------------
        mltf = nomsd
        call jedetr(mltf//'.GLOB')
        call jedetr(mltf//'.LOCL')
        call jedetr(mltf//'.ADNT')
        call jedetr(mltf//'.PNTI')
        call jedetr(mltf//'.DESC')
        call jedetr(mltf//'.DIAG')
        call jedetr(mltf//'.ADRE')
        call jedetr(mltf//'.SUPN')
        call jedetr(mltf//'.PARE')
        call jedetr(mltf//'.FILS')
        call jedetr(mltf//'.FRER')
        call jedetr(mltf//'.LGSN')
        call jedetr(mltf//'.LFRN')
        call jedetr(mltf//'.NBAS')
        call jedetr(mltf//'.DEBF')
        call jedetr(mltf//'.DEFS')
        call jedetr(mltf//'.ADPI')
        call jedetr(mltf//'.ANCI')
        call jedetr(mltf//'.NBLI')
        call jedetr(mltf//'.LGBL')
        call jedetr(mltf//'.NCBL')
        call jedetr(mltf//'.DECA')
        call jedetr(mltf//'.NOUV')
        call jedetr(mltf//'.SEQU')
        call jedetr(mltf//'.RENU')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'STOCKAGE') then
!     -----------------------------------
        stock = nomsd
        call jedetr(stock//'.SCBL')
        call jedetr(stock//'.SCDI')
        call jedetr(stock//'.SCDE')
        call jedetr(stock//'.SCHC')
        call jedetr(stock//'.SCIB')
        call jedetr(stock//'.M2LC')
        call jedetr(stock//'.LC2M')
!
        call jedetr(stock//'.SMDI')
        call jedetr(stock//'.SMDE')
        call jedetr(stock//'.SMHC')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'STOC_LCIEL') then
!     -----------------------------------
        stock = nomsd
        call jedetr(stock//'.SCBL')
        call jedetr(stock//'.SCDI')
        call jedetr(stock//'.SCDE')
        call jedetr(stock//'.SCHC')
        call jedetr(stock//'.SCIB')
        call jedetr(stock//'.M2LC')
        call jedetr(stock//'.LC2M')
!
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'STOC_MORSE') then
!     -----------------------------------
        stock = nomsd
        call jedetr(stock//'.SMDI')
        call jedetr(stock//'.SMDE')
        call jedetr(stock//'.SMHC')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'NUME_DDL') then
!     -----------------------------------
        nu = nomsd
        call detrs2('NUME_EQUA', nu//'.NUME')
        call detrs2('MLTF', nu//'.MLTF')
        call detrs2('STOCKAGE', nu//'.SLCS')
        call detrs2('STOCKAGE', nu//'.SMOS')
!
        call jedetr(nu//'.DERLI')
        call jedetr(nu//'.EXISTE')
        call jedetr(nu//'.NUM2')
        call jedetr(nu//'.NUM21')
        call jedetr(nu//'.LSUIVE')
        call jedetr(nu//'.PSUIVE')
        call jedetr(nu//'.VSUIVE')
        call jedetr(nu//'.OLDN')
        call jedetr(nu//'.NEWN')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'VARI_COM') then
!     -------------------------------------
        com = nomsd
        call assde1('CHAMP', com//'.TEMP')
        call assde1('CHAMP', com//'.HYDR')
        call assde1('CHAMP', com//'.SECH')
        call assde1('CHAMP', com//'.PHAS')
        call assde1('CHAMP', com//'.EPAN')
        call assde1('CHAMP', com//'.INST')
        call assde1('CHAMP', com//'.TOUT')
!
!     ------------------------------------------------------------------
    else if ((typ2sd .eq. 'CHAMP') .or. (typ2sd .eq. 'CHAMP_GD')) then
!     -------------------------------------------------------------
!       POUR LES CARTE, CHAM_NO, CHAM_ELEM, ET RESU_ELEM :
        champ = nomsd
        call jedetr(champ//'.CELD')
        call jedetr(champ//'.CELK')
        call jedetr(champ//'.CELV')
        call jedetr(champ//'.DESC')
        call jedetr(champ//'.LIMA')
        call jedetr(champ//'.NCMP')
        call jedetr(champ//'.NOLI')
        call jedetr(champ//'.NOMA')
        call jedetr(champ//'.PTMA')
        call jedetr(champ//'.PTMS')
        call jedetr(champ//'.REFE')
        call jedetr(champ//'.RESL')
        call jedetr(champ//'.VALE')
        call jedetr(champ//'.VALV')
!
!     ------------------------------------------------------------------
    else if ((typ2sd .eq. 'MATR_ELEM') .or. (typ2sd .eq. 'VECT_ELEM')) then
!     ---------------------------------------
        matel = nomsd
        call jeexin(matel//'.RELR', iret)
        if (iret .le. 0) goto 61
        call jelira(matel//'.RELR', 'LONUTI', nbch)
        if (nbch .gt. 0) call jeveuo(matel//'.RELR', 'L', vk24=relr)
        do i = 1, nbch
            champ = relr(i) (1:19)
            call assde1('RESUELEM', champ)
        end do
61      continue
        call jedetr(matel//'.RELR')
        call jedetr(matel//'.RELC')
        call jedetr(matel//'.RERR')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'RESULTAT') then
!     -----------------------------------
        resu = nomsd
        call jedetr(resu//'.DESC')
        call jedetr(resu//'.TACH')
        call jedetr(resu//'.TAVA')
        call jedetr(resu//'.NOVA')
        call jedetr(resu//'.ORDR')
        call jedetr(resu//'.REFD')
        call jedetr(resu//'.INDI')
        call jedetr(resu//'.RSPR')
        call jedetr(resu//'.RSPC')
        call jedetr(resu//'.RSPI')
        call jedetr(resu//'.RSP8')
        call jedetr(resu//'.RS16')
        call jedetr(resu//'.RS24')
        call jedetr(resu//'.RS32')
        call jedetr(resu//'.RS80')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'LISTE_CHARGES') then
!     -----------------------------------
        resu = nomsd
        call jedetr(resu//'.LCHA')
        call jedetr(resu//'.INFC')
        call jedetr(resu//'.FCHA')
!
!     ------------------------------------------------------------------
    else if (typ2sd .eq. 'LISTR8' .or. typ2sd .eq. 'LISTIS') then
!     --------------------------------------------------------
        liste = nomsd
        call jedetr(liste//'.LPAS')
        call jedetr(liste//'.NBPA')
        call jedetr(liste//'.BINT')
        call jedetr(liste//'.VALE')
!
!     ------------------------------------------------------------------
    else
        ASSERT(.false.)
    end if
!
70  continue
    call jedema()
end subroutine
