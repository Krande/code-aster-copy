# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`TransientGeneralizedResult` --- Results container
**************************************************
"""

import numpy

import aster
from libaster import TransientGeneralizedResult

from ..Cata.Syntax import AsException
from ..Utilities import injector


@injector(TransientGeneralizedResult)
class ExtendedTransientGeneralizedResult:
    cata_sdj = "SD.sd_dyna_gene.sd_dyna_gene"

    def _check_input_inoli(self, inoli):
        if (inoli==-1) :
            print("Nonlinearity index not specified, by default the first nonlinearity will be considered.")
            inoli = 1
        nbnoli = self._nb_nonl()
        if nbnoli == 0 :
            raise AsException("Linear calculation, no information can be retrieved.")
        if( inoli <= 0) or (inoli > nbnoli):
            raise AsException("The nonlinearity index should be a comprised between 1 and %d, the total number of nonlinearities."%(nbnoli))
        return inoli

    def _nb_nonl(self):
        desc = self.sdj.DESC.get()
        nbnoli = desc[2]
        return nbnoli

    def _print_vint_description(self, inoli):
        nltype = self._type_nonl()[inoli-1].strip()
        vintDescription = {'DIS_CHOC'      : ['F_NORMAL', 'F_TANGE1', 'F_TANGE1',
                                              'DXLOC_N1', 'DYLOC_N1', 'DZLOC_N1',
                                              'DXLOC_N2', 'DYLOC_N2', 'DZLOC_N2',
                                              'V_NORMAL', 'V_TANGE1', 'V_TANGE1',
                                              'IND_ADHE', 'VINT_FR1', 'VINT_FR2',
                                              'VINT_FR3', 'VINT_FR4', 'VINT_FR5',
                                              'VINT_FR6', 'VINT_FR7'],
                          'FLAMBAGE'      : ['F_NORMAL',
                                              'DXLOC_N1', 'DYLOC_N1', 'DZLOC_N1',
                                              'DXLOC_N2', 'DYLOC_N2', 'DZLOC_N2',
                                              'V_NORMAL', 'ENFO_PLA', 'RIGI_P_F',
                                              'ENFO_MAX'],
                           'ANTI_SISM'     : ['F_AXIAL',
                                              'DXLOC_N1', 'DYLOC_N1', 'DZLOC_N1',
                                              'DXLOC_N2', 'DYLOC_N2', 'DZLOC_N2',
                                              'V_NORMAL'],
                           'DIS_VISC'      : ['DXLOC_N1', 'DYLOC_N1', 'DZLOC_N1',
                                              'DXLOC_N2', 'DYLOC_N2', 'DZLOC_N2',
                                              'V_AXIALE', 'F_AXIALE', 'V_AXIALE',
                                              'D_AXIALE', 'PUISSANC'],
                           'DIS_ECRO_TRAC' : ['DXLOC_N1', 'DYLOC_N1', 'DZLOC_N1',
                                              'DXLOC_N2', 'DYLOC_N2', 'DZLOC_N2',
                                              'FXLOC',    'FYLOC',    'FZLOC',
                                              'DXLOC',    'DYLOC',    'DZLOC',
                                              'PUISSANC', 'PCUM',     'DXPLOC',
                                              'DYPLOC',   'DZPLOC'],
                           'ROTOR_FISS'    : ['PHI_DEGR', 'F_TANGE1', 'F_TANGE2'],
                           'COUPLAGE_EDYOS': [],
                           'RELA_EFFO_DEPL': ['DCMP_N1 ', 'FCMP_LOC', 'IND_NONZ'] ,
                           'RELA_EFFO_VITE': ['VCMP_N1 ', 'FCMP_LOC', 'IND_NONZ']  }
        print("\n" + "-"*104)
        print("Information regarding the saved internal variables for %s non linearity (index=%d)"%(nltype, inoli))
        print("-"*104)
        vintDesc = [v.center(10) for v in vintDescription[nltype]]
        indices  = [str(i+1).center(10) for i in range(len(vintDesc))]
        sep = " | "

        nblines = len(indices) // 8
        if 8*nblines<len(indices) : nblines = nblines + 1
        for i in range(nblines-1):
            print(sep.join(indices [i*8:(i+1)*8]))
            print(sep.join(vintDesc[i*8:(i+1)*8]))
            print("-"*104)
        print(sep.join(indices [8*(nblines-1):]))
        print(sep.join(vintDesc[8*(nblines-1):]))
        print("-"*104)
        return vintDesc

    def _type_nonl(self):
        Int2StrTypes = {1 : 'DIS_CHOC',
                        2 : 'FLAMBAGE',
                        3 : 'ANTI_SISM',
                        4 : 'DIS_VISC',
                        5 : 'DIS_ECRO_TRAC',
                        6 : 'ROTOR_FISS',
                        7 : 'COUPLAGE_EDYOS',
                        8 : 'RELA_EFFO_DEPL',
                        9 : 'RELA_EFFO_VITE'}

        nltypes = self.sdj.sd_nl.TYPE.get()
        return [Int2StrTypes[nltypes[i]] for i in range(len(nltypes))]

    def FORCE_NORMALE (self, inoli=-1):
        """
        Returns a 1D numpy array giving the evolution of the normal force at the archived instants"""

        inoli = self._check_input_inoli(inoli)

        nltypes = self._type_nonl()
        if not(nltypes[inoli-1] in ('DIS_CHOC', 'FLAMBAGE')) :
            dummy = self.INFO_NONL()
            raise AsException("The chosen nonlinearity index (%d) does not correspond to a DIS_CHOC or FLAMBAGE nonlinearity\nThese are the only nonlinearities that save the local normal force."%(inoli))


        vint = self.VARI_INTERNE(inoli, describe=False)
        #The normal force is saved in the first position (ind=0) of the internal variables for DIS_CHOC and FLAMBAGE nonlinearities
        return vint[:,0]

    def INFO_NONL(self):
        """
        Prints out information about the considered non linearities, returns a 2D python list (list in list) with
        the retrieved information"""

        nbnoli  = self._nb_nonl()
        if nbnoli == 0 :
            print("Linear calculation, no nonlinearities used or can be printed.")
            return None

        nltypes = self._type_nonl()
        inti    = self.sdj.sd_nl.INTI.get()

        print("-"*104)
        print("%sInformation regarding the considered non linearities"%(' '*15))
        print("-"*104)
        #      12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234
        #      ooo-----ooo+++++++++++++++++ooo---------ooo+++++++++ooo---------ooo+++++++++ooo-------------------------
        print(" |  IND  |       TYPE        |    NO1    |    NO2    |   SST1    |    SST2   |           TITLE          |")
        print("-"*104)
        Output =[None]*nbnoli
        for i in range(nbnoli):
            title, no1, no2, sst1, sst2 = inti[i*5:i*5+5]
            if len( no2.strip()) == 0 : no2  = '-'
            if len(sst1.strip()) == 0 : sst1 = '-'
            if len(sst2.strip()) == 0 : sst2 = '-'
            title = title.strip().center(25)
            no1   = no1.strip().center(9)
            no2   = no2.strip().center(9)
            sst1  = sst1.strip().center(9)
            sst2  = sst2.strip().center(9)
            sep = ' | '
            print("%s%s%s%s%s%s%s%s%s%s%s%s%s%s"%(sep,str(i+1).center(5),sep,nltypes[i].center(17),sep,no1,sep,no2,sep,sst1,sep,sst2,sep,title))
            add = [nltypes[i]]+list(inti[(i-1)*5:(i-1)*5+5])
            Output[i] = add
        print("-"*104)
        return Output

    def LIST_ARCH (self):
        """
        Returns a python list of all archived instants
        """

        disc = self.sdj.DISC.get()
        return list(disc)

    def VARI_INTERNE (self, inoli=-1, describe=True):
        """
        Returns a 2D numpy array of all internal variables for a given non linearity of index <inoli>
        """

        inoli = self._check_input_inoli(inoli)
        i = inoli-1

        vindx  = self.sdj.sd_nl.VIND.get()
        nbvint = vindx[-1]-1    # number of internal variables saved for all nonlinearities : record length of VINT

        vint    = self.sdj.sd_nl.VINT.get()
        nbsaves = len(vint) // nbvint

        start  = vindx[i  ]-1
        finish = vindx[i+1]-1
        outputLength = (finish-start)*nbsaves

        cntr = 0
        output = [0.]*outputLength
        for iord in range(nbsaves):
            for i in range(start, finish):
                output[cntr] = vint[iord*(nbvint)+i]
                cntr += 1

        output = numpy.reshape(output,(nbsaves, finish-start))

        if describe:
            dummy = self._print_vint_description(inoli)

        return output
