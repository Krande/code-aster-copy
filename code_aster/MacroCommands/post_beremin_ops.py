# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
"""
Compute:
  - Weibull stress
  - Weibull stress**m
  - Probability of failure

  - Strain model -> Stress to consider (6 components: XX, YY, ZZ, XY, XZ, YZ)
      - "PETIT"      -> "SIEF"
      - "PETIT_REAC" -> "SIEF"
      - "GDEF_LOG"   -> "VARI_ELGA"
      - "SIMO_MIEHE" -> To retrieve Kirchhof stresses ? Not obvious ...
                        Normally, we can do without
                        Use of HHO models instead ?

Integral must be computed on the initial configuration

_biblio : CR-T66-2017-132, Lorentz
"""
from math import exp
import os
import numpy as np
from libaster import EntityType

from ..Cata.Syntax import _F
from ..CodeCommands import CALC_CHAM_ELEM, CALC_CHAMP, CREA_TABLE, DEFI_GROUP
from ..Objects import ExternalVariableTraits, FieldOnCellsReal, NonLinearResult
from ..Utilities import logger, no_new_attributes

DEBUG = bool(os.environ.get("DEBUG_POST_BEREMIN", ""))


class PostBeremin:
    """
    Beremin post-processor

    Arguments:
        resultat (NonLinearResult): Mechanical code_aster result
        group_ma (str): Mesh cells group on which Beremin post-treatment
            is carried out
            Only 1 mesh cells group authorized
        deformation (str): "PETIT", "PETIT_REAC", "GDEF_LOG"
        filtre_sigm (str): "SIGM_ELGA" or "SIGM_ELMOY".
        coef_mult (float): As explained in doc aster u4.81.22 (weibull) p.15

    Returns:
        Table:

            - Weibull stress
            - Weibull stress**m
            - Failure probability

    if SIGM_MAXI in command file, MED file containing:

        - where the plasticity is active : maximum of major principal
          stress
        - where the plasticity is not active : 0
          Result of type ELGA (apply ElgaFieldToSurface filter in
          Paraview to see it)
    """

    # data
    _result = _zone = _zone_ids = _stress_option = _strain_type = None
    _intvar_idx = _stress_idx = None
    _use_hist = _use_function = None
    _coef_mult = None
    _weib_params = None
    # result to compute weibull stress
    _reswb = _rsieq = _rout = None
    _group_name = "__plast__"
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, result: NonLinearResult, strain_type: str, stress_option: str) -> None:
        """Initialization and definition of main arguments.

        Args:
            result (*Result*): Result object of the calculation.
            strain_type (str): Type of strain.
            stress_option (str): Type of stress to be used or calculated.
        """
        self._result = result
        assert stress_option in ("SIGM_ELGA", "SIGM_ELMOY"), f"unknown value: {stress_option}"
        self._stress_option = stress_option
        self._strain_type = strain_type

    def set_zone(self, group: str) -> None:
        """Define the location where the calculation will be done.

        Args:
            group (str): Mesh group of cells to be used.
        """
        self._zone = group
        self._zone_ids = self._result.getMesh().getCells(group)

    def set_indexes(self, intvar_idx: list[int] = None, stress_idx: list[int] = None) -> None:
        """Define the indexes used to extract the relevant components depending
        of the behaviour.

        Args:
            intvar_idx (list[int], optional): Indexes to access to EPSPEQ and
                INDIPLAS internal variables.
            stress_idx (list[int], optional): Indexes to extract log stresses
                tensor components.
        """
        if intvar_idx:
            self._intvar_idx = intvar_idx
        if stress_idx:
            assert len(stress_idx) in (4, 6), stress_idx
            self._stress_idx = stress_idx

    def use_history(self, value: bool) -> None:
        """Enable the use of history or just the current timestep.

        Args:
            value (bool): *True* to use the history since the beginning,
              or *False* to use values on the current timestep.
        """
        self._use_hist = value

    def set_coef_mult(self, value: float = 1.0) -> None:
        """Define the multiplicative factor to take symmetric into account.

        Args:
            value (float): Multiplicative factor (see documentation).
        """
        self._coef_mult = value

    def set_weibull_parameters(self, params):
        """Define Weibull parameters.

        Args:
            params (dict): Keywords for WEIBULL.
        """
        # FIXME may vary per groups on cells, but only one group for the moment
        self._use_function = not isinstance(params["SIGM_REFE"], (int, float))
        assert not self._use_function or "SIGM_CNV" in params, "SIGM_CNV is required!"
        self._weib_params = params.copy()

    def setup_calc_result(self):
        """Setup the result object to be used to compute the Beremin stress."""
        if self._strain_type == "GDEF_LOG":
            result = self._result
            indexes = result.getIndexes()
            # store VARI_ELGA + log stresses (from VARI_ELGA) as SIEF_ELGA
            dim = result.getMesh().getDimension()
            varn = [f"V{i}" for i in self._stress_idx]
            cmps = ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
            if dim == 2:
                cmps = ("SIXX", "SIYY", "SIZZ", "SIXY")
            dconv = dict(zip(varn, cmps))

            rvga = self._reswb = NonLinearResult()
            rvga.allocate(len(indexes))
            rvga.userName = "rvga____"

            for idx in indexes:
                chvari = result.getField("VARI_ELGA", idx)
                rvga.setField(chvari.asPhysicalQuantity("SIEF_R", dconv), "SIEF_ELGA", idx)
                # FIXME rvga.setMaterialField(result.getMaterialField(idx), idx)
                rvga.setField(result.getField("COMPORTEMENT", idx), "COMPORTEMENT", idx)
                rvga.setTime(result.getTime(idx), idx)
        else:
            self._reswb = self._result

    def add_plast_group(self, idx):
        """Add a temporary group where plastic strain exists.

        Args:
            idx (int): Storage index to be extracted.
        """
        epseq, _ = self.get_internal_variables(idx)
        cells = self.get_plast_cells(epseq, epseq.max())
        mesh = self._result.getMesh()
        mesh.setGroupOfCells(self._group_name, cells)

    def get_plast_cells(self, epseq, epsmax):
        """Return the cells ids where EPSPEQ is greater than a threshold.

        Args:
            epseq (*ComponentOnCells*): EPSPEQ component (on all cells).
            epsmax (float):

        Returns:
            list[int]: List of cells ids.
        """
        eps_thr = self._weib_params["SEUIL_EPSP_CUMU"]
        mini = 0.0 if epsmax < eps_thr else eps_thr
        cells = epseq.filterByValues(mini, epsmax, strict_maxi=False)
        return list(set(cells).intersection(self._zone_ids))

    def get_internal_variables(self, idx):
        """Return internal variables.

        Args:
            idx (int): storage index.

        Returns:
            tuple: *ComponentOnCells* objects for EPSPEQ and INDIPLAS.
        """
        v1 = f"V{self._intvar_idx[0]}"
        v2 = f"V{self._intvar_idx[1]}"
        chvari = self._result.getField("VARI_ELGA", idx)
        chvari = chvari.toSimpleFieldOnCells()
        epseq = chvari.getComponentValues(v1)
        indip = chvari.getComponentValues(v2)
        return epseq, indip

    def get_major_stress(self, idx):
        """Return the major principal stress.

        Args:
            idx (int): storage index.

        Returns:
            *ComponentOnCells*: Major principal stress.
        """
        if not self._rsieq:
            self._rsieq = self.compute_sieq()
        sieq = self._rsieq.getField("SIEQ_ELGA", idx)
        sieq = sieq.toSimpleFieldOnCells()
        return sieq.PRIN_3

    def compute_sieq(self):
        """Compute the Beremin stress."""
        logger.info("starting computation...")
        # compute fields on all timesteps, add INST/NUME_ORDRE to be limited
        self.add_plast_group(self._reswb.getLastIndex())
        input = self._reswb
        if self._stress_option == "SIGM_ELMOY":
            input = CALC_CHAMP(
                RESULTAT=self._reswb, GROUP_MA=self._group_name, CONTRAINTE="SIMY_ELGA"
            )
            for idx in input.getIndexes():
                simy = input.getField("SIMY_ELGA", idx)
                input.setField(simy, "SIEF_ELGA", idx)

        rsieq = CALC_CHAMP(RESULTAT=input, GROUP_MA=self._group_name, CRITERES="SIEQ_ELGA")
        self.clean_group()
        if DEBUG:
            for i in rsieq.getIndexes():
                print("NEW: rsieq:", i, np.sum(rsieq.getField("SIEQ_ELGA", i).getValues()))

        return rsieq

    def clean_group(self):
        """Remove temporary added group."""
        mesh = self._result.getMesh()
        DEFI_GROUP(reuse=mesh, MAILLAGE=mesh, DETR_GROUP_MA=_F(NOM=self._group_name))

    def apply_stress_correction(self, sig1, idx, cells_ids):
        """Compute the major principal stress, eventually apply SIGM_CNV/SIGM_REFE
        correction and threshold by SIGM_SEUIL.

        *The values are changed in place.*

        Args:
            sig1 (*ComponentOnCells*): PRIN_3 component of SIEQ_ELGA, changed in place.
            idx (int): Timestep index
            cells_ids (list[int]): Cells list.
        """
        logger.info("adjusting sigma1 at index %d...", idx)

        if self._use_function:
            # PRIN_3 / SIGM_REFE(TEMP) * SIGM_CNV
            temp = self.get_temperature_field(idx)
            temp.restrict(cells_ids)
            sigref = temp.apply(self._weib_params["SIGM_REFE"])
            assert sigref.abs().min() > 0.0
            sig1 /= sigref
            sig1 *= self._weib_params["SIGM_CNV"]

        if self._weib_params["SIGM_SEUIL"] > 0.0:
            # apply threshold
            valsig = sig1.getValues() - self._weib_params["SIGM_SEUIL"]
            valsig = np.where(valsig < 0.0, 0.0, valsig)
            sig1.setValues(valsig)

    def apply_threshold(self, sig1, epseq, indiplas):
        """Apply plastic strain threshold onto major stress.

        NB: *ComponentOnCells* changed in place.

        Args:
            sig1 (*ComponentOnCells*): Major stress.
            epspeq (*ComponentOnCells*): Plastic strain.
            indiplas (*ComponentOnCells*): Plastic indicator.
        """
        valthr = epseq.onSupportOf(sig1)
        indip = indiplas.onSupportOf(sig1)
        epsthr_ar = valthr.getValues() - self._weib_params["SEUIL_EPSP_CUMU"]
        epsthr_ar = np.where(epsthr_ar < 0.0, 0.0, 1.0)
        valthr.setValues(epsthr_ar)

        # * (1 if indip != 0 else 0)
        def toint_sign(array):
            return np.sign(np.array(array, dtype=int))

        indip = indip.applyOnArray(toint_sign)
        valthr *= indip
        sig1 *= valthr

    def get_temperature_field(self, idx):
        """Return the temperature field.

        Args:
            idx (idx): storage index.

        Returns:
            *ComponentOnCells*: Temperature values.
        """
        field = None
        for extvar, where in self._result.getMaterialField().getExtStateVariablesOnMeshEntities():
            varname = ExternalVariableTraits.getExternVarTypeStr(extvar.getType())
            if varname != "TEMP":
                continue
            assert where == EntityType.AllMeshEntitiesType
            field = extvar.getField()
            if not field:
                evol = extvar.getEvolutionParameter()
                assert evol, "expecting an evolution or a field, what else?!"
                field = evol.getTransientResult().getField("TEMP", idx)
            break
        assert field, "TEMP not found in VARC"
        fed = self._result.getModel().getFiniteElementDescriptor()
        field = field.toFieldOnCells(fed, "ELGA").toSimpleFieldOnCells()
        return field.TEMP

    def main(self):
        """Compute the Beremin stress."""
        result = self._reswb
        model = result.getModel()
        params = result.getAccessParameters()
        logger.info("extracting integration scheme...")
        coor_elga = CALC_CHAM_ELEM(MODELE=model, OPTION="COOR_ELGA")
        coor_elga = coor_elga.toSimpleFieldOnCells()

        logger.info("starting computation...")
        previous = None
        table = TableBeremin()

        id_store = 0
        for idx, time in zip(params["NUME_ORDRE"], params["INST"]):
            epseq, indip = self.get_internal_variables(idx)
            epseq_r = epseq.copy()
            epseq_r.restrict(self._zone_ids)
            maxepsp = epseq_r.max()
            if maxepsp <= 0.0:
                if id_store == 0:
                    # add initial state
                    table.append(0, self._result.getTime(0), self._zone, 0.0, 0.0, 0.0)
                    id_store += 1
                continue

            logger.info("getting plastified zone...")
            cells_ids = self.get_plast_cells(epseq, maxepsp)
            logger.info("computing major stress...")
            sig1 = self.get_major_stress(idx)
            sig1.restrict(cells_ids)
            self.apply_stress_correction(sig1, idx, cells_ids)
            self.apply_threshold(sig1, epseq, indip)

            if self._use_hist == "OUI":
                logger.info("computing maxi at index %d...", idx)
                if previous is None:
                    previous = sig1 * 0.0
                sig1 = sig1.maximum(previous)
                previous = sig1
            if DEBUG:
                print("NEW: sigmax:", id_store, idx, sig1.sum())

            sig1 **= self._weib_params["M"]
            self.store_sigm_maxi(id_store, time, sig1, model)

            # create a new ComponentOnCells object since it is restricted on a different group
            strwb, strwb_pm, proba = self.compute_table_values(sig1, coor_elga.W)

            table.append(id_store, time, self._zone, strwb, strwb_pm, proba)
            id_store += 1

        if id_store == 0:
            # no plastification
            table.append(id_store, 0.0, self._zone, 0.0, 0.0, 0.0)

        return table

    def compute_table_values(self, sig1, weight):
        """Compute the values to be added into the result table.

        Args:
            sig1 (*ComponentOnCells*): Major stress ^ M.
            weight (*ComponentOnCells*): Weight of each integration point,
                **changed in place**.

        Returns:
            tuple: Values for each column.
        """
        coef_volu = self._coef_mult / self._weib_params["VOLU_REFE"]
        pow_m = self._weib_params["M"]
        inv_m = 1 / pow_m
        sigma_u = self._weib_params["SIGM_CNV" if self._use_function else "SIGM_REFE"]

        logger.info("computing table values...")
        if DEBUG:
            print("NEW: sigpm:", sig1.sum())

        weight = weight.onSupportOf(sig1)
        # colonne INTE_SIXX
        intsig1pm = (sig1 * weight).sum()

        # intfin(INTE_SIXX) = (COEF_MULT/VOLU_REFE*INTE_SIXX)**(1/bere_m)
        sigma_weibull = (coef_volu * intsig1pm) ** inv_m

        # sigwm(SIGMA_WEIBULL) = SIGMA_WEIBULL**bere_m
        sigma_weibullpm = sigma_weibull**pow_m

        # probaw(SIGMA_WEIBULL) = 1-exp(-SIGMA_WEIBULL**bere_m/sigma_u**bere_m)
        #  avec sigma_u=SIMG_CNV ou SIGM_REFE
        proba_weibull = 1.0 - exp(-((sigma_weibull / sigma_u) ** pow_m))

        return sigma_weibull, sigma_weibullpm, proba_weibull

    def store_sigm_maxi(self, idx, time, sig1, model):
        """Store the major stress ^M into the output result.

        Args:
            idx (int): Storage index.
            time (float): Time value.
            sig1 (*ComponentOnCells*): Major stress ^ M.
            model (*Model*): Model object.
        """
        if not self._rout:
            return
        if idx == 0:
            self._rout.allocate(self._rsieq.getNumberOfIndexes())
        sigmax = FieldOnCellsReal(model, "ELGA", "SIEF_R")
        sigmax.setValues(0.0)
        sfield = sigmax.toSimpleFieldOnCells()
        sixx = sfield.SIXX
        sixx += sig1.expand()
        sfield.setComponentValues(sixx)
        fed = model.getFiniteElementDescriptor()
        sigmax = sfield.toFieldOnCells(fed, "RAPH_MECA", "PCONTPR")
        self._rout.setField(sigmax, "SIEF_ELGA", idx)
        self._rout.setModel(model, idx)
        self._rout.setTime(time, idx)


def post_beremin_ops(self, RESULTAT, GROUP_MA, DEFORMATION, FILTRE_SIGM, **args):
    """Main of POST_BEREMIN command."""
    post = PostBeremin(result=RESULTAT, strain_type=DEFORMATION, stress_option=FILTRE_SIGM)
    post.set_zone(GROUP_MA)
    post.set_indexes(intvar_idx=args["LIST_NUME_VARI"])
    if DEFORMATION == "GDEF_LOG":
        post.set_indexes(stress_idx=args["LIST_NUME_SIEF"])
    post.use_history(args["HIST_MAXI"] == "OUI")
    post.set_coef_mult(args["COEF_MULT"])
    post.set_weibull_parameters(args["WEIBULL"])

    post.setup_calc_result()
    if args.get("SIGM_MAXI"):
        self._rout = NonLinearResult()
        self.register(self._rout, args["SIGM_MAXI"])

    table = post.main()
    sdtable = table.create_table()

    return sdtable


class TableBeremin:
    """Helper object to build the result table."""

    types = "IRKRRR"
    cols = ["NUME_ORDRE", "INST", "GROUP_MA", "SIGMA_WEIBULL", "SIGMA_WEIBULL**M", "PROBA_WEIBULL"]
    data = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self.data = {}
        for col in self.cols:
            self.data.setdefault(col, [])

    def append(self, *values):
        """Append the values of a row (columns order must match the definition).

        Args:
            values (list[misc]): Values of the row for each parameter.
        """
        assert len(values) == len(self.cols)
        for col, value in zip(self.cols, values):
            self.data[col].append(value)

    def create_table(self):
        """Create and return the Table datastructure.

        Returns:
            Table/table_sdaster: Table datastructure.
        """
        list_kws = []
        for col, typ in zip(self.cols, self.types):
            list_kws.append({"PARA": col, "LISTE_" + typ: self.data[col]})
        tab = CREA_TABLE(LISTE=list_kws)
        return tab
