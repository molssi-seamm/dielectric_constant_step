# -*- coding: utf-8 -*-

"""Non-graphical part of the Dielectric Constant step in a SEAMM flowchart

The static dielectric constant from the fluctuations of the dipole of the cell in
equilibrium MD with fixed charges and Ewald or PPPM electrostatics:
ε = 1 + (⟨M²⟩ − ⟨M⟩²) / (3 V ε0 kB T). For an electrolyte, M is the dipole of the
neutral molecules: the ions carry a current, not a bounded dipole. The analysis is
in :mod:`seamm_md_analysis.electric`.

The MD is run before this step -- in a Loop for independent runs -- or read from
files; the step analyzes it.
"""

from collections import Counter
import importlib.resources
import logging
from pathlib import Path
import warnings

import numpy as np
from seamm_md_analysis import (
    Estimate,
    component_check,
    isotropy,
    open_md,
)
from seamm_md_analysis.electric import (
    _factor,
    cell_dipole,
    classify_molecules,
    dielectric_constant,
    dipole_correlation,
)
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
from tabulate import tabulate

import dielectric_constant_step
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
# MDAnalysis warns about optional packages (netCDF4) in the step's output
logging.getLogger("MDAnalysis").setLevel(logging.ERROR)
job = printing.getPrinter()
printer = printing.getPrinter("Dielectric Constant")

# Add this module's properties to the standard properties
path = importlib.resources.files("dielectric_constant_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: The kspace styles with conducting boundaries, for which the formula holds
EWALD_STYLES = ("ewald", "pppm")
#: e·Å in debye
DEBYE = 4.80320471
#: The most frames used for the mean molecular dipole
MAX_DIPOLE_FRAMES = 200


def hill_formula(symbols):
    """The formula of a molecule in Hill order."""
    counts = Counter(symbols)
    if "C" in counts:
        order = ["C"] + (["H"] if "H" in counts else [])
        order += sorted(e for e in counts if e not in ("C", "H"))
    else:
        order = sorted(counts)
    return "".join(e + (str(counts[e]) if counts[e] > 1 else "") for e in order)


class DielectricConstant(seamm.Node):
    """
    The non-graphical part of a Dielectric Constant step in a flowchart.

    Attributes
    ----------
    parameters : DielectricConstantParameters
        The control parameters for Dielectric Constant.
    """

    def __init__(
        self,
        flowchart=None,
        title="Dielectric Constant",
        extension=None,
        logger=logger,
    ):
        logger.debug(f"Creating Dielectric Constant {self}")
        super().__init__(
            flowchart=flowchart,
            title="Dielectric Constant",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = dielectric_constant_step.metadata
        self.parameters = dielectric_constant_step.DielectricConstantParameters()
        self._reset()

    @property
    def version(self):
        """The semantic version of this module."""
        return dielectric_constant_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return dielectric_constant_step.__git_revision__

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.paths = []  # The md.h5md file of each run

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "dielectric-constant-step"
        parser = seamm_util.getParser()

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )
        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do."""
        if P is None:
            P = self.parameters.values_to_dict()

        what = "the static dielectric constant from the fluctuations of the dipole of "
        if P["dipole of"] == "automatic":
            what += "the molecules (the neutral molecules, if there are ions)"
        else:
            what += P["dipole of"]

        source = P["md source"]
        text = f"Calculate {what}, from "
        if source == read_files:
            text += f"the MD in {P['md files']}."
        else:
            text += f"{source} earlier in the flowchart."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Dielectric Constant step."""
        next_node = super().run(printer)
        self._reset()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration()

        self.references.cite(
            raw=self._bibliography["Neumann1983"],
            alias="Neumann1983",
            module="dielectric_constant_step",
            level=1,
            note="The dielectric constant from the fluctuations of the dipole.",
        )

        self.paths = md_inputs(self, P, name="all/position")
        printer.job(
            __(
                f"Analyzing {len(self.paths)} MD run"
                f"{'s' if len(self.paths) > 1 else ''}: "
                + ", ".join(str(p) for p in self.paths),
                indent=self.indent + 4 * " ",
            )
        )
        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        self.analyze(P=P, configuration=configuration)
        return next_node

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    def _read(self, path, configuration, molecules, charges):
        """One run: the dipole of the cell in each frame, the interval between
        frames (ps), the mean volume (Å³) and temperature (K), the mean dipole of a
        molecule (e·Å), the kspace style, and whether the run is isotropic."""
        atoms = np.concatenate(molecules)
        starts = np.cumsum([0] + [len(m) for m in molecules[:-1]])
        q = charges[atoms]
        with open_md(path) as md:
            if not md.periodic:
                raise RuntimeError(f"The MD in {path} is not periodic.")
            times = md.first["particles/all/position/time"]
            unit = times.attrs.get("unit", "ps")
            if isinstance(unit, bytes):
                unit = unit.decode()
            dt = (
                Q_(float(times[1] - times[0]), unit).m_as("ps")
                if times.shape[0] > 1
                else None
            )
            T = float(md.observable("temperature", units="K").values.mean())
            kspace = md.parameters.get("kspace_style")
            isotropic = isotropy(md)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = md.universe(configuration)
            n_frames = len(u.trajectory)
            every = max(1, n_frames // MAX_DIPOLE_FRAMES)
            M = []
            V = []
            mu = []
            for i, ts in enumerate(u.trajectory):
                M.append(cell_dipole(ts.positions, atoms, q))
                V.append(ts.volume)
                if i % every == 0:
                    # The positions are unwrapped, so each molecule is whole
                    dipoles = np.add.reduceat(q[:, None] * ts.positions[atoms], starts)
                    mu.append(np.linalg.norm(dipoles, axis=1).mean())
        return {
            "M": np.array(M),
            "dt": dt,
            "V": float(np.mean(V)),
            "T": T,
            "mu": float(np.mean(mu)),
            "kspace": kspace,
            "isotropic": isotropic,
        }

    @staticmethod
    def _sample(M, V, T, n_molecules, mu):
        """ε, its diagonal components, the Kirkwood factor and |⟨M⟩| (e·Å) of one
        run or block."""
        mean = M.mean(axis=0)
        variance = M.var(axis=0)
        factor = _factor(V, T)
        total = float(variance.sum())
        return {
            "epsilon": 1.0 + factor * total,
            "components": 1.0 + 3.0 * factor * variance,
            "g_K": total / (n_molecules * mu**2) if mu > 0 else float("nan"),
            "M": float(np.linalg.norm(mean)),
        }

    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the dielectric constant."""
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")
        molecules = configuration.find_molecules(as_indices=True)
        symbols = configuration.atoms.symbols

        with open_md(self.paths[0]) as md:
            charges = md.charges("all")
        if charges is None:
            raise RuntimeError(
                f"'{self.title}': the MD has no fixed charges, as with a "
                "machine-learned forcefield, so there is no dipole to analyze."
            )
        charges = np.asarray(charges, dtype=float)

        classes = classify_molecules(charges, molecules)
        mode = P["dipole of"]
        neutral_only = mode == "the neutral molecules" or (
            mode == "automatic" and len(classes["ions"]) > 0
        )
        chosen = classes["neutral"] if neutral_only else list(range(len(molecules)))
        chosen = [m for m in chosen if m not in classes["other"]]
        if len(chosen) == 0:
            raise RuntimeError(f"'{self.title}': there are no neutral molecules.")

        notes = []
        names = Counter(
            hill_formula([symbols[i] for i in molecules[m]]) for m in chosen
        )
        text = "The dipole of " + ", ".join(f"{n} {name}" for name, n in names.items())
        if neutral_only and classes["ions"]:
            ions = Counter(
                hill_formula([symbols[i] for i in molecules[m]])
                for m in classes["ions"]
            )
            text += (
                "; the ions ("
                + ", ".join(f"{n} {name}" for name, n in ions.items())
                + ") are left out: they carry a current, not a bounded dipole"
            )
        printer.normal(__(text + ".", indent=self.indent + 4 * " "))
        if classes["other"]:
            notes.append(
                f"{len(classes['other'])} molecules have a total charge that is "
                "not an integer, which is a problem with the charges; they are "
                "left out."
            )

        chosen_molecules = [np.array(molecules[m]) for m in chosen]
        runs = [
            self._read(p, configuration, chosen_molecules, charges) for p in self.paths
        ]
        nruns = len(runs)
        n_molecules = len(chosen)

        # The independent samples: the runs, or blocks of a single run
        if nruns > 1:
            source = "runs"
            samples = [
                self._sample(r["M"], r["V"], r["T"], n_molecules, r["mu"]) for r in runs
            ]
        else:
            source = "blocks"
            run = runs[0]
            n_blocks = max(1, min(int(P["number of blocks"]), len(run["M"]) // 2))
            samples = [
                self._sample(block, run["V"], run["T"], n_molecules, run["mu"])
                for block in np.array_split(run["M"], n_blocks)
            ]

        def estimate(key, scale=1.0):
            values = [s[key] * scale for s in samples]
            if len(values) < 2:
                return Estimate(float(values[0]), None, 1, source=source)
            return Estimate.from_samples(values, source=source)

        # ε from the whole of each run, the error from the samples
        fits = [dielectric_constant(r["M"], r["V"], r["T"], n_blocks=1) for r in runs]
        epsilon = estimate("epsilon")
        epsilon.value = float(np.mean([f["epsilon"] for f in fits]))
        g_K = estimate("g_K")
        M_mean = estimate("M", DEBYE)
        mu = Estimate.from_samples([r["mu"] * DEBYE for r in runs])
        components = np.array([s["components"] for s in samples])
        check = component_check(components, names=["ε_xx", "ε_yy", "ε_zz"])
        isotropic, why = runs[0]["isotropic"]

        dt = runs[0]["dt"]
        tau = None
        try:
            acf = dipole_correlation(runs[0]["M"], dt or 1.0)
            tau = acf["tau"]
        except Exception as e:
            acf = None
            logger.warning(f"No autocorrelation of the dipole: {e}")
        n_frames = sum(len(r["M"]) for r in runs)

        table = {"Property": [], "Value": [], "± (95%)": [], "Units": []}

        def row(label, value, units=""):
            if isinstance(value, Estimate):
                v, e = value.formatted()
            else:
                v, e = ("" if value is None else f"{value:.3g}"), ""
            table["Property"].append(label)
            table["Value"].append(v)
            table["± (95%)"].append(e)
            table["Units"].append(units)

        row("ε", epsilon)
        if not isotropic:
            for name, i in (("ε_xx", 0), ("ε_yy", 1), ("ε_zz", 2)):
                row(
                    name,
                    (
                        Estimate.from_samples(components[:, i], source=source)
                        if len(components) > 1
                        else float(components[0, i])
                    ),
                )
        row("Molecular dipole μ", mu, "D")
        row("Kirkwood factor g_K", g_K)
        row("|⟨M⟩| of the cell", M_mean, "D")
        row("Dipole correlation time", tau, "ps")
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "decimal", "decimal", "left"),
        )
        length = len(text.splitlines()[0])
        printer.normal(
            __(
                "\n" + "Static dielectric constant".center(length) + "\n" + text + "\n",
                indent=8 * " ",
                wrap=False,
                dedent=False,
            )
        )
        lines = [
            f"From {n_frames} frames of {nruns} run{'s' if nruns > 1 else ''}"
            + (f", every {dt:.3g} ps" if dt else "")
            + f", for {n_molecules} molecules. "
            + epsilon.describe_error()
        ]
        if isotropic:
            lines.append(check.text)
        else:
            lines.append(
                f"The run is not treated as isotropic ({why}), so the diagonal "
                "components are given; ε is a third of their sum."
            )
        kspace = runs[0]["kspace"]
        lines.append(
            "The fluctuation formula assumes Ewald or PPPM electrostatics with "
            "conducting (tin-foil) boundaries"
            + ("." if kspace is None else f"; the runs used {kspace}.")
        )
        if kspace is not None and not kspace.split()[0].startswith(EWALD_STYLES):
            notes.append(
                f"The runs used kspace_style '{kspace}', so the value may be wrong."
            )
        if dt is not None and tau is not None:
            length_ps = len(runs[0]["M"]) * dt
            if length_ps < 50 * tau:
                notes.append(
                    f"Each run is {length_ps:.0f} ps, only {length_ps / tau:.0f} times "
                    f"the correlation time of the dipole ({tau:.1f} ps). The "
                    "dielectric constant converges slowly: use longer or more runs."
                )
        for line in lines + notes:
            printer.normal(__(line, indent=8 * " "))
        printer.normal("")
        self._print_samples(samples, source)

        state = report_state(self, printer, self.paths)

        self._plot(fits, runs, acf, dt)

        self.model = self._find_model()
        data = {
            "epsilon": epsilon.value,
            "dipole moment": mu.value,
            "Kirkwood factor": g_K.value,
            "nruns": nruns,
            "n_frames": n_frames,
        }
        if epsilon.stderr is not None:
            data["epsilon,stderr"] = epsilon.stderr
        if g_K.stderr is not None:
            data["Kirkwood factor,stderr"] = g_K.stderr
        if tau is not None:
            data["dipole correlation time"] = tau
        if state is not None:
            data.update(state_data(state))
        self.store_results(configuration=configuration, data=data)

    def _print_samples(self, samples, source):
        """The values of each run or block, in the detailed output."""
        what = "Run" if source == "runs" else "Block"
        table = {what: [], "ε": [], "ε_xx": [], "ε_yy": [], "ε_zz": [], "g_K": []}
        for i, s in enumerate(samples, start=1):
            table[what].append(i)
            table["ε"].append(f"{s['epsilon']:.2f}")
            for name, value in zip(("ε_xx", "ε_yy", "ε_zz"), s["components"]):
                table[name].append(f"{value:.2f}")
            table["g_K"].append(f"{s['g_K']:.3f}")
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("center", "decimal", "decimal", "decimal", "decimal", "decimal"),
        )
        length = len(text.splitlines()[0])
        title = f"By {what.lower()}".center(length)
        printer.verbose(
            __("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)
        )

    def _plot(self, fits, runs, acf, dt):
        """The running estimate of ε, and the autocorrelation of the dipole."""
        figure = self.create_figure(
            module_path=("seamm",),
            template="line.graph_template",
            title="Dielectric constant",
        )
        plot = figure.add_plot("running")
        x_axis = plot.add_axis("x", label="Time (ps)")
        y_axis = plot.add_axis("y", label="ε", anchor=x_axis)
        x_axis.anchor = y_axis
        for i, (fit, run) in enumerate(zip(fits, runs), start=1):
            counts, values = fit["running"]
            step = run["dt"] or 1.0
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=f"run {i}",
                x=[c * step for c in counts],
                xlabel="t",
                xunits="ps",
                y=list(values),
                ylabel="ε",
                yunits="",
            )
        if acf is not None:
            plot2 = figure.add_plot("acf")
            x2 = plot2.add_axis("x", label="Time (ps)")
            y2 = plot2.add_axis("y", label="<M(0)M(t)>/<M²>", anchor=x2)
            x2.anchor = y2
            plot2.add_trace(
                x_axis=x2,
                y_axis=y2,
                name="ACF",
                x=acf["times"],
                xlabel="t",
                xunits="ps",
                y=acf["acf"],
                ylabel="ACF",
                yunits="",
                color="red",
            )
            figure.grid_plots("running - acf")
        else:
            figure.grid_plots("running")
        path = Path(self.directory) / "dielectric_constant.graph"
        figure.dump(path)
        if "html" in self.options and self.options["html"]:
            figure.template = "line.html_template"
            figure.dump(path.with_suffix(".html"))

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    # ------------------------------------------------------------------
    # What the MD must provide
    # ------------------------------------------------------------------
    def md_requirements(self, P=None):
        """What the analysis needs from the MD: positions every 0.5 ps or so, in a
        periodic cell. It analyzes fluctuations, not the dynamics."""
        interval = (
            P["frame interval"] if P is not None else self._static("frame interval")
        )
        if interval is None or isinstance(interval, str):
            interval = Q_(0.5, "ps")
        elif not hasattr(interval, "m_as"):
            interval = Q_(float(interval), "ps")
        return Requirements(
            needs=[Need("all/position", max_time=interval, why="the dipole")],
            periodic=True,
            dynamics=False,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and is
        set up for the analysis."""
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small.
        """
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the fluctuation formula is for a periodic cell."
            )
        return [f"'{self.title}': {w}" for w in size_warnings(configuration)]
