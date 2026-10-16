# -*- coding: utf-8 -*-

"""Non-graphical part of the Dielectric Constant step in a SEAMM flowchart

The static dielectric constant from the fluctuations of the dipole of the cell in
equilibrium MD with fixed charges and Ewald or PPPM electrostatics:
ε = 1 + (⟨M²⟩ − ⟨M⟩²) / (3 V ε0 kB T). For an electrolyte, M is the dipole of the
neutral molecules: the ions carry a current, not a bounded dipole. The analysis is
in :mod:`seamm_md_analysis.electric`.
"""

from collections import Counter
import importlib.resources
import logging
from pathlib import Path
import sys
import traceback
import warnings

import numpy as np
from seamm_md_analysis import open_md
from seamm_md_analysis.electric import (
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
    request_sampling,
    run_own_md,
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
    subflowchart : seamm.Flowchart
        The MD to run, when this step runs its own MD.

    parameters : DielectricConstantParameters
        The control parameters for Dielectric Constant.
    """

    def __init__(
        self,
        flowchart=None,
        title="Dielectric Constant",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        logger.debug(f"Creating Dielectric Constant {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Dielectric Constant", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Dielectric Constant",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = dielectric_constant_step.metadata
        self.parameters = dielectric_constant_step.DielectricConstantParameters()
        self._file_handler = None
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

        self.subflowchart._parser = self.flowchart._parser

        parser_exists = parser.exists(parser_name)
        super().create_parser(name=parser_name)
        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

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

        source = P.get("md source", run_own_md)
        if source == run_own_md:
            text = f"Calculate {what}, from {P['nruns']} runs of the MD below.\n\n"
        else:
            text = f"Calculate {what}, from "
            if source == read_files:
                text += f"the MD in {P['md files']}.\n\n"
            else:
                text += f"{source} earlier in the flowchart.\n\n"
            short = True

        if not short:
            self.subflowchart.root_directory = self.flowchart.root_directory
            node = self.subflowchart.get_node("1").next()
            while node is not None:
                text += __(node.description_text()).__str__()
                text += "\n"
                node = node.next()

        return self.header + "\n" + __(text, indent=4 * " ").__str__()

    def run(self):
        """Run a Dielectric Constant step."""
        next_node = super().run(printer)
        self._reset()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        _, configuration = self.get_system_configuration()
        if (
            P["md source"] == run_own_md
            and configuration is not None
            and configuration.periodicity != 3
        ):
            raise RuntimeError(
                "The Dielectric Constant step needs a periodic system: the "
                "fluctuation formula is for a periodic cell."
            )

        self.references.cite(
            raw=self._bibliography["Neumann1983"],
            alias="Neumann1983",
            module="dielectric_constant_step",
            level=1,
            note="The dielectric constant from the fluctuations of the dipole.",
        )

        handlers = self._quiet_handlers()
        try:
            if P["md source"] == run_own_md:
                self._run_own_md(P)
            else:
                self.paths = md_inputs(self, P, name="all/position")
                printer.job(
                    __(
                        f"Analyzing {len(self.paths)} MD run"
                        f"{'s' if len(self.paths) > 1 else ''}: "
                        + ", ".join(str(p) for p in self.paths),
                        indent=self.indent + 4 * " ",
                    )
                )
        finally:
            self._restore_handlers(handlers)

        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        self.analyze(P=P, configuration=configuration)
        return next_node

    def _quiet_handlers(self):
        """Send only job-level output to job.out and the screen while MD runs."""
        result = []
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ) or isinstance(handler, logging.StreamHandler):
                result.append((handler, handler.level))
                handler.setLevel(printing.JOB)
        return result

    @staticmethod
    def _restore_handlers(handlers):
        for handler, level in handlers:
            handler.setLevel(level)

    def _run_own_md(self, P):
        """Run the MD in the sub-flowchart several times, keeping each run's file."""
        changes, missing = request_sampling(self.subflowchart, self.md_requirements(P))
        for change in changes:
            printer.job(__(f"Note: {change}.", indent=self.indent + 4 * " "))
        if len(missing) > 0:
            raise RuntimeError(
                "The MD in the sub-flowchart cannot provide what is needed: "
                + " ".join(missing)
            )

        first_node = self.subflowchart.get_node("1").next()
        if first_node is None:
            raise RuntimeError(
                "The sub-flowchart of the Dielectric Constant step is empty."
            )
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()
        self.subflowchart.executor = self.flowchart.executor

        nruns = P["nruns"]
        fmt_run = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            run_dir = Path(self.directory) / f"run_{run:{fmt_run}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt_run}}"))

            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.paths.append(self._read_run(run, run_dir))

        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None

    def _read_run(self, run, run_dir):
        """The md.h5md of one run of the sub-flowchart with the positions."""
        paths = []
        for path in sorted(run_dir.glob("**/md.h5md")):
            with open_md(path) as md:
                if md.has_particles("position"):
                    paths.append(path)
        if len(paths) == 0:
            raise RuntimeError(f"There are no positions for run {run} in {run_dir}.")
        return paths[-1]

    # ------------------------------------------------------------------
    # The analysis
    # ------------------------------------------------------------------
    def _read(self, path, configuration, atoms, charges):
        """The dipole in each frame, the interval between frames (ps), and the
        mean volume (Å³) and temperature (K) of one run."""
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
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = md.universe(configuration)
            M = []
            V = []
            for ts in u.trajectory:
                M.append(cell_dipole(ts.positions, atoms, charges))
                V.append(ts.volume)
        return np.array(M), dt, float(np.mean(V)), T

    def analyze(self, indent="", P=None, configuration=None, **kwargs):
        """Analyze, print, plot and store the dielectric constant."""
        if configuration is None:
            raise RuntimeError("There is no structure to analyze.")
        molecules = configuration.find_molecules(as_indices=True)
        symbols = configuration.atoms.symbols

        with open_md(self.paths[0]) as md:
            charges = md.charges("all")
            kspace = md.parameters.get("kspace_style")
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
        if len(chosen) == 0:
            raise RuntimeError(f"'{self.title}': there are no neutral molecules.")
        atoms = np.array(sorted(i for m in chosen for i in molecules[m]))

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
            atoms = np.array(
                [
                    i
                    for i in atoms
                    if all(i not in molecules[m] for m in classes["other"])
                ]
            )
        if kspace is None:
            notes.append(
                "The MD does not record its electrostatics. The formula assumes Ewald "
                "or PPPM summation with conducting boundaries."
            )
        elif not kspace.split()[0].startswith(EWALD_STYLES):
            notes.append(
                f"The MD used kspace_style '{kspace}'. The formula holds for Ewald or "
                "PPPM summation with conducting boundaries, so the value may be "
                "wrong."
            )

        runs = [self._read(p, configuration, atoms, charges[atoms]) for p in self.paths]
        nruns = len(runs)
        n_blocks = P["number of blocks"] if nruns == 1 else 1
        if nruns == 1 and P["md source"] == run_own_md:
            n_blocks = 1
        fits = [dielectric_constant(M, V, T, n_blocks=n_blocks) for M, _, V, T in runs]
        values = np.array([f["epsilon"] for f in fits])
        epsilon = float(values.mean())
        if nruns > 1:
            stderr = float(values.std(ddof=1) / np.sqrt(nruns))
            how = f"the standard error over the {nruns} runs"
        else:
            stderr = fits[0]["stderr"]
            how = (
                f"from {fits[0]['n_blocks']} blocks of the run"
                if stderr is not None
                else "no error estimate: use several runs or blocks"
            )
        T = float(np.mean([r[3] for r in runs]))
        V = float(np.mean([r[2] for r in runs]))
        dt = runs[0][1]
        n_frames = sum(len(r[0]) for r in runs)

        tau = None
        try:
            acf = dipole_correlation(runs[0][0], dt or 1.0)
            tau = acf["tau"]
        except Exception as e:
            acf = None
            logger.warning(f"No autocorrelation of the dipole: {e}")

        table = {
            "Property": ["ε", "Dipole correlation time", "Temperature", "Volume"],
            "Value": [
                f"{epsilon:.1f}",
                "" if tau is None else f"{tau:.2f}",
                f"{T:.2f}",
                f"{V:.1f}",
            ],
            "±": ["" if stderr is None else f"{stderr:.1f}", "", "", ""],
            "Units": ["", "ps", "K", "Å^3"],
        }
        text = tabulate(
            table, headers="keys", tablefmt="simple_outline", disable_numparse=True
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
        printer.normal(
            __(
                f"From {n_frames} frames of {nruns} run{'s' if nruns > 1 else ''}"
                + (f", every {dt:.3g} ps" if dt else "")
                + f"; the error is {how}.",
                indent=8 * " ",
            )
        )
        if dt is not None and tau is not None:
            length_ps = len(runs[0][0]) * dt
            if length_ps < 50 * tau:
                notes.append(
                    f"Each run is {length_ps:.0f} ps, only {length_ps / tau:.0f} times "
                    f"the correlation time of the dipole ({tau:.1f} ps). The "
                    "dielectric constant converges slowly: use longer or more runs."
                )
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        printer.normal("")

        self._plot(fits, runs, acf, dt)

        self.model = self._find_model()
        data = {
            "epsilon": epsilon,
            "T": T,
            "V": V,
            "nruns": nruns,
            "n_frames": n_frames,
        }
        if stderr is not None:
            data["epsilon,stderr"] = stderr
        if tau is not None:
            data["dipole correlation time"] = tau
        self.store_results(configuration=configuration, data=data)

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
            step = run[1] or 1.0
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
        periodic cell."""
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
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed."""
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze."""
        if self._static("md source") == run_own_md:
            return
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the fluctuation formula is for a periodic cell."
            )

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
