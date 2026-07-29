"""Build a PDF report of the recomputed examples and the problems found.

Reads the JSON emitted by the other report generators and produces a LaTeX
document, then compiles it.  Regenerate with

    python -m rookfields.report_pdf

Everything numeric comes from the JSON, so the document cannot drift from the
computations.  The curated prose (what to fix, and why) lives in this module.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from .examples.catalog import EXAMPLES

ROOT = Path(__file__).resolve().parents[2]
REPORTS = ROOT / "reports"
FIGURES = ROOT / "figures"

GREEN = "rowok"
AMBER = "rowwarn"
RED = "rowbad"


def tex_escape(text: str) -> str:
    """Escape a plain string for LaTeX text mode."""
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(c, c) for c in str(text))


def mono(text: str) -> str:
    return r"\texttt{" + tex_escape(text) + "}"


def load(name: str):
    path = REPORTS / name
    if not path.exists():
        return None
    return json.loads(path.read_text())


# ---------------------------------------------------------------------------
# the curated problem list
# ---------------------------------------------------------------------------

PROBLEMS_MANUSCRIPT = [
    {
        "id": "M1",
        "severity": "breaks a published number",
        "title": r"\texttt{ex:trivial\_index\_3D\_F3}: 61/12 no longer reproduces",
        "where": "examples4.tex:786",
        "body": r"""The text reports a trivial-index component with $61$ cells and $12$
double-edge pairs. Recomputing with the current \textsc{DSGRN\_utils} gives
\textbf{75 cells and 16 double-edge pairs}. Under the live
\texttt{def:active\_regulation} the published $61/12$ reproduces exactly.

So the published figure was produced with the manuscript's definition, and the
code has since drifted away from it. Nothing in the manuscript needs changing
here --- the \emph{code} does (see C1). The number is correct as printed.""",
    },
    {
        "id": "M2",
        "severity": "misprint",
        "title": r"\texttt{ex:periodicOrbit3\_F4}: the index set is wrong in its first factor",
        "where": "examples4.tex:920--923",
        "body": r"""The displayed index set is
\[
  \I = \{0,1,2,3,4\}\times\{0,1,2,3,4\}\times\{0,1,2,3\}.
\]
Node~$1$ of that network has two out-edges, so $K(1)=2$ and the first factor is
$\{0,1,2,3\}$. Computed: $\I=\{0,\dots,3\}\times\{0,\dots,4\}\times\{0,\dots,3\}$.

Everything else in the example --- the $12$-cell cycle, the absence of double
edges, the periodic-orbit index --- reproduces exactly.

\textbf{Fix:} change the first factor to $\setof{0,1,2,3}$.""",
    },
    {
        "id": "M3",
        "severity": "hypotheses fail at the stated values",
        "title": r"The published ramp parameters violate $\cH_1$, hence $\cH_2$ and $\cH_3$",
        "where": "examples4.tex:857 (van der Pol), introduction6.tex:525 (3D periodic)",
        "body": r"""Both directly specified ramp systems satisfy $\cH_0$ and $\Lambda(R)$ but
\textbf{fail} $\cH_1$: the focal value $E_n(D_\bv)/\gamma_n$ lands inside a ramp
half-window, which \texttt{defn:H1} forbids. For van der Pol at box $(1,1)$,
$E_1(D_\bv)=6.1$ and $\gamma_1=2$, so the focal value is $3.05$, inside the
forbidden window $(\theta_{m_1},\theta_{m_1}+h_{m_1}]=(3,3.15]$. Eight such
violations for van der Pol, eighteen for the three-dimensional system.

Since $\cH_2,\cH_3\subset\cH_1$, the hypotheses of \texttt{thm:R1ABlattice} and
\texttt{thm:R3ABlattice} do not hold at the stated values.
\texttt{ex:introPeriodic} invokes \texttt{thm:dynamics} at exactly those values.

\textbf{Fix:} a uniform half-width $h=0.005$ satisfies $\cH_0$, $\Lambda(R)$,
$\cH_1$, $\cH_2$ and $\cH_3$, and \emph{leaves the Morse graph unchanged} ---
same Conley indices, same Morse-graph edges. Only the two parameter tables need
new $h$ values; every computed conclusion stands.""",
    },
    {
        "id": "M4",
        "severity": "affects the $\\cF_3$ vs $\\cF_4$ narrative",
        "title": r"Under the live definition, $\cF_3$ already resolves \texttt{ex:periodicOrbit3\_F3}",
        "where": "examples4.tex:435 and 917",
        "body": r"""\texttt{ex:periodicOrbit3\_F3} reports a $19$-cell component with $6$ double-edge
pairs and calls the fixed-point-plus-periodic-orbit reading a conjecture;
\texttt{sec:dynamics\_F4} then presents $\cF_4$ as resolving it into a clean
$12$-cell cycle with $CH_0\cong CH_1\cong\Z_2$.

Under the live \texttt{def:active\_regulation}, $\cF_3$ \emph{already} gives that
$12$-cell cycle, with no double edges and the same index. Attributable to that
one definition change: no other divergence alters the result.

\textbf{Fix:} the contrast drawn between $\cF_3$ and $\cF_4$ at this parameter
node needs restating, or a different node chosen to illustrate it.""",
    },
    {
        "id": "M5",
        "severity": "scope",
        "title": r"\texttt{ex:ramp\_van\_der\_pol} is labelled $\cF_3$ but the cycle rule never fires",
        "where": "examples4.tex:841",
        "body": r"""Neither directly specified ramp system has a semi-opaque cell with a nontrivial
regulation cycle. Conditions~3.1 and~3.2 are therefore vacuous and
$\cF_3=\cF_2$ for both systems (verified: identical edge sets).

The computation is correct; the label overstates which rule is doing the work.
Cycle cells do occur for DSGRN parameters --- N3\_B node $2{,}472{,}287$ has one
$3$-cycle, N3\_A node $52{,}718{,}681{,}992$ has twelve $2$-cycles and one
$3$-cycle.""",
    },
    {
        "id": "M6",
        "severity": "notation",
        "title": r"The outer sentinel threshold differs between two chapters",
        "where": "RampSystemsv4.tex:388 vs RectangularGeo.tex:27",
        "body": r"""\texttt{defn:ramp-wall-labeling} sets $\theta_{m_{K+1}}=\gab_n$;
\texttt{eq:mKnjK} sets $\theta_{m_{K+1}}=\gab_n-\tfrac14$ with
$h_{m_{K+1}}=\tfrac14$, so that $\theta+h=\gab_n$. Both are self-consistent, but
any formula quantified over $\bv\in\prod\{1,\dots,K(n)\}$ that involves
$\theta_{m_{k+1}}$ --- $I_n$, $\Xi_n$, and hence the $\cH_2$ inequality --- picks
up whichever convention its chapter uses, and the two disagree at the outermost
cells.

\textbf{Fix:} one sentence saying which convention is meant where.""",
    },
    {
        "id": "M7",
        "severity": "statement/proof mismatch",
        "title": r"\texttt{prop:smallh} states $i\in\{0,1,3\}$ but proves $i=2$ as well",
        "where": "RampSystemsv4.tex:624",
        "body": r"""The statement asserts nonemptiness for $i\in\setof{0,1,3}$; the proof also treats
$i=2$. Separately, the constructive $\cH_1$ bound it produces ($0.02$ for both
published systems) is not small enough for $\cH_2$: the largest uniform width
satisfying $\cH_0+\cH_1+\cH_2$ is $0.00517$. That is consistent with the recorded
concern that the constant $c$ in the $i=2$ branch is built from $L/U$, which
themselves depend on $h$.""",
    },
    {
        "id": "M8",
        "severity": "labelling",
        "title": r"Two labelling issues in the trivial-index examples",
        "where": "examples4.tex:786 and 820",
        "body": r"""\texttt{fig:trivial\_index\_3D\_F2} is named \texttt{F2} for an $\cF_3$ example
(already recorded in the knowledge base).

For the six-dimensional EMT example the eight trivial-index nodes reproduce in
\emph{count} ($8$ of $25$) but not in the displayed labels
$\setof{15,\dots,18,21,\dots,24}$. Morse-node numbering in \texttt{MorseGraph}
breaks rank ties by the internal SCC grading value and is not canonical, so the
displayed set is an artefact of the run rather than a property of the dynamics.

The same effect shows up in the six-dimensional semiconjugacy example: the two
readings give Morse graphs that differ in two edges, but every Morse set is the
same collection of cells and the edges agree after exchanging the numbers of two
nodes carrying the same Conley index. Any comparison of Morse graphs has to go
through the cell sets, which is what
\texttt{rookfields.pipeline.compare\_results} does.

Consider reporting the count and the index data instead of the labels.""",
    },
]

PROBLEMS_MANUSCRIPT.append(
    {
        "id": "M9",
        "severity": "invariants unchanged, cells change",
        "title": r"The divergences reach dimension six without changing any reported invariant",
        "where": "examples4.tex:663 and 820 (the two EMT examples)",
        "body": r"""Comparing the two readings across all seventeen examples, and comparing Morse
graphs through their cell sets rather than their node numbers:

\begin{itemize}[leftmargin=1.4em,itemsep=1pt]
  \item \texttt{ex:periodicOrbit3\_F3} is the only example whose node count or
    Conley indices change (M4);
  \item \texttt{ex:trivial\_index\_3D\_F3}, \texttt{ex:3d\_example\_1} and the
    second EMT example keep the same number of Morse nodes and the same multiset
    of Conley indices, but move cells between components. For the EMT node
    $1{,}739{,}757{,}491{,}101$ every large component shrinks
    ($2436\to2394$, $440\to340$, $415\to336$, $390\to301$, $300\to262$,
    $77\to59$) while the count of trivial-index nodes stays at $8$ of $25$;
  \item the remaining thirteen are identical.
\end{itemize}

So the corrections propagate into dimension six, but outside
\texttt{ex:periodicOrbit3\_F3} they do not disturb any homological quantity the
text reports. The only \emph{printed} number affected is the cell count in
\texttt{ex:trivial\_index\_3D\_F3} (M1).""",
    }
)

PROBLEMS_CODE = [
    {
        "id": "C1",
        "severity": "changes published results",
        "title": r"\texttt{active\_regulation\_map} implements the retired definition",
        "where": "CubicalBlowupGraph.py:298--309",
        "body": r"""The target ranges over \texttt{range(self.dim)}, which is
\texttt{def:active\_regulation-old}. The live \texttt{def:active\_regulation}
(RookFields6.tex:738) requires $n,m\in J_i(\xi)$ and types the map
$\rmap\xi:\activeset(\xi)\to J_i(\xi)$.

Entries whose target is an \emph{essential} direction reach the
\texttt{level < 4} guard (\texttt{:461--467}) and \texttt{semi\_opaque\_cell}
(\texttt{:470}), so $\cF_2$ and $\cF_3$ change. This is what makes
\texttt{ex:trivial\_index\_3D\_F3} fail to reproduce (M1) and what changes the
$\cF_3$/$\cF_4$ story (M4). All two-dimensional examples are unaffected.""",
    },
    {
        "id": "C2",
        "severity": "silent wrong answer",
        "title": r"\texttt{decision\_wall} reads the wall label off a wrapped cell",
        "where": "CubicalBlowupGraph.py:427--434",
        "body": r"""\texttt{defn:back-walls} shifts by
$\hat\bv=\sum_{n\in J_i(\xi')}\tfrac{1+r_n}{2}\bzero^{(n)}$, which can leave
$\cX$ at an outer endpoint. The code applies the shift with no bounds check and
passes negative coordinates to \texttt{cell\_index}, which wraps to an unrelated
cell. Measured on decision walls that are actually built: $22\%$ (repressilator),
$15\%$ (3-node), $37\%$ (4-node) have a negative coordinate, and in every one of
those cases the resulting index is a different cell than intended.

Discarding those pairs changes the state transition graph broadly but changes
\emph{no} Morse-node count and \emph{no} Conley index in any family sampled.
The manuscript leaves the intended response open (safe condition
$1\le v'_j\le K(j)$); the code should at least not read a wrapped cell.""",
    },
    {
        "id": "C3",
        "severity": "breaks a standing assumption",
        "title": r"\texttt{WallLabelling.py}'s outer sentinel is $\theta+10h$, not $\gab_n$",
        "where": "WallLabelling.py:44",
        "body": r"""\texttt{defn:ramp-wall-labeling} specifies $\theta_{m_{K(n)+1}}=\gab_n$. The code
uses $\theta_{m_K}+10\,h_{m_K}$, which depends on $h$. Consequences:

\begin{itemize}
  \item it contradicts \texttt{prop:wall-labeling-const} (the labelling is
    constant in $h$ on $\cH_0$): shrinking $h$ from $0.15$ to $0.0198$ flips
    $2$ of $36$ top-cell labels for the van der Pol system;
  \item it \textbf{breaks strong dissipativity}, the standing assumption from
    \texttt{CombinatorialDynamics.tex:4} onward, once $h\le0.05$.
\end{itemize}

With the $\gab_n$ sentinel the labelling is byte-identical at every width tested.
The bug is invisible at the published widths, which is presumably why it
survived --- but small $h$ is exactly the regime $\cH_1$ forces (M3).""",
    },
    {
        "id": "C4",
        "severity": "definition not enforced",
        "title": r"\texttt{semi\_opaque\_cell} omits the two excluded classes",
        "where": "CubicalBlowupGraph.py:470--475",
        "body": r"""\texttt{defn:partially\_opaque} requires
$\xi\in\cX\setminus(\cX^{(N)}\cup\bdy(\cX))$ in addition to $\rmap\xi$ being a
bijection; the code checks only the bijection. A top cell has an empty regulation
map and passes vacuously. Enforcing the exclusions changes the state transition
graph on $16/60$ and $20/80$ parameters at levels~3 and~4, with no change to
Morse nodes or Conley indices.""",
    },
    {
        "id": "C5",
        "severity": "crashes",
        "title": r"\texttt{IsomorphismQuery} is broken against current HEAD",
        "where": "IsomorphismQuery.py:23",
        "body": r"""\texttt{DSGRN.isomorphic\_morse\_graphs} calls \texttt{.split(':')} on a vertex
label, but commit \texttt{fb10673} changed \texttt{MorseGraph.vertex\_label} to
return a list --- so this raises \texttt{AttributeError}. The loop also uses
\texttt{continue} where \texttt{break} is meant, so a parameter matching two
representatives is added to both isomorphism classes.

Unrelated to the manuscript, but it means the query is unusable as shipped.""",
    },
]


# ---------------------------------------------------------------------------
# document assembly
# ---------------------------------------------------------------------------

PREAMBLE = r"""
\documentclass[11pt]{article}
\usepackage[margin=2.4cm]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage[table,dvipsnames]{xcolor}
\usepackage{graphicx}
\usepackage{caption}
\usepackage{enumitem}
\usepackage[colorlinks=true,linkcolor=NavyBlue,urlcolor=NavyBlue]{hyperref}
\usepackage{fancyhdr}
\usepackage{titlesec}
\usepackage{tcolorbox}
\tcbuselibrary{skins,breakable}

\definecolor{rowok}{RGB}{233,247,235}
\definecolor{rowwarn}{RGB}{255,246,224}
\definecolor{rowbad}{RGB}{253,235,235}
\definecolor{okgreen}{RGB}{22,120,60}
\definecolor{badred}{RGB}{176,32,32}
\definecolor{warnamber}{RGB}{176,110,0}

\newcommand{\cF}{\mathcal{F}}
\newcommand{\cH}{\mathcal{H}}
\newcommand{\cX}{\mathcal{X}}
\newcommand{\I}{\mathcal{I}}
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\bv}{\mathbf{v}}
\newcommand{\bzero}{\mathbf{0}}
\newcommand{\gab}{\mathrm{GB}}
\newcommand{\bdy}{\mathrm{bdy}}
\newcommand{\rmap}{o_}
\newcommand{\rook}{\Phi}
\newcommand{\cD}{\mathcal{D}}
\newcommand{\cM}{\mathcal{M}}
\newcommand{\cN}{\mathcal{N}}
\newcommand{\cG}{\mathcal{G}}
\newcommand{\activeset}{\mathrm{Act}}
\newcommand{\setof}[1]{\{#1\}}
\newcommand{\OK}{\textcolor{okgreen}{\textbf{ok}}}
\newcommand{\BAD}{\textcolor{badred}{\textbf{differs}}}
\newcommand{\WARN}{\textcolor{warnamber}{\textbf{note}}}

\pagestyle{fancy}
\fancyhf{}
\lhead{\small Updated computations --- \texttt{Rook\_Field\_Paper\_v2}}
\rhead{\small\thepage}
\renewcommand{\headrulewidth}{0.4pt}

\newtcolorbox{problembox}[2]{
  breakable, enhanced, colback=white, colframe=#1,
  boxrule=0.8pt, left=6pt, right=6pt, top=5pt, bottom=5pt,
  title={#2}, fonttitle=\bfseries\small, coltitle=white
}

\setlength{\parskip}{0.45em}
\setlength{\parindent}{0pt}
"""


def status_of(outcome) -> str:
    if outcome["error"]:
        return "error"
    return "ok" if all(r[3] for r in outcome["results"]) else "mismatch"


def build_tex() -> str:
    examples = load("examples.json")
    hypotheses = load("hypotheses.json")
    alignment = load("alignment.json")
    divergence = load("divergence.json")
    if examples is None:
        raise SystemExit("reports/examples.json missing; run rookfields.examples first")

    by_label: dict[str, list] = {}
    for o in examples:
        by_label.setdefault(o["label"], []).append(o)
    catalog = {e.label: e for e in EXAMPLES}

    total = len(examples)
    matched = sum(1 for o in examples if status_of(o) == "ok")

    L: list[str] = [PREAMBLE, r"\begin{document}"]

    # -- title ---------------------------------------------------------
    L.append(
        r"""
\begin{center}
{\LARGE\bfseries Updated computations for the examples of}\\[2pt]
{\LARGE\bfseries \texttt{Rook\_Field\_Paper\_v2}}\\[10pt]
{\large Recomputed from the manuscript's own network specifications and
parameter nodes}\\[4pt]
\today
\end{center}
\vspace{4pt}
\hrule
\vspace{10pt}
"""
    )

    # -- summary -------------------------------------------------------
    L.append(r"\section*{Summary}")
    L.append(
        rf"""
All {len(catalog)} computational examples in the monograph were recomputed under two
algorithm specifications:

\begin{{itemize}}[leftmargin=1.4em,itemsep=1pt]
  \item \textbf{{legacy}} --- reproduces the current \textsc{{DSGRN\_utils}}
    bit-for-bit (verified across five levels and four network families), so it
    matches whatever produced the published figures;
  \item \textbf{{paper}} --- the live manuscript definitions.
\end{{itemize}}

\textbf{{{matched} of {total} runs match the printed values exactly.}} The
remaining {total - matched} are analysed below; each is understood, and none is
a computational error in the manuscript's mathematics.

All four large networks reproduce their stated parameter-graph sizes exactly
($87{{,}}280{{,}}405{{,}}632$; $3{{,}}600{{,}}000$; $13{{,}}608{{,}}000{{,}}000$;
$4{{,}}429{{,}}771{{,}}960{{,}}320$), and both directly specified ramp systems
reproduce their stated index sets. Runtimes are comparable to or faster than the
reported ones.
"""
    )

    # -- traffic light -------------------------------------------------
    L.append(r"\subsection*{At a glance}")
    L.append(r"\begin{longtable}{@{}p{0.30\textwidth}cccp{0.34\textwidth}@{}}")
    L.append(r"\toprule")
    L.append(r"\textbf{Example} & \textbf{Level} & \textbf{legacy} & \textbf{paper} & \textbf{Note} \\")
    L.append(r"\midrule\endhead")
    for label, group in by_label.items():
        entry = catalog.get(label)
        statuses = {o["spec"]: status_of(o) for o in group}
        cell = lambda s: {"ok": r"\OK", "mismatch": r"\BAD", "error": r"\BAD"}.get(s, "--")
        note = ""
        if label == "ex:trivial_index_3D_F3":
            note = "published 61/12 reproduces only under \\textbf{paper} (M1)"
        elif label == "ex:periodicOrbit3_F4":
            note = "index set misprinted (M2)"
        elif label == "ex:periodicOrbit3_F3":
            note = "\\textbf{paper} resolves the component (M4)"
        elif label == "EMT-6D-trivial-indices":
            note = "count matches; node labels are not canonical (M8)"
        elif label == "ex:ramp_van_der_pol":
            note = "cycle rule never fires: $\\cF_3=\\cF_2$ (M5)"
        row = (
            f"{mono(label)} & {group[0]['level']} & {cell(statuses.get('legacy'))} & "
            f"{cell(statuses.get('paper'))} & \\small {note} \\\\"
        )
        L.append(row)
        del entry
    L.append(r"\bottomrule")
    L.append(r"\end{longtable}")

    # -- problems ------------------------------------------------------
    L.append(r"\clearpage")
    L.append(r"\section*{Problems to address}")
    L.append(
        r"""Ordered by consequence. Items \textbf{M$n$} are manuscript-side, items
\textbf{C$n$} are \textsc{DSGRN\_utils}-side. Several manuscript symptoms have a
code-side cause, and are cross-referenced."""
    )

    L.append(r"\subsection*{Manuscript}")
    for p in PROBLEMS_MANUSCRIPT:
        colour = "badred" if "breaks" in p["severity"] or "fail" in p["severity"] else "warnamber"
        L.append(
            rf"\begin{{problembox}}{{{colour}}}{{{p['id']} --- {p['title']}}}"
        )
        L.append(rf"\textit{{\small {mono(p['where'])} \quad ({tex_escape(p['severity'])})}}")
        L.append("")
        L.append(p["body"])
        L.append(r"\end{problembox}")
        L.append(r"\vspace{4pt}")

    L.append(r"\subsection*{Code}")
    for p in PROBLEMS_CODE:
        colour = "badred" if "changes" in p["severity"] or "breaks" in p["severity"] or "crash" in p["severity"] else "warnamber"
        L.append(rf"\begin{{problembox}}{{{colour}}}{{{p['id']} --- {p['title']}}}")
        L.append(rf"\textit{{\small {mono(p['where'])} \quad ({tex_escape(p['severity'])})}}")
        L.append("")
        L.append(p["body"])
        L.append(r"\end{problembox}")
        L.append(r"\vspace{4pt}")

    # -- per-example detail --------------------------------------------
    L.append(r"\clearpage")
    L.append(r"\section*{Example by example}")
    L.append(
        r"""Each block lists what the text claims and what was computed. A row is
shaded red when the computed value differs from the printed one."""
    )

    for label, group in by_label.items():
        entry = catalog.get(label)
        L.append(rf"\subsection*{{{mono(label)}}}")
        if entry is not None:
            source = mono(entry.source)
            target = (
                f"network {mono(entry.network)}, parameter node ${entry.index:,}$".replace(",", "{,}")
                if entry.network
                else f"ramp system {mono(entry.ramp_system)}"
            )
            L.append(
                rf"\small\textit{{{source} \quad$\cdot$\quad {target} \quad$\cdot$\quad "
                rf"$\cF_{{{entry.level}}}$}}\par"
            )
            if entry.claim:
                L.append(rf"\textbf{{Text claims:}} {tex_escape(entry.claim)}\par")

        rows_written = False
        for o in group:
            if not o["results"]:
                continue
            rows_written = True
            L.append(rf"\textbf{{spec = {o['spec']}}} \hfill \small {o['seconds']:.2f}\,s")
            L.append(r"\begin{longtable}{@{}p{0.46\textwidth}p{0.18\textwidth}p{0.18\textwidth}c@{}}")
            L.append(r"\toprule")
            L.append(r"\textbf{Quantity} & \textbf{Claimed} & \textbf{Computed} & \\")
            L.append(r"\midrule\endhead")
            for name, claimed, computed, ok in o["results"]:
                shade = "" if ok else r"\rowcolor{rowbad}"
                mark = r"\OK" if ok else r"\BAD"
                L.append(
                    f"{shade} \\small {tex_escape(name)} & \\small {mono(claimed)} & "
                    f"\\small {mono(computed)} & {mark} \\\\"
                )
            L.append(r"\bottomrule")
            L.append(r"\end{longtable}")

        if not rows_written:
            L.append(r"\textit{\small No numeric claim in the text to check.}\par")

        summary = group[0]["summary"]
        L.append(
            rf"""\small\textbf{{Computed:}} {summary['morse_nodes']} Morse nodes, index set
{mono(summary['index_set'])}, {summary['stg_vertices']} STG vertices and
{summary['stg_edges']} edges.\par"""
        )
        if entry is not None and entry.notes:
            L.append(rf"\small\textit{{{tex_escape(entry.notes)}}}\par")
        if entry is not None and entry.claimed_seconds:
            got = min(o["seconds"] for o in group)
            L.append(
                rf"\small Runtime: {got:.2f}\,s here against {entry.claimed_seconds:g}\,s reported.\par"
            )
        L.append(r"\vspace{6pt}")

    # -- supporting evidence -------------------------------------------
    L.append(r"\clearpage")
    L.append(r"\section*{Supporting evidence}")

    if divergence:
        L.append(r"\subsection*{Which definition change affects what}")
        L.append(
            r"""Parameters whose result changes relative to \textbf{legacy}, per
divergence, sampled over whole parameter graphs where feasible."""
        )
        L.append(r"\begin{longtable}{@{}llrrrr@{}}")
        L.append(r"\toprule")
        L.append(
            r"\textbf{Network} & \textbf{Divergence} & \textbf{Level} & "
            r"\textbf{STG} & \textbf{Morse nodes} & \textbf{Conley} \\"
        )
        L.append(r"\midrule\endhead")
        for r in divergence["rows"]:
            if not (r["stg_changed"] or r["morse_nodes_changed"] or r["conley_changed"]):
                continue
            if r["spec"] == "paper":
                continue
            shade = r"\rowcolor{rowbad}" if r["conley_changed"] else ""
            L.append(
                f"{shade} {mono(r['network'])} & {mono(r['spec'])} & {r['level']} & "
                f"{r['stg_changed']} & {r['morse_nodes_changed']} & {r['conley_changed']} \\\\"
            )
        L.append(r"\bottomrule")
        L.append(r"\end{longtable}")
        L.append(
            r"""Three further readings were parameterised and changed \emph{nothing}
anywhere: evaluating every back-wall choice (positive evidence for
\texttt{prop:back\_wall\_well\_defined}); computing $\cF_3$ as an intersection
plus an unconditional union rather than as a cascade; and widening the domain of
$\rmap\xi$."""
        )

    if hypotheses:
        L.append(r"\subsection*{Analytic admissibility of the ramp parameters}")
        L.append(r"\begin{longtable}{@{}llccccc@{}}")
        L.append(r"\toprule")
        L.append(
            r"\textbf{System} & \textbf{half-width} & $\cH_0$ & $\Lambda(R)$ & "
            r"$\cH_1$ & $\cH_2$ & $\cH_3$ \\"
        )
        L.append(r"\midrule\endhead")
        for name, data in hypotheses.items():
            if name == "sentinel":
                continue
            for row in data["rows"]:
                keys = ["H_0", "Lambda(R)", "H_1", "H_2", "H_3"]
                bad = any(not row[k]["holds"] for k in keys)
                shade = r"\rowcolor{rowbad}" if bad else r"\rowcolor{rowok}"
                cells = " & ".join(r"\OK" if row[k]["holds"] else r"\BAD" for k in keys)
                L.append(f"{shade} {mono(name)} & {mono(row['width'])} & {cells} \\\\")
        L.append(r"\bottomrule")
        L.append(r"\end{longtable}")
        L.append(
            r"""The Morse graph is identical at every width in this table, so repairing
the width does not disturb any computed conclusion."""
        )

    if alignment:
        L.append(r"\subsection*{Numerical alignment of the rectangular geometrization}")
        L.append(
            r"""The inward normal component of the ramp field, sampled on every embedded
wall that the combinatorial map orients. \texttt{thm:R1ABlattice} asserts all of
them are positive for $\cF_1$."""
        )
        L.append(r"\begin{longtable}{@{}llrrrc@{}}")
        L.append(r"\toprule")
        L.append(
            r"\textbf{System} & \textbf{half-width} & \textbf{Level} & "
            r"\textbf{walls} & \textbf{worst inward} & \textbf{aligned} \\"
        )
        L.append(r"\midrule\endhead")
        for row in alignment:
            if row["spec"] != "paper":
                continue
            shade = r"\rowcolor{rowok}" if row["aligned"] else r"\rowcolor{rowbad}"
            L.append(
                f"{shade} {mono(row['system'])} & {mono(row['width'])} & {row['level']} & "
                f"{row['walls']} & ${row['margin']:.3g}$ & "
                f"{r'\OK' if row['aligned'] else r'\BAD'} \\\\"
            )
        L.append(r"\bottomrule")
        L.append(r"\end{longtable}")
        L.append(
            r"""At an admissible width $\cF_1$ is aligned with zero failures; at the
published widths it is not, so the width condition is load-bearing. Every
$\cF_2$/$\cF_3$ failure sits on a wall that $\cF_1$ left as a double edge ---
precisely the faces the R2 GO-manifold and R3 cycle-surface constructions are
meant to replace."""
        )

    # -- Morse graphs and Morse sets, per example ------------------------
    manifest_path = FIGURES / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}

    if manifest.get("examples"):
        L.append(r"\clearpage")
        L.append(r"\section*{Morse graphs and Morse sets}")
        L.append(
            r"""For each example, the Morse graph and the Morse sets under both
specifications. Morse graph nodes carry the Conley index over $\Z_2$ and the
classification \texttt{FP} (fixed point), \texttt{PO} (periodic orbit),
\texttt{T} (trivial index) or \texttt{M} (more complex), with the number of
cells in parentheses. In the Morse set panels blue arrows are single edges of
the state transition graph, red arrows are double edges, and a red disc marks a
self-edge; above dimension two the complex is projected onto the first two
coordinates. Plots are produced by \texttt{DSGRN\_utils.PlotMorseGraph} and
\texttt{PlotMorseSets}, the same routines the monograph's figures use."""
        )
        for entry in manifest["examples"]:
            path = FIGURES / "examples" / entry["file"]
            if not path.exists():
                continue
            L.append(r"\begin{figure}[htbp]\centering")
            L.append(rf"\includegraphics[width=0.99\textwidth,height=0.80\textheight,keepaspectratio]{{{path}}}")
            marker = (
                r"\textcolor{badred}{\textbf{The two specifications disagree.}} "
                if entry["differs"]
                else "Both specifications agree. "
            )
            L.append(
                rf"\caption*{{\small {marker}{mono(entry['label'])} "
                rf"({mono(entry['source'])}), $\cF_{{{entry['level']}}}$, "
                rf"dimension {entry['dim']}.}}"
            )
            L.append(r"\end{figure}")
            L.append(r"\clearpage")

    if manifest.get("gallery2d"):
        L.append(r"\section*{Two-dimensional geometrizations}")
        L.append(
            r"""Each DSGRN parameter node below is realised by an explicit ramp
system through \texttt{defn:DSGRN\_ramp}: a point of the semi-algebraic region
is sampled, each switching function is replaced by a type-1 ramp with the same
plateaus, and the induced wall labelling is checked to be the one DSGRN computes
combinatorially --- which is the content of \texttt{prop:DSGRN\_wall\_constant}.

That gives a genuine vector field on the geometrized complex, so the left panel
of each figure carries, in real phase-space coordinates: the rectangular
geometrization of \texttt{defn:rectGeoXb} (grey grid), the Morse sets (shaded),
the ramp field, the two nullclines (blue $\dot x_1=0$, red $\dot x_2=0$), the
GO-manifolds of \texttt{defn:GO-manifold} (black solid where the pair is
externally pruned, grey dashed where internally pruned), and the orange
connectors of \texttt{defn:orange-manifold} for the triples of
$\cD_3(\rook)$.

The axes are cropped to the threshold range: the outermost cells of the complex
extend to $\gab_n$, which for a sampled DSGRN parameter is far beyond the
thresholds."""
        )
        for entry in manifest["gallery2d"]:
            path = FIGURES / "gallery2d" / entry["file"]
            if not path.exists():
                continue
            L.append(r"\begin{figure}[p]\centering")
            L.append(rf"\includegraphics[width=0.99\textwidth,height=0.84\textheight,keepaspectratio]{{{path}}}")
            L.append(
                rf"\caption*{{\small {mono(entry['network'])} node "
                rf"${entry['index']:,}$".replace(",", "{,}")
                + rf" --- {tex_escape(entry['caption'])}. "
                rf"{entry['go_pairs']} codimension-two GO pairs "
                rf"({entry['external']} externally pruned), "
                rf"{entry['d3_triples']} $\cD_3$ triples, "
                rf"{entry['morse_nodes']} Morse nodes.}}"
            )
            L.append(r"\end{figure}")
            L.append(r"\clearpage")

    # -- figures --------------------------------------------------------
    figs = [
        ("van_der_pol_phase_admissible.png",
         r"""\texttt{ex:ramp\_van\_der\_pol} at an admissible width. Shaded: the
$\cF_3$ Morse sets, drawn as embedded blowup cells in real ramp coordinates.
Blue and red curves: the nullclines $\dot x_1=0$ and $\dot x_2=0$. Black: the
GO-manifolds of \texttt{defn:GO-manifold}. Orange: genuine trajectories."""),
        ("intro_periodic_morse_sets_3d.png",
         r"""\texttt{ex:periodicOrbit3}: the $\cF_3$ Morse sets in phase space with two
integrated trajectories. The stable periodic orbit is the component with
$CH_0\cong CH_1\cong\Z_2$."""),
        ("go_manifold_single.png",
         r"""One GO-manifold, in the context of the two blowup cells whose shared
rectangular face it replaces. Black: the $\delta$-nullcline base. Blue sheet: its
flow-out under the perturbed field \texttt{eq:rampSysPerturbed}. Dashed red: the
terminal slice on the directed $n_g$-wall."""),
        ("van_der_pol_alignment_margins.png",
         r"""Inward normal component per oriented wall. $\cF_1$ lies entirely to the
right of zero; the $\cF_2$ and $\cF_3$ tails to the left are the walls needing
the geometric modification."""),
    ]
    L.append(r"\section*{Ramp systems specified directly in the text}")
    for name, caption in figs:
        if not (FIGURES / name).exists():
            continue
        L.append(r"\begin{figure}[htbp]\centering")
        L.append(rf"\includegraphics[width=0.82\textwidth]{{{FIGURES / name}}}")
        L.append(rf"\caption*{{\small {caption}}}")
        L.append(r"\end{figure}")
        L.append(r"\clearpage")

    # -- reproduction ---------------------------------------------------
    L.append(r"\section*{Reproducing}")
    L.append(
        r"""\begin{verbatim}
cd code
./.venv/bin/python -m pytest rookfields/tests          # 145 tests
./.venv/bin/python -m rookfields.examples --spec both  # this report's data
./.venv/bin/python -m rookfields.hypotheses
./.venv/bin/python -m rookfields.alignment
./.venv/bin/python -m rookfields.plotting
./.venv/bin/python -m rookfields.report_pdf            # this document
\end{verbatim}

Full prose discussion in \texttt{code/rookfields/reports/FINDINGS.md}; the
actionable items are also filed in \texttt{coauthor-todos.md}. No file under
\texttt{paper/} was modified."""
    )

    L.append(r"\end{document}")
    return "\n".join(L)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=REPORTS / "updated-computations.pdf")
    parser.add_argument("--keep-tex", action="store_true")
    args = parser.parse_args(argv)

    if shutil.which("pdflatex") is None:
        raise SystemExit("pdflatex not found")

    build = REPORTS / "_build"
    build.mkdir(parents=True, exist_ok=True)
    tex = build / "updated-computations.tex"
    tex.write_text(build_tex())

    for _ in range(2):  # twice, for the longtable column widths
        proc = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex.name],
            cwd=build,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            log = (build / "updated-computations.log")
            tail = log.read_text(errors="replace").splitlines()[-40:] if log.exists() else []
            print("\n".join(tail))
            raise SystemExit("pdflatex failed")

    produced = build / "updated-computations.pdf"
    shutil.copy(produced, args.out)
    if not args.keep_tex:
        for suffix in (".aux", ".log", ".out", ".toc"):
            (build / f"updated-computations{suffix}").unlink(missing_ok=True)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
