"""Build a static, cautious report for the prespecified response analysis."""

from html import escape
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def write_analysis_report(output, scores, included, results):
    output = Path(output)
    combination = included[included.arm == "durvalumab_sbrt"].copy()
    stats = {row["signature"]: row for row in results}
    names = list(stats)
    colors = {0: "#63758a", 1: "#187d78"}

    fig, axes = plt.subplots(1, len(names), figsize=(6.1 * len(names), 5.1), squeeze=False)
    for axis, name in zip(axes[0], names):
        values = scores.loc[combination.sample_id, name].to_numpy(dtype=float)
        labels = combination.mpr.to_numpy(dtype=int)
        for group in (0, 1):
            group_values = values[labels == group]
            jitter = np.linspace(-0.09, 0.09, len(group_values)) if len(group_values) > 1 else np.zeros(len(group_values))
            axis.scatter(group + 1 + jitter, group_values, s=58, color=colors[group],
                         edgecolor="white", linewidth=0.8, zorder=3, label=f"{'MPR' if group else 'No MPR'} (n={len(group_values)})")
            if len(group_values):
                median = float(np.median(group_values))
                axis.plot([group + 0.78, group + 1.22], [median, median], color="#14283d", linewidth=2.3, zorder=4)
        result = stats[name]
        axis.set_title(name, loc="left", fontsize=13, fontweight="bold", color="#173047")
        axis.set_xticks([1, 2], ["No MPR", "MPR"])
        axis.set_ylabel("Locked signature score")
        axis.grid(axis="y", color="#e5ebf0", linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(0.02, 0.98,
                  f"Rank probability: {result['rank_probability']:.3f}\nExact p: {result['p_exact_two_sided']:.4g} · Holm p: {result['p_holm']:.4g}",
                  transform=axis.transAxes, va="top", ha="left", fontsize=9,
                  bbox={"boxstyle": "round,pad=.4", "facecolor": "white", "edgecolor": "#dce5ec"})
        axis.legend(frameon=False, loc="best", fontsize=9)
    fig.suptitle("Pretreatment signature scores by MPR in the durvalumab + SBRT arm", x=0.06, ha="left",
                 fontsize=15, color="#173047", fontweight="bold")
    fig.text(0.06, 0.015, "Points are individual samples; bars mark group medians. Association only; this is not a prediction or treatment-benefit estimate.",
             fontsize=9, color="#607186")
    fig.tight_layout(rect=(0, 0.05, 1, 0.91))
    fig.savefig(output / "score_by_mpr.png", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    rows = []
    for result in results:
        rows.append(
            "<tr>"
            f"<th scope=\"row\">{escape(result['signature'])} <span>{escape(result['role'])}</span></th>"
            f"<td>{result['n_mpr']}</td><td>{result['n_non_mpr']}</td>"
            f"<td>{result['median_mpr']:.4g}</td><td>{result['median_non_mpr']:.4g}</td>"
            f"<td>{result['rank_probability']:.3f} ({result['bootstrap_ci_low']:.3f}–{result['bootstrap_ci_high']:.3f})</td>"
            f"<td>{result['p_exact_two_sided']:.4g}</td><td>{result['p_holm']:.4g}</td></tr>"
        )
    n_mpr = int((combination.mpr == 1).sum())
    n_no_mpr = int((combination.mpr == 0).sum())
    html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="Exploratory fixed-signature association analysis in a small NSCLC radiotherapy cohort.">
<title>Response association analysis</title>
<style>
:root{{--ink:#173047;--muted:#607186;--paper:#f3f6f8;--line:#dce5ec;--teal:#187d78;--white:#fff}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}}
header{{padding:42px max(22px,calc((100vw - 1060px)/2));background:linear-gradient(120deg,#10263d,#1d5364);color:#fff}}
header p{{color:#d8e8ef;margin:8px 0 0}}main{{max-width:1060px;margin:26px auto;padding:0 20px}}section{{background:#fff;border:1px solid var(--line);border-radius:14px;padding:22px;margin-bottom:16px}}
h1{{font-size:clamp(26px,4vw,39px);line-height:1.1;margin:7px 0}}h2{{font-size:20px;margin:0 0 12px}}.eyebrow{{text-transform:uppercase;letter-spacing:.12em;font-size:11px;font-weight:750;color:#8ed0c9;margin:0}}
.cards{{display:flex;gap:12px;flex-wrap:wrap;margin:16px 0}}.card{{background:#f3f6f8;border:1px solid var(--line);border-radius:10px;padding:13px 16px;min-width:170px}}.card strong{{display:block;font-size:22px}}.muted{{color:var(--muted)}}
img{{max-width:100%;height:auto;border:1px solid var(--line);border-radius:10px}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%;min-width:760px}}th,td{{text-align:left;border-bottom:1px solid var(--line);padding:11px 10px;vertical-align:top}}th span{{display:inline-block;color:var(--muted);font-weight:500;font-size:12px}}
.callout{{background:#fff7e8;border-left:4px solid #d78b45;padding:13px 16px;border-radius:5px}}a{{color:var(--teal)}}footer{{color:var(--muted);padding:8px 0 30px}}nav a{{color:white;margin-right:16px}}
</style></head><body><header><p class="eyebrow">Prespecified fixed-signature analysis · exploratory</p><h1>Signature scores and major pathologic response</h1>
<p>Pretreatment GSE253564 profiles in the durvalumab + SBRT arm, compared by verified MPR labels.</p></header>
<main><section><h2>Samples analyzed</h2><div class="cards"><div class="card"><span class="muted">MPR</span><strong>{n_mpr}</strong></div><div class="card"><span class="muted">No MPR</span><strong>{n_no_mpr}</strong></div></div>
<p class="muted">Outcome labels are read from the supplied verified metadata file. Scores are computed from the locked gene lists and weights before labels are used.</p></section>
<section><h2>Score distributions</h2><img src="score_by_mpr.png" alt="Pretreatment signature scores for MPR and non-MPR groups"><p class="muted">Each point is one sample; horizontal marks show group medians. Signature panels use separate vertical scales.</p></section>
<section><h2>Association estimates</h2><div class="table-wrap"><table><thead><tr><th>Signature</th><th>MPR n</th><th>No MPR n</th><th>Median MPR</th><th>Median no MPR</th><th>Rank probability (95% bootstrap interval)</th><th>Exact p</th><th>Holm p</th></tr></thead><tbody>{''.join(rows)}</tbody></table></div>
<p>Rank probability is the chance that a randomly selected MPR sample has a higher score than a randomly selected no-MPR sample, plus half the tie probability. The two-sided p-value comes from an exact label permutation test. Holm adjustment covers the two prespecified signatures.</p></section>
<section class="callout"><strong>Interpretation limit.</strong> This is a small, retrospective cohort association analysis. The intervals are descriptive and unstable at this sample size. The analysis does not fit a classifier, estimate clinical prediction performance, establish a causal treatment benefit, or validate a biomarker.</section>
<section><h2>Reproducibility files</h2><ul><li><a href="statistics.csv">Statistical results</a></li><li><a href="scores.csv">Per-sample signature scores</a></li><li><a href="sample_flow.csv">Inclusion and exclusion audit</a></li><li><a href="sensitivity.csv">Prespecified score sensitivity checks</a></li><li><a href="gene_qc.csv">Gene-level expression QC</a></li><li><a href="manifest.json">Input checksums and run settings</a></li></ul><p><a href="../exploration/index.html">Back to the data explorer</a></p></section>
<footer>Exploratory research output; not for clinical decision-making.</footer></main></body></html>'''
    (output / "analysis_report.html").write_text(html)
