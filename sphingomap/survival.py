"""Kaplan-Meier / log-rank and Cox models for score-based patient strata."""
from __future__ import annotations

import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test


def prepare(score: pd.Series, clinical: pd.DataFrame, time_col: str = "OS.time",
            event_col: str = "OS") -> pd.DataFrame:
    df = pd.concat([score.rename("score"), clinical[[time_col, event_col]]], axis=1, join="inner")
    df = df.dropna()
    df = df[df[time_col] > 0]
    df["high"] = df["score"] >= df["score"].median()
    return df.rename(columns={time_col: "time", event_col: "event"})


def median_split_logrank(df: pd.DataFrame) -> dict:
    hi, lo = df[df.high], df[~df.high]
    res = logrank_test(hi.time, lo.time, event_observed_A=hi.event, event_observed_B=lo.event)
    cph = CoxPHFitter().fit(df[["time", "event", "score"]], "time", "event")
    s = cph.summary.loc["score"]
    return {
        "n": len(df), "events": int(df.event.sum()),
        "logrank_p": float(res.p_value),
        "cox_hr_per_unit": float(s["exp(coef)"]),
        "cox_hr_ci_low": float(s["exp(coef) lower 95%"]),
        "cox_hr_ci_high": float(s["exp(coef) upper 95%"]),
        "cox_p": float(s["p"]),
    }


def km_curves(df: pd.DataFrame) -> dict[str, KaplanMeierFitter]:
    fits = {}
    for label, sub in (("High", df[df.high]), ("Low", df[~df.high])):
        fits[label] = KaplanMeierFitter(label=f"{label} (n={len(sub)})").fit(sub.time, sub.event)
    return fits
