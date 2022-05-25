#!/usr/bin/env python3
from __future__ import annotations

from typing import Optional, Any, Dict
from pathlib import Path
import json
from dataclasses import dataclass

import appdirs

_CONFIG_PATH = Path(appdirs.user_config_dir(appname="TaxI", appauthor="iTaxoTools"))
if not _CONFIG_PATH.exists():
    _CONFIG_PATH.mkdir(parents=True)


def get_config(name: str) -> Optional[Dict[str, Any]]:
    config_path = _CONFIG_PATH / name
    if not config_path.exists():
        return None
    try:
        with open(config_path) as config_file:
            return json.load(config_file)
    except json.JSONDecodeError:
        return None
    except FileNotFoundError:
        return None


DEFAULT_SCORES_DICT = {
    "gap penalty": -8,
    "gap extend penalty": -1,
    "end gap penalty": -1,
    "end gap extend penalty": -1,
    "match score": 1,
    "mismatch score": -1,
}


@dataclass
class AlignmentScores:
    gap_penalty: int
    gap_extend_penalty: int
    end_gap_penalty: int
    end_gap_extend_penalty: int
    match_score: int
    mismatch_score: int

    @classmethod
    def _from_scores_dict(cls, scores_dict: Dict[str, int]) -> AlignmentScores:
        return cls(
            gap_penalty=scores_dict["gap penalty"],
            gap_extend_penalty=scores_dict["gap extend penalty"],
            end_gap_penalty=scores_dict["end gap penalty"],
            end_gap_extend_penalty=scores_dict["end gap extend penalty"],
            match_score=scores_dict["match score"],
            mismatch_score=scores_dict["mismatch score"],
        )

    def _to_scores_dict(self) -> Dict[str, int]:
        return {
            "gap penalty": self.gap_penalty,
            "gap extend penalty": self.gap_extend_penalty,
            "end gap penalty": self.end_gap_penalty,
            "end gap extend penalty": self.end_gap_extend_penalty,
            "match score": self.match_score,
            "mismatch score": self.mismatch_score,
        }


@dataclass
class Config:
    alignment_scores: AlignmentScores

    @classmethod
    def load(cls) -> Config:
        scores = get_scores()
        return cls(alignment_scores=AlignmentScores._from_scores_dict(scores))

    def save(self) -> None:
        with open(_CONFIG_PATH / "scores.json", mode="w") as scores_file:
            json.dump(self.alignment_scores._to_scores_dict(), scores_file)


def get_scores() -> Dict[str, int]:
    scores = get_config("scores.json")
    if not scores:
        with open(_CONFIG_PATH / "scores.json", mode="w") as scores_file:
            json.dump(DEFAULT_SCORES_DICT, scores_file)
        return DEFAULT_SCORES_DICT
    for score_name, value in scores.items():
        if not isinstance(value, int):
            raise ValueError(
                f"The value for '{score_name}' in data/scores.tab is not a number"
            )
    return scores
