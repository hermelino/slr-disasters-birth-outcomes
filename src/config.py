"""Carregador de configuração YAML com resolução de variáveis de ambiente."""

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml


@dataclass
class PubMedConfig:
    email: str = ""
    api_key: str = ""
    tool_name: str = "disaster_gestational_review"


@dataclass
class ScopusConfig:
    api_key: str = ""
    institutional_token: str = ""


@dataclass
class BvsConfig:
    base_url: str = "https://pesquisa.bvsalud.org/portal/"


@dataclass
class ApiConfig:
    pubmed: PubMedConfig = field(default_factory=PubMedConfig)
    scopus: ScopusConfig = field(default_factory=ScopusConfig)
    bvs: BvsConfig = field(default_factory=BvsConfig)


@dataclass
class DateRange:
    start: str = "1990/01/01"
    end: str = "2026/02/21"


@dataclass
class SearchConfig:
    date_range: DateRange = field(default_factory=DateRange)
    languages: List[str] = field(default_factory=list)
    max_results_per_db: int = 10000


@dataclass
class DedupConfig:
    fuzzy_threshold: int = 90
    title_normalization: bool = True


@dataclass
class OutputConfig:
    directory: str = "output"
    formats: List[str] = field(default_factory=lambda: ["csv", "ris"])
    timestamp: bool = True


@dataclass
class Config:
    api: ApiConfig = field(default_factory=ApiConfig)
    search: SearchConfig = field(default_factory=SearchConfig)
    dedup: DedupConfig = field(default_factory=DedupConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


ENV_VAR_PATTERN = re.compile(r"\$\{(\w+)\}")


def _resolve_env_vars(value: str) -> str:
    """Substitui referências ${VAR_NAME} pelo valor da variável de ambiente."""
    def replacer(match):
        var_name = match.group(1)
        return os.environ.get(var_name, "")
    return ENV_VAR_PATTERN.sub(replacer, value)


def _resolve_dict(d: dict) -> dict:
    """Resolve variáveis de ambiente recursivamente em um dicionário."""
    resolved = {}
    for k, v in d.items():
        if isinstance(v, str):
            resolved[k] = _resolve_env_vars(v)
        elif isinstance(v, dict):
            resolved[k] = _resolve_dict(v)
        elif isinstance(v, list):
            resolved[k] = [
                _resolve_env_vars(item) if isinstance(item, str) else item
                for item in v
            ]
        else:
            resolved[k] = v
    return resolved


def load_config(config_path: str = "config.yaml") -> Config:
    """Carrega e valida a configuração a partir de um arquivo YAML."""
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Arquivo de configuração não encontrado: {config_path}")

    with open(path, "r", encoding="utf-8") as f:
        raw = yaml.safe_load(f) or {}

    raw = _resolve_dict(raw)

    config = Config()

    # API
    api_raw = raw.get("api", {})
    pm = api_raw.get("pubmed", {})
    config.api.pubmed = PubMedConfig(
        email=pm.get("email", ""),
        api_key=pm.get("api_key", ""),
        tool_name=pm.get("tool_name", "disaster_gestational_review"),
    )
    sc = api_raw.get("scopus", {})
    config.api.scopus = ScopusConfig(
        api_key=sc.get("api_key", ""),
        institutional_token=sc.get("institutional_token", ""),
    )
    bvs = api_raw.get("bvs", {})
    config.api.bvs = BvsConfig(
        base_url=bvs.get("base_url", "https://pesquisa.bvsalud.org/portal/"),
    )

    # Search
    search_raw = raw.get("search", {})
    dr = search_raw.get("date_range", {})
    config.search = SearchConfig(
        date_range=DateRange(
            start=dr.get("start", "1990/01/01"),
            end=dr.get("end", "2026/02/21"),
        ),
        languages=search_raw.get("languages", []),
        max_results_per_db=search_raw.get("max_results_per_db", 10000),
    )

    # Dedup
    dedup_raw = raw.get("dedup", {})
    config.dedup = DedupConfig(
        fuzzy_threshold=dedup_raw.get("fuzzy_threshold", 90),
        title_normalization=dedup_raw.get("title_normalization", True),
    )

    # Output
    out_raw = raw.get("output", {})
    config.output = OutputConfig(
        directory=out_raw.get("directory", "output"),
        formats=out_raw.get("formats", ["csv", "ris"]),
        timestamp=out_raw.get("timestamp", True),
    )

    # Validações
    if not config.api.pubmed.email:
        raise ValueError(
            "PubMed requer um email (api.pubmed.email). "
            "Defina no config.yaml ou via variável de ambiente."
        )

    return config
