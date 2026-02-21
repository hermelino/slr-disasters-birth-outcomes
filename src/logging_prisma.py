"""Logging estruturado e contadores PRISMA para rastreabilidade da revisão."""

import json
import logging
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List


@dataclass
class DatabaseStats:
    """Estatísticas de uma base de dados individual."""
    database_name: str = ""
    search_query: str = ""
    search_date: str = ""
    records_identified: int = 0
    records_retrieved: int = 0


@dataclass
class PrismaLog:
    """Registro completo para diagrama PRISMA 2020."""
    execution_start: str = ""
    execution_end: str = ""
    databases: List[DatabaseStats] = field(default_factory=list)
    total_identified: int = 0
    total_retrieved: int = 0
    duplicates_removed_doi: int = 0
    duplicates_removed_fuzzy: int = 0
    total_duplicates_removed: int = 0
    records_after_dedup: int = 0
    output_csv: str = ""
    output_ris: str = ""


def setup_logging(log_file: str, verbose: bool = False) -> None:
    """Configura logging com handler de arquivo e console."""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)

    # Limpar handlers existentes
    root_logger.handlers.clear()

    # Handler de arquivo (detalhado)
    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler.setFormatter(file_fmt)
    root_logger.addHandler(file_handler)

    # Handler de console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    console_handler.setFormatter(console_fmt)
    root_logger.addHandler(console_handler)


def save_prisma_log(prisma: PrismaLog, output_path: str) -> str:
    """Salva o log PRISMA como JSON."""
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(asdict(prisma), f, ensure_ascii=False, indent=2)
    logging.getLogger(__name__).info("PRISMA log salvo: %s", output_path)
    return output_path


def print_summary(prisma: PrismaLog) -> None:
    """Imprime resumo formatado no console."""
    print("\n" + "=" * 60)
    print("  RESUMO DA EXTRAÇÃO BIBLIOGRÁFICA")
    print("=" * 60)

    for db in prisma.databases:
        print(f"\n  {db.database_name}:")
        print(f"    Identificados: {db.records_identified}")
        print(f"    Recuperados:   {db.records_retrieved}")

    print(f"\n  {'─' * 40}")
    print(f"  Total identificados:     {prisma.total_identified}")
    print(f"  Total recuperados:       {prisma.total_retrieved}")
    print(f"  Duplicatas (DOI):        {prisma.duplicates_removed_doi}")
    print(f"  Duplicatas (fuzzy):      {prisma.duplicates_removed_fuzzy}")
    print(f"  Total duplicatas:        {prisma.total_duplicates_removed}")
    print(f"  Registros únicos:        {prisma.records_after_dedup}")

    print(f"\n  Arquivos gerados:")
    if prisma.output_csv:
        print(f"    CSV: {prisma.output_csv}")
    if prisma.output_ris:
        print(f"    RIS: {prisma.output_ris}")

    print("=" * 60 + "\n")
