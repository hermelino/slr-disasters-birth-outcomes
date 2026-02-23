# Estratégia de busca — Revisão sistemática

## Tema

Efeitos de desastres ambientais em indicadores gestacionais, neonatais e infantis.

A literatura está dispersa entre epidemiologia, saúde pública, economia da saúde e ciências ambientais.

## Bases de dados

### Primárias (utilizadas no pipeline)

| Base | Justificativa |
|---|---|
| PubMed/MEDLINE | Indispensável para estudos epidemiológicos e clínicos sobre mortalidade materna, infantil e prematuridade associadas a desastres |
| Web of Science | Boa cobertura multidisciplinar, útil para cruzar com literatura econômica e ambiental |
| Scopus | Ampla cobertura em ciências da saúde e ambientais, com boas ferramentas de análise bibliométrica |
| LILACS (BVS) | Fundamental para literatura brasileira e latino-americana, incluindo o contexto do SUS |

### Complementares (busca manual) — *não implementado*

| Base | Justificativa |
|---|---|
| EMBASE | Forte em farmacologia e saúde pública, complementa o PubMed |
| SciELO | Acesso aberto com boa cobertura de periódicos brasileiros e regionais |
| REPIDISCA (PAHO/OPAS) | Repositório especializado em saúde ambiental nas Américas |
| Global Health (CAB Abstracts) | Foco em saúde em países de baixa e média renda |
| Disaster Lit (NLM) | Repositório especializado em literatura sobre desastres e saúde pública |
| ProQuest Dissertations | Teses e dissertações sobre o tema |

### Periódicos-chave para monitorar — *não implementado*

- **Saúde-desastres:** Disaster Medicine and Public Health Preparedness, Prehospital and Disaster Medicine, Environmental Health Perspectives
- **Saúde materna e infantil:** The Lancet, Paediatric and Perinatal Epidemiology, BMC Pregnancy and Childbirth, Maternal and Child Health Journal
- **Epidemiologia ambiental:** Environmental Research, American Journal of Epidemiology, International Journal of Environmental Research and Public Health

### Literatura cinzenta — *não implementado*

Repositórios da OMS, UNICEF, UNFPA, Ministério da Saúde (DATASUS) e Defesa Civil/CEMADEN. Muito do que existe sobre o Brasil está em relatórios técnicos, não em periódicos indexados.

---

## Lógica de busca

```
(Bloco 1: Desastres) AND (Bloco 2: Desfechos gestacionais e neonatais)
```

Filtro temporal: 2000–2025

## Bloco 1 — Desastres ambientais

disaster, natural disaster, flood, earthquake, hurricane, typhoon, tornado, wildfire, tsunami, drought, cyclone, chemical spill, oil spill, volcanic eruption, landslide, mudslide, nuclear accident, radioactive hazard release, dam failure, dam collapse, mining disaster, environmental disaster, climate disaster, mass casualty

## Bloco 2 — Desfechos gestacionais e neonatais

**Gestacionais:** pregnancy outcome, pregnancy complication, premature birth, preterm birth, preterm delivery, low birth weight (LBW), small for gestational age (SGA), gestational age, spontaneous abortion, stillbirth, preeclampsia, pre-eclampsia, eclampsia, fetal death, birth weight, maternal mortality, miscarriage, gestational diabetes, gestational hypertension, intrauterine growth restriction (IUGR), premature rupture, birth outcome, adverse birth, pregnancy loss

**Neonatais:** infant mortality, neonatal mortality, perinatal mortality, neonatal death, infant death, Apgar score, congenital abnormalities, congenital anomalies, congenital malformations, birth defect, NICU, neonatal intensive care, neonatal outcome, fetal growth retardation, infant newborn diseases

---

## Adaptação por base

| Base | Campo de busca | Vocabulário controlado |
|---|---|---|
| PubMed | MeSH Terms + Title/Abstract [tiab] | MeSH (com explosion hierárquico) |
| Scopus | TITLE-ABS-KEY (título, abstract, keywords) | Não |
| Web of Science | TS — Topic (título, abstract, keywords) | Não |
| BVS/LILACS | DeCS/MeSH (pt/en/es) | DeCS |

---

## Resultados da busca (PRISMA)

Data da execução: 2026-02-22/23

### Identificação

| Base | Registros identificados |
|---|---|
| PubMed | 3.116 |
| Scopus | 1.439 |
| Web of Science | 518 |
| BVS/LILACS | 21 |
| **Total** | **5.094** |

### Deduplicação

| Etapa | Registros |
|---|---|
| Total identificados | 5.094 |
| Duplicatas por DOI | 679 |
| Duplicatas por similaridade (fuzzy) | 69 |
| **Total duplicatas removidas** | **748** |
| Registros após deduplicação | 4.346 |
| Removidos manualmente (dados fora do padrão) | 1 |
| **Registros para triagem** | **4.345** |

### Fluxo PRISMA

```
Identificação (n = 5.094)
    PubMed:  3.116
    Scopus:  1.439
    WoS:       518
    BVS:        21
        |
        v
Duplicatas removidas (n = 748)
    DOI exato:        679
    Título similar:    69
        |
        v
Registros após deduplicação (n = 4.346)
        |
        v
Removidos manualmente (n = 1)
    Dados fora do padrão (abstract com texto completo do artigo)
        |
        v
Registros para triagem (n = 4.345)
```

---

## Colunas do arquivo CSV de resultados

| Coluna | Descrição |
|---|---|
| source_db | Base de origem (pubmed, scopus, wos, bvs) |
| source_id | Identificador na base de origem (PMID, EID, UT, LILACS ID) |
| doi | Digital Object Identifier |
| title | Título do artigo |
| authors | Autores (separados por ponto e vírgula) |
| journal | Nome do periódico |
| year | Ano de publicação |
| volume | Volume do periódico |
| issue | Número/fascículo |
| pages | Páginas |
| issn | ISSN do periódico |
| abstract | Resumo do artigo |
| keywords | Palavras-chave dos autores (separadas por ponto e vírgula) |
| mesh_terms | Descritores MeSH indexados (separados por ponto e vírgula) |
| matched_disasters_mesh | Termos MeSH do artigo que correspondem ao Bloco 1 (desastres) |
| matched_gestational_mesh | Termos MeSH do artigo que correspondem aos desfechos gestacionais |
| matched_neonatal_mesh | Termos MeSH do artigo que correspondem aos desfechos neonatais |
| language | Idioma da publicação |
| publication_type | Tipo de publicação (journal article, review, etc.) |
| url | Link para o artigo |

---

## Estrutura do projeto

```
des_ambientais_e_ind_gestacionais/
|
|-- main.py                          # Ponto de entrada do pipeline
|-- config.yaml                      # Configuração (API keys, filtros, output)
|-- config.example.yaml              # Exemplo de configuração (sem credenciais)
|-- requirements.txt                 # Dependências Python
|
|-- README.md                         # Este arquivo (estratégia de busca)
|-- ACESSO_BASES.md                  # Instruções de acesso a cada base de dados
|
|-- src/
|   |-- models.py                    # BibRecord — modelo normalizado de registro bibliográfico
|   |-- config.py                    # Carregamento e validação do config.yaml
|   |-- logging_prisma.py            # Logging, estatísticas PRISMA e resumo final
|   |
|   |-- extractors/
|   |   |-- base.py                  # BaseExtractor — classe abstrata (build_query, search, fetch_records)
|   |   |-- pubmed_extractor.py      # PubMed via E-utilities (API). MeSH matching e explosion
|   |   |-- scopus_extractor.py      # Scopus via pybliometrics (API) ou import manual (CSV/RIS)
|   |   |-- wos_extractor.py         # Web of Science via Starter API ou import manual (RIS/Tab-delimited)
|   |   |-- bvs_extractor.py         # BVS/LILACS — somente import manual (RIS)
|   |
|   |-- exporters/
|   |   |-- csv_exporter.py          # Exportação CSV (UTF-8-BOM, separador ponto e vírgula)
|   |   |-- ris_exporter.py          # Exportação RIS (formato de intercâmbio bibliográfico)
|   |
|   |-- dedup/
|       |-- deduplicator.py          # Deduplicação em 2 fases: DOI exato + similaridade fuzzy de título
|
|-- output/                          # Resultados gerados pelo pipeline
|   |-- <base>_<timestamp>/          # Exportação individual por execução de base
|   |-- <timestamp>/                 # Execução combinada (todas as bases + dedup)
|       |-- <base>_<ts>.csv/.ris     # Registros individuais por base
|       |-- results_<ts>.csv/.ris    # Registros combinados pós-deduplicação
|       |-- duplicates_<ts>.csv      # Registros duplicados removidos
|       |-- search_log_<ts>.log      # Log completo da execução
|
|-- tests/                           # Testes
```
