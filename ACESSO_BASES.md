# Guia de Acesso - Bases de Dados Bibliograficas

Instrucoes de acesso e configuracao para cada base de dados utilizada
no pipeline de revisao sistematica.

---

## 1. PubMed/MEDLINE

**Tipo de acesso:** API automatica (NCBI E-utilities)
**Autenticacao:** Email obrigatorio + API key opcional

### Configuracao

1. Preencha seu email em `config.yaml`:
   ```yaml
   api:
     pubmed:
       email: "seu.email@universidade.edu.br"
   ```

2. (Opcional) Obtenha uma API key para aumentar o rate limit de 3 para 10 req/s:
   - Acesse: https://www.ncbi.nlm.nih.gov/account/settings/
   - Crie uma conta NCBI (gratuita)
   - Gere a API key em "API Key Management"
   - Configure em `config.yaml`:
     ```yaml
     api:
       pubmed:
         api_key: "sua_api_key_aqui"
     ```

### Execucao

```bash
python main.py --databases pubmed
```

O pipeline executa a busca e recupera os registros automaticamente.

---

## 2. Scopus

**Tipo de acesso:** Manual (proxy institucional) ou API automatica
**Autenticacao:** Proxy CAPES ou API key Elsevier

### Opcao A - Exportacao manual via proxy CAPES (recomendado)

1. Acesse o Scopus pelo proxy institucional:
   ```
   https://www-scopus-com.ez11.periodicos.capes.gov.br/search/form.uri?display=advanced
   ```

2. Gere a query de busca:
   ```bash
   python main.py --dry-run --databases scopus
   ```
   A query sera salva em `output/scopus_query.txt`.

3. No portal Scopus:
   - Preencha duas linhas de busca (Search within: Article title, Abstract, Keywords):
     - Linha 1: colar LINHA 1 do `scopus_query.txt` (termos de desastres)
     - Operador: AND
     - Linha 2: colar LINHA 2 do `scopus_query.txt` (termos de desfechos)
   - Clique em "+ Add date range" e defina 1990-2026
   - Clique em "Search"
   - Exporte: Export -> CSV (com Abstract) ou RIS
   - Selecione todos os campos (Citation information + Abstract + Keywords)

4. Importe o arquivo exportado:
   ```bash
   python main.py --import-scopus caminho/do/arquivo.csv
   ```

**Formatos suportados:** CSV (.csv) e RIS (.ris)

### Opcao B - API automatica (Elsevier/pybliometrics)

1. Obtenha uma API key gratuita (uso academico):
   - Acesse: https://dev.elsevier.com/apikey/manage
   - Crie uma conta Elsevier Developer
   - Registre uma nova API key

2. Configure em `config.yaml`:
   ```yaml
   api:
     scopus:
       api_key: "sua_api_key_aqui"
       institutional_token: ""  # opcional, para acesso completo
   ```

3. Execute:
   ```bash
   python main.py --databases scopus
   ```

**Nota:** Via API sem assinatura institucional, abstract e autores podem ficar vazios. O pipeline tenta enriquecer via AbstractRetrieval (mais lento). A Opcao A (manual) garante dados completos.

---

## 3. Web of Science

**Tipo de acesso:** Manual (proxy institucional) ou API automatica
**Autenticacao:** Proxy CAPES ou API key Clarivate

### Opcao A - Exportacao manual via proxy CAPES (recomendado)

1. Acesse o WoS pelo proxy institucional:
   ```
   https://www-webofscience-com.ez11.periodicos.capes.gov.br/wos/woscc/smart-search
   ```

2. Gere a query de busca:
   ```bash
   python main.py --dry-run --databases wos
   ```
   A query sera salva em `output/wos_query.txt`.

3. No portal WoS:
   - Abra "Advanced Search"
   - Cole a query do arquivo `wos_query.txt`
   - Execute a busca
   - Clique em "Export" -> selecione formato:
     - **RIS** (recomendado, compativel com Zotero/Mendeley)
     - Ou **Tab-delimited** (Full Record)
   - Selecione "Full Record" para incluir abstract
   - Salve o arquivo

4. Importe o arquivo exportado:
   ```bash
   python main.py --import-wos caminho/do/arquivo.ris
   ```

### Opcao B - API automatica (Clarivate Starter API)

1. Obtenha uma API key:
   - Acesse: https://developer.clarivate.com/apis/wos-starter
   - Registre-se para a Starter API (gratuita, limite 50K registros/ano)

2. Configure em `config.yaml`:
   ```yaml
   api:
     wos:
       api_key: "sua_api_key_aqui"
   ```

3. Execute:
   ```bash
   python main.py --databases wos
   ```

**Nota:** A API Starter nao retorna abstract. Para abstract completo,
use a Opcao A (exportacao manual).

---

## 4. BVS/LILACS

**Tipo de acesso:** Manual (exportacao pelo portal)
**Autenticacao:** Nenhuma

O portal BVS utiliza protecao CDN (Bunny Shield) que bloqueia acesso programatico. A busca deve ser feita manualmente no navegador.

### Passo a passo

1. Gere a query de busca:
   ```bash
   python main.py --dry-run --databases bvs
   ```
   A query sera salva em `output/bvs_query.txt`.

2. Acesse o portal BVS:
   ```
   https://pesquisa.bvsalud.org/portal/
   ```

3. No portal:
   - Cole a query na busca avancada
   - Aplique filtro: Base de dados = LILACS
   - Verifique o total de resultados
   - Clique em "Exportar" (icone de download)
   - Formato: "RIS (para Reference Manager, Zotero, etc.)"
   - Salve o arquivo .ris

4. Importe o arquivo exportado:
   ```bash
   python main.py --import-bvs caminho/do/arquivo.ris
   ```

**Formatos suportados:** RIS (.ris) e CSV (.csv)

---

## Pipeline completo

### Execucao com todas as bases automaticas

```bash
python main.py --databases pubmed,scopus
```

### Importando bases manuais junto com automaticas

```bash
python main.py --databases pubmed --import-scopus output/scopus_20260222_213251/scopus_export.ris --import-bvs output/bvs_20260222_214312/bvs_export.ris --import-wos output/wos_20260222_213650/wos_export.ris

python main.py --databases pubmed --import-scopus scopus_export.csv --import-bvs bvs_export.ris --import-wos wos_export.ris
```

### Importando todas as bases de RIS (sem busca online)

```bash
python main.py --databases "" --import-pubmed output/pubmed_20260222_212228/pubmed_20260222_212228.ris --import-scopus output/scopus_20260222_213251/scopus_export.ris --import-bvs output/bvs_20260222_214312/bvs_export.ris --import-wos output/wos_20260222_213650/wos_export.ris
```

### Apenas gerar queries (sem executar buscas)

```bash
python main.py --dry-run
```

Isso gera os arquivos `scopus_query.txt`, `wos_query.txt` e `bvs_query.txt` no diretorio de saida para uso manual nas respectivas interfaces web.

### Exportar sem deduplicacao

```bash
python main.py --databases pubmed --import-scopus scopus.csv --import-wos wos.ris --skip-dedup
```

---

## Arquivos de saida

O pipeline gera os seguintes arquivos no diretorio `output/<timestamp>/`:

| Arquivo | Descricao |
|---|---|
| `pubmed_<ts>.csv/.xlsx/.ris` | Registros individuais do PubMed |
| `scopus_<ts>.csv/.xlsx/.ris` | Registros individuais do Scopus |
| `wos_<ts>.csv/.xlsx/.ris` | Registros individuais do WoS |
| `bvs_<ts>.csv/.xlsx/.ris` | Registros individuais da BVS |
| `results_<ts>.csv/.xlsx/.ris` | Registros combinados (pos-deduplicacao) |
| `duplicates_<ts>.csv` | Registros duplicados removidos |
| `prisma_log_<ts>.json` | Log PRISMA da revisao sistematica |
| `scopus_query.txt` | Query para busca manual no Scopus |
| `wos_query.txt` | Query para busca manual no WoS |
| `bvs_query.txt` | Query para busca manual na BVS |

---

## Resumo de requisitos por base

| Base | Acesso | Requisito | Custo |
|---|---|---|---|
| PubMed | API automatica | Email (+ API key opcional) | Gratuito |
| Scopus | Manual (proxy) | Acesso CAPES/institucional | Gratuito via CAPES |
| Scopus | API automatica | API key Elsevier | Gratuito (academico) |
| Web of Science | Manual (proxy) | Acesso CAPES/institucional | Gratuito via CAPES |
| Web of Science | API automatica | API key Clarivate | Gratuito (Starter) |
| BVS/LILACS | Manual (portal) | Nenhum | Gratuito |
