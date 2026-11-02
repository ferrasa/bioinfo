# Aula prática UEPG — RNA-Seq em *Arabidopsis thaliana*

![Ubuntu](https://img.shields.io/badge/Ubuntu-24.04%2B-E95420?logo=ubuntu&logoColor=white)
![RNA--Seq](https://img.shields.io/badge/RNA--Seq-Pipeline-087E8B)
![Arabidopsis](https://img.shields.io/badge/Arabidopsis-TAIR10-238636)
![STAR](https://img.shields.io/badge/STAR-Alinhamento-2563EB)
![featureCounts](https://img.shields.io/badge/featureCounts-Quantificação-9A6700)
![DESeq2](https://img.shields.io/badge/DESeq2-Expressão_diferencial-8250DF)


**BioProject NCBI:** [PRJNA904764](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904764)  

**Tema:** resposta transcricional ao estresse hídrico em tecido aéreo (*shoot*) de *Arabidopsis*.  

**Pergunta científica:** quais genes diferem em expressão entre plantas irrigadas (*Watered_Day1*) e plantas submetidas ao tratamento de seca (*Drought_Day1*) no mesmo ponto de coleta (*Day1*)?  

**Pipeline:** SRA Toolkit → FastQC/MultiQC → STAR → featureCounts → DESeq2.  

**Sistema operacional:** Linux Ubuntu.  

**Desenho experimental:** seis runs RNA-Seq single-end: três controles e três tratados.  

**Dados verificados:** `META/SraRunTable.csv`, exportado do NCBI SRA Run Selector e incluído neste material; `META/samples.tsv` já está preenchido para os scripts.

## Amostras selecionadas (seleção fechada)

Estas seis bibliotecas foram identificadas diretamente no arquivo `SraRunTable.csv` do NCBI. Os nomes dos grupos seguem o campo `treatment` do arquivo. Os valores em MB são **decimais**, derivados de `Bytes / 1.000.000`, relativos aos dados SRA informados pelo NCBI; não representam o espaço necessário após conversão para FASTQ ou alinhamento.

| Grupo | Run (SRA) | Sample Name / Library Name | Réplica | SRA (MB) | AvgSpotLen (nt) |
|---|---|---|---:|---:|---:|
| Controle — `Watered_Day1` | `SRR22397138` | `shoot_Watered_Day1_rep1` | 1 | 88,6 | 143 |
| Controle — `Watered_Day1` | `SRR22397137` | `shoot_Watered_Day1_rep2` | 2 | 72,1 | 137 |
| Controle — `Watered_Day1` | `SRR22397135` | `shoot_Watered_Day1_rep3` | 3 | 82,4 | 144 |
| Seca — `Drought_Day1` | `SRR22397141` | `shoot_Drought_Day1_rep1` | 1 | 94,2 | 146 |
| Seca — `Drought_Day1` | `SRR22397140` | `shoot_Drought_Day1_rep2` | 2 | 99,3 | 138 |
| Seca — `Drought_Day1` | `SRR22397139` | `shoot_Drought_Day1_rep3` | 3 | 94,0 | 135 |

**Volume total no SRA:** controle 243,2 MB + tratamento 287,6 MB = **530,8 MB** (aproximadamente 506,2 MiB). As seis bibliotecas somam cerca de 1,22 bilhão de bases sequenciadas segundo a coluna `Bases`. O número real de reads de cada run deverá ser obtido a partir dos FASTQs ou do log da conversão; `AvgSpotLen` é a média por run, não o tamanho fixo de todas as reads.

**Critérios usados para reduzir possíveis fatores de confusão:** todas as seis bibliotecas estão descritas no CSV como `RNA-Seq`, `SINGLE`, `TRANSCRIPTOMIC`, organismo `Arabidopsis`, tecido `shoot`, ecótipo `Col-0`, estádio `vegetative` e equipamento `Illumina NovaSeq 6000`. Os tratamentos diferem, mas o dia de coleta é `Day1` nos dois grupos. As réplicas `rep1–rep3` são tratadas como biológicas **conforme a identificação experimental no depósito**; sua independência deve ser confirmada no método do estudo original antes de formular conclusões científicas.

**Comparação DGE:** `tratamento` (Drought_Day1) **versus** `controle` (Watered_Day1). Portanto, um `log2FoldChange` positivo indica maior expressão no grupo sob seca e um valor negativo, maior expressão no grupo irrigado.

## Estrutura do projeto

```text
PRJNA904764/
├── README.md
├── FASTQ/                 # FASTQ compactados das seis bibliotecas
├── BAM/                   # BAM alinhados, índices e relatórios STAR
├── COUNT/                 # featureCounts e matriz bruta por gene
├── DGE/                   # resultados de expressão diferencial
├── REF/                   # genoma TAIR10, anotação GTF, índice STAR
├── META/
│   ├── SraRunTable.csv    # CSV original do NCBI (139 runs; para conferência)
│   └── samples.tsv        # seleção de seis runs pronta para executar
├── QC/                    # FastQC e MultiQC
├── LOGS/                  # registros das execuções
├── SRA/                   # download local temporário
├── TMP/                   # arquivos temporários
└── scripts/               # programas
```

Todos os comandos abaixo devem ser executados na raiz `PRJNA904764/`. A pasta principal do pacote ZIP se chama `Curso_RNAseq_Arabidopsis_UEPG`; após descompactar, você pode renomeá-la para `PRJNA904764` (o nome não altera o funcionamento dos scripts).

## Etapa 0 — Instalar ferramentas e preparar o ambiente

> [!NOTE]
> Os comandos pressupõem Linux Ubuntu e ferramentas instaladas. Configure `THREADS=2` em computadores com recursos limitados.

Em Linux com `micromamba` instalado (pode-se usar `mamba` ou `conda` com a mesma lista de pacotes):

```bash
micromamba create -y -n rnaseq-uepg -c conda-forge -c bioconda \\
  sra-tools star subread samtools fastqc multiqc \\
  r-base bioconductor-deseq2 python curl
micromamba activate rnaseq-uepg
prefetch --version
fasterq-dump --version
STAR --version
featureCounts -v
samtools --version | head -1
Rscript -e 'packageVersion("DESeq2")'
```

Na pasta extraída do pacote:

```bash
# Alternativamente, renomeie a pasta extraida "Curso_RNAseq_Arabidopsis_UEPG" para "PRJNA904764".
cd PRJNA904764
mkdir -p FASTQ BAM COUNT DGE REF QC META LOGS SRA TMP
export THREADS=2
cat META/samples.tsv
```

Para máquinas com poucos recursos, comece com **2 threads** e reserve **10–15 GB de espaço livre** para arquivos SRA, FASTQ, BAM e temporários. **4 GB de RAM** são um ponto de partida a verificar, com **6–8 GB** preferíveis; os requisitos reais dependem do índice STAR, do ambiente R e do sistema. Em máquinas muito limitadas, uma alternativa didática para alinhamento é o HISAT2 (não implementado neste pacote). Instale os pacotes antes da aula para não consumir o tempo de laboratório.

## Etapa 1 — Conferir a seleção e os metadados

O arquivo `META/samples.tsv` está **pronto** no pacote; você não precisa descobrir os SRRs ou criar a matriz manualmente. Este é o seu conteúdo (separado por tabulações):

```tsv
run condition replicate
SRR22397138 controle  1
SRR22397137 controle  2
SRR22397135 controle  3
SRR22397141 tratamento  1
SRR22397140 tratamento  2
SRR22397139 tratamento  3
```

As colunas significam:

- `run`: accession SRR, usado para baixar e identificar cada FASTQ e BAM;

- `condition`: fator do DESeq2, com valores **exatamente** `controle` ou `tratamento` (essas grafias devem ser mantidas, pois o script R as espera);

- `replicate`: número da réplica dentro de cada grupo.

Para **reconstruir e conferir** o arquivo a partir do CSV enviado, execute:

```bash
python3 scripts/00_metadados.py --csv META/SraRunTable.csv --list
python3 scripts/00_metadados.py \\
  --csv META/SraRunTable.csv \\
  --control shoot_Watered_Day1 \\
  --treated shoot_Drought_Day1 \\
  --out META/samples.tsv
cat META/samples.tsv
```

O script `00_metadados.py` aceita tanto o formato `SraRunTable.csv` (colunas `Library Name`, `Assay Type`, `Organism` e `Bytes`) quanto o formato `SraRunInfo.csv` (colunas `LibraryName`, `LibraryStrategy`, `ScientificName` e `size_MB`), verifica os critérios básicos e compara os fatores biológicos presentes. Ele exige três réplicas numeradas 1, 2 e 3 por grupo, sem ambiguidades de SRR.

Também é possível baixar uma tabela atualizada do NCBI, **opcionalmente**, sem baixar reads:

```bash
curl -fL --retry 3 \\
  'https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA904764' \\
  -o META/runinfo.csv
python3 scripts/00_metadados.py --csv META/runinfo.csv --list
```

Se desejar repetir a seleção no `runinfo.csv`, use os mesmos argumentos `--control` e `--treated` apresentados anteriormente (substituindo somente o caminho CSV). Se algum campo não aparecer, utilize o arquivo `META/SraRunTable.csv` que acompanha o material. Para conferir no portal: [NCBI SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA904764).

**Discussão:** por que um run (SRR) não é necessariamente sinônimo de réplica biológica? Por que preservar o mesmo tecido, genótipo e ponto de coleta é importante para atribuir diferenças ao tratamento?

## Etapa 2 — Baixar raw reads com SRA Toolkit

> [!WARNING]
> O `fasterq-dump` utiliza espaço temporário adicional. Reserve pelo menos 10–15 GB livres e baixe somente os seis SRRs selecionados.

```bash
bash scripts/01_download.sh
ls -lh FASTQ/
```

Para cada SRR, `prefetch` baixa o SRA em `SRA/`, `fasterq-dump` converte para FASTQ single-end, `gzip -1` compacta o resultado. O script remove a cópia `.sra` **após** verificar a integridade do FASTQ compactado. O limite de download por run é 250 MB. Cada amostra é processada separadamente para economizar disco.

**Atenção ao armazenamento:** os 530,8 MB somados na tabela representam os arquivos **SRA compactados** informados pelo NCBI; `fasterq-dump` cria dados descompactados e arquivos temporários possivelmente muito maiores. Não baixe 139 experimentos de uma vez. Se um run não for single-end, interrompa e atualize a seleção; o script é intencionalmente específico para single-end.

**Discussão:** explique as quatro linhas do formato FASTQ e a distinção entre sequência e scores Phred.

## Etapa 3 — Baixar referência e indexar com STAR

> [!IMPORTANT]
> Use genoma FASTA e anotação GTF da mesma montagem e release (TAIR10/Ensembl Plants 63). Misturar referências pode comprometer o alinhamento e a contagem.

```bash
bash scripts/02_referencia.sh
ls -lh REF/
ls REF/STAR/
```

Fonte: Ensembl Plants release **63**, genoma **TAIR10**, com `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa` e `Arabidopsis_thaliana.TAIR10.63.gtf` da **mesma release**. Os comandos aplicam `--genomeSAindexNbases 12` e `--genomeSAsparseD 2` para diminuir o índice de um genoma pequeno. `--sjdbOverhang 99` é o valor inicial fixado no script para esta aula, mas os metadados indicam `AvgSpotLen` entre 135 e 146 nt. Como exercício, discuta o valor apropriado para o comprimento efetivo das reads (recomendação geral: comprimento das reads − 1) e, se necessário, ajuste o parâmetro antes de construir o índice. Não confunda o comprimento médio do run com o comprimento exato de cada read.

Se o download do arquivo GTF retornar HTTP 404, verifique o diretório oficial da release e ajuste o nome em **todos** os scripts. Não misture versões de referência e anotação. Se a máquina não tiver memória suficiente para indexar com STAR, discuta a troca por HISAT2.

## Etapa 4 — Controle de qualidade

```bash
bash scripts/03_qc.sh
# Abrir localmente QC/FASTQC/SRR*_fastqc.html
# Relatório consolidado: QC/multiqc_report.html
```

Inspecionar: distribuição de qualidade por base, conteúdo GC, adaptadores, duplicação e comprimento de reads. **Não execute trimming automaticamente sem evidência.** Se houver adaptadores/baixa qualidade, discuta e acrescente `fastp` ou `cutadapt` como extensão opcional.

## Etapa 5 — Alinhamento

```bash
bash scripts/align.sh
ls -lh BAM/*.bam
cat BAM/SRR*_Log.final.out | grep -E 'Uniquely mapped reads %|Number of input reads'
```

Alinha **uma amostra por vez**, com duas threads, lê `.fastq.gz` com `zcat`, escreve BAM ordenado por posição genômica e cria `.bai` com `samtools index`. Todos os BAM devem existir antes de contar. Confira o mapeamento único e as reads não alinhadas.

**Discussão:** diferença entre `SAM`, `BAM`, ordenação por coordenadas e alinhamento com splice junctions.

## Etapa 6 — Quantificação gene-level

> [!IMPORTANT]
> Verifique se a biblioteca é *stranded* antes de definir `STRAND`. O valor `0` assume ausência de informação de orientação.

```bash
# Padrão do exercício: não stranded. Ajustar somente após conferir a biblioteca.
export STRAND=0
bash scripts/count.sh
head -n 4 COUNT/count_matrix.tsv
cat COUNT/featureCounts.txt.summary
```

`featureCounts -t exon -g gene_id -s 0` agrega alinhamentos que atingem éxons de um gene. O script usa o formato específico do conjunto SINGLE: **não utiliza `-p`**. Caso o protocolo seja stranded, a direção correta deve ser conferida; alterar `STRAND=1` ou `STRAND=2` pode mudar drasticamente as contagens. Verifique `Assigned`, `Unassigned_NoFeatures` e `Unassigned_Ambiguity` na saída.

O script `make_matrix.py` lê cabeçalhos, associa BAM por accession real, e cria automaticamente a matriz **genes x 6 amostras**, sem supor posições fixas como coluna 13 ou 18.

**Discussão:** counts são contagens brutas, não TPM/FPKM, e o DESeq2 exige contagens brutas.

## Etapa 7 — Expressão diferencial

> [!CAUTION]
> A baixa profundidade deste conjunto limita o poder estatístico. Não interprete genes sem significância como prova de ausência de resposta à seca.

```bash
bash scripts/run_deseq2.sh
ls -lh DGE/
cat DGE/summary.txt
head -n 10 DGE/DGE_significant.tsv
```

Modelo estatístico: `design = \~ condition` e contraste **tratamento / controle**, isto é, **Drought_Day1 / Watered_Day1**. Portanto:

- `log2FoldChange > 0`: maior expressão em **Drought_Day1** (seca).

- `log2FoldChange < 0`: maior expressão em **Watered_Day1** (irrigado).

- `padj`: p-valor ajustado para múltiplos testes, por Benjamini-Hochberg.

- Lista didática `DGE_significant.tsv`: `padj < 0,05` e `abs(log2FoldChange) >= 1` (critério de filtragem posterior; o teste de hipótese padrão não exige limiar de fold-change).

Saídas:

| Arquivo | Interpretação |
|---|---|
| `DGE_all.tsv` | Todos os genes testados, com `baseMean`, `log2FoldChange`, erro padrão, estatística, `pvalue` e `padj` |
| `DGE_significant.tsv` | Genes que passaram nos limiares didáticos |
| `normalized_counts.tsv` | Contagens normalizadas para inspeção (não são a entrada do DESeq2) |
| `library_sizes.tsv` | Total de reads atribuídas por amostra |
| `summary.txt` | Quantidade de genes, resumo estatístico e versões dos pacotes |
| `dispersion.png` | Estimativas de dispersão |
| `MA.png` | Relação expressão média vs fold-change |
| `PCA.png` | Variabilidade entre as seis amostras |
| `sample_distances.png` | Distâncias entre amostras após transformação VST |
| `pvalues.png` | Distribuição dos p-values |
| `volcano.png` | Efeito (log2FC) vs significância ajustada |

**Limitação científica:** como o conjunto foi escolhido para ser pequeno, não interpretar ausência de significância como ausência de resposta biológica. Investigue amostras discrepantes no PCA e taxas baixas de mapeamento, em vez de concluir automaticamente que o experimento falhou.

## Etapa 8 — Questões

1. Explique por que o formato FASTQ tem quatro linhas por read e como os scores de qualidade são representados.

2. Descreva como o STAR identifica reads que atravessam junções éxon-éxon.

3. Calcule o percentual de reads unicamente mapeadas de cada amostra e identifique anomalias.

4. Explique por que `-t exon -g gene_id` é apropriado e por que `-t gene -g gene_name` pode produzir contagens problemáticas.

5. Verifique quantos genes passaram no pré-filtro; justifique o critério de 10 counts em 3 amostras.

6. Interprete os PC1 e PC2 no gráfico PCA. As réplicas biológicas de uma mesma condição formam grupos?

7. Diferencie `pvalue` de `padj` e explique o controle da taxa de falsas descobertas (FDR).

8. Se `log2FoldChange=2`, qual o fator de aumento no grupo tratado? E se `log2FoldChange=-1`?

9. Pesquise a função de dois genes com `padj < 0,05` e discuta sua possível relação com respostas à seca, citando TAIR/UniProt/artigos.

10. Compare a interpretação com a análise completa do artigo original e discuta os efeitos do menor número de reads.

## Referências

1. [NCBI — BioProject PRJNA904764](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA904764)

2. [NCBI — SRA Run Selector (PRJNA904764)](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA904764); metadados em `META/SraRunTable.csv`

3. [NCBI — Run SRR22397139](https://trace.ncbi.nlm.nih.gov/Traces/index.html?run=SRR22397139)

4. [Artigo — Gonzalez et al., eLife 2024](https://elifesciences.org/articles/84747)

5. [NCBI — Downloading SRA data](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/)

6. [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

7. [Ensembl Plants — Arabidopsis TAIR10](https://jun2026-plants.ensembl.org/Arabidopsis_thaliana/Info/Index)

8. [DESeq2 — vignette oficial](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
