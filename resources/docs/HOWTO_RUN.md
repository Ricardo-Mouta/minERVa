# Como Executar o Pipeline minERVa (Estratégia Atual)

Este guia descreve a execução do minERVa utilizando a infraestrutura atual, que tira proveito do download automático de genomas, bancos de dados locais temporários e submissão de trabalhos para clusters HPC com SLURM.

## Pré-requisitos
Certifique-se de que o **Snakemake** e o **plugin SLURM** estejam instalados no seu ambiente conda. 
Caso não tenha um ambiente, você pode criar um a partir do arquivo fornecido:
```bash
conda env create -f envs/snakemake.yaml
conda activate snakemake-base
```

## Configuração do Pipeline
Antes de rodar, verifique e atualize os parâmetros no arquivo `config.yaml`.
- **`genomes_list`**: O arquivo de texto contendo os Accessions do NCBI (ex: `genomes.txt`) a serem analisados, um por linha.
- **`gbk_dir`**: Diretório contendo os arquivos no formato GenBank.
- **`targets`**: Lista dos alvos virais a serem buscados (ex: `pol`, `gag`, `LTR`).
- **`blast_params`**: As opções de busca do BLAST.

## Comportamento do Pipeline
Na configuração atual, o minERVa otimiza recursos e evita sobrecargas:
- **Download Local:** O pipeline baixa automaticamente cada genoma usando a ferramenta oficial `ncbi-datasets-cli`.
- **Tolerância a Falhas:** Caso a conexão com a API do NCBI sofra quedas intermitentes, o sistema realiza até 3 tentativas automáticas de download (`retries: 3`).
- **Economia de Disco:** Os arquivos FASTA dos genomas e os diretórios do banco BLAST são gerados apenas pelo tempo necessário (`temp()`). Eles são automaticamente apagados assim que todas as buscas contra o genoma correspondente finalizam.

## Rodando no SLURM
Para submeter o fluxo completo em um cluster SLURM de maneira automatizada, onde o Snakemake irá gerenciar as dependências e despachar as regras (`download_genome`, `makeblastdb`, `local_blastn`, etc) como *jobs* individuais:

```bash
snakemake --use-conda --executor slurm --jobs 50
```

Se o seu servidor exigir parâmetros específicos como conta ou partição:
```bash
snakemake --use-conda --executor slurm --default-resources slurm_account=SUA_CONTA slurm_partition=SUA_PARTICAO --jobs 50
```

> **Dica:** O uso da flag `--rerun-incomplete` é sempre recomendado para descartar execuções anteriores que possam ter sido interrompidas de forma abrupta:
> ```bash
> snakemake --use-conda --executor slurm --jobs 50 --rerun-incomplete
> ```
