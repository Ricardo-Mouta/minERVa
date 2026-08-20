#!/usr/bin/env python3

# Ricardo Mouta - 23/01/2025
# A ideia geral desse script é automatizar a captura da diversidade de um ERV dentro do genoma de seus hospedeiros.
# Com esses dados em mãos, será possível fazer análises downstream usando os outptus desse script. 
# Futuramente eu conectarei essas outras análises a esse script para poder ter realmente um programa de análise de ERVs. 
# Por enquanto estou nomeando ele 'minERVa', ou 'miniature ERV analyser'.
# Como usar esse script: python3 minERVa_021.py -genome 'file_name' -query 'file_name' -further arguments

import os, argparse, csv, re
from datetime import date

################################################################################################################################################################
                                                        # ARGUMENTOS DO PROGRAMA #
################################################################################################################################################################

# Aqui eu crio a variável que vai conter todos os argumentos que serão passados no início do programa.
parser = argparse.ArgumentParser(
    description= "minERVa is a program that mines and performs preliminary analyses on ERV data, with a focus on recent ERVs that retain genic regions. The following arguments can be used to suit minERVa to your needs:",
    epilog = '''This program requires an UNIX command line and requires the installation of:\n 
1. BLAST+: https://www.ncbi.nlm.nih.gov/books/NBK569861/ \n
2. MAFFT: https://mafft.cbrc.jp/alignment/software/ \n
3. IQ-TREE3: https://iqtree.github.io/ \n
4. FIGTREE: https://github.com/rambaut/figtree/releases \n
5. Python packages: os, argparse, csv, re, datetime. \n'''
    )

# E a partir daqui eu vou detalhar cada um dos argumentos que o programa reconhece.
## 1º: o genoma ou genomas a serem utilizados.
parser.add_argument('-genome', '-g', '-subject',
    type = str,
    required = True,
    nargs = '+',
    help = 'path to the input file(s) containing the genomes to be mined in .FASTA format. It is possible to input multiple genomes as a single concatenated file or in different files by specifying each path separated by a space.'
    )

## 2º: uma flag que representa se os genomas já foram previamente transformados em uma base de dados para buscas do BLAST. Se essa flag for ativada, a etapa de criar a base de dados é pulada, basicamente.
parser.add_argument('-dbready', '-ready',
    action = 'store_true',
    help = 'Signals that the file provided by the genome argument is already a database for the BLAST+ program',
    )

## 3º: o query a ser utilizado. Esse programa requer apenas que uma sequência de LTR do ERV a ser analisado seja inserido aqui.
parser.add_argument('-query','-q',
    type = str,
    required = True,
    nargs = '+',
    help = 'path to the input file containing the LTR sequence to be used as query in BLAST searches.'
    )

## 4º: a cobertura mínima a ser capturada durante a busca por BLAST. Talvez seja necessário modificar esse valor de acordo com o ERV. O valor base é 50.
parser.add_argument('-hspcoverage','-coverage',
    type = int,
    required = False,
    default = 50,
    help = 'a numerical value between 0 and 100 that specifies the threshold for High-scoring pairs coverage to be used in the BLAST search. Default is 50')

## 5º: a identidade mínima a ser capturada durante a busca por BLAST. Não utilizo a princípio, mas pode ser necessário para outros ERVs.
parser.add_argument('-identity','-id',
    type = int,
    required = False,
    default = 0,
    help = 'a numerical value between 0 and 100 that specifies the threshold for identity to be used in the BLAST search. Default is that this threshold is not used')

## 6º: a distância em nucleotídeos a ser considerada como o 'contexto genômico'. Praticamente é uma margem de segurança para garantir que o ERV inteiro será capturado ao redor do resultado do BLAST.
parser.add_argument('-context_range','-contextr',
    type = int,
    required = False,
    default = 20000,
    help = 'a numerical value representing the number of base pairs around the hit to remove duplicate results. Default is 20000 base pairs.')

## 7º: uma flag para que a base de dados do BLAST feita por esse programa seja mantida após o fim da análise.
parser.add_argument('-dbkeep', '-keep',
    action = 'store_true',
    help = 'Signals that the database prepared by this program should not be deleted by the end of the program run. This argument is ignored if the option -dbready is active.',
    )

## 8º: a distância em nucleotídeos a ser considerada como a 'margem de erro' do BLAST. Praticamente será utilizada para compensar possíveis imprecisões na captura das LTRs pelo BLAST.
parser.add_argument('-LTR_range','-LTRr',
    type = int,
    required = False,
    default = 500,
    help = 'a numerical value representing the number of base pairs around the hit to ensure a perfect trim of the ERV limits. Default is 500 base pairs.')

## 9º: Um apelido para a corrida. Será usado como como um prefixo para os arquivos e nome da pasta em que eles serão reunidos ao fim da corrida.
parser.add_argument('-prefix', '-run', '-nick',
    type = str,
    required = False,
    default = 'ERV',
    help = 'A string that contains a prefix or nickname for a single minERVa run.')

## 10º: a flag to signal that the BLAST search should be performed remotelly.
parser.add_argument('-remote',
    action = 'store_true',
    help = 'Signals that the the BLAST search should be performed remotelly. If this argument is given, the query must be given locally, but instead of a file, the -genome argument must be fed with a genome accession number'
    )

## jogar todos esses argumentos numa variável simples.
args = parser.parse_args()

################################################################################################################################################################
                                                            # NOMEAÇÃO DE ARQUIVOS #
################################################################################################################################################################
# LOGISTICA: para adicionar o dia de hoje aos arquivos gerados.
hoje = date.today().strftime('%Y_%m_%d')

# Nomear os outputs para otimizar o chamamento deles ao longo do script
BLAST_results = f"{hoje}_blt_{args.prefix}_results.txt" # capturar o nome do arquivo de resultado do BLAST
BLAST_location_LTR = f"{hoje}_bla_{args.prefix}_loc_LTR.txt" # capturar o nome do arquivo de localização de LTRs
BLAST_location_genic = f'{hoje}_bla_{args.prefix}_loc_blast_genic.txt' # capturar o nome do arquivo de localização de cópias gênicas
FASTA_regions_LTR = f"{hoje}_ltr_{args.prefix}_regions_LTR.fasta" # capturar o nome do arquivo com as sequências de LTRs solo
FASTA_regions_genic = f"{hoje}_gen_{args.prefix}_regions_genic.fasta" # capturar o nome do arquivo com as sequências gênicas
FASTA_sig_LTR = f"{hoje}_ltr_{args.prefix}_signature_LTR.fasta" # capturar o nome dos FASTAs cortados pelas assinaturas
FASTA_sig_genic = f"{hoje}_gen_{args.prefix}_signature_genic.fasta" # capturar o nome dos FASTAs cortados pelas assinaturas
FASTA_aln = f"{hoje}_aln_{args.prefix}_ref.fasta" # capturar o nome do arquivo fasta alinhado
RESULT_arv = f"{hoje}_arv_{args.prefix}_ref" # capturar o nome geral para os resultados do IQ-TREE
RESULT_arv_anotada = f'{hoje}_arv_{args.prefix}_ref_annotated.treefile' # capturar o arquivo da árvore anotada
RESULT_directory = f"{hoje}_{args.prefix}_results" # capturar o nome da pasta onde ficarão os resultados

################################################################################################################################################################
                                                            # BUSCA POR BLAST #
################################################################################################################################################################
# 1. DECISÃO: Essa estrutura de decisão é responsável por realizar a busca por BLAST. Ela segue a seguinte lógica:
## 2. BLAST: SE a base de dados estiver pronta, o programa automaticamente faz um blast usando essa base de dados.
if args.remote == True: # SE o genoma estiver remoto
    if args.identity == 0: # BLAST cobertura 
        os.system('blastn -db refseq_genomes -query {} -remote -entrez_query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {}'.format(*args.query, *args.genome, BLAST_results, args.hspcoverage))
    else: # BLAST com critério de similaridade
        os.system('blastn -db refseq_genomes -query {} -remote -entrez_query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {} -perc_identity {}'.format(*args.query, *args.genome, BLAST_results, args.hspcoverage, args.identity))
elif args.dbready == True: # SE a base ded dados está pronta localmente
    if args.identity == 0: # BLAST cobertura
        command = 'blastn -db {} -query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {}'.format(*args.genome,*args.query,BLAST_results, args.hspcoverage)
        os.system(command)
    else: # BLAST similaridade
        command = 'blastn -db {} -query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {} -perc_identity {}'.format(*args.genome,*args.query,BLAST_results, args.hspcoverage,args.identity)
        os.system(command)

## 3. DECISÃO: SE a base de dados não estiver pronta, ele vai produzir uma, e em seguida fazer um BLAST contra ela:
else: 
    
    ### 3.1.1. CRIAÇÃO DE BANCO DE DADOS: SE houver mais de um arquivo de genoma, concatena os arquivos e cria uma base de dados com o arquivo concatenado.
    if len(args.genome)>=2: 
        command = 'cat {} > genomas.fasta'.format(' '.join(args.genome)) 
        os.system(command) ### concatenar arquivos
        command = 'makeblastdb -in genomas.fasta -dbtype nucl -parse_seqids -out db_{} -title db_{}'.format(args.prefix, args.prefix) 
        os.system(command) ### base de dados (com o arquivo concatenado)
    
    ### 3.1.2. CRIAÇÃO DE BANCO DE DADOS: SE houver apenas um genoma, prepara o banco de dados diretamente.
    else: 
        os.system('makeblastdb -in {} -dbtype nucl -parse_seqids -out db_{} -title db_{}'.format(*args.genome,args.prefix, args.prefix))  ### base de dados (com input do usuário)

    ### 3.2 BLAST: Usa o banco de dados preparado para realizar a busca por BLAST
    if args.identity == 0: ### BLAST cobertura
        os.system('blastn -db db_{} -query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {}'.format(args.prefix, *args.query, BLAST_results,args.hspcoverage)) 
    else: ### BLAST similaridade
        os.system('blastn -db db_{} -query {} -out {} -outfmt "6 sseqid sstart send sstrand" -qcov_hsp_perc {} -perc_identity {}'.format(args.prefix, *args.query, BLAST_results,args.hspcoverage,args.identity)) 

################################################################################################################################################################
                                                            # RECUPERAÇÃO DE FASTAS #
################################################################################################################################################################
# Nessa etapa, o resultado do BLAST é utilizado para obter as sequências FASTA correspondentes ao hit e mais uma margem de erro ao redor.
# 1. INPUT: Ler o resultado do BLAST dentro do Python
with open(f'{BLAST_results}') as input:
    blast = csv.reader(input.readlines(), delimiter = '\t')

# 2. VARIÁVEIS
possivelLTR = 'nada'
anterior = 'nada'
inicio = 0
fim = 0
fita = 'nada'

# 3. LOOPING: Produzir um arquivo reformatado a partir do resultado do BLAST:
## 3.1. LOOPING: criar e manter aberto o arquivo
locLTR = open(f'{BLAST_location_LTR}','w')

with open(f'{BLAST_location_genic}','w') as loc:

## 3.2. LOOPING: para cada linha do resultado do BLAST (BLAST_results.txt)
    for line in sorted(blast):
## 3.3. FILTRAGEM: completar sequências que possuam duas LTRs.
        if anterior == str(line[0]) and inicio - args.context_range  < int(line[1]) < fim + args.context_range: # SE este hit está na mesma sequência que o hit anterior E o início da sequência está dentro da margem de contexto viral estabelecida (default = 20.000 nt)? 
            if str(line[3]) == 'plus' or int(line[1]) < int(line[2]): # SE a cópia está na fita direta
                if inicio > int(line[1]): # SE a posição inicial deste hit for a 5' da posição inicial do hit anterior
                    if  int(line[1]) > args.LTR_range: # SE a posição inicial for maior que margem de erro
                        inicio =  int(line[1]) - args.LTR_range # subtrair margem de erro dela
                    else: # SE a posição inicial for menor que margem de erro
                        inicio = 1 # substituir por 1
                if fim < int(line[2]): # SE a posição final deste hit for a 3' da posição final do hit anterior
                    fim = int(line[2]) + args.LTR_range # acrescentar margem de erro à posição final
            else:
                if inicio > int(line[2]):
                    if  int(line[2]) > args.LTR_range: # SE a posição inicial for maior que margem de erro
                        inicio =  int(line[1]) - args.LTR_range # subtrair margem de erro dela
                    else: # SE a posição inicial for menor que margem de erro
                        inicio = 1 # substituir por 1
                if fim < int(line[1]):
                    fim = int(line[1]) + args.LTR_range # acrescentar margem de erro à posição final
            possivelLTR = ''
## 3.4. Escrita dos arquivos
        else:            
            if possivelLTR == '' and anterior != 'nada':
                loc.write(f'{anterior} {inicio}-{fim} {fita}\n') # escrever a nova linha do arquivo genico
            elif anterior != 'nada':
                locLTR.write(possivelLTR) # escrever a nova linha do arquivo de LTR solo
# 3.5. LOOPING: resultados do BLAST são rearranjados para o formato novo
            # 3.5.1. REFORMATAÇÃO: captura a linha para resultados na fita líder
            anterior = str(line[0])
            if str(line[3]) == 'plus' or int(line[1]) < int(line[2]):
                if  int(line[1]) > args.LTR_range: # SE a posição inicial for maior que margem de erro
                    inicio =  int(line[1]) - args.LTR_range # subtrair margem de erro dela
                else: # SE a posição inicial for menor que margem de erro
                    inicio = 1 # substituir por 1
                fim =  int(line[2]) + args.LTR_range # acrescentar margem de erro à posição final
                fita = line[3]
            
            # 3.5.2. REFORMATAÇÃO: captura a linha para resultados na fita reversa
            else:
                if int(line[2]) > args.LTR_range: # SE a posição inicial for maior que margem de erro
                    inicio = int(line[2]) - args.LTR_range # subtrair margem de erro dela
                else: # SE a posição inicial for menor que margem de erro
                    inicio = 1 # substituir por 1
                fim = int(line[1]) + args.LTR_range # adicionar margem de erro na posição final
                fita = line[3]
            possivelLTR = (f'{anterior} {inicio}-{fim} {fita}\n') # escrever a nova linha do arquivo novo
# 3.6 Toques finais: para escrever a última linha
    if possivelLTR == '' and anterior != 'nada':
        loc.write(f'{anterior} {inicio}-{fim} {fita}\n') # escrever a nova linha do arquivo genico
    elif anterior != 'nada':
        locLTR.write(possivelLTR) # escrever a nova linha do arquivo de LTR solo

## 3.6. FECHAMENTO: fechar o arquivo reformatado (loc_blast.txt) pronto.
locLTR.close()
loc.close() 
input.close()

## 3.7. Sumário dos resultados para ser mostrado na linha de comando:
with open(f'{BLAST_location_genic}','r') as loc:
    lines = len(loc.readlines())
    print('Total Number of copies containing genic regions:', lines)

with open(f'{BLAST_location_LTR}','r') as locLTR:
    lines = len(locLTR.readlines())
    print('Total Number of solo LTR copies:', lines)

locLTR.close()
loc.close()

## 3.8. FASTA: Usar os arquivos reformatados (loc_blast_genic.txt e loc_blast_all.txt) pra retornar a sequência com o contexto genômico em formato FASTA.
os.system(f'blastdbcmd -db db_{args.prefix} -entry_batch {BLAST_location_LTR} -out {FASTA_regions_LTR} -line_length 99999')
os.system(f'blastdbcmd -db db_{args.prefix} -entry_batch {BLAST_location_genic} -out {FASTA_regions_genic} -line_length 99999')

# O arquivo '{prefix}_regions_genic.fasta' contém o contexto genômico dos hits encontrados pelo BLAST com duas LTRs no mesmo alcance no formato FASTA.
# O arquivo '{prefix}_regions_LTR.fasta' contém o contexto genômico dos hits encontrados pelo BLAST com uma LTR solo no formato FASTA.

################################################################################################################################################################
                                                            # DELIMITAÇÃO DE SEQUÊNCIAS #
################################################################################################################################################################

# INPUT: abrir o arquivo de contexto genômico e salvar seu texto no terminal
copies = [f'{FASTA_regions_LTR}',f'{FASTA_regions_genic}']

for type in copies:
    with open('{}'.format(type),'r+') as f:
        ervloc = f.readlines()

        ervnovo = []

        PATTERNS = [
        (re.compile(r"[ :\.-]"), "_"),              # Substitui caracteres especiais
        (re.compile(r",.+$|isolate_"), ""),     # Remove texto após vírgula, número de acesso e isolate
        (re.compile(r"chromosome_"), "c"),          # Abrevia "chromosome"
    ]
        f.close()

    # REFORMATAÇÃO: remoção de caracteres especiais e expressões comuns dos headers dos FASTAS.
    for line in ervloc:
        if line.startswith('>'):
            for encontrar, substituir in PATTERNS:
                line = encontrar.sub(substituir, line)
        ervnovo.append(line)

    ervloc = ()

    # OUTPUT: reescrever o arquivo de contexto genômico (erv_regions.fasta) com as modificações
    with open(f'{type}','w') as f:
        f.writelines(ervnovo) # escrever o arquivo de fato.
    f.close()

# LOOPING: definir os limites dos ERVs através das assinaturas.
# VARIAVEIS
header = ' '
count = 0
seqatual = ' '
seqtotal = ' '
gen_headers = [] # TESTE: figtree anotado

# OUTPUT: escrever o arquivo trimado contendo apenas os ERVs (quando for possível)
with open(FASTA_regions_genic, 'r') as genicos:
    gen_lines = genicos.readlines()
    genicos.close()

with open(f'{FASTA_sig_genic}', 'w') as assinaturado:
     
    # LOOPING: a cada linha do fasta original
    for line in gen_lines:
        
    # DECISÃO: se a linha for um header, capturá-lo. 
        if '>' in line:
            header = line.rstrip()
            gen_headers.append(header) # TESTE: figtree anotado
            # DECISÃO: se a linha conter a sequência, analisar e reescrever as linhas
        else:
            seqatual = line
                
                # DECISÃO: se as assinaturas estiverem na posição correta, reescrever a linha trimando o fasta e adicionando a assinatura ao header.
            if seqatual[(args.LTR_range-4):args.LTR_range] == seqatual[-(args.LTR_range+1):-(args.LTR_range-3)] and seqatual[(args.LTR_range-4):args.LTR_range] != '' and seqatual[-(args.LTR_range+1):-(args.LTR_range-3)] != '':                    
                count += 1
                seqtotal = f'{header}_{seqatual[496:500]}\n{seqatual[496:-497]}\n'
            elif seqatual[(args.LTR_range-6):args.LTR_range] == seqatual[-(args.LTR_range+1):-(args.LTR_range-5)] and seqatual[(args.LTR_range-6):args.LTR_range] != '' and seqatual[-(args.LTR_range+1):-(args.LTR_range-5)] != '':
                count += 1
                seqtotal = f'{header}_{seqatual[494:500]}\n{seqatual[494:-495]}\n'
                # DECISÃO: se a assinatura não for encontrada, simplesmente adicionar '_NS' ao fim do header.
            else:
                seqtotal = f'{header}_NS\n{seqatual}\n'
                
                # OUTPUT: adicionar a nova linha ao arquivo novo.
            assinaturado.write(seqtotal)
        # OUTPUT: fechar o arquivo reescrito

    assinaturado.close()
# OUTPUT:  mostrar quantas assinaturas foram encontradas
print(f'Signatures: {count} found in {len(gen_lines)/2:.0f} input sequences.')

################################################################################################################################################################
                                                            # FILOGENIA INICIAL #
################################################################################################################################################################
# FILOGENIA: série de comandos UNIX para chamar programas que irão: 
## 1. alinhar os resultados e as referências,
## 2. inferir uma árvore filogenética

os.system(f'mafft --auto --thread -1 {FASTA_sig_genic} > {FASTA_aln}')
os.system(f'mafft --thread -1 --add {FASTA_aln} --keeplength aln_retroviridae_references.fasta > {RESULT_arv}.fasta')
os.system(f'iqtree3 -s ./{RESULT_arv}.fasta -m MFP -bb 10000 -alrt 1000 -T AUTO --undo --prefix {RESULT_arv}')
os.system(f'figtree ./{RESULT_arv}.treefile')

## FILOGENIA: visualizar a árvore filogenética através do Figtree. Por enquanto ainda estou tendo problemas quanto à abertura dessa árvore pelo figtree, mas estou trabalhando para melhorar.
# TESTE: figtree anotado
'''taxanexus = "'\n'".join(gen_headers)

with open(f'{RESULT_arv}.treefile','r') as tree:
    nexus = tree.read()
    tree_annotation = f"""#NEXUS
#NEXUS
begin taxa;
	dimensions ntax={27+len(gen_headers)};
	taxlabels
    'Alpharetrovirus_avileu|Gallus_gallus'[$!color=#4ecb8d]
    'Betaretrovirus_murmamtum|Mus_musculus'[$!color=#ff9d3a]
    'Betaretrovirus_ovijaa|Ovis_aries'[$!color=#ff9d3a]
	'Deltaretrovirus_priTlym1|Chlorocebus_tantalus'[$!color=#f9e858]
	'Deltaretrovirus_priTlym2|Homo_sapiens'[$!color=#f9e858]
	'Deltaretrovirus_priTlym3|Cercocebus_torquatus'[$!color=#f9e858]
	'Deltaretrovirus_bovleu|Bos_taurus'[$!color=#f9e858]
	'Epsilonretrovirus_waldersar|Sander_vitreus'[$!color=#000dff]
	'Gammaretrovirus_gibleu|Hylobates_sp'[$!color=#ff73b6]
	'Gammaretrovirus_woomonsar|Lagothrix_sp'[$!color=#ff73b6]
	'Gammaretrovirus_felleu|Felis_catus'[$!color=#ff73b6]
	'Gammaretrovirus_koa|Phascolarctos_cinereus'[$!color=#ff73b6]
	'Gammaretrovirus_murleu|Mus_musculus'[$!color=#ff73b6]
	'Lentivirus_capartenc|Capra_hircus'[$!color=#c701ff]
	'Lentivirus_felimdef|Felis_catus'[$!color=#c701ff]
	'Lentivirus_humimdef1|Homo_sapiens'[$!color=#c701ff]
	'Lentivirus_equinfane|Equus_caballus'[$!color=#c701ff]
	'Lentivirus_humimdef2|Homo_sapiens'[$!color=#c701ff]
	'Lentivirus_pum|Puma_concolor'[$!color=#c701ff]
	'Lentivirus_simimdef|Chlorocebus_aethiops'[$!color=#c701ff]
	'Equispumavirus_equcab|Equus_caballus'[$!color=#d83034]
	'Bovispumavirus_bostau|Bos_taurus'[$!color=#d83034]
	'Felispumavirus_felcat|Felis_catus'[$!color=#d83034]
	'Prosimiispumavirus_otocra|Otolemur_crassicaudatus'[$!color=#d83034]
	'Simiispumavirus_atespp|Ateles_sp'[$!color=#d83034]
	'Simiispumavirus_gorgorgor|Gorilla_gorilla_gorilla'[$!color=#d83034]
	'Simiispumavirus_sapxan|Sapajus_xanthosternos'[$!color=#d83034]
    {taxanexus}
;
end;

Begin Trees;
    Tree {args.prefix} = [&R] {nexus}
End;

Begin FigTree;
   Set	appearance.colorAttribute="User selection";
   Set	Legend.isShown=true;
   Set	Legend.attributeName="rvgroup";
   Set	NodeLabels.isShown="true";
   Set	NodeLabels.displayAttribute="label";
   Set	NodeLabels.fontSize=10;
   Set	NodeLabels.fontName "Arial";
   set	trees.order="true";
   set	trees.orderType="increasing";
   set  trees.rootingType="Midpoint";


End;"""

tree.close()

with open(f'{RESULT_arv_anotada}', 'w') as trenot:
    trenot.write(tree_annotation)
    trenot.close()
'''

################################################################################################################################################################
                                                            # ARQUIVOS TEMPORÁRIOS #
################################################################################################################################################################
# LOGISTICA: Criar diretório para os resultados
os.system(f"mkdir {RESULT_directory}")
os.system(f'mkdir {RESULT_arv}')

# LOGISTICA: Jogar todos os arquivos gerados para esse diretório
os.system(f'mv {BLAST_results} {BLAST_location_LTR} {BLAST_location_genic} {FASTA_regions_genic} {FASTA_regions_LTR} {FASTA_aln} {FASTA_sig_genic} {RESULT_directory}/')
os.system(f'mv {RESULT_arv}.* {RESULT_arv}/')
os.system(f'mv {RESULT_arv}/ {RESULT_directory}/{RESULT_arv}/')
if args.dbready == False:
    os.system(f'mv db_{args.prefix}* {RESULT_directory}')


# Ideia para anotação:
#* Aceitar novo input:
# * annotation OU anot
#### múltiplas sequências como input
#
#* Delimitação de sequência
## Retira múltiplos dos resultados para encontrar a região genômica 
### Vou ter que aumentar o limite de vezes que a repretição pode ocorrer
## OU usa o resultado da versão atual como base de dados direto 
#
#* Anotação propriamente dita
## Faz um novo BLAST usando as regiões genômicas como database, com parâmetros relaxados (25% identidade)
#### Cria FASTAS para cada gene separadamente
#### Cria tabela com cada cópia e o conteúdo dela 
# formato: 
# Copy  location    GAG gag_length  POL pol_length  ENV env_length)
# {header da sequência/subject do blast}    posição de início-fim   Y/N comprimento no subject  Y/N comprimento no subject  Y/N comprimento no subject

