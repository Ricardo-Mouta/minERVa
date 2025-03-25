# Ricardo Mouta - 23/01/2025
# Essa é a minha primeira tentativa de fazer um script de Python. 
# A ideia é simples: automatizar o uso de alguns comandos do BLAST 
# para sair com um alinhamento dos resultados. Com o tempo, quero 
# ir acresecntando mais opções a esse script para que ele seja 
# mais facilmente modificável.

# Como usar esse script: python3 BLAST_RM -genome 'file_name' -query 'file_name' -further arguments

import os, argparse, csv

# Os argumentos usados aqui são o genoma e o query para o BLAST. 
parser = argparse.ArgumentParser(
    description= "minERVa is a program that mines and performs preliminary analyses on ERV data, with a focus on recent ERVs that retain genic regions",
    epilog = 'This program requires an UNIX command line and requires the installation of the BLAST+ program. Instructions for installation are available at the "BLAST® Command Line Applications User Manual"\n.'
    )

# Dessa forma, estou tentando possibilitar que múltiplos genomas sejam adicionados de uma vez. 
parser.add_argument('-genome', '-g',
    type = str,
    required = True,
    nargs = '+',
    help = 'path to the input file(s) containing the genomes to be mined in .FASTA format. It is possible to input multiple genomes as a single concatenated file or in different files by specifying each path separated by a space.'
    )

parser.add_argument('-dbready', '-ready',
    action = 'store_true',
    help = 'Signals that the file provided by the genome argument is already a database for the BLAST+ program',
    )

parser.add_argument('-query','-q',
    type = str,
    required = True,
    nargs = '+',
    help = 'path to the input file(s) containing the sequences to be used as query in BLAST searches. It is possible to input multiple queries as a single concatenated file or in different files by specifying each path separated by a space.'
    )

parser.add_argument('-hspcoverage','-coverage',
    type = int,
    required = False,
    default = 50,
    help = 'a numerical value between 0 and 100 that specifies the threshold for High-scoring pairs coverage to be used in the BLAST search. Default is 50')

parser.add_argument('-identity','-id',
    type = int,
    required = False,
    default = 0,
    help = 'a numerical value between 0 and 100 that specifies the threshold for identity to be used in the BLAST search. Default is that this threshold is not used')

parser.add_argument('-context_range','-range',
    type = int,
    required = False,
    default = 8000,
    help = 'a numerical value representing the number of base pairs to extend the hit in order to extract the entire ERV from its genome and to remove duplicate results. Default is 8000 base pairs.')

args = parser.parse_args()

if args.dbready == False:
    parser.add_argument('-dbkeep', '-keep',
        action = 'store_true',
        help = 'Signals that the database prepared by this program should not be deleted by the end of the program run. This argument is ignored if the option -dbready is active.',
        )
    keep = parser.parse_args()

print('O(s) genoma(s) usado(s) foi(ram) {}'.format(' '.join(args.genome)))
print('Você usou {} arquivo(s) para formar sua base de dados.'.format(len(args.genome)))
print('O query usado foi {}'.format(*args.query))
print('A cobertura  usada foi {}'.format(args.hspcoverage))
print('A identidade usada foi {}'.format(args.identity))

# Essa estrutura de decisão é responsável por realizar a busca por BLAST. Ela segue a seguinte lógica:
if args.dbready == True: ## SE for usado o argumento de que o banco de dados está pronto, o programa automaticamente faz um blast usando esse banco de dados.
    if args.identity == 0: # e não houver critério de similaridade
        command = 'blastn -db {} -query {} -out BLAST_results.txt -outfmt "6 sseqid sstart send" -qcov_hsp_perc {}'.format(*args.genome,*args.query,args.hspcoverage)
        os.system(command)
    else: # e houver critério de similaridade
        command = 'blastn -db {} -query {} -out BLAST_results.txt -outfmt "6 sseqid sstart send" -qcov_hsp_perc {} -perc_identity {}'.format(*args.genome,*args.query,args.hspcoverage,args.identity)
        os.system(command)
else: ## SE o banco de dados não estiver pronto:
    if len(args.genome)>=2: ### SE houver mais de um arquivo de genoma,
        command = 'cat {} > genomas.fasta'.format(' '.join(args.genome)) 
        os.system(command) ### o programa concatena esses arquivos
        command = 'makeblastdb -in genomas.fasta -dbtype nucl -parse_seqids -out dbtemp -title dbtemp' 
        os.system(command) ### e em seguida prepara o banco de dados com os arquivos concatenados.
    else: ### SE houver apenas um genoma
        command = 'makeblastdb -in {} -dbtype nucl -parse_seqids -out dbtemp -title dbtemp'.format(*args.genome)
        os.system(command)  ### o programa prepara o banco de dados diretamente
    if args.identity == 0: ### Em seguida, o programa vai usar o banco de dados preparado para realizar a busca por BLAST
        command = 'blastn -db dbtemp -query {} -out BLAST_results.txt -outfmt "6 sseqid sstart send" -qcov_hsp_perc {}'.format(*args.query,args.hspcoverage)
        os.system(command)
    else:
        command = 'blastn -db dbtemp -query {} -out BLAST_results.txt -outfmt "6 sseqid sstart send" -qcov_hsp_perc {} -perc_identity {}'.format(*args.query,args.hspcoverage,args.identity)
        os.system(command)
    args.genome = 'dbtemp'

# Agora vou ler o resultado do BLAST dentro do Python
with open('BLAST_results.txt') as input:
    blast = csv.reader(input.readlines(), delimiter = '\t')

anterior = 'nada'
inicio = 0
# E vou escrever um novo arquivo formatado para ser usado como input do blastdbcmd:
with open('loc_blast.txt','w') as loc:
    for line in blast: # a cada linha do output do BLAST
        if anterior == line[0] and inicio -args.context_range < int(line[1]) < inicio +args.context_range:
            print('resultado redundante removido')
        else:
            anterior = line[0]
            inicio = int(line[1])
            if int(line[1]) >  int(line[2]): # SE está na fita reversa
                line2 = int(line[1]) + args.context_range # adicionar 8000 na posição final
                if int(line[2]) > args.context_range: # SE a posição inicial for maior que 8000
                    line1 = int(line[2]) - args.context_range # subtrair 8000 dela
                else: # SE a posição inicial for menor que 8000
                    line1 = 1 # substituir por 1
                loc.write('{} {}-{} minus\n'.format(line[0],line1,line2)) # escrever a nova linha do arquivo novo
            else: # SE está na fita líder
                line2 =  int(line[2]) + args.context_range # acrescerntar 8000 a posição final
                if  int(line[1]) > args.context_range: # SE a posição inicial for maior que 8000
                    line1 =  int(line[1]) - args.context_range # subtrair 8000 dela
                else: # SE a posição inicial for menor que 8000
                    line1 = 1 # substituir por 1
                loc.write('{} {}-{} plus\n'.format(line[0],line1,line2)) # escrever a nova linha do arquivo novo
            

loc.close() # fechar o arquivo para escrita

'''
# para ler o arquivo final
loc = open("loc_blast.txt", "r")
print(loc.read())
'''

# Depois que o arquivo estiver pronto, mandar ele pro blasdbcmd pra retornar a sequência com o contexto genômico em formato FASTA.
command = 'blastdbcmd -db ./{} -entry_batch loc_blast.txt -out erv_regions.txt'.format(''.join(args.genome))
os.system(command)

'''
k = open('erv_regions.txt','r')
print(k.read())
k.close()
'''

if keep.dbkeep == False:
    import re
    def purge(dir, pattern):
        for f in os.listdir(dir):
            if re.search(pattern, f):
                os.remove(os.path.join(dir, f))
    purge('.','dbtemp*+')

# BLAST com LTR. Preciso de:
## 1. argumento para passar o arquivo fasta com a LTR
## 2. comando do BLAST
## 3. 
'''
O próximo passo a partir daqui seria a anotação das cópias, que por enquanto estou fazendo pelo Geneious.
Quando eu aprender como fazer isso por linha de comando eu posso adicionar esse passo aqui. 
São alguns objetivos para essa etapa:
# 1. Anotar as regiões do genoma de enFeLV - bem ou mal daria pra fazer com o BLAST
# 2. Capturar a assinatura genômica de cada cópia e adicionar ao header do FASTA correspondente
# 3. Checar se a assinatura genômica está igual no início e fim das LTRs - mensagem de erro + orientação de fazer manualmente caso isso não seja possível
# 4. Cortar as cópias de enFeLV no início e fim das LTRs

Em seguida, quero fazer uma filogenia preliminar:
# 1. Alinhar sequências
# 2. Trimar
# 3. Fazer a filogenia
# 4. Abrir o arquivo final
'''