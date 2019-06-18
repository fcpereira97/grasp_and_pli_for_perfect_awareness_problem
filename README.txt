Este diretório contém arquivos de implementação e de soluções
do projeto final da disciplina MO824 - Tópicos em Otimização Combinatória,
ministrada pelo professor Dr. Fábio Usberti.

Discente: Felipe de Carvalho Pereira (230214)

Repositório do projeto no GitHub:
https://github.com/fcpereira97/grasp_and_pli_for_perfect_awareness_problem.git

Descrição dos subdiretórios e arquivos:

* codigos-fonte/

	- grasp-pa.cpp

		Implementação na linguagem C++ de heurística baseada em
		GRASP para o Perfect Awareness Problem. Foram implementadas três diferentes
		estratégias para a fase de construção: standard, random plus greedy e 
		sampled greedy. Para a fase de busca local, duas estratégias do tipo
		first improvemente foram implementadas, embora apenas uma dela tenha
		sido utilizada nos experimentos. Os detalhes técnicos serão abordados
		no relatório final do projeto.


	- pli-pa.cpp

		Implementação na linguagem C++ de um modelo PLI para o Perfect Awareness Problem.
		A implementação faz uso de bibliotecas do resolvedor comercial Gurobi para solucionar
		a instância fornecida de maneira exata. A descrição do modelo estará presente no
		relatório do projeto.

* executaveis/

	- grasp-pa-executavel

		Executavel produzido a partir da implementação grasp-pa.cpp. O comando de
		execução deve obedecer o seguinte formato:

		./grasp-pa-executavel A B C D time_limit input_path output_path

		O argumento A é um inteiro correspondente à estratégia que será utilizada
		na fase de construção (1 para standard, 2 para random plus greedy e 3 para o
		sampled greedy).

		O argumento B é um inteiro que corresponde à quantidade de vértices que
		serão inseridos na solução parcial a casa passo durante a fase de construção.

		O argumento C é um inteiro que corresponde à estratégia de busca local que será
		utilizada (1 para a busca local nº 1 e 2 para a busca local nº 2).

		O argumento D é um número real que corresponde a porcentagem da vizinhança que
		será avaliada na busca local nº 2, caso esta seja escolhida.

		O argumento time_limit é um inteiro que correspondem ao tempo de máximo de execução,
		em segundos.

		Os argumentos input_path e output_path são strings que correspondem,
		respectivamente, aos caminhos do arquivo de entrada (instância) e arquivo
		de saída.

		O arquivo de entrada deve estar no formato:

		N M
		v_i v_j
		.
		.
		.
		v_k v_l

		onde N é o número de vértices, M é o número de arestas e as M linhas seguintes
		correspondem a cada um dos pares de vértices que descrevem cada aresta do grafo.

		O arquivo de saída conterá:
			+ Estratégia de construção utilizada
			+ Estratégia de busca local utilizada
			+ Instância
			+ Número de iterações realizadas
			+ Média dos valores das soluções obtidas a cada iteração
			+ Valor da melhor solução obtida
			+ Número de rodadas da propagação da melhor solução obtida
			+ Quantidade de nós com status "conhecedor" ao final da propagação da melhor solução obtida
			+ Nº da iteração em que a melhor solução foi obtida

		No console, um log será impresso com estas informações,
		o valor de cada nova solução incumbent e também a descrição da melhor
		solução obtida, isto é, os vértices do seed set.


	- pli-pa-executavel

		Executavel produzido a partir da implementação pli-pa.cpp. O comando de
		execução deve obedecer o seguinte formato:

		./pli-pa-executavel input_path

		O argumento input_path é uma string que contém o caminho do arquivo
		de entrada (instância).

		No console, um log padrão do Gurobi será impresso com as informações
		sobre o processo de execução, além dos limitantes primais e duais obtidos
		e da solução obtida.

		Serão produzidos dois arquivos de saída:
			+ "model.lp" contendo o programa PLI da instância executada
			+ "model" contendo a descrição da solução obtida

* outputs/

	- grasp-pa-outputs.csv

		Compilado de outputs das execuções das herísticas sobre as instâncias. As colunas representam:

		Instancia - instância utilizada
		Construcao - estratégia de construção utilizada
		Busca local - estratégia de busca local utilizada
		N iteracoes - número de iterações realizadas
		Media solucoes - média dos valores das soluções obtidas a cada iteração
		Melhor solucao - valor da melhor solução obtida
		N rodadas melhor sol - número de rodadas da propagação da melhor solução obtida
		N conhecedores - quantidade de nós com status "conhecedor" ao final da propagação da melhor solução obtida
		Iteracao melhor sol - nº da iteração em que a melhor solução foi obtida

	- pli-pa-karate-model.lp
	- pli-pa-karate-solution.sol
	- pli-pa-karate-output.txt

		Arquivos de saída da execução do modelo PLI para a instância karate.gml
		O primeiro arquivo correspondem ao modelo PLI, o segundo à solução
		obtida e o último ao output padrão do Gurobi.

	- pli-pa-jazz-model.lp
	- pli-pa-jazz-solution.sol
	- pli-pa-jazz-output.txt

		Arquivos de saída da execução do modelo PLI para a instância jazz.net
		O primeiro arquivo correspondem ao modelo PLI, o segundo à solução
		obtida e o último ao output padrão do Gurobi.

* logs/
	
	- log_grasp_standard.txt
	- log_grasp_sampled_greedy.txt
	- log_grasp_random_plus_greedy.txt

		Arquivos de log obtidos do console das execuções de todas as instâncias
		para o a heurística GRASP implementada, variando o método de construção.
		O primeiro arquivo contém o log da execução do GRASP com a estratégia 
		standard, o segundo com a estratégia sampled greedy e o terceiro com
		random plus greedy.

* instancias/

	- descricao_instancias.csv

		Descrição das instâncias utilizadas nas execuções. O arquivo
		corresponde a uma tabela com os seguintes dados:

		+ Nome da instância
		+ Número de vértices
		+ Número de arestas
		+ Desidade do grafo

		Os arquivos das instâncias podem ser consultados no repositório do projeto
		no GitHub.