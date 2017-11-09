#include <iostream>
#include <string>
#include <fstream>
#include <iostream> 
#include <regex>
#include <iterator>

#include <locale> 

// codigo terminado para realizar o parseTree do arquivo recebido

// observações sobre
/*

	Qtd '(' = numero de nos internos
	nos internos = folhas -1
	nos internos * 2 + 1 = total nos
*/

using namespace std;

const int MAXN = 100010; // numero maximo  de elementos que um arquivo newick pode conter
int numTotalNos =0, numNosInternos =0, numNosFolhas = 0, numElemPut =0, tamanhoArrays =0; // dados sobre nos da arvore

// array para salvar caracteriscas da arvore
// seu tamnho sera determinado dentro da função parsingTree 
string * species;
int * pai;
int * filhoEsquerda;
int * filhoDireita;
float * comprimentoRamo;
float * comprimentoRamoRaiz;




int fileLines =0;
int filePutLines =0;

string fileLine; // salva cada caracter do arquivo, em um vetor de char
string filePut[MAXN];
void parsingTree();
void preencherVetor(int * vetor, const int SIZE, const int VALOR);
void preencherVetor(string * vetor, const int SIZE, const string VALOR);

int main()
{

	ifstream infile;
	infile.open("newick_input.tree");

	if(!infile){

		cout<<"Error opening newick file..." <<endl;
	}

	else{

		char c;
  		while (infile.get(c)) {
  			fileLine +=c;
  			fileLines++;

  		}
  		

	}



	infile.close();


	infile.open("put.list");
	if(!infile){

		cout<<"Error opening put file..." <<endl;
	}

	else{

		
	    while (! infile.eof() ) //enquanto end of file for false continua
	    {
	      getline (infile,filePut[filePutLines]); // como foi aberto em modo texto(padrão)
	                             //e não binário(ios::bin) pega cada linha
	      filePutLines++;
	    }


	}

	infile.close();



	parsingTree();

	return 0;
}

void preencherVetor(int * vetor, const int SIZE, const int VALOR){

	for (int i = 0; i < SIZE; ++i)
	{
		vetor[i] = VALOR;
	}

}

void preencherVetor(float * vetor, const int SIZE, const int VALOR){

	for (int i = 0; i < SIZE; ++i)
	{
		vetor[i] = VALOR;
	}

}

void preencherVetor(string * vetor, const int SIZE, const string VALOR){

	for (int i = 0; i < SIZE; ++i)
	{
		vetor[i] = VALOR;
	}
}

void parsingTree(){

	int aParen = 0, fParen =0; // parentes abrindo e fechando '(' ')'
	int quantElementosFile=0, posPai = -1, posPaiFile = -1;
	int comma=0;
	string folha =" ", interno =" ", elementoAtual=" ", elemenAnterior=" ", elemPai = " "; // salva o atual e o ultimo elemento
	string filhoLeft=" ", filhoRight= " ", comprimeRamoLeft ="", comprimeRamoRight = "";
	int auxiliarNumNos =0, auxiliarGeral =0, auxilarPreencherVetor =0; // usado para fazer as trocas de elementos no vetor
	int indexFilhoLeft =-1, indexFilhoRight =-1, indexStartReplace = -1;
	bool alphabeticModeOn = false; // 
	// regex 
	std::smatch m;
  	std::regex e ("\\([^(^)]+\\)");
  	std::regex folhas("\\([A-z_]+|,[A-z_]+"); // achar todas as folhas e separar no vetor
  	std::regex internos("\\)[A-z_]+");
	quantElementosFile = fileLines; // qnts elementos o arquivo tem
	
	
	// primeira varredura apenas para verificar inconsistencias
	for(int i = 0; i < quantElementosFile; i++){ // faz uma varredura no arquivo
		elementoAtual = fileLine[i];
		if(elementoAtual == "(") aParen++;
		if(elementoAtual == ")") fParen++;
		if(elementoAtual == ",") comma++;
	}

	if(aParen != fParen){

		cout<< "Arquivo inconsistente, parentes não balanceados" <<endl;
		throw;
	}

	numNosInternos = aParen; // nos internos = numero de parenteses abertos
	numTotalNos = numNosInternos * 2 +1; // total de nos na arvore
	numNosFolhas = numNosInternos +1; // nos folhas

	tamanhoArrays = numTotalNos + (filePutLines * 2) + 1;
	// vetores q irao salvar infor sobre arvore
	species = new string[tamanhoArrays];
	pai = new int[tamanhoArrays];
	filhoEsquerda = new int[tamanhoArrays];
	filhoDireita = new int[tamanhoArrays];
	comprimentoRamo = new float[tamanhoArrays];
	comprimentoRamoRaiz = new float[tamanhoArrays];

	preencherVetor(species, tamanhoArrays, "#");
	preencherVetor(pai, tamanhoArrays, -2);
	preencherVetor(filhoEsquerda, tamanhoArrays, -2);
	preencherVetor(filhoDireita, tamanhoArrays, -2);
	preencherVetor(comprimentoRamo, tamanhoArrays, 0);
	preencherVetor(comprimentoRamoRaiz, tamanhoArrays, 0);


  	// preencher vetor com todas as species
	// usando o regex para pegar todos os folhas
	string copyNewick = fileLine;
	while (std::regex_search (copyNewick,m,folhas)) {

		folha = "";
	    for (int i=0; i<m.size(); ++i) {

	    	auxiliarGeral = m.position(i) +1; // posicão do match
	    	folha = copyNewick[auxiliarGeral];
	    	alphabeticModeOn = true;


	    	while(alphabeticModeOn) {
	    		auxiliarGeral++;
	    		elementoAtual = copyNewick[auxiliarGeral];


	    		// conferir se o proximo elemento ainda é valido
	    		std::locale loc;
				for (std::string::iterator it=elementoAtual.begin(); it!=elementoAtual.end(); ++it)
				{
				    
				    if ( (std::isalpha(*it,loc)) or (elementoAtual=="_") or (elementoAtual=="-") ){

				    	folha += elementoAtual;
				      
					}
					else{

						alphabeticModeOn = false; // terminou a folha
					}
				}
	    		
	    	}
  		}

  		species[auxilarPreencherVetor++] = folha;

	    copyNewick = m.suffix().str();

  	}

  	// preencher vetor com todas as species
	// usando o regex para pegar todos os nos internos
	copyNewick = fileLine;
	auxilarPreencherVetor = numNosFolhas + (filePutLines * 2) + 1;
	while (std::regex_search (copyNewick,m,internos)) {

		interno = "";
	    for (int i=0; i<m.size(); ++i) {

	    	auxiliarGeral = m.position(i) +1; // posicão do match
	    	interno = copyNewick[auxiliarGeral];
	    	alphabeticModeOn = true;


	    	while(alphabeticModeOn) {
	    		auxiliarGeral++;
	    		elementoAtual = copyNewick[auxiliarGeral];


	    		// conferir se o proximo elemento ainda é valido
	    		std::locale loc;
				for (std::string::iterator it=elementoAtual.begin(); it!=elementoAtual.end(); ++it)
				{
				    
				    if ( (std::isalpha(*it,loc)) or (elementoAtual=="_") or (elementoAtual=="-") ){

				    	interno += elementoAtual;
				      
					}
					else{

						alphabeticModeOn = false; // terminou a interno
					}
				}
	    		
	    	}
  		}

  		species[auxilarPreencherVetor++] = interno;

	    copyNewick = m.suffix().str();

  	}

  	// preencher novos put
  	string auxiliarPut[2], auxiliar, put;

  	for (int linePut = 0; linePut < filePutLines; linePut++)
  	{

  		auxiliar = filePut[linePut];

  		put = ""; 
  		auxiliarGeral = 0;
  		alphabeticModeOn = false;
 		
	    for (int elemenIndex = 0; elemenIndex < auxiliar.length(); elemenIndex++)
	    {	

	       if ( isspace(auxiliar[elemenIndex]) and alphabeticModeOn) 
	       {
	       		auxiliarPut[auxiliarGeral++] = put;
	       		
	           	put = "";
	           	alphabeticModeOn = false;
	       }

	       else{

	       		if ( !isspace(auxiliar[elemenIndex]) ){ 

	       			alphabeticModeOn = true;
	        		put += auxiliar[elemenIndex];
	        	}
	        	
	       }
	 
	    }

	    if(put != ""){

	    	auxiliarPut[auxiliarGeral] = put;
	    }

	    // insert no array especies

	    species[numNosFolhas + linePut] = auxiliarPut[0];


	    for (int index = 0; index < tamanhoArrays; index++)
	    {
	    	if(species[index] == auxiliarPut[1]){

	    		pai[numNosFolhas + linePut ] = index;
	    		break;
	    	}
	    }
	    
  	}

  	pai[tamanhoArrays-1] = -1; // no raiz não tem um pai
  	
  	/*

  	for (int i = 0; i < tamanhoArrays; i++)
  	{
  		cout<<species[i] <<" ";
  	}
  	cout<<endl;

  	for (int i = 0; i < tamanhoArrays; i++)
  	{
  		cout<<pai[i] <<" ";
  	}
  	cout<<endl;
  	*/

  	




	// logica se da no principio de achar todos os nos folhas pares, em cada loop, dai verificamos o seu devido pai
	// e os "eliminamos" da arvore, criando novos filhos folhas.
	// Para isso, estamos usando a biblioteca Redex, para achar os matchs e fazer o replace em seguida.
	// links: http://www.cplusplus.com/reference/regex/regex_search/
	//		  http://www.cplusplus.com/reference/regex/match_results/position/
	//        http://www.cplusplus.com/reference/regex/regex_replace/
	
	
	//regex logica
	// enquanto tivermos nos para buscar, vamos tirar as folhas
	// sobrara no final apenas o pai raiz

	while(auxiliarNumNos < numTotalNos -1){

		bool filhoLeftAtivo = false, filhoRighttAtivo = false;
		bool ramoLeftAtivo = false, ramoRightAtivo = false, filhoDone = false;
		filhoLeft = "";
		filhoRight = "";
		comprimeRamoLeft = "";
		comprimeRamoRight = "";
		std::regex_search ( fileLine, m, e );
    	
    	elementoAtual = fileLine[m.position(0)]; // primeiro paranteses dos nos folhas achados
    	
    	auxiliarGeral = m.position(0);

    	indexStartReplace = m.position(0);
    	// achar entao o pai dos elementos(folhas)
    	while(elementoAtual !=")"){ // enquanto não achar o parentese entao que fecha tais folhas

    		// passar entao em todos os elementos dentro do '()'


    		std::locale loc;
			for (std::string::iterator it=elementoAtual.begin(); it!=elementoAtual.end(); ++it)
			{
			    
			    if ( (std::isalpha(*it,loc)) or (elementoAtual=="_") or (elementoAtual=="-") ){

			    	if( ! filhoRighttAtivo){

    					filhoLeft += elementoAtual;
    					filhoLeftAtivo = true;
    				
	    			}
	    			else{

	    				filhoRight += elementoAtual;
	    				
	    			}

	    			filhoDone = false;
			      
				}
				else{ 

					// se terminou de achar todo o nome do filho.
					if(! filhoDone){
						// mudar para filho direito, e setar comprimento ramo esquerdo ativo
						if( filhoLeftAtivo){

							filhoRighttAtivo = true;
							filhoLeftAtivo = false;

							ramoLeftAtivo = true;
							ramoRightAtivo = false;
						}
						else if( filhoRighttAtivo){

							filhoRighttAtivo = false;
							filhoLeftAtivo = true;

							ramoRightAtivo = true;
							ramoLeftAtivo = false;
						}

						filhoDone = true;
					}

					// comprimento ramo
					// pegar o valor do comprimento do ramo
					// pegar caracter por caracter, e no final transformar em float
					if ( (std::isdigit(*it,loc)) or (elementoAtual==".") ){

						if(ramoLeftAtivo){

							comprimeRamoLeft += elementoAtual;

						}

						else if (ramoRightAtivo){

							comprimeRamoRight += elementoAtual;

						}

					}
					
				}
				
			}
    	
    		// proximo elemento
    		auxiliarGeral++;
    		elementoAtual = fileLine[auxiliarGeral];

    	}

    	posPaiFile = auxiliarGeral + 1;
    	elemPai = fileLine[posPaiFile]; // pega a posição do pai daqueles filhos, que esta logo apos fechar o parenteses
  		
  		// achar o index entao dos filhos tirados e do pai

    	for(int i=0; i<tamanhoArrays; i++){

    		if(species[i] == elemPai){

    			posPai = i;

    			if( (indexFilhoLeft != -1) and (indexFilhoRight != -1) ) break; // parar se ja achou indexes
    		}


    		else if(species[i]==filhoRight){

    			
    			indexFilhoRight = i;
    			

    			if( (indexFilhoLeft != -1) and (posPai != -1) ) break; 
    		}

    		else if(species[i]==filhoLeft){

    			
    			indexFilhoLeft = i;
    			

    			if( (indexFilhoRight != -1) and (posPai != -1) ) break;
    		}

    	}

    	// preencher vetores

    	pai[indexFilhoLeft] = posPai;
    	pai[indexFilhoRight] = posPai;
    	filhoEsquerda[posPai] = indexFilhoLeft;
    	filhoDireita[posPai] = indexFilhoRight;

    	// comprimento do ramo
    	try{
	    	comprimentoRamo[indexFilhoRight] = ::atof(comprimeRamoRight.c_str());
	    	comprimentoRamo[indexFilhoLeft] = ::atof(comprimeRamoLeft.c_str());
    	}catch(exception e){

    	}

	  	// fazer o replace dos elementos retirado
	  	for (int i = indexStartReplace ; i < posPaiFile; i++)
	  	{
	  		fileLine[i] = ' ';
	  		
	  	}

	  	posPai = -1;
	  	// reset variaveis
  		filhoRight = "";
  		filhoLeft = "";
  		comprimeRamoLeft = "";
  		comprimeRamoRight = "";
  		indexFilhoRight = -1;
  		indexFilhoLeft = -1;
		auxiliarNumNos = auxiliarNumNos + 2; // ou seja, foi retirado 2 filhos

	}

	// Calcular comprimento do ramo ate a raiz
	// usando busca em profundidade
	bool done = false;
	bool folhaDone = false;
	bool rightModeOn = false;
	float compRamo = 0;
	int posRamo = tamanhoArrays -1;

	while(not done){

		// primeiramente, faz uma busca profunda, pela esquerda(mas na vdd tanto faz), e busca um no folha
		// com isso, sabemos a profundidade de todos os outros folhas, restando então apenas os nos internos
		// essa regra se aplica apenas para arvores filogeneticas


		while(not folhaDone){

			if(filhoDireita[posRamo] < 0){ // ou seja, não tem filho(folha)
				folhaDone = true;

				// temos então o comprimento de todos os folhas da arvore
				// atualizar de todas as folhas então
				for (int i = 0; i < numNosFolhas + filePutLines; i++)
				{
					comprimentoRamoRaiz[i] = compRamo;
				}

				posRamo = pai[posRamo]; // volta entao a posição ramo 1 posição, pois chegou no limite da arvore(folha)
				compRamo -= comprimentoRamo[filhoDireita[posRamo]]; // volta entao o comp ramo para o anterior
				
				rightModeOn = true; // loop esta para o lado direito
				break;
			}

			compRamo += comprimentoRamo[filhoDireita[posRamo]]; // ou seja, o comprimento daquele filho ate  
																// o pai atual do loop

			posRamo = 	filhoDireita[posRamo]; // proximo filho a direita												
			
		}
		

		// fazer a busca em profundidade agr para os nos internos


		// se os dois filhos da raiz, ja tiverem seus comprimentos achados,
		// entao significa q a busca em profundidade foi concluida
		if( (comprimentoRamoRaiz[filhoDireita[tamanhoArrays-1]] > 0) and 
			(comprimentoRamoRaiz[filhoEsquerda[tamanhoArrays-1]] > 0)   ) {

			done = true;
			break;
		}


		// // se busca esta para o lado direito.....
		// if(rightModeOn){

		// cheka se elemento atual ainda tem filho 
		if( (filhoDireita[posRamo] >=0 or filhoEsquerda[posRamo] >=0)  ){

			// se tiver filho da direita e o comp dele ainda n foi calculado
			if(filhoDireita[posRamo] >=0 and comprimentoRamoRaiz[filhoDireita[posRamo]] <= 0 ){

				compRamo += comprimentoRamo[filhoDireita[posRamo]]; // ou seja, o comprimento daquele filho ate  
																// o pai atual do loop
				// nova posRamo é entao aquele filho da direita
				posRamo = filhoDireita[posRamo];
				rightModeOn = true; // loop foi para esquerda
			}

			 // se tiver filho da esquerda e o comp dele ainda n foi calculado
			else if(filhoEsquerda[posRamo] >=0 and comprimentoRamoRaiz[filhoEsquerda[posRamo]] <= 0 ){

				compRamo += comprimentoRamo[filhoEsquerda[posRamo]]; 
				// nova posRamo é entao aquele filho da direita
				posRamo = filhoEsquerda[posRamo];
				rightModeOn = false; // loop foi pra esquerda
			}

			// ou seja, aquela sub arvore esta concluida
			else{


				comprimentoRamoRaiz[posRamo] = compRamo;

				posRamo = pai[posRamo]; // volta entao a posição ramo 1 posição, pois chegou no limite da arvore(folha)
				
				if (rightModeOn) // se estavamos para direita
				{
					compRamo -= comprimentoRamo[filhoDireita[posRamo]];
				}
				else{

					compRamo -= comprimentoRamo[filhoEsquerda[posRamo]];
				}
				

			}
	
		}


		
		

	}

	for(int i=0; i<tamanhoArrays; i++){

		cout<<species[i]<<" ";
	}
	cout<<endl;

	for(int i=0; i<tamanhoArrays; i++){

		cout<<pai[i]<<" ";
	}
	cout<<endl;


	for(int i=0; i<tamanhoArrays; i++){

		cout<<filhoEsquerda[i]<<" ";
	}


	cout<<endl;

	for(int i=0; i<tamanhoArrays; i++){

		cout<<filhoDireita[i]<<" ";
	}
	cout<<endl;


	for(int i=0; i<tamanhoArrays; i++){

		cout<<comprimentoRamo[i]<<" ";
	}
	cout<<endl;

	for(int i=0; i<tamanhoArrays; i++){

		cout<<comprimentoRamoRaiz[i]<<" ";
	}
	cout<<endl;


	

}




