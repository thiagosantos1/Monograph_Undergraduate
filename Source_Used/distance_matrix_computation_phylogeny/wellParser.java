//import java.io.BufferedReader;
//import java.io.File;
//import java.io.FileNotFoundException;
//import java.io.FileReader;
//import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.StringTokenizer;
import java.io.*;

public class wellParser {
	private int MAX_FILHOS = 2;
	private int TAM_Z = 5000;
	private int pos = 0;
	private int lp = 0;
	private int rp = TAM_Z - 1;
	private int nfol = 0;
	private int idx_ni = 0; // indice para acrescentar no interno
	private int ind_insercao = 0;
	private int n_insercoes = 0;
	private String newick = null;
	private int z[] = null;    // configuração
	private float zb[] = null; // banching
	private int zq[] = null; // qtd de filhos na sub-árvore
	private int ze[] = null; // qtd de especies na sub-árvore
	private int nz[] = null; // new z array
	private float nz_br[] = null; // banching
	private int nz_qf[] = null; // qtd de filhos na sub-árvore
	private int nz_qe[] = null; // qtd de especies na sub-árvore
	private int nz_p[] = null; // pai
	private int nz_f1[] = null; // filho 1
	private int nz_f2[] = null; // filho 2
	private HashMap<String, Integer> nomes;
	private HashMap<Integer, String> idsNomes;
	private HashMap<String, Integer> meus_nomes;
	private HashMap<Integer, String> meus_idsNomes;
	private int ids;
	private int dummy_ids;
	private String vazio;
	private String dummy_name;

	public static void main(String[] args) {
		
		File folder = new File("./trees");
    	File[] listOfFiles = folder.listFiles();
    	String arquivo;
    	int ntrees = listOfFiles.length;
		
		wellParser pf = new wellParser();
		
		arquivo = "./trees/" + listOfFiles[0].getName();
		pf.parseArvore(arquivo);
	//	pf.parseNovas();
		pf.Gera_Arquivo(arquivo, false, ntrees);
		System.out.println(arquivo + " ");
		
		for (int i = 1; i < listOfFiles.length; i++) {
			arquivo = "./trees/" + listOfFiles[i].getName();
			pf.parseArvore(arquivo);
		//	pf.parseNovas();
			pf.Gera_Arquivo(arquivo, true, ntrees);
			System.out.println(arquivo + " ");
		}
	}	
	
	public void parseArvore(String arquivo) {
		try {
		//	File f = new File("newick_input.tree");
			File f = new File(arquivo);
			BufferedReader br = new BufferedReader(new FileReader(f));
			String linha = br.readLine();
			
			parseNewick(linha);
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void parseNewick(String linha) {
		// inicializa parse
		this.pos = 0;
		this.newick = linha;
		this.z = new int[this.TAM_Z];
		this.zb = new float[this.TAM_Z];
		this.zq = new int[this.TAM_Z];
		this.ze = new int[this.TAM_Z];
		this.ids = 0;
		this.dummy_ids = 0;
		this.nomes = new HashMap<String, Integer>();
		this.idsNomes = new HashMap<Integer, String>();
		this.vazio = "";
		this.dummy_name = "X";
		this.lp = 0;
		this.rp = TAM_Z - 1;
		this.nfol = 0;
		this.idx_ni = 0; // indice para acrescentar no interno
		this.ind_insercao = 0;
		this.n_insercoes = 0;
		this.meus_nomes = new HashMap<String, Integer>();
		this.meus_idsNomes = new HashMap<Integer, String>();
		this.nz = new int[this.TAM_Z];
		this.nz_br = new float[this.TAM_Z];
		this.nz_qf = new int[this.TAM_Z];
		this.nz_qe = new int[this.TAM_Z];
		this.nz_p = new int[this.TAM_Z];
		this.nz_f1 = new int[this.TAM_Z];
		this.nz_f2 = new int[this.TAM_Z];

		// inicializa vetor z
		for (int i=0; i<this.TAM_Z; i++) {
		   z[i] = -2;
		   nz[i] = -2;
		}
		
		// le sub-árvore
		int idRaiz = leSubarvore();
		
		this.idx_ni = this.nfol;
		imprime(idRaiz);
		
		// ponto virgula final
		pos++;
	}
	
	private void imprime(int idRaiz) {
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", i);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.print("----");
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", z[i]);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", zq[i]);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%.5f ", zb[i]);
		//System.out.println();

		//System.out.println( toNewick(idRaiz) );		
		
		//Set<String> keys = this.nomes.keySet();
		//for (String nome : keys) {
		//	Integer id = this.nomes.get(nome);
		//	System.out.printf("%d -> [%s]\n", id, nome);
		//}
		
		//System.out.println( toPreOrder(idRaiz, -1) );
		String guarda1 = toPreOrder(idRaiz, -1);
		
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", nz[i]);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", nz_p[i]);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", nz_f1[i]);
		//System.out.println();
		//for (int i=0; i<this.TAM_Z; i++) System.out.printf("%d ", nz_f2[i]);
		//System.out.println();
		
		//Set<String> meus_keys = this.meus_nomes.keySet();
		//for (String nome : meus_keys) {
		//	Integer id = this.meus_nomes.get(nome);
		//	System.out.printf("%d -> [%s]\n", id, nome);
		//}
		//System.out.println( Meu_toNewick(TAM_Z - 1) );
		String guarda2 = Meu_toNewick(TAM_Z - 1);
	}

	private int leSubarvore() {
		if (Character.isLetter(this.newick.charAt(this.pos))) {
			Object folha[] = leNo();
			
			// Dados da folha
			int idFolha      = Integer.parseInt(folha[0].toString());
			float brFolha    = Float.parseFloat(folha[1].toString());
			
			// Branching da folha
			this.zb[idFolha] = brFolha;
			
			// Qtd de nós aqui
			//this.zq[idFolha] = 1;
			this.zq[idFolha] = 0; // conta altura
			this.ze[idFolha] = 1; // conta especies
			
			// Atribui como sendo raiz, por default
			this.z[idFolha] = -1;
			
			// Conta quantas folhas tem a arvores
			this.nfol++;
			
			return idFolha;
		}
		else if (this.newick.charAt(this.pos) == '(') {
			int filhos[] = new int[MAX_FILHOS];
			int posFi = 0;
			int zqParcial = 0;
			int zeParcial = 0;
			int maior = 0;
			
			do {
				this.pos++; // salta o caractere separador
				
				if (posFi > 1) {
					System.out.println("Not a binary tree: more than 2 children.");
					//continue;
					System.exit(0);
				}
				filhos[posFi] = leSubarvore();
			//	zqParcial += this.zq[filhos[posFi]];
				zeParcial += this.ze[filhos[posFi]];
				if (this.zq[filhos[posFi]] > maior) maior = this.zq[filhos[posFi]];
				posFi++;
			} 
			while (newick.charAt(this.pos) == ',');
			
			if (posFi == 1) {
					System.out.println("Not a binary tree: only one child.");
					System.exit(0);
			}
			
			zqParcial += maior;
			
			// salta o parentese fechando
			this.pos++; 

			// dados da raiz da sub-arvore
			Object raiz[] = leNo(); 
			
			int idRaiz = Integer.parseInt(raiz[0].toString());
			float brRaiz = Float.parseFloat(raiz[1].toString());
			
			// Guarda o branching 
			this.zb[idRaiz] = brRaiz;
			
			// Guarda a quantidade de nós nesta sub-arvore
			this.zq[idRaiz] = zqParcial;
			this.zq[idRaiz]++;   // conta o no (interno) corrente tambem

			// Guarda a quantidade de especies nesta sub-arvore
			this.ze[idRaiz] = zeParcial;

			for (int j=0; j<posFi; j++) {
				this.z[filhos[j]] = idRaiz;
			}
			
			// Por default, cada raiz é considerada raiz da árvore
			// Se ela nao for raiz, quando chegar no parse do seu pai, o -1 será substituido
			// pelo id do pai
			this.z[idRaiz] = -1;

			return idRaiz;
		}
		
		return -1;
	}

	private Object[] leNo() {
		int ini = this.pos;
		Object no[] = new Object[2];
		
		while (pos < this.newick.length()-1 && this.newick.charAt(pos) != ':') this.pos++;

		// Nao tem rotulo
		if (ini == pos) {
			this.vazio += " ";
			no[0] = idNo(this.vazio);
			no[0] = idNo(this.dummy_name + this.dummy_ids++);
		}
		else {
			no[0] = idNo(this.newick.substring(ini, this.pos));
		}

		// salta os dois pontos
		this.pos++;
		ini = this.pos;
		
		while (this.pos < this.newick.length()-1 && (Character.isDigit(newick.charAt(pos)) || newick.charAt(pos) == '.')) this.pos++;
		
		if (ini == pos) {
			no[1] = 0;
		}
		else {
			no[1] = this.newick.substring(ini, this.pos);
		}
		
		return no;
	}
	
	public void parseNovas() {
		try {
			File f = new File("novasEspecies.list");
			BufferedReader br = new BufferedReader(new FileReader(f));
			String linha = null;
			
			this.ind_insercao = this.nfol;
			while ((linha = br.readLine()) != null) {
				StringTokenizer st = new StringTokenizer(linha, "\t");
				
				String parte1 = st.nextToken();
				String parte2 = st.nextToken();
				
				StringTokenizer esps = new StringTokenizer(parte1, ",");

				if (esps.countTokens() == 2) {
					String token1 = esps.nextToken();
					String token2 = esps.nextToken();
					int id1 = Meu_idNo(parte2, 0);
					int id2 = Meu_idNo(token1, 1);
					int id3 = Meu_idNo(token2, 1);
			
			        this.nz[id2] = id1;
			        this.nz[id3] = id1;
			
					//System.out.println("Esp1: " + idNo(esps.nextToken()) + ", Esp2: " + idNo(esps.nextToken()) + " -> " + idNo(parte2));
					// -fica- System.out.printf("%d %d %d %d\n", 2, idNo(parte2), idNo(token1), idNo(token2));
					// -fica- System.out.printf("%d %d %d %d\n", 2, id1, id2, id3);
					// Meu_idNo(nome1, 1) no children
					// Conta quantos nos tem a arvores
			        this.nfol += 2;
			        this.n_insercoes += 2;
				}
				else {
					String token1 = esps.nextToken();
					int id1 = Meu_idNo(parte2, 0);
					int id2 = Meu_idNo(token1, 1);
			
			        this.nz[id2] = id1;
			        this.nz_f1[id2] = -2;
					this.nz_f2[id2] = -2;
					
					//System.out.println("Esp: " + idNo(esps.nextToken()) + " -> " + idNo(parte2));
					// -fica- System.out.printf("%d %d %d -1\n", 1, idNo(parte2), idNo(token1));
					// -fica- System.out.printf("%d %d %d -1\n", 1, Meu_idNo(parte2,0), Meu_idNo(token1,1));
			        // Conta quantos nos tem a arvores
			        this.nfol++;
			        this.n_insercoes++;
				}
			}
			
			br.close();
			//System.out.printf("\n %d \n", this.nfol);
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private int idNo(String nome) {
		if (this.nomes.get(nome) != null) {
			return (int) this.nomes.get(nome);
		}
		else {
			int id = this.ids++;
			this.nomes.put(nome, id);
			this.idsNomes.put(id, nome);
			return id;
		}
	}
		
	private int Meu_idNo(String nome, int tipo_no) {
		if (this.meus_nomes.get(nome) != null) {
			return (int) this.meus_nomes.get(nome);
		}
		else {
		    int id;
		    if (tipo_no == 0) { // no interno
		    	id = this.rp--;
		    } else { // no folha
		    	id = this.lp++;
		    }
			this.meus_nomes.put(nome, id);
			this.meus_idsNomes.put(id, nome);
			return id;
		}
	}
	
	private String toNewick(int idRaiz) {
		ArrayList<Integer> children = getChildren(idRaiz);
		
		if (children.size() == 0) { /* Não tem filhos */ 
			return this.idsNomes.get(idRaiz).trim() + ":" + this.zb[idRaiz];
		}
		else { /* Tem filhos */
			String saida = "(" + toNewick(children.get(0));
			
			for (int i=1; i<children.size(); i++)
				saida = saida + "," + toNewick(children.get(i));
			
			saida += ")";
			
			String nome = this.idsNomes.get(idRaiz).trim(); 
			if (this.zb[idRaiz] > 0 || nome.length() > 0) saida += nome + ":" + this.zb[idRaiz];
			
			return saida;
		}
	}
	
	private String Meu_toNewick(int idRaiz) {
		
		if (this.nz_f1[idRaiz] < 0) { /* Não tem filhos */ 
			return this.meus_idsNomes.get(idRaiz).trim() + ":" + this.nz_br[idRaiz];
		}
		else { /* Tem filhos */
			String saida = "(" + Meu_toNewick(this.nz_f1[idRaiz]);
			
			saida = saida + "," + Meu_toNewick(this.nz_f2[idRaiz]);
			
			saida += ")";
			
			String nome = this.meus_idsNomes.get(idRaiz).trim(); 
			if (this.nz_br[idRaiz] > 0 || nome.length() > 0) saida += nome + ":" + this.nz_br[idRaiz];
			
			return saida;
		}
	}
	
	private ArrayList<Integer> getChildren(int idRaiz) {
		ArrayList<Integer> children = new ArrayList<Integer>();
		
		for (int i=0; i<this.TAM_Z; i++) if (this.z[i] == idRaiz) children.add(i);
		
		return children;
	}
	
	private String toPreOrder(int idRaiz, int pai) {
		ArrayList<Integer> children = getChildren(idRaiz);
		
		if (children.size() == 0) { /* Não tem filhos */ 
			String nome1 = this.idsNomes.get(idRaiz).trim();
			int meu_id = Meu_idNo(nome1, 1);
			
			this.nz[meu_id] = meu_id;
			this.nz_br[meu_id] = this.zb[idRaiz];
			this.nz_qf[meu_id] = this.zq[idRaiz];
			this.nz_qe[meu_id] = this.ze[idRaiz];
			this.nz_p[meu_id] = pai;
			this.nz_f1[meu_id] = -2;
			this.nz_f2[meu_id] = -2;
			return this.idsNomes.get(idRaiz).trim() + ":" + this.zb[idRaiz];
		}
		else { /* Tem filhos */
		    String saida = "";
		    String nome0 = this.idsNomes.get(idRaiz).trim();
			int meu_id = Meu_idNo(nome0, 0);
			
			this.nz[meu_id] = meu_id;
			this.nz_br[meu_id] = this.zb[idRaiz];
			this.nz_qf[meu_id] = this.zq[idRaiz];
			this.nz_qe[meu_id] = this.ze[idRaiz];
		    this.nz_p[meu_id] = pai;
			if (this.zb[idRaiz] > 0 || nome0.length() > 0) {
			   saida += nome0 + ":" + this.zb[idRaiz];
		    }
		    int f1 = children.get(0);
		    int f1_id = 0;
		    String nomef1 = this.idsNomes.get(f1).trim();
		    ArrayList<Integer> granchildren = getChildren(f1);
		    if (granchildren.size() == 0) { /* Não tem filhos */ 
			   f1_id = Meu_idNo(nomef1, 1);
			} else {
			   f1_id = Meu_idNo(nomef1, 0);
			}
		    this.nz_f1[meu_id] = f1_id;
			saida += "(" + toPreOrder(f1,meu_id);
			
			for (int i=1; i<children.size(); i++) {
				int f2 = children.get(i);
				int f2_id = 0;
				String nomef2 = this.idsNomes.get(f2).trim();
                ArrayList<Integer> granchild = getChildren(f2);
		        if (granchild.size() == 0) { /* Não tem filhos */ 
			   		f2_id = Meu_idNo(nomef2, 1);
				} else {
			   		f2_id = Meu_idNo(nomef2, 0);
				}
		    	this.nz_f2[meu_id] = f2_id;
				saida = saida + "," + toPreOrder(f2,meu_id);
			}
			saida += ")";
			return saida;
		}
	}
	
	public void Gera_Arquivo(String arquivo, boolean append, int ntrees) {
		int nulo = -1;
	// Stream to write file
		FileOutputStream fout;		

		try
		{
		    // Open an output stream
		    fout = new FileOutputStream ("wellParser.out", append);
		    String linha = "";

			int valor = 0;
            int nnos = 2 * this.nfol;
            int dif = TAM_Z - nnos;
            
            this.idx_ni = nnos - this.idx_ni;
            if (append) {
            	linha = linha + arquivo + "\n" + nnos + " " + this.idx_ni + "\n";
            } else {
            	linha = linha + ntrees + " " + nnos + "\n" + arquivo + "\n" + nnos + " " + this.idx_ni + "\n";
            }
            linha = linha + this.ind_insercao + " " + this.n_insercoes++ + "\n";
            // vetor principal
            for(int i=0; i<this.nfol; i++) {
                valor = nz[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
                valor = nz[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            linha += "\n";
            // vetor com simbolos
            for(int i=0; i<this.nfol; i++) {
            	if (nz[i] >= 0) {
            		linha = linha + this.meus_idsNomes.get(i).trim() + " ";
            	} else {
            		linha = linha + "# ";
            	}
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
            	if (nz[i] >= 0) {
            		linha = linha + this.meus_idsNomes.get(i).trim() + " ";
            	} else {
            		linha = linha + "# ";
            	}
            }
            linha += "\n";
            // vetor de distancia nos ramos
            for(int i=0; i<this.nfol; i++) {
            	linha = linha + nz_br[i] + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
            	linha = linha + nz_br[i] + " ";
            }
            linha += "\n";
            // vetor de quantidade de filhos
            for(int i=0; i<this.nfol; i++) {
            	linha = linha + nz_qf[i] + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
            	linha = linha + nz_qf[i] + " ";
            }
            linha += "\n";
            // vetor de quantidade de especies
            for(int i=0; i<this.nfol; i++) {
            	linha = linha + nz_qe[i] + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
            	linha = linha + nz_qe[i] + " ";
            }
            linha += "\n";
            // pai
            for(int i=0; i<this.nfol; i++) {
                valor = nz_p[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
                valor = nz_p[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            linha += "\n";
            // filho1
            for(int i=0; i<this.nfol; i++) {
                valor = nz_f1[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
                valor = nz_f1[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            linha += "\n";
            // filho2
            for(int i=0; i<this.nfol; i++) {
                valor = nz_f2[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
            for(int i=TAM_Z-this.nfol; i<TAM_Z; i++) {
                valor = nz_f2[i];
            	if (valor > nnos-1) valor = valor - dif; 
            	linha = linha + valor + " ";
            }
//            linha += "\n";

		    // Print a line of texto
		    new PrintStream(fout).println (linha);

		    // Close our output stream
		    fout.close();		
		}
		// Catches any error conditions
		catch (IOException e)
		{
			System.err.println ("Unable to write to file");
			System.exit(-1);
		}
	
	}
}