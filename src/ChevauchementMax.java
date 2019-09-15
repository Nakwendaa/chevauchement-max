import java.lang.Math;

public class ChevauchementMax {
	
	// Variable de score match, mismatch et indel pour un alignement de chevauchement maximal. Utilisés par toutes les méthodes.
	private static int match, mismatch, indel;
	
	// Imprimer une table de programmation dynamique
	private static void printMatrix(int[][] tableProgDyn) {
		for(int i = 0; i < tableProgDyn.length; i++) {
			for(int j = 0; j < tableProgDyn[i].length; j++) {
				if (j + 1 < tableProgDyn[i].length && tableProgDyn[i][j+1] < 0 && tableProgDyn[i][j+1] > - 10) {
					System.out.print(Integer.toString(tableProgDyn[i][j]) + "   ");
				}
				else if (j + 1 < tableProgDyn[i].length && tableProgDyn[i][j+1] <= -10) {
					System.out.print(Integer.toString(tableProgDyn[i][j]) + "  ");
				}
				else if (j + 1 < tableProgDyn[i].length && tableProgDyn[i][j+1] >= 10) {
					System.out.print(Integer.toString(tableProgDyn[i][j]) + "   ");
				}
				else {
					System.out.print(Integer.toString(tableProgDyn[i][j]) + "    ");
				}
			}
			System.out.println("");
		}
	}
	
	// On veut ici calculer à partir de deux séquences x et y la table de programmation dynamique associé pour leur chevauchement et les
	// scores match, mismatch et indel.
	private static int[][] getMatriceChevauchement(String x, String y) {
		
		int score1, score2, score3;	
		int[][] matrix = new int[x.length() + 1][y.length() + 1]; // On initialise la table de programmation dynamique et toutes ses cases sont à 0.
																// V(0,j) et V(i,0), pour tout i et j, sont donc initialisés à 0.
		for(int i = 1; i < matrix.length; i++) { 
			for(int j = 1; j < matrix[0].length; j++) { // Cette boucle permet alors de remplir la table de programmation dynamique ligne par ligne
													  // en utilisant la formule de récurrence donnée à la question 3.
				score1 = matrix[i - 1][j] + indel;
				score2 = matrix[i][j - 1] + indel;
				if (x.charAt(i - 1) == y.charAt(j - 1)) {
					score3 = matrix[i - 1][j - 1] + match;
				}
				else {
					score3 = matrix[i - 1][j - 1] + mismatch;
				}
				matrix[i][j] = Math.max(score1, Math.max(score2, score3));
			}
		}
		return matrix;
	}
	
	// On veut ici calculer à partir de deux séquences x et y la table de programmation dynamique associé pour le meilleur alignement
	// d'un facteur de x avec y en entier.
	private static int[][] getMatriceFacteur(String x, String y) {
		
		int score1, score2, score3;
		int[][] matrix = new int[x.length() + 1][y.length() + 1]; // On initialise la table de programmation dynamique et toutes ses cases sont à 0.
																// D(0,j) et D(i,0), pour tout i et j, sont donc initialisés à 0.
		for (int j = 0; j < matrix[0].length; j++) { // On met la ligne i = 0 à jour avec notre nouvelle condition initiale.
			matrix[0][j] = j * indel;
		}
		for(int i = 1; i < matrix.length; i++) { 
			for(int j = 1; j < matrix[0].length; j++) { // Cette boucle permet alors de remplir la table de programmation dynamique ligne par ligne
													  // en utilisant la formule de récurrence donnée à la question 3.
				score1 = matrix[i - 1][j] + indel;
				score2 = matrix[i][j - 1] + indel;
				if (x.charAt(i - 1) == y.charAt(j - 1)) {
					score3 = matrix[i - 1][j - 1] + match;
				}
				else {
					score3 = matrix[i - 1][j - 1] + mismatch;
				}
				matrix[i][j] = Math.max(score1, Math.max(score2, score3));
			}
		}
		return matrix;
	}
	
	// On veut obtenir la séquence de la protéine codée par le gene X

	
	// On retrouve la case optimale dans une matrice dans le cas d'un chevauchement.
	private static int[] caseOptimaleChevauchement(int[][] matrix) {
		
		int derniereLigne = matrix.length - 1;
		int derniereColonne = matrix[0].length - 1;
		
		int max1 = matrix[derniereLigne][0];
		int max2 = matrix[0][derniereColonne];
		
		int[] index1 = new int[2];
		int[] index2 = new int[2];
		
		for (int i = 1; i < matrix[derniereLigne].length; i++) {
			if (max1 < matrix[derniereLigne][i]) {
				max1 = matrix[derniereLigne][i];
				index1[0] = derniereLigne;
				index1[1] = i;
			}
		}
		
		for (int j = 1; j < matrix.length; j++) {
			if (max2 < matrix[j][derniereColonne]) {
				max2 = matrix[j][derniereColonne];
				index2[0] = j;
				index2[1] = derniereColonne;
			}
		}
		
		if (max1 > max2) {
			return index1;
		}
		else {
			return index2;
		}
	}

	// On retrouve la case optimale dans une matrice dans le cas d'un chevauchement.
	private static int[] caseOptimaleFacteur(int[][] matrix) {
		
		int derniereColonne = matrix[0].length - 1; // Longueur de Y
		int max = matrix[0][derniereColonne]; // Valeur à la case ligne 0, derniere colonne
		int[] index = new int[2]; // Tableau qui contient (i,j) la case où se trouve la solution optimale
		
		for (int i = 1; i < matrix.length; i++) {
			if (max < matrix[i][derniereColonne]) {
				max = matrix[i][derniereColonne];
				index[0] = i;
			}
		}
		index[1] = derniereColonne;
		return index;
	}

	
	// Retrouver le chevauchement maximal d'une chaîne x et y
	private static String[] algoChevauchementMax(String x, String y) {
		int[][] matrix = getMatriceChevauchement(x, y);
		int[] index = caseOptimaleChevauchement(matrix);
		int i = index[0];
		int j = index[1];

		String seq1 = x.substring(i, x.length()); // La case optimale (i,j) nous donne l'information de où s'arrête le chevauchement, et donc l'alignement.
		String seq2 = y.substring(j, y.length()); // On récupère donc le reste de chaque chaîne pour mettre en forme l'alignement.
		
		// On ajoute les gaps à la fin de la séquence ayant la plus petite longueur. En effet, cela signifie que son suffixe a été aligné
		// et donc que cette séquence a besoin d'indels
		while(seq1.length() > seq2.length()) {
			seq2 = seq2.concat("-");
		}
		while(seq1.length() < seq2.length()) {
			seq1 = seq1.concat("-");
		}
		
		// On remonte dans la matrice/table jusqu'à atteidre une ligne i = 0 ou une colonne j = 0
		while (i > 0 && j > 0) {
			if ((matrix[i][j] == matrix[i - 1][j - 1] + match) && x.substring(i - 1, i).equals(y.substring(j - 1, j))
				|| matrix[i][j] == matrix[i - 1][j - 1] + mismatch && !x.substring(i - 1, i).equals(y.substring(j - 1, j))) {
				seq1 = x.substring(i - 1, i).concat(seq1);
				seq2 = y.substring(j - 1, j).concat(seq2);
				i -= 1;
				j -= 1;
			}
			
			else if (matrix[i][j] == matrix[i][j - 1] + indel) {
				seq1 = "-".concat(seq1);
				seq2 = y.substring(j - 1, j).concat(seq2);
				j -= 1;
			}
			else if (matrix[i][j] == matrix[i - 1][j] + indel) {
				seq1 = x.substring(i - 1, i).concat(seq1);
				seq2 = "-".concat(seq2);
				i -= 1;
			}
		}
		
		// On ajoute les caractères non-insérés de la séquence qui n'a pas été alignée.
		seq1 = x.substring(0, i).concat(seq1);
		seq2 = y.substring(0, j).concat(seq2);
		
		// Puis on ajoute les gaps au début de la séquence ayant la plus petite longueur 		
		while(seq1.length() > seq2.length()) {
			seq2 = "-".concat(seq2);
		}
		while(seq1.length() < seq2.length()) {
			seq1 = "-".concat(seq1);
		}
		
		return new String[]{seq1, seq2, Integer.toString(matrix[index[0]][index[1]])};	
	}

	
	// Retrouver le facteur de x qui a un alignement maximal avec y en entier. Retourune un tableau de string contenant les chaines x et y alignées et le score de l'alignement.
	private static String[] algoFacteurMax(String x, String y) {

		int[][] matrix = getMatriceFacteur(x, y);
		int[] index = caseOptimaleFacteur(matrix);
		int i = index[0];
		int j = index[1];
		
		String seq1 = x.substring(i, x.length()); // Séquence de x restante non alignée
		String seq2 = "";
		
		while (seq2.length() < seq1.length()) { // On a aligné Y sur son entièreté, donc le seul cas possible où on aura besoin de tracer des indels
			seq2 = seq2.concat("-");			// sera à la fin de la chaîne de Y.
		}
		
		// On remonte jusqu'à une case optimale qui se trouve sur la ligne j = 0;
		
		while (j > 0) {
			if (i != 0 && ((matrix[i][j] == matrix[i - 1][j - 1] + match) && x.substring(i - 1, i).equals(y.substring(j - 1, j))
					|| matrix[i][j] == matrix[i - 1][j - 1] + mismatch && !x.substring(i - 1, i).equals(y.substring(j - 1, j)))) {
				seq1 = x.substring(i - 1, i).concat(seq1);
				seq2 = y.substring(j - 1, j).concat(seq2);
				i -= 1;
				j -= 1;
			}
			
			else if (matrix[i][j] == matrix[i][j - 1] + indel) {
				seq1 = "-".concat(seq1);
				seq2 = y.substring(j - 1, j).concat(seq2);
				j -= 1;
			}
			else if (i != 0 && matrix[i][j] == matrix[i - 1][j] + indel) {
				seq1 = x.substring(i - 1, i).concat(seq1);
				seq2 = "-".concat(seq2);
				i -= 1;
			}
		}
		
		seq1 = x.substring(0, i).concat(seq1);
		while(seq1.length() > seq2.length()) { // On rajoute des indels seulement Y puisque Y a été aligné en entier de manière optimale.
			seq2 = "-".concat(seq2);			
		}
		
		return new String[]{seq1, seq2, Integer.toString(matrix[index[0]][index[1]])};
	}
	
	
	// Retourne la séquence de la protéine codée par le gene X
	private static String getProteine(String x) {
		
		int xLen = x.length();
		String proteine, seqActuel, meilleureSeq; // proteine contient la séquence de proteine du geneX pour un frame donné.
												 // seqActuel contient la séquence de protéine actuellement décompté pour un frame donné.
												 // meilleureSeq contiendra la meilleure séquence de protéine (c-à-d la plus longue) pour un frame donné.
												 //
		String[] proteinesCandidates = new String[3];
		int check = 1; // Variable de check qui permet de contrôler le début d'une séquence ou sa fin. Si check = 1, c'est qu'il n'y a pas
						// de séquence de protéine en cours de comptage;
		
		for (int j = 0; j < 3; j++) { // Boucle permettant de répéter l'opération d'analyse pour 3 frames distincts.
			proteine = "";
			seqActuel = "";
			meilleureSeq = "";
			check = 1;
			for (int i = 0; i < (xLen - j)/ 3; i++) { // Boucle qui permet d'analyser un frame donné: traduire une séquence de bases en une séquence d'acides aminés,
				String prot = "";					  // et on en déduit la meilleur séquence de protéine
													  // prot contiendra la traduction d'un triplet en acide aminé.
				switch(x.substring(i*3 + j , i*3 + 3 + j)) {
					case "ATT": case "ATC": case "ATA": case "AUU": case "AUC": case "AUA":
						prot = "I"; break;
					case "CTT": case "CTC": case "CTA": case "CTG": case "TTA": case "TTG":
					case "CUU": case "CUC": case "CUA": case "CUG": case "UUA": case "UUG":
						prot = "L"; break;
					case "GTT": case "GTC": case "GTA": case "GTG":
					case "GUU": case "GUC": case "GUA": case "GUG":
						prot = "V"; break;
					case "TTT": case "TTC": case "UUU": case "UUC":
						prot = "F"; break;
					case "ATG": case "AUG":
						prot = "M"; break;
					case "TGT": case "TGC": case "UGU": case "UGC":
						prot = "C"; break;
					case "GCT": case "GCC": case "GCA": case "GCG": case "GCU":
						prot = "A"; break;
					case "GGT": case "GGC": case "GGA": case "GGG": case "GGU":
						prot = "G"; break;
					case "CCT": case "CCC": case "CCA": case "CCG":  case "CCU":
						prot = "P"; break;
					case "ACT": case "ACC": case "ACA": case "ACG": case "ACU":
						prot = "T"; break;
					case "TCT": case "TCC": case "TCA": case "TCG": case "AGT": case "AGC":  
					case "UCU": case "UCC": case "UCA": case "UCG": case "AGU":
						prot = "S"; break;
					case "TAT": case "TAC":  case "UAU": case "UAC":
						prot = "Y"; break;
					case "TGG": case "UGG":
						prot = "W"; break;
					case "CAA": case "CAG":
						prot = "Q"; break;
					case "AAT": case "AAC": case "AAU":
						prot = "N"; break;
					case "CAT": case "CAC": case "CAU":
						prot = "H"; break;
					case "GAA": case "GAG":
						prot = "E"; break;
					case "GAT": case "GAC": case "GAU":
						prot = "D"; break;
					case "AAA": case "AAG":
						prot = "K"; break;
					case "CGT": case "CGC": case "CGA": case "CGG": case "AGA": case "AGG":
						prot = "R"; break;
					case "TAA": case "TAG": case "TGA": case "UAA": case "UAG": case "UGA":
						prot = "<STOP>"; break;
					default:
						prot = "ERREUR"; break;
				}
				if (check == 1 && prot.equals("M")) {  // Si on voit un M, début de séquence protéinique
					seqActuel = "M";
					check = 0;
				}
				else if (check == 0 && prot.equals("<STOP>")) { // Si on décompte une séquence protéinique, mais qu'un codon STOP apparait
					if (seqActuel.length() > meilleureSeq.length()) { // Si notre séquence de protéine décomptée jusqu'au codon stop est plus grande que les
						meilleureSeq = seqActuel;					  // séquences de protéine décomptées jusque là....
					}
					seqActuel = "";
					check = 1; // Fin de séquence de protéine dans le frame
				}
				else if (check == 0) {
					seqActuel = seqActuel.concat(prot); // Sinon si on est actuellement dans une séquence protéinique, on rajoute l'acide aminé à la séquence de protéine.
				}

				proteine = proteine.concat( "" + prot);
			}
			// FIN de boucle (for; i < (xLen - j)/ 3) qui correspond à l'analyse d'un frame du geneX.
			
			if (check == 0 && seqActuel.length() > meilleureSeq.length()) { // Un début de séquence protéinique existe jusqu'à la fin du gène X sans codon stop.
				meilleureSeq = seqActuel;
			}
			
			// On imprime les informations
			
			if (j == 0) {
				System.out.println("Frame 1: ");
			}
			else if (j == 1) {
				System.out.println("Frame 2: ");
			}
			else {
				System.out.println("Frame 3: ");
			}
			proteinesCandidates[j] = meilleureSeq;
			System.out.println(proteine);
			System.out.println("Séquence protéinique candidate: " + meilleureSeq);
			System.out.println();
		}
		// FIN de boucle (for; j < 3)
		
		// On retourne la plus longue séquence protéinique parmi les meilleures séquences protéiniques de chaque frame.
		
		int max = proteinesCandidates[0].length();
		int index = 0;
		for (int i = 1; i < proteinesCandidates.length; i++) {
			if (max < proteinesCandidates[i].length()) {
				index = i;
				max = proteinesCandidates[i].length();
			}
		}
		return proteinesCandidates[index];
	}
	
	
	
	public static void main(String[] args) {
		
		// Score de similarité pour la formule de récurrence. Valable pour tous les algorithmes utilisés.
		match = 1;
		mismatch = -1;
		indel = -2;
		
		// EXERCICE 1
		System.out.println("########################################################################################################################################################################");
		System.out.println("EXERCICE 1:\n");
		
		String x = "GATACGTCACGTGCACGG";
		String y = "ACGCATT";
		int[][] tableProgDyn = getMatriceChevauchement(x, y);
		System.out.println("Table de programmation dynamique de X et Y pour le chevauchement maximal:\n");
		printMatrix(tableProgDyn);
		
		System.out.println();
		String [] resultats = algoChevauchementMax(x,y);
		System.out.println("X: ".concat(resultats[0]).concat("\nY: ").concat(resultats[1]).concat("\nScore: ").concat(resultats[2]).concat("\n"));
		
		
		
		// EXERCICE 2
		System.out.println("########################################################################################################################################################################");
		System.out.println("EXERCICE 2:\n");
		String[] frags = {"AACTCTCTACTGCTTTCCCC", "CTACTGCTTTCCCCGCCGGA", "CTTTCCCCGCCGGAACCTTCAC", "TAAATTACAACTCTCTACTA" };
		for (int i = 0; i < frags.length; i++) {
			for (int j = i + 1; j < frags.length; j++) {
				resultats = algoChevauchementMax(frags[i], frags[j]);
				System.out.println("S".concat(Integer.toString(i+1)).concat(": ").concat(resultats[0]).concat("\nS").concat(Integer.toString(j+1)).concat(": ").concat(resultats[1]).concat("\nScore: ").concat(resultats[2]).concat("\n"));
			}
		}
		
		
		
		// EXERCICE 3
		System.out.println("########################################################################################################################################################################");
		System.out.println("EXERCICE 3:\n");
		
		String z = "CTCTCCTACAGAGCTTAAATTACAACTCTCTACTGCTTTCCCCGCCGGAACCTTCACACCAGTCACACGT" +
			"ATGTCTCAGAGCAACCGGGAGCTGGTGGTCGACTTTCTCTCCTACAAGCTTTCCCAGAAAGGATACAGCT" +
			"GGAGTCAGTTTAGTGATGTCGAAGAGAATAGGACTGAGGCCCCAGAAGAAACTGAAGCAGAGAGGGAGAC" +
			"CCCCAGTGCCATCAATGGCAACCCATCCTGGCACCTGGCGGATAGCCCGGCCGTGAATGGAGCCACTGGC" +
			"CACAGCAGCAGTTTGGATGCGCGGGAGGTGATTCCCATGGCAGCAGTGAAGCAAGCGCTGAGAGAGGCAG" +
			"GCGATGAGTTTGAACTGCGGTACCGGAGAGCGTTCAGTGATCTAACATCCCAGCTTCACATAACCCCAGG" +
			"GACCGCGTATCAGAGCTTTGAGCAGGTAGTGAATGAACTCTTTCGGGATGGAGTAAACTGGGGTCGCATC" +
			"GTGGCCTTTTTCTCCTTTGGCGGGGCACTGTGCGTGGAAAGCGTAGACAAGGAGATGCAGGTATTGGTGA" +
			"GTCGGATTGCAAGTTGGATGGCCACCTATCTGAATGACCACCTAGAGCCTTGGATCCAGGAGAACGGCGG" +
			"CTGGGACACTTTTGTGGATCTCTACGGGAACAATGCAGCAGCCGAGAGCCGGAAAGGCCAGGAGCGCTTC" +
			"AACCGCTGGTTCCTGACGGGCATGACTGTGGCTGGTGTGGTTCTGCTGGGCTCACTCTTCAGTCGGAAGU" +
			"AGGTCGTTGGTATGGTATGAGTGTAAGTAAGAAAAAAAAAAAAAAAAAAAAAA";
		
		String proteine = getProteine(z);
		System.out.println("Séquence protéinique probable:");
		System.out.println(proteine);
		System.out.println("Longueur de la protéine X: " + proteine.length());
		System.out.println();
		
		
		
		
		// EXERCICE 4
		System.out.println("########################################################################################################################################################################");
		System.out.println("EXERCICE 4:\n");
		System.out.println("Table de programmation dynamique de X et Y pour le meilleur alignement entre un facteur de X et Y:");
		System.out.println();
		resultats = algoFacteurMax(x,y);
		tableProgDyn = getMatriceFacteur(x,y);
		printMatrix(tableProgDyn);
		System.out.println();
		System.out.println("Alignement:");
		System.out.println("X: ".concat(resultats[0]).concat("\nY: ").concat(resultats[1]).concat("\nScore: ").concat(resultats[2]).concat("\n"));
		System.out.println("X' = ACGTCACG" );
	}

}
