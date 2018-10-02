import java.awt.*;
/**
 * A Gel of a biological experiment
 * @author Guillaume Arras (Guillaume.Arras@curie.fr)
 * @version 15/05/2009
 */
public class Gel {
	//Parameters
	private Spot [] legel;
	private String name;
	
	//Constructor
	public Gel(){
		legel = null;
		name = null;
	}

	public Gel(String name){
		legel = null;
		this.name = name;
	}
	
	public Gel(String name, Spot[] legel){
		this.name = name;
		this.legel = legel;
	}
	
	//Methods
	/**
	 * Add a spot to the gel
	 * @paral s the spot to add from the list
	 */
	public void addSpot(Spot s){
		
		if(legel != null){
			int len = legel.length;
			Spot [] ajout = new Spot[len+1];
			for(int i = 0 ; i < len ; i++){
				ajout[i] = legel[i];
			}
			ajout[len] = s;
			legel = ajout;
		}else{
			legel = new Spot[1];
			legel[0] = s;
		}
	}
	
	public void removeSpot(int ID_SPOT){

		if( (legel != null) && (legel.length > 1) ){
			
			Spot [] copie = new Spot[legel.length-1];
			int indC = 0;
			for(int i = 0 ; i < legel.length ; i++){
				if(legel[i].getIdSpot() != ID_SPOT){
					copie[indC] = legel[i];
					indC++;
				}
			}

			legel = copie;
			
		}else{
			legel = null;
		}

	}
	
	public void removeSpot(Spot s){
		if(legel != null){
			int i = 0;
			while(i < legel.length){
				if(legel[i].equalsto(s)){
					removeSpot(i);
					break;
				}
				i++;
			}
		}
	}
	
	public void removeSpot(String s){
		if(legel != null){
			int i = 0;
			while(i < legel.length){
				if(legel[i].getName() ==s){
					removeSpot(i);
					break;
				}
				i++;
			}
		}
	}
	
	public String toString(){
		String ret = "";
		if(legel != null){
				ret = "Nombre de spots annotes sur le gel:\t"+legel.length+"\n";
				for(int i = 0 ; i < legel.length ; i++){
				ret = ret+legel[i].toString()+legel[i].isSelected();
			}
		}else{
			ret = "gel non annote\n";
		}
		return ret;
	}
	
	public Spot isOnGel(int dx, int dy, double scale){
		Spot search = null;
		for(int i = 0; i < legel.length ; i++){
			if( legel[i].isAround(dx, dy , scale)){
				search = legel[i];
			}
		}
		return search;
	}
	
	public String getName(){ return name; }
	
	public void drawSpots(Graphics g){
		for(int i = 0 ; i < legel.length ; i++){
			legel[i].paint(g);
		}
	}
	
	public boolean isnotnull(){
		if(legel == null){
			return false;
		}
		return true;
	}
	
	/**
	 * Method which tells the list of all spots that are available on a gel
	 * @return
	 */
	public String [] getSpots(){
		if(legel == null){
			return null;
		}
		int taille = legel.length;
		String [] listeSpot = new String[taille];
		
		for(int i = 0; i < taille ; i++){
			listeSpot[i] = legel[i].getName();
		}
		return listeSpot;
	}
	
	public boolean [] getStatusSpots(){
		if(legel == null){
			return null;
		}
		int taille = legel.length;
		boolean [] statSpot = new boolean[taille];
		for(int i = 0; i < taille ; i++){
			statSpot[i] = legel[i].isSelected();
		}
		return statSpot;
	}
	
	public void removeSelection(){
		for(int i = 0 ; i < legel.length ; i++){
			legel[i].deselect();
		}
	}
	
	public Spot setSelectionSpot(int ID_SPOT){
		int t=0;
		for(int i = 0 ; i < legel.length ; i++){
			legel[i].deselect();
			if( legel[i].getIdSpot() == ID_SPOT){
				legel[i].select();
				t = i;
			}
		}
		return legel[t];
	}
	
	public String searchMajProt(String s){
		String ret = null;
		for(int i = 0 ; i < legel.length ; i++){
			if(legel[i].getProt().compareTo(s) == 0 ){
				if(ret == null){ ret = ""; }
				ret = ret + legel[i].getName() + "\n";
			}
		}
		return ret;
	}
	
	public Spot getSpotSelected(){
		for(int i = 0 ; i < legel.length ; i++){
			if( legel[i].isSelected() ){
				return legel[i];
			}
		}
		return null;
	}
	
	public Spot getSpot(int ID_SPOT){
		for(int i = 0 ; i < legel.length ; i++){
			if( legel[i].getIdSpot() == ID_SPOT){
				return legel[i];
			}
		}
		return null;
	}
	
	public void updateSpot(int ID_SPOT, Spot s){
		for(int i = 0; i < legel.length ; i++){
			if(legel[i].getIdSpot() == ID_SPOT){
				System.out.println("Old Spot="+legel[i]);
				System.out.println("New Spot="+s);
				legel[i] = s;
				legel[i].setIdSpot(ID_SPOT);
				legel[i].select();
			}
		}
	}
	
	public Spot [] getAllSPots(){
		return legel;
	}
	
	public void setName(String s){ name = s; }
	
	public static void main(String [] args){
		Spot test = new Spot();
		Gel testg = new Gel();
		System.out.println(testg);
		testg.addSpot(test);
		System.out.println(testg);
		testg.addSpot(test);
		System.out.println(testg);
		testg.removeSpot(test);
		System.out.println(test.equalsto(new Spot(0,0,"spot Fictif",45.0008,75.9865)));
		System.out.println(testg);
		testg.removeSpot(test);
		System.out.println(testg);
	}
}
