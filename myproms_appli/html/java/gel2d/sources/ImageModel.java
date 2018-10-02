/**
 * Class reprensting the model of data that are used by the Applet.
 * This class handle the connection with the myProMS database to
 * print correctly the informations on each gels
 * @author Guillaume Arras (Guillaume.Arras@curie.fr)
 * @version 14/05/2009
 */
public class ImageModel{
	// Parameters
	private Gel mongel;
	
	// Constructor
	public ImageModel(String gelname){
		mongel = new Gel(gelname);
	}
	
	// Methods
	public String getImageName(){
		return mongel.getName();	
	}
	
	public void addSpot(Spot s){
		mongel.addSpot(s);
	}
	
	public Gel getGel(){ return mongel; }
	
	public void removeSelection(){
		mongel.removeSelection();
	}
	
	public Spot setSelectionedSpot(int ID_SPOT){
		return mongel.setSelectionSpot(ID_SPOT);
	}
	
	public void removeSelectedSpot(int ID_SPOT){
		mongel.removeSpot(ID_SPOT);
	}
	
	public String searchMajProt(String s){
		return mongel.searchMajProt(s);
	}
	
	public Spot getSpotSelected(){
		return mongel.getSpotSelected();
	}
	
	public Spot getSpot(int ID_SPOT){
		return mongel.getSpot(ID_SPOT);
	}
	
	public String [] getSpotNames(){
		return mongel.getSpots();
	}
	
	/**
	 * Get the list of samples not already associated to a spot so as to give a choice to
	 * the user.
	 * A String is sent to java via AJAX procedure. The samples are separated with the
	 * '##' and the information of each sample (ie sample name and sample id) is seperated
	 * by '@@'
	 * @return
	 */
	public String[][] getlistfiles(String s){
		String[] firstSplit = s.split("@@");
		int nSamp = firstSplit.length;
		String [][] all;
		
		/* No sample was found not associated */
		if(nSamp < 2){
			all = new String[2][2];
			all[0][0] = "New Sample";
			all[0][1] = null;
			all[1][0] = "NONE";
			all[1][1] = null;
		}else{
			/* Just one sample has been found not associated to a spot */
			if( nSamp == 2){
				
				all = new String[3][2];
				all[0][0] = "New Sample";
				all[0][1] = null;
				all[1][0] = firstSplit[1];
				all[1][1] = firstSplit[0];
				all[2][0] = "NONE";
				all[2][1] = null;
			}else{
				String[] fs = s.split("##");
				all = new String[fs.length+2][2];
				all[0][0] = "New Sample";
				all[0][1] = null;
				for(int i = 0 ; i < fs.length ; i++){
					String[] description = fs[i].split("@@");
					all[i+1][0] = description[1];
					all[i+1][1] = description[0];
				}
				all[fs.length+1][0] = "NONE";
				all[fs.length+1][1] = null;
			}
		}
		return all;
	}
	
	public String getName(){ return mongel.getName(); }
	
	public void updateSpot(Spot newer){
		mongel.updateSpot(newer.getIdSpot(), newer);
	}

}
