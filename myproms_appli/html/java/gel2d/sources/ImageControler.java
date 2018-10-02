import java.awt.event.*;
import javax.swing.event.*;
import javax.swing.*;
import netscape.javascript.*;// To compile it -> neet to add plugin.jar to your classpath
/**
 * Class reprensting the controler of data that are used by the Applet (MVC structure).
 * This class gets what the ImageApplet sends and updates the ImageModel and the ImageView
 * if it is necessary. It handles the mouse movements.
 * @author Guillaume Arras (Guillaume.Arras@curie.fr)
 * @version 14/05/2009
 * @see ImageApplet
 * @see ImageModel
 */

public class ImageControler implements ActionListener, MouseMotionListener, MouseListener, ListSelectionListener{
	// Parameters
	private ImageView vue;
	private ImageModel modele;
	private JSObject window;// Get all the javascript methods from the applet...
	private boolean editable;
	private boolean isAllowed;
	private int mx;
	private int my;
	private boolean isNotAdding=true;

	// Constructor
	public ImageControler(ImageView v, ImageModel m, boolean b){
		vue = v;
		modele = m;
		editable = b;
		isAllowed=false;
	}
	
	// Methods
	public void setWindow(JSObject window){
		this.window=window;
	}
	
	public void mouseDragged(MouseEvent e){
	}
	
	public void mouseClicked(MouseEvent e){
		// This test avoid to create a spot outside the image...
		if(editable){// Check if the user can create spots...
			if( (e.getClickCount() == 2) 
					&& (e.getX() <= (int)((double)(vue.getImageWidth())*vue.getScale())) 
					&& (e.getY() <= (int)((double)(vue.getImageHeight())*vue.getScale()))){
				// Create a new Spot
				mx=e.getX();
				my=e.getY();
				window.eval("getFreeSamples(\"Create\");");
			}
		}
	}
	/**
	 * Launched where a click is performed --> to get where the click has been done
	 * @param e
	 */
	public void mousePressed(MouseEvent e) {
		//this.setCursor(new Cursor(Cursor.MOVE_CURSOR));
		if(modele.getGel().isnotnull()){
			int x = e.getX();
			int y = e.getY();
			Spot s = modele.getGel().isOnGel(x, y, vue.getScale());
			if (s != null){
				if(s.isSelected()){
					s.deselect();
					window.eval("updateFrame(\"select\",0)");
					modele.removeSelection();
					vue.unselectSpot();
					vue.repaint();
				}else{
					s.select();
					modele.setSelectionedSpot(s.getIdSpot());
					vue.setSpotSelected(s.getIdSpot());
					vue.repaint();
				}
			}
		}
	}
	
	/**
	 * Handles the position of the mouse on the gel
	 */
	public void mouseMoved(MouseEvent e) {

		if(modele.getGel().isnotnull()){// if there is at least one spot defined in the gel
			/*
			 * Get the position of the mouse on the area of the picture
			 */
			int x = e.getX();
			int y = e.getY();
			Spot s = modele.getGel().isOnGel(x, y, vue.getScale());

			if( s != null ){
				/*
				 * A spot area has been found. It needs to be updated so as to print information about the spot
				 * like his isoelectric point,...
				 */
				vue.setTextPopup(s);// Information to print
				vue.showPopup(x,y);// Show the popup generated from the spot itself
			}else{
				vue.setVisiblePopup(false);// Hide the Popup information
			}	
		}
	}
	
	public void mouseEntered(MouseEvent e) {}

	public void mouseExited(MouseEvent e) {}

	public void mouseReleased(MouseEvent e) {}
	
	public void editOption(String mysamp){
		// If no spot selected --> Print an information that no spot has been selected
		if(modele.getSpotSelected() == null){
			vue.showDialogEdit();
		}else{
			// Show the dialog that is going to be used for the edition
			Spot s = modele.getSpotSelected();
						
			if(s.getSampleID() < 0){// No association
				vue.createDialog(s,modele.getlistfiles(mysamp),"NONE",-1);
			}else{
				vue.createDialog(s,modele.getlistfiles(mysamp),s.getSampleName(),s.getSampleID());
			}
		}		
	}
	
	/**
	 * Launch the dialog to create a new spot in the database
	 * @param mysamp String send via JS to get all the samples that are not associated
	 */
	public void createOption(String mysamp){
		vue.createDialog(mx,my,modele.getlistfiles(mysamp));
	}

	public void actionPerformed(ActionEvent e) {
		
		if(!isAllowed){
			isAllowed=true;
		}
		
		if((e.getActionCommand()).compareTo("Center on spot") == 0 ){
			Spot s = modele.getSpotSelected();
			if(s != null){
    			vue.scrollsTo(s.getX(), s.getY());
    			vue.repaint();
    		}else{
    			vue.showDialogEdit();
    		}
		}
		
		//Selection of a zoom factor
		if((e.getActionCommand()).compareTo("comboBoxChanged") == 0 ){
			JComboBox cb = (JComboBox)e.getSource();
	        String factorZoom = (String)cb.getSelectedItem();
	        if( factorZoom.compareTo("Fit to page") == 0){
	        	vue.setScaleFactor(1.00);
	        }else{
	        	String mult = factorZoom.split(" ")[1];
	        	vue.setScaleFactor(Double.parseDouble(mult));
	        }
	        vue.repaint();
		}
		
		if( (e.getActionCommand()).compareTo("send") == 0){
			if(vue.checkDialog()){// if fields of spot are not empty
				Spot ajout = vue.getDialogInfo();
				ajout.setSampleName(vue.getSampleName());
				ajout.setSampleID(vue.getSampleId());
				String action="&spotName="+ajout.getName()+"&x_pos="+ajout.getX()+"&y_pos="+ajout.getY()+"&pi="+ajout.getPi()+"&pw="+ajout.getPw()+"&intensity="+ajout.getIntensity()+"&externalid="+ajout.getExternalid()+"&sampleid="+ajout.getSampleID()+"&sampleName="+ajout.getSampleName();
				window.eval("manageSpot(\""+action+"\",\"createSpot\")");
			}
		}
		
		if( (e.getActionCommand()).compareTo("search") == 0 ){
			if(vue.getSearch().replaceAll(" ", "").isEmpty()){
				vue.setResearchResult(vue.getSearch().replaceAll(" ", ""));
			}else{
				int alias=0,description=0;
	    		if(vue.isASelected()){
	    			alias=1;
	    		}
	    		if(vue.isPSelected()){
	    			description=1;
	    		}
	    		if(alias == 0 && description == 0){
	    			vue.setResearchResult(null);
	    		}else{
	    			window.eval("searchProtein('"+vue.getSearch()+"',"+alias+","+description+")");
	    		}
			}
    	}
		
		if( (e.getActionCommand()).compareTo("Delete spot") == 0 ){
			if(modele.getGel().getSpots() != null){
				if(vue.getSpotSelected() > -1){
					Spot spotSelected = modele.getSpotSelected();
					if( vue.showDialogDelete(spotSelected.getName()) == 0){
						window.eval("deleteSpot("+spotSelected.getIdSpot()+","+spotSelected.getSampleID()+")");
					}
				}else{
					vue.showDialogEdit();
				}
			}else{
				vue.showDialogEdit();
			}
    	}
		
		if( (e.getActionCommand()).compareTo("Edit spot") == 0 ){
			window.eval("getFreeSamples(\"Edit\")");
		}
		
		if( (e.getActionCommand()).compareTo("update") == 0 ){
			
			/* Traiter le cas de l'edition d'un spot avec Javascript*/
			if(vue.checkDialog()){
				Spot older = modele.getSpotSelected();
				Spot update = vue.getDialogInfo();
				String action="&idspot="+older.getIdSpot()+"&spotName="+update.getName()+"&x_pos="+update.getX()+"&y_pos="+update.getY()+"&pi="+update.getPi()+"&pw="+update.getPw()+"&intensity="+update.getIntensity()+"&externalid="+update.getExternalid()+"&sampleid="+vue.getSampleId()+"&sampleName="+vue.getSampleName();
				if(vue.isNewSample()){
					action+="&isnew=1";
				}else{
					action+="&isnew=0";
				}
				window.eval("manageSpot(\""+action+"\",\"editSpot\")");
				window.eval("updateFrame(\"edit\","+older.getIdSpot()+")");
			}
		}
		
	}
	
	/**
	 * Method that allow the manipulation of the JList. That method is called when an
	 * element is selected or deselected
	 */
	public void valueChanged(ListSelectionEvent e){
		System.out.println("valueChanged()");
		if( (!e.getValueIsAdjusting()) && (modele.getGel().getSpots() != null) && (isNotAdding)){// Check if it is a selection event
			int spotSelected = vue.getSpotSelected();
			if(spotSelected > -1){   /* When there is a deselection done by clicking on a spot directly,
									  * it is needed to avoid to do all the following steps
									  */
				Spot s = modele.setSelectionedSpot(spotSelected);
				window.eval("updateFrame(\"select\","+s.getIdSpot()+")");
				//vue.setSpotSelected(spotSelected);
				if(vue.isNotOnVue(s.getX(),s.getY())){
					vue.scrollsTo(s.getX(),s.getY());
				}
				vue.repaint();
			}
		}
	}
	
	public void addSpot(Spot s){
		isNotAdding=false;
		modele.addSpot(s);
		vue.updateSpotList(s.getName(),s.getIdSpot(),true);
		if(isAllowed){
			modele.setSelectionedSpot(s.getIdSpot());
			vue.setSpotSelected(s.getIdSpot());
		}
		vue.repaint();
		isNotAdding=true;
	}
	
	public void updateSpot(Spot s, String action){
		if(action.compareTo("createSpot")==0){
			System.out.println("updateSpot: cas create");
			window.eval("updateFrame(\"create\","+s.getIdSpot()+")");
			addSpot(s);
		}
		if(action.compareTo("editSpot")==0){
			System.out.println("updateSpot: cas edit");
			window.eval("updateFrame(\"update\","+s.getIdSpot()+")");
			String oldname=modele.getSpot(s.getIdSpot()).getName();
			modele.updateSpot(s);
			if(s.getName().compareTo(oldname) != 0){
				vue.updateSpotList(oldname,s.getIdSpot(), false);
				vue.updateSpotList(s.getName(),s.getIdSpot(), true);
				vue.setSpotSelected(s.getIdSpot());
			}
			vue.repaint();
		}
	}
	
	public void displaySearch(String result){
		vue.setResearchResult(result);
	}
	
	public void deleteSpot(int ID_SPOT){
		Spot s = modele.setSelectionedSpot(ID_SPOT);
		modele.removeSelectedSpot(ID_SPOT);
		window.eval("updateFrame(\"delete\","+ID_SPOT+")");
		vue.updateSpotList(s.getName(),ID_SPOT,false);
		vue.repaint();
	}
	
	public void printSearch(String s){
		if(s==null){
			vue.setResearchResult("Protein not found");
		}else{
			if(s.isEmpty()){
				vue.setResearchResult("No string entered");
			}else{
				String[] idspots = s.split("@@");
				String toPrint = "Protein found in:\n";
				modele.removeSelection();
				vue.unselectSpot();
				for(int i = 0 ; i < idspots.length ; i++ ){
					Spot spot = modele.getSpot(Integer.parseInt(idspots[i]));
					spot.hasProtein();
					toPrint = toPrint + spot.getName() + "\n";
				}
				vue.setResearchResult(toPrint);
				vue.repaint();
			}
		}
	}

}
