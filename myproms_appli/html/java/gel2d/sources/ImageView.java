import java.awt.*;
import java.io.IOException;
import java.net.*;
import javax.imageio.ImageIO;
import javax.swing.*;
import netscape.javascript.*;
/**
 * This class is a component of myProMS software and display
 * a gel image in an html page with its annotation whether it
 * was annotated or not by ImageMaster.
 * @author Guillaume Arras (Guillaume.Arras@curie.fr)
 * @version 12/05/2009
 */

public class ImageView extends JApplet {
	// Parameters
	private static final long serialVersionUID = 42L;// just a static parameter
	private ImageModel modele;// MVC -> Model
	private ImageControler controleur;// MVC -> Controler
	private ImagePanel zoneimagep;// Panel to print the image
	private DefaultListModel spotlistmodel;//
	private int[] idSpotList=null;
	private JList spotlist;
	private JScrollPane imageScroller;
	private SpotForm createSpot;
	private String header ="No current spot";
	private JTextField mysearch;
	private JTextArea research_result;
	private JCheckBox alias;
	private JCheckBox description;
	private JComboBox zoomList;
	private boolean editable;
	private JScrollPane listScroller;
	
	public ImageView(){
		super();
	}
	
	// Methods	
	public void init () {
		
		modele = new ImageModel(getParameter("image_file"));
		
		// Is it possible to edit a spot or not ?
		int userPermission = Integer.parseInt(getParameter("editable"));
		if( userPermission == 0){
			editable = false;
		}else{
			editable = true;
		}
		controleur = new ImageControler(this,modele,editable);
		
		// Creation of the menu of the applet
		//mylayout.setVgap(30);//Increase the space between buttons
		/* The left panel*/
		JPanel leftuppanel = new JPanel(new GridLayout(2,1));
		JPanel leftpanel = new JPanel(new GridLayout(2,1));
		JPanel research = new JPanel(new GridLayout(2,1));
		/*The buttons*/
		JPanel zoneBoutton = new JPanel(new GridLayout(5,1));//6 lines - 1 column
		JPanel recherche = new JPanel(new GridLayout(2,1));
		
		JButton bdelete = new JButton("Delete spot");
		bdelete.addActionListener(controleur);
		bdelete.setEnabled(editable);
		
		JButton edit = new JButton("Edit spot");
		edit.addActionListener(controleur);

		JButton center = new JButton("Center on spot");
		center.addActionListener(controleur);
		
		if(!editable){
			bdelete.setSelected(false);
			edit.setSelected(false);
		}
		
		
		// Le JMenu du zoon
		
		String[] t = new String[7];
		t[0] = "Fit to page";
		for(int i = 1 ; i < 7 ; i++){
			t[i] = "X "+(int)(Math.pow(2.00, (double)(i)));
		}
		
		//JPanel zoom = new JPanel(new GridLayout(2,1));
		
		//Zoom comboBox element
		zoomList = new JComboBox(t);
		zoomList.setSelectedIndex(0);
		zoomList.addActionListener(controleur);
		JTextField zoomlab = new JTextField("Select zoom factor");
		zoomlab.setBackground(zoneBoutton.getBackground());
		zoomlab.setFont(zoomlab.getFont().deriveFont(Font.BOLD));
		zoomlab.setEditable(false);
		//zoneBoutton.add(zoom);
		
		zoneBoutton.add(zoomlab);
		zoneBoutton.add(zoomList);
		zoneBoutton.add(center);
		zoneBoutton.add(bdelete);
		zoneBoutton.add(edit);
		
		/*The search area*/
		JLabel enter = new JLabel("Search for protein:");
		JButton search = new JButton("SEARCH");
		search.setActionCommand("search");
		recherche.add(enter);
		search.addActionListener(controleur);
		mysearch = new JTextField(15);
		recherche.add(mysearch);
		
		alias = new JCheckBox("name in myProMS");
		alias.setSelected(true);
		description = new JCheckBox("description");
		description.setSelected(true);
		
		JPanel result = new JPanel(new GridLayout(3,1));
		result.add(alias);
		result.add(description);
		result.add(search);
		
		research.add(recherche);
		research.add(result);
		
		/*The search print*/
		research_result = new JTextArea(5,5);
		research_result.setEditable(false);
		research_result.setAutoscrolls(true);

		/*Finalisation of the left panel*/
		leftuppanel.add(zoneBoutton);
		leftuppanel.add(research);
		
		leftpanel.add(leftuppanel);
		leftpanel.add(new JScrollPane(research_result));
		
		Container conteneur = getContentPane();//get pane to print the elements !
		conteneur.setLayout(new BorderLayout());//set the Layout of the Major Paint

		conteneur.add(leftpanel, BorderLayout.WEST);
		JLabel title = new JLabel("Image loaded: "+modele.getName(),JLabel.CENTER);
		title.setFont(title.getFont().deriveFont(Font.BOLD));
		conteneur.add(title,BorderLayout.NORTH);

		try{
			System.out.println(getParameter("gel2D_file"));
			zoneimagep = new ImagePanel(ImageIO.read(new URL(getParameter("gel2D_file"))),new Dimension(getX(),getY()));
		}catch (IOException ex) {
			System.out.println("Loading of "+modele.getImageName()+" image file didn't succeed");
		}
		imageScroller = new JScrollPane(zoneimagep,
            	JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            	JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		imageScroller.setAutoscrolls(true);
		imageScroller.setPreferredSize(new Dimension(getX(),getY()));
		conteneur.add(imageScroller, BorderLayout.CENTER);
		
		zoneimagep.addMouseListener(controleur);
		zoneimagep.addMouseMotionListener(controleur);
		zoneimagep.setGel(modele.getGel());

		spotlistmodel = new DefaultListModel();
		spotlistmodel.addElement(header);

		spotlist = new JList(spotlistmodel);
		spotlist.setFixedCellWidth(110);
		spotlist.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		spotlist.addListSelectionListener(controleur);
		spotlist.setVisibleRowCount(-1);
		listScroller = new JScrollPane(spotlist,JScrollPane.VERTICAL_SCROLLBAR_ALWAYS,
            	JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		listScroller.getViewport().setView(spotlist);
		listScroller.setAutoscrolls(true);
		
		zoneimagep.setScaleFactor(1.00);
		conteneur.add(listScroller, BorderLayout.EAST);
		setVisible(true);
		JSObject window = JSObject.getWindow(this);
		window.eval("getSpotList()");
		controleur.setWindow(window);
	}
	
	public void setGel(Gel g){
		zoneimagep.setGel(g);
	}
	
	public double getScale(){ return zoneimagep.getScale(); }
	
	public void setTextPopup(Spot s){
		zoneimagep.setTextPopup(s);
	}

	public void setVisiblePopup(boolean b){
		zoneimagep.setVisiblePopup(b);
	}
	
	public void showPopup(int a , int b){
		zoneimagep.showPopup(a,b);
	}
	
	/**
	 * Method to get a Frame to allow the creation of a JDialog (need a frame to be set)
	 * @return
	 */
	public Frame findParentFrame(){
		Component c = getParent();
		
		while(true){
			if(c instanceof Frame)
				return (Frame)c;
			
			c = c.getParent();
			}	
	}
	
	/**
	 * Method that create a dialog when double clicking has occured somewhere on the gel.
	 * @param x
	 * @param y
	 */
	public void createDialog(int x, int y, String[][] listfiles){
		int dx = (int)(x/zoneimagep.getScale());// Update the x coordinate from the mouse
		int dy = (int)(y/zoneimagep.getScale());// Update the y coordinate from the mouse
		createSpot = new SpotForm(findParentFrame(),"Fill the spot information",dx,dy,listfiles);
		createSpot.addActionListener(controleur);
		createSpot.setVisible(true);
	}
	
	public void createDialog(Spot s, String[][] listfiles, String NAME, int ID_SAMPLE){
		createSpot = new SpotForm(findParentFrame(),s,listfiles,NAME,ID_SAMPLE);
		createSpot.addActionListener(controleur);
		createSpot.setVisible(true);
	}
	
	/**
	 * Create a spot with the information entered in the dialog
	 * @return
	 */
	public Spot getDialogInfo(){
		Spot nouveau = createSpot.getSpot();
		createSpot.dispose();
		return nouveau; 
		}
	
	public boolean isNewSample(){
		return createSpot.isNewSample();
	}
	
	public boolean checkDialog(){ return createSpot.checkFormat(); }
	
	/**
	 * Update the JList by adding or removing the spot name from the list
	 * @param s
	 * @param b if true, the name is added, if not, it is removed
	 */
	public void updateSpotList(String s, int ID_SPOT, boolean b){
		if(b){
			//Check if it is the first spot to add to the list
			if(spotlistmodel.contains(header)){
				spotlistmodel.remove(0);
			}
			//Add the name to the list
			if(idSpotList == null){
				idSpotList = new int[1];
				idSpotList[0] = ID_SPOT;
				spotlistmodel.addElement(s);
			}else{
				int pos=-1,tag=0;				
				int len = idSpotList.length;
				int [] ajout = new int[len+1];
				for(int i = 0 ; i < idSpotList.length ; i++){
					/*Test if the new string is lexicographically smaller*/
					if((spotlistmodel.getElementAt(i).toString().compareToIgnoreCase(s) > 0) && (pos<0)){
						spotlistmodel.add(i,s);
						pos=i;
						ajout[tag]=ID_SPOT;
						tag++;
					}
					ajout[tag]=idSpotList[i];
					tag++;
				}
				if(pos<0){
					spotlistmodel.addElement(s);
					ajout[len] = ID_SPOT;
				}
				idSpotList = ajout;
			}
			
			/*//Part that works
			spotlistmodel.addElement(s);
			if(idSpotList != null){
				int len = idSpotList.length;
				int [] ajout = new int[len+1];
				for(int i = 0 ; i < len ; i++){
					ajout[i] = idSpotList[i];
				}
				ajout[len] = ID_SPOT;
				idSpotList = ajout;
			}else{
				idSpotList = new int[1];
				idSpotList[0] = ID_SPOT;
			}
			//spotlist.setSelectedIndex(idSpotList.length-1);
			//setSpotSelected(ID_SPOT);*/
			
		}else{
			if(idSpotList != null && idSpotList.length > 1){
				int [] copie = new int[idSpotList.length-1];
				int indC = 0;
				int posElement=0;
				for(int i = 0 ; i < idSpotList.length ; i++){
					if(idSpotList[i] != ID_SPOT){
						copie[indC] = idSpotList[i];
						indC++;
					}else{
						posElement=i;
					}
				}
				idSpotList = copie;
				spotlistmodel.remove(posElement);
			}else{
				spotlistmodel.removeAllElements();
				spotlistmodel.addElement(header);
				idSpotList=null;
			}
		}
	}
	
	public void updateSpotList(String name, int ID_SPOT){
		for(int i = 0 ; i < idSpotList.length ; i++){
			if(idSpotList[i] == ID_SPOT){
				spotlistmodel.set(i, name);
			}
		}
	}
	
	/**
	 * Get the spot selected in the SpotList
	 * @return ID_SPOT ID_SPOT of the selectioned spot
	 */
	public int getSpotSelected(){
		//System.out.println("ImageView getSpotSelected()");
		if(spotlist.getSelectedIndex() != -1){
			int pos = spotlist.getSelectedIndex();
			//System.out.println("pos="+pos+"\tidSpotList[pos]="+idSpotList[pos]);
			return idSpotList[pos];
		}
		return -1;
	}
	
	public void setSpotSelected(int ID_SPOT){
		//System.out.println("setSpotSelected() ID_SPOT="+ID_SPOT);
		for(int i = 0 ; i < idSpotList.length ; i++){
			//System.out.println("i="+i+"\tidSpotList[i]="+idSpotList[i]);
			if(ID_SPOT == idSpotList[i]){
				spotlist.setSelectedIndex(i);
				break;
			}
		}
		//System.out.println("End setSpotSelected()");
	}
	
	public void unselectSpot(){
		spotlist.clearSelection();
		//spotlist.setSelectedIndices(new int[] {});
	}
	
	public int getImageWidth(){
		return (int)(zoneimagep.getImageWidth());
	}
	
	public int getImageHeight(){
		return (int)(zoneimagep.getImageHeight());
	}
	
	/**
	 * If a spot is selected on the image,
	 * the JViewport is updated so as to be focused on the spot
	 * @param a X coordinate of a spot
	 * @param b Y coordinate of a spot
	 */
	public void scrollsTo(int a, int b){
		int x,y;
		Rectangle visible = zoneimagep.getVisibleRect();// Get the visible rectangle of the JViewport
		int maxY = zoneimagep.getHeight();
		int maxX = zoneimagep.getWidth();
		x = (int)(a*zoneimagep.getScale()) - visible.width/2;//get the x position of the spot according to scale
		y = (int)(b*zoneimagep.getScale()) - visible.height/2;//get the y position of the spot according to scale
		// Checks if the coordinates are inside the picture
		if( x < 0 ){
			x = 0;
		}
		if( y < 0){
			y = 0;
		}
		if( y + visible.height > maxY){
			y = maxY - visible.height;
		}
		if( x + visible.width > maxX ){
			x = maxX - visible.width;
		}
		//System.out.println("ImageView_Scrollsto -> x: "+x+" y: "+y+" a: "+a+" b: "+b);

		imageScroller.getViewport().setViewPosition(new Point(x,y));
		imageScroller.validate();
	}
	
	/**
	 * Inform the user that no spot has been seleceted and that it is not possible to perform
	 * the task wanted (edit the information of the spot)
	 */
	public void showDialogEdit(){
		JOptionPane.showMessageDialog(findParentFrame(), "No spot selected",
				"alert", JOptionPane.ERROR_MESSAGE); 
	}
	
	public int showDialogDelete(String name){
		
		Object[] options = {"Of course","No way!"};
		int n = JOptionPane.showOptionDialog(findParentFrame(),
				"Would you like to remove "+name+" from gel?",
				"Delete confirmation",
				JOptionPane.YES_NO_OPTION,
				JOptionPane.WARNING_MESSAGE,
				null,     //do not use a custom Icon
				options,  //the titles of buttons
				options[1]); //default button title
		return n;
	}
	
	public String getSearch(){ return mysearch.getText(); }
	
	public void setScaleFactor(double d){ zoneimagep.setScaleFactor(d); }
	
	public String getSampleName(){ return createSpot.getSampleName(); }
	
	public int getSampleId(){ return createSpot.getSampleId(); }
	
	public boolean getNoneChecked() { return createSpot.getNoneChecked(); }
	
	public boolean isNotOnVue(int x, int y){
		double d = zoneimagep.getScale();
		Rectangle rec = zoneimagep.getVisibleRect();
		if( x*d > rec.x && x*d < rec.x + rec.width && y*d > rec.y && y*d < rec.y + rec.height ){
			return false;
		}
		return true;
	}
	
	public boolean isASelected(){ return alias.isSelected(); }
	
	public boolean isPSelected(){ return description.isSelected(); }
	
	public void setResearchResult(String s){
		research_result.setText(s);
	}
	
	/* Methods to communicate with */
	
	public void sendToControler(int di){
		int i = zoomList.getSelectedIndex();
		zoomList.setSelectedIndex(i+di);
	}
	
	public void addSpot(int ID_SPOT, String name, int x, int y, double pi, double pw, 
			double intensity, String externalid, int ID_SAMPLE, String sampleName, String proteinMaj){
		Spot s = new Spot(ID_SPOT,name,x,y,pi,pw,intensity,externalid,ID_SAMPLE,sampleName,proteinMaj);
		controleur.addSpot(s);
	}
	
	public void addSpot(String info, String action){
		Spot s = new Spot(info);
		controleur.updateSpot(s,action);
	}
	
	public void sendFreeSamples(String s, String action){
		if(s.compareTo("@#NONE#@")==0){
			s="";
		}
		if(action.compareTo("Edit") == 0){
			controleur.editOption(s);
		}
		if(action.compareTo("Create") == 0){
			controleur.createOption(s);
		}
	}
	
	public void deleteSpot(int ID_SPOT){
		controleur.deleteSpot(ID_SPOT);
	}
	
	public void printSearch(String spotList){
		if(spotList.isEmpty() || spotList.compareTo("@@NONE@@")==0){
			controleur.printSearch(null);
		}else{
			controleur.printSearch(spotList);
		}
	}
	
	public void updateSpot(String spot, String action){
		System.out.println("Spot non céé: spot="+spot+" action="+action);
		Spot s = new Spot(spot);
		System.out.println("Spot céé: spot="+spot+" action="+action);
		controleur.updateSpot(s,action);
	}

}