import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
/**
 * Class that handles the edit of a spot. Send information to the controler
 * that creates a new Spot intense.
 * @author Guillaume ARRAS (Guillaume.Arras@curie.fr)
 */

public class SpotForm extends JDialog{
	// Parameters
	private static final long serialVersionUID = 42L;// just a static parameter
	private JLabel pil;
	private JLabel pwl;
	private JLabel namel;
	private JTextField namet;
	private JTextField pit;
	private JTextField pwt;
	private JLabel pwll;
	private JLabel intensl;
	private JTextField intenst;
	private JLabel externl;
	private JTextField externt;
	private JButton send;
	private JButton cancel;
	private int x;
	private int y;
	private String myActionCommand;
	private JComboBox listefiles;
	private String[][] sampleNameandIDs;
	private boolean newchecked;
	private String chosen;
	private boolean nonechecked;
	
	// Constructor
	/**
	 * For the creation of a new spot
	 */
	public SpotForm(Frame owner, String s, int x, int y, String[][] listf){
		super(owner,s);
		newchecked = false;
		nonechecked = false;
		setLayout(new GridLayout(7,2));// 7 rows - 2 columns
		namel = new JLabel("Spot name:",JLabel.RIGHT);
		namet = new JTextField(15);
		sampleNameandIDs = listf;
		listefiles = new JComboBox(getListNames());
		pil = new JLabel("Isoelectric point (pI):",JLabel.RIGHT);
		pit = new JTextField(10);
		pwl = new JLabel("Molecular weight:",JLabel.RIGHT);
		pwll = new JLabel("kDa",JLabel.LEFT);
		pwt = new JTextField(10);
		JPanel weight=new JPanel(new GridLayout(1,2));
		weight.add(pwt);
		weight.add(pwll);
		externl = new JLabel("External identifier:",JLabel.RIGHT);
		intensl = new JLabel("Intensity:",JLabel.RIGHT);
		externt = new JTextField(10);
		intenst = new JTextField(10);
		send = new JButton("Add spot");
		cancel = new JButton("Cancel");
		cancel.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				dispose();
			}
		});
		myActionCommand = "send";
		this.x = x;
		this.y = y;
		
		System.out.println("MaximumRowCount="+getListNames().length);
		listefiles.setSelectedIndex(getListNames().length - 1);
		
		add(namel);
		add(namet);
		add(pil);
		add(pit);
		add(pwl);
		add(weight);
		add(intensl);
		add(intenst);
		add(externl);
		add(externt);
		add(new JLabel("Sample",JLabel.RIGHT));
		add(listefiles);
		add(send);
		add(cancel);
		
		Dimension maxS = pil.getPreferredSize();
		setPreferredSize(new Dimension(maxS.width*2+10,maxS.height*10));
		setResizable(false);
		pack();
	}
	
	/**
	 * For the edition of a spot
	 * @param owner
	 * @param s
	 */
	public SpotForm(Frame owner, Spot s, String[][] listf, String NAME, int ID_SAMPLE){
		super(owner,"Modify spot information");
		newchecked = false;
		nonechecked = false;
		setLayout(new GridLayout(7,2));// 6 rows - 2 columns
		namel = new JLabel("Spot name:",JLabel.RIGHT);
		namet = new JTextField(15);
		namet.setText(s.getName());
		sampleNameandIDs = listf;
		listefiles = new JComboBox(getListNames());
		if(NAME.compareTo("NONE") != 0){
			addListName(NAME,ID_SAMPLE);
		}
		if(NAME.compareTo("NONE") == 0){
			int n = listf.length-1;
			listefiles.setSelectedIndex(n);
		}
		pil = new JLabel("Isoelectric point (pI):",JLabel.RIGHT);
		pit = new JTextField(10);
		pit.setText(s.getPi());
		pwl = new JLabel("Molecular weight:",JLabel.RIGHT);
		pwll = new JLabel("kDa",JLabel.LEFT);
		pwt = new JTextField(10);
		pwt.setText(s.getPw());
		JPanel weight=new JPanel(new GridLayout(1,2));
		weight.add(pwt);
		weight.add(pwll);
		externl = new JLabel("External identifier:",JLabel.RIGHT);
		intensl = new JLabel("Intensity:",JLabel.RIGHT);
		externt = new JTextField(10);
		externt.setText(s.getExternalid());
		intenst = new JTextField(10);
		intenst.setText(s.getIntensity());
		cancel = new JButton("Cancel");
		cancel.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e){
				dispose();
			}
		});

		send = new JButton("Save");
		myActionCommand = "update";
			
		x = s.getX();
		y = s.getY();
		
		add(namel);
		add(namet);
		add(pil);
		add(pit);
		add(pwl);
		add(weight);
		add(intensl);
		add(intenst);
		add(externl);
		add(externt);
		add(new JLabel("Sample:",JLabel.RIGHT));
		add(listefiles);//Update to put the name of the actual spot
		add(send);
		add(cancel);
		
		Dimension maxS = pil.getPreferredSize();
		setPreferredSize(new Dimension(maxS.width*2+10,maxS.height*10));
		setResizable(false);
		pack();
	}
	
	// Methods
	public void addActionListener(ActionListener a){
		send.setActionCommand(myActionCommand);
		send.addActionListener(a);
	}
	
	public String toString(){ return "test toString de SpotForm"; }
	
	public boolean checkFormat(){
		/*if(pit.getText().isEmpty() && pwt.getText().isEmpty()){
			return false;
		}*/
		return true;
	}
	
	public Spot getSpot(){
		String value = (String)listefiles.getSelectedItem();
		double pi;
		double pw;
		if(pit.getText().isEmpty()){
			pi=-1;
		}else{
			pi=Double.parseDouble(pit.getText().replace(',','.'));
		}
		if(pwt.getText().isEmpty()){
			pw=-1;
		}else{
			pw=Double.parseDouble(pwt.getText().replace(',','.'));
		}
		
		if( value.compareTo("NONE") == 0 ){
			nonechecked = true;
		}
		if( value.compareTo("New Sample") == 0 ){
			newchecked = true;
			chosen = JOptionPane.showInputDialog("Please type new Sample name:");
			if(chosen == null || chosen.compareTo("") == 0){
				newchecked = false;
				nonechecked = true;
				chosen = "NONE";
			}
		}
		if(intenst.getText().compareTo("")==0){
			return new Spot(x,y,namet.getText(),
					pi,
					pw,
					externt.getText()
					);
		}
		
		return new Spot(x,y,namet.getText(),
				pi,
				pw,
				Double.parseDouble(intenst.getText().replace(',','.')),
				externt.getText()
				);
	}
	
	public boolean isNewSample(){ return newchecked; }
	
	private String[] getListNames(){
		String[][] liste = sampleNameandIDs;
		String[] s = new String[liste.length];
		for(int i = 0 ; i < liste.length ; i++){
			s[i] = liste[i][0];
		}
		return s;
	}
	
	public int getListIdSample(int pos){
		return Integer.parseInt(sampleNameandIDs[pos][1]);
	}
	
	public String getSampleName(){
		if(chosen==null){
			int pos = listefiles.getSelectedIndex();
			return sampleNameandIDs[pos][0];
		}
		return chosen;
	}
	
	public int getSampleId(){
		int pos = listefiles.getSelectedIndex();
		if(sampleNameandIDs[pos][1] == null){
			return -1;
		}
		return Integer.parseInt(sampleNameandIDs[pos][1]);
	}
	
	private void addListName(String name, int id){
		int taille = sampleNameandIDs.length;
		String[][] repro = new String[taille+1][2];
		repro[0][0] = sampleNameandIDs[0][0];
		repro[0][1] = sampleNameandIDs[0][1];
		repro[1][0] = name;
		repro[1][1] = ""+id;
		for(int i = 1 ; i < taille ; i++){
			repro[i+1][0] = sampleNameandIDs[i][0];
			repro[i+1][1] = sampleNameandIDs[i][1];
		}
		sampleNameandIDs = repro;
		//System.out.println("SpotForm.addListName : "+name);
		listefiles.insertItemAt(name, 1);
		listefiles.setSelectedIndex(1);
	}
	
	public boolean getNoneChecked(){ return nonechecked; }

}
