import javax.swing.*;
import java.awt.*;
import java.util.*;
import java.awt.image.*;
import java.awt.Graphics2D;

/**
 * Class that handle the graphic composant of the Panel
 * @author Guillaume ARRAS (Guillaume.Arras@curie.fr)
 *
 */
public class ImagePanel extends JPanel implements Observer{
	// Parameters
	private static final long serialVersionUID = 42L;
	private BufferedImage image;
	private double scale=1.0;
	private double scalefactor;
	private Dimension size;
	private JPopupMenu myPopup;
	private JLabel nameS;
	private JLabel piS;
	private JLabel pwS;
	private JLabel exID;
	private JLabel iten;
	private JLabel protMaj;
	private Gel legel;
	
	// Constructor
	public ImagePanel(BufferedImage image, Dimension size){
		super();
		this.size = size;
		setLayout(new FlowLayout(FlowLayout.CENTER));
		setBackground(Color.GRAY);
		setMinimumSize(this.size);
		setMaximumSize(new Dimension(2*image.getHeight(),2*image.getWidth()));
		this.image = image;
		//initScale(image.getWidth(),image.getHeight());
		//addMouseListener(this);
		//addMouseMotionListener(this);
		myPopup = new JPopupMenu();
		nameS = new JLabel();
		piS = new JLabel();
		pwS = new JLabel();
		iten = new JLabel();
		exID = new JLabel();
		protMaj = new JLabel();
		myPopup.add(nameS);
		myPopup.add(piS);
		myPopup.add(pwS);
		myPopup.add(iten);
		myPopup.add(exID);
		myPopup.add(protMaj);
		legel = null;
		}
	
	// Methods
	public void update(Observable obs, Object obj) {}
	
	public void paintComponent(Graphics g){
		super.paintComponents(g);
		Rectangle test = this.getVisibleRect();
		this.initScale(test);
		// Let the scroll pane know to update itself and its scrollbars.
		setPreferredSize(new Dimension((int)((double)image.getWidth()*scale),
				(int)((double)image.getHeight()*scale)));
		// Paint in gray if it is necessary
		//doLayout();
		revalidate();
		//System.out.println("Debut2\n"+this.getHeight()+"\n"+this.getWidth()+"\n"+this.getX()+
		//		"\n"+this.getY()+"\n"+this.getPreferredSize()+"\n"+this.getVisibleRect());
		
		g.setColor(Color.GRAY);
        g.fillRect(0, 0, getWidth(), getHeight());
		
		Graphics2D g2 = (Graphics2D) g;

		g2.scale(scale, scale);
		g2.drawImage(image, 0, 0, this);
		if(legel.isnotnull()){
			legel.drawSpots(g);
		}
	}
	
	public int getImageHeight(){ return image.getHeight(); }
	
	public int getImageWidth(){ return image.getWidth(); }
	
	public double getScale(){ return scale; }
	
	public void setTextPopup(Spot s){
		nameS.setText("  Name: "+s.getName()+"  ");
		piS.setText("  pI: "+s.getPi()+"  ");
		if(s.getPw().compareTo("")==0){
			pwS.setText("  Molecular weight:   ");
		}else{
			pwS.setText("  Molecular weight: "+s.getPw()+" kDa  ");
		}
		if(s.getProt() == null){
			protMaj.setText("  Top protein:   ");
		}else{
			protMaj.setText("  Top protein: "+s.getProt()+"  ");
		}
		iten.setText("  Intensity: "+s.getIntensity()+"  ");
		if(s.getExternalid() == null){
			exID.setText("  External identifier:   ");
		}else{
			exID.setText("  External identifier: "+s.getExternalid()+"  ");
		}

	}
	
	/**
	 * Give the size of the JLabel that is the biggest
	 */
	public int getMaxLengthWidth(){
		int maximum = nameS.getWidth();
		if( piS.getWidth() > maximum){ maximum = piS.getWidth(); }
		if( pwS.getWidth() > maximum){ maximum = pwS.getWidth(); }
		if( protMaj.getWidth() > maximum){ maximum = protMaj.getWidth(); }
		if( iten.getWidth() > maximum){ maximum = iten.getWidth(); }
		if( exID.getWidth() > maximum){ maximum = exID.getWidth(); }
		return maximum;
	}
	
	public int getMaxLengthHeight(){
		return (piS.getHeight()+pwS.getHeight()+nameS.getHeight()+protMaj.getHeight()+iten.getHeight()+exID.getHeight());
	}
	
	public void setVisiblePopup(boolean b){
		myPopup.setVisible(b);
	}
	
	public void showPopup(int a, int b){
		int maxW = getMaxLengthWidth();
		int maxH = getMaxLengthHeight();
		Rectangle visible = getVisibleRect();
		
		if(a + maxW > visible.getWidth()){
			a = a - 15 - maxW;
		}else{
			a = a + 15;
		}
		if(b + maxH > visible.getHeight()){
			b = b - 15 - maxH;
		}
		myPopup.show(this, a, b);
	}
	
	public void setGel(Gel g){
		legel = g;
	}
	
	/**
	 * Method to initiate the Fit to page value for the zoom
	 * invocated by the constructor
	 */
	public void initScale(Rectangle visibleSize){
		int width = image.getWidth();
		int height = image.getHeight();
		
		double division=1.00; 
		
		if((double)(height) / visibleSize.getHeight() > 1.00){
			division = visibleSize.getHeight() / (double)(height);
		}
		
		if( (double)(width)*division - visibleSize.getWidth() > 0.00 ){
			division = visibleSize.getWidth() / (double)(width);
		}

		scale = division*Math.sqrt(scalefactor);
	}
	
	public void setScaleFactor(double d){ scalefactor = d; }
}