import java.awt.*;
/**
 * Class that gets all the information you need about a spot
 * @author Guillaume ARRAS (Guillaume.Arras@curie.fr)
 * @version 15/05/2009
 */
public class Spot {
	//Parameters
	private int x;//x localisation in the gel
	private int y;//y localisation in the gel
	private String name;// name of the given spot
	private double pi;//Isoelectric point of the spot
	private double pw;//weight of the spot
	private String proteinMaj;//The most present protein of the spot
	private boolean isSelected;
	private static final int spotArea = 20;//which pixels should be included in the search of an environment
	private int crosssize=7;
	private int ID_SPOT;//Id representation of a spot in myProMS database
	private double intensity=-1.0;
	private String externalid;
	private int ID_SAMPLE;
	private String sampleName;
	private boolean spotColor=false;
	
	//Constructors
	public Spot(){// By default, a fake spot to test the class
		x = 50;
		y = 36;
		name = "spot test";
		pi = 45.0008;
		pw = 75.9865;
		isSelected = false;
		proteinMaj = "NONE";
	}
	
	/**
	 * Information sent by javascript to create a new Spot
	 * @param info
	 */
	public Spot(String info){
		System.out.println(info);
		String[] jsInfo=info.split("##");
		this.ID_SPOT = Integer.parseInt(jsInfo[0]);
		System.out.println("ID_SPOT="+ID_SPOT);
		this.name = jsInfo[1];
		System.out.println("name="+name);
		this.x = Integer.parseInt(jsInfo[2]);
		System.out.println("x="+x);
		this.y = Integer.parseInt(jsInfo[3]);
		System.out.println("y="+y);
		if(jsInfo[4].isEmpty()){
			this.pi = -1;
		}else{
			this.pi = Double.parseDouble(jsInfo[4]);
		}
		System.out.println("pi="+pi);
		if(jsInfo[5].isEmpty()){
			this.pw = -1;
		}else{
			this.pw = Double.parseDouble(jsInfo[5]);
		}
		System.out.println("pw="+pw);
		if(jsInfo[6].isEmpty()){
			this.intensity = -1;
		}else{
			this.intensity = Double.parseDouble(jsInfo[6]);
		}
		System.out.println("intensity="+intensity);
		if(jsInfo[7].isEmpty()){
			this.externalid = null;
		}else{
			this.externalid = jsInfo[7];
		}
		System.out.println("externalid="+externalid);
		if(jsInfo[8].isEmpty()){
			this.proteinMaj=null;
		}else{
			this.proteinMaj=jsInfo[8];
		}
		System.out.println("proteinMaj="+proteinMaj);
		if(Integer.parseInt(jsInfo[9])<0){
			this.ID_SAMPLE = -1;
			this.sampleName = null;
		}else{
			this.ID_SAMPLE = Integer.parseInt(jsInfo[9]);
			this.sampleName = jsInfo[10];
		}
		System.out.println("ID_SAMPLE="+ID_SAMPLE);
		System.out.println("sampleName="+sampleName);
	}
	
	public Spot(int ID_SPOT, String name, int x, int y, double pi, double pw, 
			double intensity, String externalid, int ID_SAMPLE, String sampleName, String proteinMaj){
		this.ID_SPOT=ID_SPOT;
		this.ID_SAMPLE=ID_SAMPLE;
		this.name=name;
		this.x=x;
		this.y=y;
		this.pi=pi;
		this.pw=pw;
		this.intensity=intensity;
		this.proteinMaj=proteinMaj;
		this.externalid=externalid;
		this.sampleName=sampleName;
	}

	
	public Spot(int x, int y, String name, double pi, double pw){
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.externalid = null;
	}
	
	public Spot(int idspot,int x, int y, String name, double pi, double pw){
		ID_SPOT = idspot;
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.externalid = null;
	}
	
	public Spot(int idspot,int x, int y, String name, double pi, double pw,
			double intensity, String externalid){
		ID_SPOT = idspot;
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.intensity = intensity;
		this.externalid = externalid;
	}
	
	public Spot(int idspot,int x, int y, String name, double pi, double pw, String externalid){
		ID_SPOT = idspot;
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.externalid = externalid;
	}
	
	public Spot(int x, int y, String name, double pi, double pw, 
			double intensity, String externalid){
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.intensity = intensity;
		this.externalid = externalid;
	}
	
	public Spot(int x, int y, String name, double pi, double pw, String externalid){
		this.x = x;
		this.y = y;
		this.name = name;
		this.pi = pi;
		this.pw = pw;
		isSelected = false;
		proteinMaj = "NONE";
		this.externalid = externalid;
	}
	
	//Method
	public void paint(Graphics g){
		if(spotColor){//The spot contains a search protein
			g.setColor(Color.RED);
			for(int i = 2*crosssize ; i < 3*crosssize ; i++){
				g.drawOval(x-i,y-i,2*i,2*i);
				g.drawOval(x-(i+1),y-(i+1),2*(i+1),2*(i+1));
				g.drawOval(x-(i+2),y-(i+2),2*(i+2),2*(i+2));
				g.drawOval(x-(i-1),y-(i-1),2*(i-1),2*(i-1));
				g.drawOval(x-(i-2),y-(i-2),2*(i-2),2*(i-2));
			}
		}
		if(isSelected){
			g.setColor(Color.RED);
			g.drawLine(x, y-crosssize, x, y+crosssize);
			g.drawLine(x-crosssize, y, x+crosssize, y);
		}else{
			g.setColor(new Color(51,153,255));
			g.drawLine(x, y-crosssize, x, y+crosssize);
			g.drawLine(x-crosssize, y, x+crosssize, y);
		}
	}
	
	public String toString(){
		return " Spot-Name : \t "+name+" \t pI: \t "+pi+" \t pW: \t "+pw+" \t X : \t"+x+"\t Y : \t"+y+"\t id_spot: \t"+ID_SPOT+"\t intensity : \t"+intensity+"\t externalid : \t"+externalid+ "\t sampleName : \t" + sampleName + "\t ID_SAMPLE : \t" + ID_SAMPLE +"\n";
	}
	
	public boolean equalsto(Spot s){
		if( (name.equals(s.getName())) && (x == s.getX()) && (y == s.getY())
				&& (Double.toString(pi).compareTo(s.getPi()) == 0) 
				&& (Double.toString(pw).compareTo(s.getPw()) == 0) ){
			return true;
		}
		return false;
	}
	
	public int getX(){ return x; }
	public int getY(){ return y; }
	public String getName(){ return name; }
	public String getPi(){
		if(pi<0.0){
			return "";
		}
		return Double.toString(pi);
	}
	public String getPw(){
		if(pw<0.0){
			return "";
		}
		return Double.toString(pw);
	}
	public String getIntensity(){
		if(intensity<0.0){
			return "";
		}
		return Double.toString(intensity);
	}
	public String getExternalid(){ return externalid; }
	public int getSampleID(){ return ID_SAMPLE; }
	public String getSampleName(){ return sampleName; }
	public void setSampleName(String s){ sampleName=s; }
	public void setSampleID(int i){ ID_SAMPLE=i; }
	
	/**
	 * Class that checks if a location on the image may be the fact of the spot or not
	 * by searching in the area of the coordinates given and the scale factor of the image
	 * @param dx X coordinate of the location searched
	 * @param dy Y coordinate of the location searched
	 * @param scale Scale factor of the image
	 * @return Whether the spot is around or not according to spotArea defined
	 */
	public boolean isAround(int dx, int dy, double scale){
		/*
		 * If the spot is around, return true
		 */
		if( (dx < (int)((x + spotArea)*scale ) ) && (dx > (int)((x - spotArea)*scale))
				&& (dy < (int) ((y + spotArea)*scale) ) && (dy > (int)((y - spotArea)*scale) ) ){
			return true;
		}
		return false;
	}
	
	public boolean isSelected(){ return isSelected; }
	
	public void select(){
		isSelected = true;
		spotColor = false;
	}
	
	public void deselect(){
		isSelected = false;
		spotColor = false;
	}
	
	public void hasProtein(){
		spotColor = true;
	}
	
	public String getProt(){ return proteinMaj; }
	
	public void setProt(String s){ 
		if(s == null){
			proteinMaj = "NONE";
		}else{
			proteinMaj = s;
		}
	}
	
	public void setIdSpot(int i){ ID_SPOT = i; }
	
	public int getIdSpot(){ return ID_SPOT; }
}
