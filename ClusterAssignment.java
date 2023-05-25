package cmsc420_s23;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class ClusterAssignment<LPoint extends LabeledPoint2D> {
	
	private SMkdTree<LPoint> tree;
	private ArrayList<LPoint> centers;
	private ArrayList<LPoint> sites;
	private LPoint startCenter1;

	public ClusterAssignment(int rebuildOffset, Rectangle2D bbox, LPoint startCenter) throws Exception { 
		//Create a new structure with an empty kd-tree
		tree = new SMkdTree<LPoint>(rebuildOffset, bbox, startCenter);
		centers = new ArrayList<LPoint>();
		sites = new ArrayList<LPoint>();
		
		//Initialize the cluster center set to consist of startCenter
		centers.add(startCenter);
		startCenter1 = startCenter;
	}
	
	//Add a new site in the structure
	public void addSite(LPoint site) throws Exception { 
		tree.insert(site);
		sites.add(site); //for listAssignments()
	}
	
	public void deleteSite(LPoint site) throws Exception { 
		//If there was a rebuild...
		if(tree.delete(site.getPoint2D()) == true) {
			//Generate contender lists for nodes of the tree
			for(LPoint point : centers) {
				tree.addCenter(point);
			}
		}
	}
	
	public void addCenter(LPoint center) throws Exception { 
		centers.add(center); //but then by this logic... we would only want centers to add center if addcenter added it to contenders
		tree.addCenter(center);
	}
	
	public int sitesSize() { 
		return tree.size();
	}
	
	public int centersSize() {
		return centers.size();
	}
	
	public void clear() { 
		tree.clear();
		//Removes all centers except startCenter
		centers.clear();
		centers.add(startCenter1);	
	}
	
	public ArrayList<String> listKdWithCenters() { 
		return tree.listWithCenters();
	}
	
	public ArrayList<String> listCenters() {
		ArrayList<String> list = new ArrayList<String>();
		
		sortByLabel(centers);

		for(LPoint c : centers) {
			String str = c.getLabel() + ": (" + c.getX() + "," + c.getY() + ")";
			list.add(str);
		}
		
		return list;	
	}
	
	public LPoint deleteCenter(LPoint center) {
		return null;
	}
	
	public ArrayList<String> listAssignments() {
		ArrayList<String> str = new ArrayList<>();
		str = tree.listAssignments();
		return str;
	}
	
	public void sortByLabel(List<LPoint> centers) {
		Comparator<LPoint> labelComparator = new Comparator<LPoint>() {
			@Override
	        public int compare(LPoint p1, LPoint p2) {
				return p1.getLabel().compareTo(p2.getLabel());
	        }
	    };
	    Collections.sort(centers, labelComparator);
	}
}
