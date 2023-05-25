package cmsc420_s23;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;


public class SMkdTree<LPoint extends LabeledPoint2D> {
	
	
	// ------------------------------------------------------------------------
	// The following class is not required, but you may find it helpful. It
	// represents the triple (site, center, squared-distance). Feel free to
	// delete or modify.
	// ------------------------------------------------------------------------
	public class AssignedPair implements Comparable<AssignedPair> {
		public LPoint site; // a site
		public LPoint center; // its assigned center
		public double distanceSq; // the squared distance between them
		
		public int compareTo(AssignedPair o) {
			
			if(this.distanceSq > o.distanceSq) {
				return 1;
			} else if (this.distanceSq < o.distanceSq) {
				return -1;
			} else { //If the distance is the same, use the site to compare
				LPoint lp1 = this.site;
				LPoint lp2 = o.site;
				
				if(lp1.getX() > lp2.getX()) {
					return 1;
				} else if(lp1.getX() < lp2.getX()) {
					return -1;
				} else {
					if(lp1.getY() > lp2.getY()) {
						return 1;
					} else if(lp1.getY() < lp2.getY()) {
						return -1;
					}else {
						return 0;
					}
				}
			}
		} // for sorting
		
		//Constructor
		public AssignedPair(LPoint site, LPoint center, double distanceSq) {
			this.site = site;
			this.center = center;
			this.distanceSq = site.getPoint2D().distanceSq(center.getPoint2D());
		}
		
	}

	// -----------------------------------------------------------------
	// Generic Node Type
	// -----------------------------------------------------------------

	private abstract class Node { // generic node type
		
		//private LinkedList<LPoint> contenders;
		
		InternalNode parent; // our parent

		// -----------------------------------------------------------------
		// Standard dictionary helpers
		// -----------------------------------------------------------------

		abstract LPoint find(Point2D q); // find point in subtree

		abstract Node insert(LPoint pt, Rectangle2D cell) throws Exception; // insert a single point
		
		abstract Node delete(Point2D pt) throws Exception; // delete point

		abstract ArrayList<String> list(); // list entries in reverse preorder

		abstract LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited); // nearest-neighbor helper

		// -----------------------------------------------------------------
		// Other utilities
		// -----------------------------------------------------------------

		abstract String debugPrint(String prefix); // print for debugging

		abstract public String toString(); // string representation
		
		abstract ArrayList<LPoint> bruteForceNN(Point2D center); // brute-force nearest neighbor
		
		abstract Node rebuildAfterInsert(LPoint pt, Rectangle2D cell); // check for rebuild
		
		abstract void buildPointList(List<LPoint> pts); // build a list of points in this subtree
		
		abstract Node rebuild(Rectangle2D cell); // rebuild subtree rooted at this node
		
		abstract int getSize(); // get number of points in subtree
		
		abstract ExternalNode leftMost(); // get leftmost descendant
		
		//New
		abstract void addCenter(LPoint center, Rectangle2D rootCell);
		
		abstract ArrayList<String> listWithCenters();
		
		abstract public String toString2(); // string representation

		abstract public ArrayList<AssignedPair> listAssignments();
		private LinkedList<LPoint> contenders;

	}

	// -----------------------------------------------------------------
	// Internal node
	// -----------------------------------------------------------------

	private class InternalNode extends Node {

		int cutDim; // the cutting dimension (0 = x, 1 = y)
		double cutVal; // the cutting value
		Node left, right; // children
		int size; // number of points
		int insertCt; // insert counter
		
		//Contenders list
		LinkedList<LPoint> contenders;

		/**
		 * Constructor
		 */
		InternalNode(int cutDim, double cutVal) {
			this.parent = null;
			this.cutDim = cutDim;
			this.cutVal = cutVal;
			this.left = new ExternalNode(null);
			this.right = new ExternalNode(null);
			this.insertCt = 0;
			this.size = 0;
			this.contenders = new LinkedList<>();
			
			//Every node will have at least one contender, even when the tree is empty
			this.contenders.add(startCenter);
		}

		// -----------------------------------------------------------------
		// Internal-node Utilities
		// -----------------------------------------------------------------

		boolean onLeft(Point2D pt) { // in the left subtree? (for Point2D)
			return pt.get(cutDim) < cutVal;
		}

		boolean onLeft(LPoint pt) { // in the left subtree? (for LPoint)
			return pt.get(cutDim) < cutVal;
		}

		Rectangle2D leftPart(Rectangle2D cell) { // left part of cell
			return cell.leftPart(cutDim, cutVal);
		}

		Rectangle2D rightPart(Rectangle2D cell) { // right part of cell
			return cell.rightPart(cutDim, cutVal);
		}

		// -----------------------------------------------------------------
		// Internal-node helpers
		// -----------------------------------------------------------------

		/**
		 * Find point in this subtree
		 */
		LPoint find(Point2D q) {
			if (onLeft(q)) return left.find(q);
			else return right.find(q);
		}

		/**
		 * Insert a single point.
		 */
		Node insert(LPoint pt, Rectangle2D cell) throws Exception{
			if (onLeft(pt)) {
				left = left.insert(pt, leftPart(cell));
				left.parent = this;
			} else {
				right = right.insert(pt, rightPart(cell));
				right.parent = this;
			}
			size += 1; // increment our size
			insertCt += 1; // increment our insert counter
			return this;
		}
		
		/**
		 * Check/rebuild if needed after insertion. (Internalnode)
		 */
		Node rebuildAfterInsert(LPoint pt, Rectangle2D cell) {
			if (insertCt > (size + rebuildOffset)/2) { // time to rebuild
				//Saved contenders
				LinkedList<LPoint> savedContenders = this.contenders;
				//Node returned from bulkcreate
				Node bulkCreateReturned = rebuild(cell);
				
				//For each one of the contenders, call addCenter on the node returned
				for(LPoint point : savedContenders) {
					bulkCreateReturned.addCenter(point, cell);
				}
				
				return bulkCreateReturned;
				
			} else if (onLeft(pt)) {
				left = left.rebuildAfterInsert(pt, leftPart(cell));
				left.parent = this;
			} else {
				right = right.rebuildAfterInsert(pt, rightPart(cell));
				right.parent = this;
			}
			return this;
		}

		/**
		 * Delete a point from this subtree. We recurse on the side
		 * containing the point. A null value signifies that our
		 * child is an external node that has lost its last point,
		 * and we respond by returning the other child (since we are
		 * no longer needed).
		 */
		Node delete(Point2D pt) throws Exception {
			if (onLeft(pt)) {
				left = left.delete(pt);
				left.parent = this;
			}
			else {
				right = right.delete(pt);
				right.parent = this;
			}
			size -= 1; // one less point
			return this;
		}

		/**
		 * Get a list of the nodes in reverse preorder.
		 *
		 * Adds current key in parentheses, followed by left and right subtrees.
		 */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			list.addAll(right.list()); // add right
			list.addAll(left.list()); // add left
			return list;
		}
		
		//listWithCenters
		ArrayList<String> listWithCenters() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString2()); // add this node
			list.addAll(right.listWithCenters()); // add right
			list.addAll(left.listWithCenters()); // add left
			return list;
		}

		/**
		 * Nearest neighbor helper. The closest point seen so far is
		 * best. We search the closer subtree first and the other 
		 * subtree if it is still viable.
		 */
		LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			// child cells
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal);
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal);
			if (onLeft(center)) { // left side is closer
				best = left.nearNeigh(center, leftPart(cell), best, visited);
				double bestDist = (best == null ? Double.POSITIVE_INFINITY : center.distanceSq(best.getPoint2D()));
				if (rightCell.distanceSq(center) <= bestDist) {
					best = right.nearNeigh(center, rightPart(cell), best, visited);
				}
			} else { // right side is closer
				best = right.nearNeigh(center, rightPart(cell), best, visited);
				double bestDist = (best == null ? Double.POSITIVE_INFINITY : center.distanceSq(best.getPoint2D()));
				if (leftCell.distanceSq(center) <= bestDist) {
					best = left.nearNeigh(center, leftPart(cell), best, visited);
				}
			}
			return best;
		}

		/**
		 * Rebuild the subtree rooted at this node. We first compute a sorted list of
		 * external nodes and then rebuild the subtree recursively.
		 */
		Node rebuild(Rectangle2D cell) {
			List<LPoint> pts = new ArrayList<LPoint>(); // list for points
			buildPointList(pts); // collect the points in this subtree
			return bulkCreate(pts, cell, this.contenders); // build a new tree
		}

		/**
		 * Builds a list of points in this subtree.
		 */
		void buildPointList(List<LPoint> pts) {
			left.buildPointList(pts); // add left
			right.buildPointList(pts); // add right
		}
		
		/**
		 * Get leftmost descendant. This is used in the iterator.
		 */
		ExternalNode leftMost() {
			return left.leftMost();
		}

		// -----------------------------------------------------------------
		// Internal-node utilities
		// -----------------------------------------------------------------

		/**
		 * Get subtree size.
		 */
		int getSize() {
			return size;
		}

		/**
		 * String representation.
		 */
		public String toString() {
			String cutIndic = (cutDim == 0 ? "x=" : "y=");
			return "(" + cutIndic + cutVal + ") " + size + ":" + insertCt;
		}
		
		//for listWithCenters()
		public String toString2() {
			String cutIndic = (cutDim == 0 ? "x=" : "y=");
			
			String contenderStr = "";
			
			//Sort contenders by alphabetical order
			sortByLabel(contenders);
			
			if(contenders.size() <= 10) {
				/*for(LPoint c : contenders) {
					contenderStr += c.getLabel();
					contenderStr += " ";
				}*/
				
				int i;
				for(i = 0; i < contenders.size() - 1; i++) {
			        String label = contenders.get(i).getLabel();
			        //if (!contenderStr.contains(label)) { // check if label is already in contenderStr
			            contenderStr += label;
			            contenderStr+= " ";
			        //}
				}
				contenderStr += contenders.get(i).getLabel();
			} else {
				for(int i = 0; i < 10; i++) {
			        String label = contenders.get(i).getLabel();
			        //if (!contenderStr.contains(label)) { // check if label is already in contenderStr
			            contenderStr += label;
			            if(i != 9) {
				            contenderStr+= " ";
			            }
			        //}
				}
				contenderStr += "...";
			}
			
			return "(" + cutIndic + cutVal + ") " + size + ":" + insertCt + " => {" + contenderStr + "}";
		}

		/**
		 * Brute-force nearest neighbor (for validation)
		 */
		ArrayList<LPoint> bruteForceNN(Point2D center) {
			ArrayList<LPoint> leftClose = left.bruteForceNN(center); // get left and right closest
			ArrayList<LPoint> rightClose = right.bruteForceNN(center);
			if (leftClose.isEmpty()) { // if either is empty, return the other
				return rightClose;
			} else if (rightClose.isEmpty()) {
				return leftClose;
			} else { // both are nonempty
				double leftDist = center.distanceSq(leftClose.get(0).getPoint2D());
				double rightDist = center.distanceSq(rightClose.get(0).getPoint2D());
				if (leftDist < rightDist) { // left is closer
					return leftClose;
				} else if (rightDist < leftDist) { // right is closer
					return rightClose;
				} else { // same distances
					leftClose.addAll(rightClose); // merge both lists
					return leftClose;
				}
			}
		}

		/**
		 * Debug print subtree.
		 *
		 * This prints the subtree for debugging.
		 */
		String debugPrint(String prefix) {
			return  right.debugPrint(prefix + "| ") + System.lineSeparator() + // right child
					prefix + toString() + System.lineSeparator() + // this node
					left.debugPrint(prefix + "| "); // left child
		}

		@Override
		void addCenter(LPoint center, Rectangle2D rootCell) {
	        double rNew = rootCell.maxDistanceSq(center.getPoint2D());
	        double rMin = rNew;
	        
	        //Go through contenders, save min value among all of these including center as rMin
	        for (LPoint point : contenders) {
	            rNew = rootCell.maxDistanceSq(point.getPoint2D());
	            if(rNew < rMin) {
	            	rMin = rNew;
	            }
	        }
	        
	        //Go through contenders a second time, compute ri value.  If it is bigger than rMin,
	        //remove that center from contenders
	        LinkedList<LPoint> removing = new LinkedList<>();
	        for(LPoint point : contenders) {
	        	rNew = rootCell.distanceSq(point.getPoint2D());
	        	
	        	//If rNew is bigger than rMin, remove this point from contenders
	        	if(rNew > rMin) {
	        		removing.add(point);
	        	}
	        }
	        contenders.removeAll(removing);
	        
	        //If center's ri is bigger than Rmin, don't add it
	        if(rootCell.distanceSq(center.getPoint2D()) <= rMin) {
	        	
	        	if(!contenders.contains(center)) {
		        	contenders.add(center);
	        	}
	        	
	        	//Recursively do it on left and right subtrees, DOUBLE CHECK THIS PART
	        	this.left.addCenter(center, rootCell.leftPart(cutDim, cutVal));
	        	this.right.addCenter(center, rootCell.rightPart(cutDim, cutVal));
	        }
		}
		
		public ArrayList<AssignedPair> listAssignments() {
			ArrayList<AssignedPair> apList = new ArrayList<>();
			apList.addAll(this.left.listAssignments());
			apList.addAll(this.right.listAssignments());
			return apList;
		}
	}

	// -----------------------------------------------------------------
	// External node
	// -----------------------------------------------------------------

	private class ExternalNode extends Node {

		LPoint point; // the point (null if no point)
		
		//Contenders list
		LinkedList<LPoint> contenders;

		/**
		 * Default constructor.
		 */
		ExternalNode(LPoint point) {
			this.point = point;
			this.parent = null;
			this.contenders = new LinkedList<>();
			
			//Every node will have at least one contender, even when the tree is empty
			this.contenders.add(startCenter);
		}

		// -----------------------------------------------------------------
		// External-node helpers
		// -----------------------------------------------------------------

		/**
		 * Find point in external node.
		 */
		LPoint find(Point2D q) {
			if (point != null && point.getPoint2D().equals(q)) {
				return point;
			} else {
				return null;
			}
		}

		/**
		 * Insert a single point.
		 */
		Node insert(LPoint pt, Rectangle2D cell) throws Exception {
			if (point == null) {
				point = pt;
				return this;
			}
			else {
				if (point.getPoint2D().equals(pt.getPoint2D())) {
					throw new Exception("Insertion of duplicate point");
				}
				ArrayList<LPoint> pts = new ArrayList<LPoint>();
				pts.add(point);
				pts.add(pt);
				return bulkCreate(pts, cell, this.contenders); // create new tree with both points
			}
		}

		/**
		 * Rebuild if needed after insert. (No effect on externals.)
		 */
		Node rebuildAfterInsert(LPoint pt, Rectangle2D cell) {
			return this;
		}

		/**
		 * Delete a point from this node.
		 */
		Node delete(Point2D pt) throws Exception {
			if (point == null || !point.getPoint2D().equals(pt)) {
				throw new Exception("Deletion of nonexistent point");
			}
			point = null;
			return this;
		}

		/**
		 * List the contents of this node
		 */
		ArrayList<String> list() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString()); // add this node
			return list;
		}
		
		//listWithCenters
		ArrayList<String> listWithCenters() {
			ArrayList<String> list = new ArrayList<String>();
			list.add(toString2()); // add this node
			return list;
		}

		/**
		 * Nearest neighbor. Returns the closest point to center. 
		 */
		LPoint nearNeigh(Point2D center, Rectangle2D cell, LPoint best, ArrayList<LPoint> visited) {
			if (visited != null && this.point != null)
				visited.add(this.point);
			return closerOf(center, best, point);
		}
		
		/**
		 * Rebuild this subtree. (Does nothing for external nodes.)
		 */
		Node rebuild(Rectangle2D cell) {
			return this;
		}

		/**
		 * Builds a list of points in this subtree.
		 */
		void buildPointList(List<LPoint> pts) {
			if (point != null) pts.add(point);
		}
		
		/**
		 * Get leftmost descendant.
		 */
		ExternalNode leftMost() {
			return this;
		}

		// -----------------------------------------------------------------
		// External-node Utilities
		// -----------------------------------------------------------------

		/**
		 * Get subtree size.
		 */
		int getSize() {
			return point == null ? 0 : 1;
		}

		/**
		 * String representation.
		 */
		public String toString() {
			String contents = new String();
			contents += "[" + (point == null ? "null" : point.toString()) + "]";
			return contents;
		}
		
		//for listWithCenters
		public String toString2() {
			String contents = new String();
			
			
			String contenderStr = "";
			
			//Sort contenders
			sortByLabel(contenders);
			
			if(contenders.size() <= 10) {
				/*for(LPoint c : contenders) {
					contenderStr += c.getLabel();
					contenderStr += " ";
				}*/
				int i;
				for(i = 0; i < contenders.size() - 1; i++) {
			        String label = contenders.get(i).getLabel();
			        //if (!contenderStr.contains(label)) { // check if label is already in contenderStr
			            contenderStr += label;
			            contenderStr+= " ";
			        //}
				}
				contenderStr += contenders.get(i).getLabel();
			} else {
				for(int i = 0; i < 10; i++) {
			        String label = contenders.get(i).getLabel();
			        //if (!contenderStr.contains(label)) { // check if label is already in contenderStr
			            contenderStr += label;
			            if(i != 9) {
				            contenderStr+= " ";
			            }
			        //}
				}
				contenderStr += "...";
			}
			

						
			contents += "[" + (point == null ? "null" : point.toString()) + "] => {" + contenderStr + "}";
			return contents;
		}

		/**
		 * Debug print.
		 */
		String debugPrint(String prefix) {
			return prefix + toString();
		}

		/**
		 * Nearest neighbor by brute force (for validation). Returns
		 * the point in this node if it is non-null.
		 */
		ArrayList<LPoint> bruteForceNN(Point2D center) {
			ArrayList<LPoint> result = new ArrayList<LPoint>();
			if (point != null) result.add(point);
			return result;
		}

		@Override
		void addCenter(LPoint center, Rectangle2D rootCell) {
	        double rNew = rootCell.maxDistanceSq(center.getPoint2D());
	        double rMin = rNew;
	        
	        //Go through contenders, save min value among all of these including center as rMin
	        for (LPoint point : contenders) {
	            rNew = rootCell.maxDistanceSq(point.getPoint2D());
	            if(rNew < rMin) {
	            	rMin = rNew;
	            }
	        }
	        
	        //Go through contenders a second time, compute ri value.  If it is bigger than rMin,
	        //remove that center from contenders
	        LinkedList<LPoint> removing = new LinkedList<>();
	        for(LPoint point : contenders) {
	        	rNew = rootCell.distanceSq(point.getPoint2D());
	        	
	        	//If rNew is bigger than rMin, remove this point from contenders
	        	if(rNew > rMin) {
	        		removing.add(point);
	        	}
	        }
	        contenders.removeAll(removing);
	        
	        //If center's ri is bigger than Rmin, don't add it
	        if(rootCell.distanceSq(center.getPoint2D()) <= rMin) {
	        	
	        	if(!contenders.contains(center)) {
		        	this.contenders.add(center);
	        	}
	        	
	        }
		}
		
		public ArrayList<AssignedPair> listAssignments() {
			ArrayList<AssignedPair> apList = new ArrayList<>();
			
			//Access list of contenders if not null
			if(this.point != null) {
				LPoint currBest = null;
				double closestD = -1;
				
				for (int i = 0; i < this.contenders.size(); i++) {
					LPoint c = this.contenders.get(i);

					if (currBest == null) {
						currBest = c;
						closestD = this.point.getPoint2D().distanceSq(c.getPoint2D());
					}
					if (currBest != null) {
						double currClosest = this.point.getPoint2D().distanceSq(currBest.getPoint2D());
						double checker = this.point.getPoint2D().distanceSq(c.getPoint2D());
						if (checker < currClosest) {
							currBest = c;
							closestD = checker;
							//The distances match, now check x and y
						} else if (checker == currClosest) {
							if(currBest.getPoint2D().getX() < c.getPoint2D().getX()) {//Pick smaller x
								closestD = checker;
							} else if(currBest.getPoint2D().getX() > c.getPoint2D().getX()) {
								currBest = c;
								closestD = checker;
							} else {
								if(currBest.getPoint2D().getY() < c.getPoint2D().getY()) {//Pick smaller x
									closestD = checker;
								} else if(currBest.getPoint2D().getY() > c.getPoint2D().getY()) {
									currBest = c;
									closestD = checker;
								}
							}
						}
					}
				}
				
				AssignedPair newPair = new AssignedPair(this.point, currBest, closestD);
				apList.add(newPair);
			}
			return apList;
		} 
	}
	
	// -----------------------------------------------------------------
	// Iterator
	// -----------------------------------------------------------------
	
	private class LPointIterator implements Iterator<LPoint> {
	    private ExternalNode next; // external node with the next point

		/**
		 * Constructor. If tree is empty, the iterator is null. Otherwise
		 * Initialize to first point in leftmost leaf.
		 */
	    public LPointIterator() {
	    	if (nPoints == 0) {
	    		next = null;
	    	} else {
		    	next = root.leftMost();
	    	}
	    }

		/**
		 * More elements remaining?
		 */
	    public boolean hasNext() {
	    	advanceToNonNull();
	        return next != null; 
	    }

		/**
		 * Return current and advance to next.
		 */
		public LPoint next() throws NoSuchElementException {
			advanceToNonNull();
			if (next == null) {
				throw new NoSuchElementException();
			}
			LPoint result = next.point; // final result
			advanceOnce();
	    	return result;
		}
		
		/**
		 * Advance to next non-null point.
		 */
		private void advanceToNonNull() {
			while (next != null && next.point == null) {
				advanceOnce();
			}
		}

		/**
		 * Advance to next (which might be null).
		 */
		private void advanceOnce() {
			Node v = next;
			InternalNode u = v.parent; // move along right chain
			while (u != null && u.right == v) {
				v = u;
				u = u.parent;
			}
			if (u == null) { // hit the root?
				next = null;
			} else {
				next = u.right.leftMost(); // leftmost leaf of right child
			}
		}
	}

	// -----------------------------------------------------------------
	// Tree private data
	// -----------------------------------------------------------------

	private final int rebuildOffset; // offset for rebuild counters
	private Rectangle2D rootCell; // the bounding box
	private Node root; // root of the tree
	private int nPoints; // number of points in the tree
	private int deleteCt; // deletion counter
	
	//New
	private LPoint startCenter;

	// -----------------------------------------------------------------
	// Tree public members
	// -----------------------------------------------------------------

	/**
	 * Creates an empty tree.
	 */
	public SMkdTree(int rebuildOffset, Rectangle2D rootCell, LPoint startCenter) {
		this.rebuildOffset = rebuildOffset;
		this.rootCell = new Rectangle2D(rootCell);

		//New
		this.startCenter = startCenter;
		
		//The initial structure is an empty External node whose contender list is a single point
		this.root = new ExternalNode(null);
		//Contenders is already initialized in our externalnode constructor
		
		clear();
	}
	
	/**
	 * Return a list of sites whose centers had changed as a result of the addition
	 * (sites that have now been assigned to the newly inserted center)
	 * 
	 * This function will be called in ClusterAssignments
	 */
	void addCenter(LPoint center) {
		root.addCenter(center, rootCell);
	}


	public ArrayList<String> listAssignments() {
		ArrayList<String> ret = new ArrayList<>();
		ArrayList<AssignedPair> pairs = this.root.listAssignments();
		Collections.sort(pairs);
		for(AssignedPair p : pairs) {
			ret.add("[" + p.site.getLabel() + "->" + p.center.getLabel() + "] distSq = " + p.distanceSq);
		}
		return ret;
	}
	
	/**
	 * Clear the tree.
	 */
	public void clear() {
		root = new ExternalNode(null);
		nPoints = 0;
		deleteCt = 0;
	}

	/**
	 * Number of points.
	 */
	public int size() {
		return nPoints;
	}

	/**
	 * Delete counter value.
	 */
	public int deleteCount() {
		return deleteCt;
	}

	/**
	 * Find an point in the tree.
	 */
	public LPoint find(Point2D q) {
		return root.find(q);
	}

	/**
	 * Insert a point.
	 */
	public void insert(LPoint pt) throws Exception {
		if (!rootCell.contains(pt.getPoint2D())) { // outside bounds?
			throw new Exception("Attempt to insert a point outside bounding box");
		}
		root = root.insert(pt, rootCell); // insert
		root = root.rebuildAfterInsert(pt, rootCell); // rebuild if needed
		root.parent = null;
		nPoints += 1; // update point count
	}
	
	/**
	 * Delete a point. Note that the point being deleted does not need to match
	 * fully. It suffices that it has enough information to satisfy the comparator.
	 * 
	 * We first check whether the point is in the tree, using find. This is
	 * necessary, because (due to fact that points lying on the splitting line
	 * be in either subtree) it is difficult for the deletion helper to check
	 * for nonexistent points.
	 */
	public boolean delete(Point2D pt) throws Exception {
		root = root.delete(pt);
		root.parent = null;
		nPoints -= 1; // one fewer point
		deleteCt += 1; // increment deletion counter
		if (deleteCt > nPoints) {
			
			root = root.rebuild(rootCell); // time to rebuild
			
			deleteCt = 0; // reset the delete counter
			
			return true; //This tells us if there was a rebuild
		}
		return false; //No rebuild
	}

	/**
	 * Get a list of entries in preorder
	 */
	public ArrayList<String> list() {
		return root.list();
	}
	
	//Listwithcenters
	public ArrayList<String> listWithCenters() {
		return root.listWithCenters();
	}

	/**
	 * Find the nearest point to center. Ties are broken 
	 * lexicographically.
	 */
	public LPoint nearestNeighbor(Point2D center) {
		LPoint best = null;
		LPoint nn = root.nearNeigh(center, rootCell, best, null);
		return nn;
	}

	/**
	 * Same as nearestNeighbor, but we return the nodes visited.
	 */
	public ArrayList<LPoint> nearestNeighborVisit(Point2D center) {
		LPoint best = null;
		ArrayList<LPoint> visited = new ArrayList<LPoint>();
		root.nearNeigh(center, rootCell, best, visited);
		Collections.sort(visited, new ByXThenY());
		return visited;
	}
	
	/**
	 * Generate an iterator.
	 */
	public LPointIterator iterator() {
		return new LPointIterator();
	}

	// -----------------------------------------------------------------
	// Tree Utilities for sorting and splitting.
	// -----------------------------------------------------------------
	
	/**
	 * Create a new subtree from a list of points. This does most 
	 * of the work of the sliding midpoint method. If there is at
	 * most point point, we create a new external node. Otherwise, 
	 * we determine the longest dimension of the current cell
	 * (breaking ties in favor of x) and make this the cutting
	 * dimension. We then compute the cell's midpoint and make this
	 * the cutting value. We sort the points by the cutting 
	 * dimension and partition them about the cutting
	 * value. Points that lie on the cutting value are placed in
	 * the right subtree. If the resulting partition is trivial, we
	 * set the cutting value to the coordinate (largest or smallest).
	 */
	Node bulkCreate(List<LPoint> pts, Rectangle2D cell, LinkedList<LPoint> contenders) {
		
//		Node u = bulkCreateHelper(pts, cell); // build a new tree
//		return u;
		if (pts.size() == 0) { // no points
			return new ExternalNode(null);
		} else if (pts.size() == 1) { // one point
			return new ExternalNode(pts.get(0));
		} else { // two or more points
			int cutDim = getWideSide(cell); // get cutting dimension

			sortPts(pts, cutDim); // sort along cutting dimension
			LPoint leftmost = pts.get(0); // leftmost/rightmost points
			LPoint rightmost = pts.get(pts.size()-1);
			if (leftmost.get(cutDim) == rightmost.get(cutDim)) {
				cutDim = 1 - cutDim; // reverse splitting direction
			}
			double cutVal = (cell.low.get(cutDim) + cell.high.get(cutDim))/2; // default cutting value
			if (rightmost.get(cutDim) < cutVal) { // all points on left?
				cutVal = rightmost.get(cutDim);
			} else if (leftmost.get(cutDim) > cutVal) { // all points on right?
				cutVal = leftmost.get(cutDim);
			}
			int splitIndex = 0; // where to split
			while (splitIndex < pts.size() && pts.get(splitIndex).get(cutDim) < cutVal) {
				splitIndex++;
			}
													// left/right cells
			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal);
			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal);
													// the new node
			InternalNode u = new InternalNode(cutDim, cutVal);
			u.left = bulkCreate(pts.subList(0, splitIndex), leftCell, u.left.contenders);
			u.right = bulkCreate(pts.subList(splitIndex,  pts.size()), rightCell, u.right.contenders);
			u.left.parent = u.right.parent = u;
			u.size = u.left.getSize() + u.right.getSize();
			
			if(contenders != null){
				for(LPoint l: contenders) {
					u.addCenter(l, cell);
				}
			}
			return u;
		}
	}
	
//	Node bulkCreateHelper(List<LPoint> pts, Rectangle2D cell) {
//		if (pts.size() == 0) { // no points
//			return new ExternalNode(null);
//		} else if (pts.size() == 1) { // one point
//			return new ExternalNode(pts.get(0));
//		} else { // two or more points
//			int cutDim = getWideSide(cell); // get cutting dimension
//
//			sortPts(pts, cutDim); // sort along cutting dimension
//			LPoint leftmost = pts.get(0); // leftmost/rightmost points
//			LPoint rightmost = pts.get(pts.size()-1);
//			if (leftmost.get(cutDim) == rightmost.get(cutDim)) {
//				cutDim = 1 - cutDim; // reverse splitting direction
//			}
//			double cutVal = (cell.low.get(cutDim) + cell.high.get(cutDim))/2; // default cutting value
//			if (rightmost.get(cutDim) < cutVal) { // all points on left?
//				cutVal = rightmost.get(cutDim);
//			} else if (leftmost.get(cutDim) > cutVal) { // all points on right?
//				cutVal = leftmost.get(cutDim);
//			}
//			int splitIndex = 0; // where to split
//			while (splitIndex < pts.size() && pts.get(splitIndex).get(cutDim) < cutVal) {
//				splitIndex++;
//			}
//													// left/right cells
//			Rectangle2D leftCell = cell.leftPart(cutDim, cutVal);
//			Rectangle2D rightCell = cell.rightPart(cutDim, cutVal);
//													// the new node
//			InternalNode u = new InternalNode(cutDim, cutVal);
//			u.left = bulkCreateHelper(pts.subList(0, splitIndex), leftCell);
//			u.right = bulkCreateHelper(pts.subList(splitIndex,  pts.size()), rightCell);
//			u.left.parent = u.right.parent = u;
//			u.size = u.left.getSize() + u.right.getSize();
//			return u;
//		}
//	}

	/**
	 * Get the wider side of a rectangle.
	 */
	private static int getWideSide(Rectangle2D cell) {
		return (cell.getWidth(0) >= cell.getWidth(1) ? 0 : 1);
	}

	/**
	 * Return a reference to the closer of two points. If the two
	 * distances are the same, it returns the point having the smaller
	 * (x,y) coordinates in lexicographical order.
	 */
	private LPoint closerOf(Point2D q, LPoint p1, LPoint p2) {
		double dist1 = (p1 == null ? Double.POSITIVE_INFINITY : q.distanceSq(p1.getPoint2D()));
		double dist2 = (p2 == null ? Double.POSITIVE_INFINITY : q.distanceSq(p2.getPoint2D()));
		if      (p1 == null || dist1 > dist2) return p2;
		else if (p2 == null || dist2 > dist1) return p1;
		else { // same distances and both are non-null
			if      (p1.getX() < p2.getX()) return p1;
			else if (p2.getX() < p1.getX()) return p2;
			else { // ...and same x-coordinates
				if (p1.getY() <= p2.getY()) return p1;
				else return p2;
			}
		}
	}


	/**
	 * Sort points according to given dimension.
	 */
	private void sortPts(List<LPoint> pts, int d) {
		if (d == 0) { // sort by the cutting dimension
			Collections.sort(pts, new ByXThenY());
		} else {
			Collections.sort(pts, new ByYThenX());
		}
	}

	/**
	 * Comparators for sorting (by X or by Y). Both are lexicographical.
	 */
	private class ByXThenY implements Comparator<LPoint> { // lexicographic (x,y)
		public int compare(LPoint pt1, LPoint pt2) {
			double x1 = pt1.getX();
			double x2 = pt2.getX();
			if (x1 < x2)
				return -1;
			else if (x1 > x2)
				return +1;
			else {
				double y1 = pt1.getY();
				double y2 = pt2.getY();
				if (y1 < y2)
					return -1;
				else if (y1 > y2)
					return +1;
				else {
					return 0;
				}
			}
		}
	}

	private class ByYThenX implements Comparator<LPoint> { // lexicographic (y,x)
		public int compare(LPoint pt1, LPoint pt2) {
			double y1 = pt1.getY();
			double y2 = pt2.getY();
			if (y1 < y2)
				return -1;
			else if (y1 > y2)
				return +1;
			else {
				double x1 = pt1.getX();
				double x2 = pt2.getX();
				if (x1 < x2)
					return -1;
				else if (x1 > x2)
					return +1;
				else {
					return 0;
				}
			}
		}
	}
	
				
	    
	public void sortByLabel(List<LPoint> contenders) {
		Comparator<LPoint> labelComparator = new Comparator<LPoint>() {
			@Override
	        public int compare(LPoint p1, LPoint p2) {
				return p1.getLabel().compareTo(p2.getLabel());
	        }
	    };
	    Collections.sort(contenders, labelComparator);
	}
	    
}
