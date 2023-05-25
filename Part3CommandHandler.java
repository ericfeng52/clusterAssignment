package cmsc420_s23;

// YOU SHOULD NOT MODIFY THIS FILE

import java.util.Scanner;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

/**
 * Command handler. This reads a single command line, processes the command (by
 * invoking the appropriate method(s) on the data structure), and returns the
 * result as a string.
 */

public class Part3CommandHandler {

	private boolean initialized; // have we initialized the structure yet?
	private int rebuildOffset; // rebuild offset for the kd-tree
	private Rectangle2D bbox; // bounding box for the point set
	private HashMap<String, Airport> airports; // airport codes seen so far
	private ClusterAssignment<Airport> clusterAssignment; // cluster-assignment structure
	private FarthestFirst<Airport> farthestFirst; // farthest-first structure

	/**
	 * Initialize command handler
	 */
	public Part3CommandHandler() {
		initialized = false;
		airports = new HashMap<String, Airport>();
		clusterAssignment = null;
		farthestFirst = null;
	}

	/**
	 * Process a single command and return the string output. Each command begins
	 * with a command followed by a list of arguments. The arguments are separated
	 * by colons (":").
	 */

	public String processCommand(String inputLine) throws Exception {
		String output = new String(); // for storing summary output
		Scanner line = new Scanner(inputLine);
		try {
			line.useDelimiter(":"); // use ":" to separate arguments
			String cmd = (line.hasNext() ? line.next() : ""); // next command
			// -----------------------------------------------------
			// INITIALIZE
			// - this command must come first in the input
			// - sets the bounding box
			// -----------------------------------------------------
			if (cmd.compareTo("initialize") == 0) {
				rebuildOffset = line.nextInt(); // save rebuild offset
				double xMin = line.nextDouble(); // bounding box
				double xMax = line.nextDouble();
				double yMin = line.nextDouble();
				double yMax = line.nextDouble();
				if (xMin > xMax || yMin > yMax) {
					throw new Exception("Error - Invalid bounding box dimensions");
				}
				if (initialized) {
					throw new Exception("Error - Structure already initialized");
				} else {
					bbox = new Rectangle2D(new Point2D(xMin, yMin), new Point2D(xMax, yMax));
					output += "initialize: rebuild-offset = " + rebuildOffset + " bounding-box = " + bbox + System.lineSeparator();
					initialized = true;
				}
			}
			// -----------------------------------------------------
			// COMMENT string
			// - comment line for the output
			// -----------------------------------------------------
			else if (cmd.compareTo("comment") == 0) {
				String message = line.next(); // read the comment
				output += "comment: " + message + System.lineSeparator();
			}
			// -----------------------------------------------------
			// ALERT string
			// - comment with output also sent to stderr
			// -----------------------------------------------------
			else if (cmd.compareTo("alert") == 0) {
				String message = line.next(); // read the comment
				output += "alert: " + message + System.lineSeparator();
				System.err.println("... alert: " + message);
			}
			// -----------------------------------------------------
			// CA-INITIALIZE
			// - initializes an empty CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-initialize") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport start = readAirport(line); // read in airport
				output += "ca-initialize(start = " + start.getCode() + "): ";
				clusterAssignment = new ClusterAssignment<Airport>(rebuildOffset, bbox, start);
				output += "successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-INITIALIZE
			// - initializes an empty FF structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-initialize") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport start = readAirport(line); // read in airport
				output += "ff-initialize(start = " + start.getCode() + "): ";
				farthestFirst = new FarthestFirst<Airport>(rebuildOffset, bbox, start);
				output += "successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-ADD-SITE code city x y
			// - insert site in CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-add-site") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport ap = readAirport(line); // read in airport
				output += "ca-add-site(" + ap.getCode() + "): ";
				clusterAssignment.addSite(ap); // add to CA structure
				output += "successful {" + ap.getString("attributes") + "}" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-ADD-SITE code city x y
			// - insert site in FF structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-add-site") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport ap = readAirport(line); // read in airport
				output += "ff-add-site(" + ap.getCode() + "): ";
				farthestFirst.addSite(ap); // add to FF structure
				output += "successful {" + ap.getString("attributes") + "}" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-ADD-CENTER code city x y
			// - insert center in CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-add-center") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport ap = readAirport(line); // read in airport
				output += "ca-add-center(" + ap.getCode() + "): ";
				clusterAssignment.addCenter(ap); // add to CA structure
				output += "successful {" + ap.getString("attributes") + "}" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-DELETE-SITE code
			// - delete site from CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-delete-site") == 0) {
				confirmInitialized(); // confirm that we are initialized
				String code = line.next(); // get parameters
				output += "ca-delete-site(" + code + "): ";
				Airport ap = airports.get(code); // look up the airport
				if (ap == null) { // no such airport?
					throw new Exception("Deletion of nonexistent airport code");
				}
				clusterAssignment.deleteSite(ap); // delete from kd-tree
				airports.remove(code); // delete from dictionary
				output += "successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-LIST-KD
			// - list the kd-tree for the CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-list-kd") == 0) {
				confirmInitialized(); // confirm that we are initialized
				ArrayList<String> list = clusterAssignment.listKdWithCenters();
				if (list == null)
					throw new Exception("Error - clusterAssignment.listKdWithCenters returned a null result");
				output += "ca-list-kd:" + System.lineSeparator();
				for (String s : list) { // output the list
					output += "  " + s + System.lineSeparator();
				}
				output += treeStructure("ca-list-kd", list); // summarize tree contents (indented)
			}
			// -----------------------------------------------------
			// FF-LIST-KD
			// - list the kd-tree for the FF structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-list-kd") == 0) {
				confirmInitialized(); // confirm that we are initialized
				ArrayList<String> list = farthestFirst.listKdWithCenters();
				if (list == null)
					throw new Exception("Error - farthestFirst.listKdWithCenters returned a null result");
				output += "ff-list-kd:" + System.lineSeparator();
				for (String s : list) { // output the list
					output += "  " + s + System.lineSeparator();
				}
				output += treeStructure("ff-list-kd", list); // summarize tree contents (indented)
			}
			// -----------------------------------------------------
			// FF-TRAVERSAL-SIZE
			// - how many points in the traversal
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-traversal-size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = farthestFirst.traversalSize(); // get the traversal's current size
				output += "ff-traversal-size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-SITES-SIZE
			// - how many sites are in CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-sites-size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = clusterAssignment.sitesSize();
				output += "ca-sites-size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-SITES-SIZE
			// - how many sites are in FF structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-sites-size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = farthestFirst.sitesSize(); 
				output += "ff-sites-size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// CA-CENTERS-SIZE
			// - how many sites are in CA structure
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-centers-size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = clusterAssignment.centersSize(); 
				output += "ca-centers-size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-TRAVERSAL-SIZE
			// - current size of the traversal
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-centers-size") == 0) {
				confirmInitialized(); // confirm that we are initialized
				int size = farthestFirst.traversalSize();
				output += "ff-centers-size: " + size + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-EXTRACT-NEXT
			// - extract the next point of the traversal
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-extract-next") == 0) {
				confirmInitialized(); // confirm that we are initialized
				Airport result = farthestFirst.extractNext();
				if (result == null) {
					output += "ff-extract-next: [No more points remaining]" + System.lineSeparator();
				} else {
					output += "ff-extract-next: " + result.getLabel() + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// FF-LIST-TRAVERSAL
			// - list the points of the traversal
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-list-traversal") == 0) {
				confirmInitialized(); // confirm that we are initialized
				ArrayList<Airport> traversal = farthestFirst.getTraversal();
				if (traversal == null)
					throw new Exception("Error - farthestFirst.getTraversal returned a null result");
				output += summarizeTraversal(traversal); // summarize the traversal
			}
			// -----------------------------------------------------
			// CA-LIST-CENTERS
			// - list all the centers
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-list-centers") == 0) {
				ArrayList<String> list = clusterAssignment.listCenters();
				if (list == null) {
					throw new Exception("clusterAssignment.listCenters returned a null result");
				}
				output += "ca-list-centers:" + System.lineSeparator();
				for (String s : list) {
					output += "  " + s + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// CA-LIST-ASSIGNMENTS
			// - list all the assignments
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-list-assignments") == 0) {
				ArrayList<String> list = clusterAssignment.listAssignments();
				if (list == null) {
					throw new Exception("clusterAssignment.listAssignments returned a null result");
				}
				output += "ca-list-assignments:" + System.lineSeparator();
				for (String s : list) {
					output += "  " + s + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// FF-LIST-CENTERS
			// - list all the centers
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-list-centers") == 0) {
				ArrayList<String> list = farthestFirst.listCenters();
				if (list == null) {
					throw new Exception("farthestFirst.listCenters returned a null result");
				}
				output += "ff-list-centers:" + System.lineSeparator();
				for (String s : list) {
					output += "  " + s + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// FF-LIST-ASSIGNMENTS
			// - list all the assignments
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-list-assignments") == 0) {
				ArrayList<String> list = farthestFirst.listAssignments();
				if (list == null) {
					throw new Exception("farthestFirst.listAssignments returned a null result");
				}
				output += "ff-list-assignments:" + System.lineSeparator();
				for (String s : list) {
					output += "  " + s + System.lineSeparator();
				}
			}
			// -----------------------------------------------------
			// CA-CLEAR
			// -----------------------------------------------------
			else if (cmd.compareTo("ca-clear") == 0) {
				confirmInitialized(); // confirm that we are initialized
				clusterAssignment.clear(); // clear the cluster assignment structure
				airports.clear(); // clear the airports map
				output += "ca-clear: successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// FF-CLEAR
			// -----------------------------------------------------
			else if (cmd.compareTo("ff-clear") == 0) {
				confirmInitialized(); // confirm that we are initialized
				farthestFirst.clear(); // clear the farthest-first structure
				airports.clear(); // clear the airports map
				output += "ff-clear: successful" + System.lineSeparator();
			}
			// -----------------------------------------------------
			// Invalid command or empty
			// -----------------------------------------------------
			else {
				if (cmd.compareTo("") == 0)
					System.err.println("Error: Empty command line (Ignored)");
				else
					System.err.println("Error: Invalid command - \"" + cmd + "\" (Ignored)");
			}
			line.close();
		} 
		// -----------------------------------------------------
		// Error processing
		// -----------------------------------------------------
		catch (Exception e) { // exception thrown?
			if (e.getMessage() == null) {
				output += "Failure due to unexpected exception (probably runtime error)" + System.lineSeparator();
			} else {
				output += "Failure due to exception: \"" + e.getMessage() + "\"" + System.lineSeparator();
			}
			e.printStackTrace(System.err); // print stack trace
		} catch (Error e) { // error occurred?
			System.err.print("Operation failed due to error: " + e.getMessage());
			e.printStackTrace(System.err);
		} finally { // always executed
			line.close(); // close the input scanner
		}
		return output; // return summary output
	}
	
	/**
	 * Read a new airport from the input scanner, check for duplicates, and
	 * add it to the map.
	 */
	Airport readAirport(Scanner line) throws Exception {
		String code = line.next(); // get parameters
		String city = line.next();
		double x = line.nextDouble();
		double y = line.nextDouble();
		Airport ap = new Airport(code, city, x, y); // create airport object
		Airport ap2 = airports.get(code);
		if (ap2 != null) { // code already exists?
			throw new Exception("Insertion of duplicate airport code");
		} else {
			airports.put(ap.getCode(), ap); // add to dictionary
		}
		return ap;
	}


	/**
	 * Confirm that the data structure has been initialized, or throw an exception.
	 */
	void confirmInitialized() throws Exception {
		if (!initialized) {
			throw new Exception("Error: First command must be 'initialize'.");
		}
	}

	/**
	 * Summarize the results of a find command.
	 */
	static String summarizeFind(Point2D q, Airport result) {
		String output = new String("find" + q + ": ");
		if (result != null) {
			output += result.getLabel() + System.lineSeparator();
		} else {
			output += "[not found]" + System.lineSeparator();
		}
		return output;
	}

	/**
	 * Summarize the results of a search command.
	 */
	static String summarizeNNSearch(Point2D q, Airport result) {
		String output = new String("nearest-neighbor" + q + ": ");
		if (result != null) {
			output += result.getLabel() + " (d=" + q.distanceSq(result.getPoint2D()) + ")" + System.lineSeparator();
		} else {
			output += "[not found]" + System.lineSeparator();
		}
		return output;
	}

	/**
	 * Summarize the results of a traversal.
	 */
	static String summarizeTraversal(ArrayList<Airport> result) {
		final int entriesPerLine = 20;
		final int entriesFirstLine = 0;
		String output = new String("list-traversal:");
		for (int i = 0; i < result.size(); i++) {
			Airport ap = result.get(i);
			if (i % entriesPerLine == entriesFirstLine)
				output += System.lineSeparator() + "   ";
			output += " " + ap.getLabel();
		}
		output += System.lineSeparator();
		return output;
	}

	/**
	 * Print the tree contents with indentation.
	 */
	static String treeStructure(String title, ArrayList<String> entries) {
		String output = title + " [tree structure]:" + System.lineSeparator();
		Iterator<String> iter = entries.iterator(); // iterator for the list
		if (iter.hasNext()) { // tree is nonempty
			output += treeStructureHelper(iter, "  "); // print everything
		}
		return output;
	}

	/**
	 * Recursive helper for treeStructure. The argument iterator specifies the next
	 * node from the preorder list to be printed, and the argument indent indicates
	 * the indentation to be performed (of the form "| | | ...").
	 */
	static String treeStructureHelper(Iterator<String> iter, String indent) {
		final String levelIndent = "| "; // the indentation for each level of the tree
		String output = "";
		if (iter.hasNext()) {
			String entry = iter.next(); // get the next entry
			Boolean isExtern = (entry.length() > 0 && entry.charAt(0) == '['); // external?
			if (isExtern) { // print external node entry
				output += indent + entry + System.lineSeparator();
			} else {
				output += treeStructureHelper(iter, indent + levelIndent); // print right subtree
				output += indent + entry + System.lineSeparator(); // print this node
				output += treeStructureHelper(iter, indent + levelIndent); // print left subtree
			}
		} else {
			System.err.println("Unexpected trailing elements in entries list"); // shouldn't get here!
		}
		return output;
	}
}