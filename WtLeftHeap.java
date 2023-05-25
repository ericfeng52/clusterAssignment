package cmsc420_s23;

import java.util.ArrayList;

// ---------------------------------------------------------------------
// Author: Dave Mount
// For: CMSC 420
// Date: Spring 2023
//
// This is a weight-balanced implementation of a (max) leftist heap data 
// structure. This is a priority queue, which in addition to supporting 
// insertion and extract-max, it supports merging of two heaps. The
// twist is that rather than using the null path length, the subtree
// with the fewer nodes is on the right side.
// ---------------------------------------------------------------------

public class WtLeftHeap<Key extends Comparable<Key>, Value> {

	// -----------------------------------------------------------------
	// Node of the heap
	// -----------------------------------------------------------------

	private class Node {
		Key key; // key (priority)
		Value value; // value (application dependent)
		Node left; // children
		Node right;
		Node parent; // parent
		int weight; // subtree weight
		Locator loc; // this node's locator

		/**
		 * Basic constructor.
		 */
		Node(Key x, Value v) {
			this.key = x;
			this.value = v;
			this.left = null;
			this.right = null;
			this.parent = null;
			this.weight = 1;
			this.loc = null;
		}
	}

	// -----------------------------------------------------------------
	// Locator - Used to locate a previously inserted item
	// -----------------------------------------------------------------

	public class Locator {
		private Node node;

		private Locator(Node node) { // basic constructor
			this.node = node;
			node.loc = this;
		}

		private Node get() { // get the associated node
			return node;
		}
		
		public Value getValue() { // just used for debugging
			return node.value;
		}
	}
	
	// -----------------------------------------------------------------
	// Private members
	// -----------------------------------------------------------------

	private Node root; // root of tree
	private int size; // number of entries

	// -----------------------------------------------------------------
	// Local utilities
	// -----------------------------------------------------------------

	/**
	 * Access a node's weight (allowing for null pointers).
	 */
	int getWeight(Node u) {
		return (u == null ? 0 : u.weight);
	}

	/**
	 * Merge helper. Merge subtrees rooted at u and v. We do
	 * not update the parent link of the new root node. This
	 * is the responsibility of the calling function.
	 */
	Node merge(Node u, Node v) {
		if (u == null)
			return v; // if one is empty, return the other
		if (v == null)
			return u;
		if (u.key.compareTo(v.key) < 0) { // swap so that u is larger
			Node t = u;
			u = v;
			v = t;
		}
		if (u.left == null) { // u must be a leaf
			u.left = v; // put v on its left
			v.parent = u; // update v's parent
		} else { // merge v on right and swap if needed
			u.right = merge(u.right, v); // recursively merge u's right subtree
			u.right.parent = u; // update right's parent
			if (u.left.weight < u.right.weight) { // not leftist?
				Node t = u.left; // swap children
				u.left = u.right;
				u.right = t;
			}
		}
		u.weight = 1 + getWeight(u.left) + getWeight(u.right); // update weight
		return u; // return the root
	}
	
	/**
	 * Swap node contents. This swaps key, value, and updates the
	 * locators.
	 */
	void swapContents(Node u, Node v) {
		Key tempKey = v.key;
		v.key = u.key;
		u.key = tempKey;
		Value tempValue = v.value;
		v.value = u.value;
		u.value = tempValue;
		Locator tempLoc = v.loc;
		v.loc = u.loc;
		u.loc = tempLoc;
		v.loc.node = v;
		u.loc.node = u;
	}

	/**
	 * Sifts a node up to the proper heap position.
	 */
	Node siftUp(Node u) {
		Node p = u.parent;
		if (p != null && u.key.compareTo(p.key) > 0) { // out of order?
			swapContents(u, p); // swap node contents
			return siftUp(p); // repeat with parent
		} else {
			return u; // done
		}
	}

	/**
	 * Sifts a node down to the proper heap position.
	 */
	Node siftDown(Node u) {
		Node ul = u.left;
		Node ur = u.right;
		if (ul == null || ur != null && ur.key.compareTo(ul.key) > 0) {
			ul = ur; // make ul the larger non-null child
		}
		if (ul != null && u.key.compareTo(ul.key) < 0) { // out of order?
			swapContents(u, ul); // swap node contents
			return siftDown(ul); // repeat
		} else {
			return u; // done
		}
	}

	/**
	 * List the nodes of subtree in reverse preorder. By "reverse" we mean that the
	 * right subtree comes before the left.
	 */
	ArrayList<String> reversePreorder(Node u) {
		ArrayList<String> list = new ArrayList<String>();
		if (u == null) {
			list.add("[]"); // null link indicator
		} else {
			list.add("(" + u.key + ", " + u.value + ") [" + u.weight + "]"); // add this node
			list.addAll(reversePreorder(u.right)); // process right
			list.addAll(reversePreorder(u.left)); // process left
		}
		return list;
	}
	
	/**
	 * Check that a locator is correct.
	 */
	void validateLocator(Locator loc) throws Exception {
		if (loc == null) 
			throw new Exception("Invalid locator");
		Node u = loc.node;
		if (u == null) 
			throw new Exception("Invalid locator");
		while (u.parent != null) {
			u = u.parent;
		}
		if (u != root) 
			throw new Exception("Invalid locator");
	}


	// -----------------------------------------------------------------
	// Public members
	// -----------------------------------------------------------------

	/**
	 * Clear the entire structure.
	 */
	public void clear() {
		root = null;
		size = 0;
	}

	/**
	 * Constructor.
	 */
	public WtLeftHeap() {
		clear();
	}

	/**
	 * Insert key-value pair.
	 */
	public Locator insert(Key x, Value v) {
		Node u = new Node(x, v); // create single node
		Locator loc = new Locator(u);
		root = merge(root, u); // merge it with the root
		size += 1;
		return loc;
	}

	/**
	 * Merge with another heap. This clears the other heap.
	 */
	public void mergeWith(WtLeftHeap<Key, Value> h2) {
		if (h2 == null || h2 == this) return; // trivial - ignore
		root = merge(root, h2.root); // merge the two trees
		root.parent = null;
		size = size + h2.size;
		h2.clear();
	}

	/**
	 * Extract the maximum item from the heap and return its value.
	 */
	public Value extract() throws Exception {
		if (size == 0) {
			throw new Exception("Extract from empty heap");
		}
		Value result = root.value; // final result
		root.loc.node = null; // wreck locator (for safety)
		root = merge(root.left, root.right); // merge subtrees
		if (root != null) root.parent = null;
		size -= 1;
		return result;
	}

	/**
	 * Update key. This modifies the key value of a node (given by its
	 * locator). If the new value is strictly smaller, it sifts the 
	 * new key down. If larger, it sifts it up.
	 */
	public void updateKey(Locator loc, Key x) throws Exception {
		validateLocator(loc); // check that locator is valid
		Node u = loc.get();
		if (x.compareTo(u.key) < 0) { // decrease?
			u.key = x; // modify and sift down
			u = siftDown(u);
		} else if (x.compareTo(u.key) > 0) { // increase?
			u.key = x; // modify and sift up
			u = siftUp(u);
		}
		// if (DEBUG) debugPrint("After update key:");
	}

	/**
	 * Return the maximum key (without extraction)
	 */
	public Key peekKey() {
		return (root == null ? null : root.key);
	}

	/**
	 * Return the value associated with the maximum key (without extraction).
	 */
	public Value peekValue() {
		return (root == null ? null : root.value);
	}

	/**
	 * Size of the heap.
	 */
	public int size() {
		return size;
	}

	/**
	 * Get a list of entries in reverse (right-left) preorder.
	 */
	public ArrayList<String> list() {
		ArrayList<String> list = new ArrayList<String>();
		list.addAll(reversePreorder(root));
		return list;
	}
}
