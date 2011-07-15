#ifndef genomic_Tree_h
#define genomic_Tree_h

#include <iostream>
#include <stdexcept>

using namespace std;

namespace tree {
	
	template <typename T>
	struct Node
	{
		T value;
		Node<T> *left, *right;
		bool leaf;
		Node() {}
		Node(T nodeValue)
		: value(nodeValue), leaf(true), left(NULL), right(NULL) {}
		Node(T nodeValue, Node* leftChild, Node* rightChild)
		: value(nodeValue), left(leftChild), right(rightChild), leaf(false) {}
	};
	
	// Tree derived classes have same interface, different mechanism
	// Use Public Function Calls Virtual Protected idiom
	template <typename T, typename Key>
	class Tree
	{
	protected:
		Node<T>* root;
		size_t count;
		virtual bool _find(Node<T>* subroot, const Key& key, T& item) const = 0;
		virtual Node<T>* _insert(Node<T>* subroot, const T& item) = 0;
		virtual Node<T>* _remove(Node<T>* subroot , const Key& key, Node<T>*& item) = 0;
		virtual void _print(Node<T>* subroot, int level) const = 0;
	private:
		void __clear(Node<T>* subroot) {
			if (subroot == NULL) return;
			__clear(subroot->left);
			__clear(subroot->right);
			delete subroot;
		}
	public:
		Tree() : root(NULL), count(0) {}
		virtual ~Tree() {
			clear();
		}
		void clear() {
			__clear(root);
			root = NULL;
			count = 0;
		}
		void insert(const T& item) {
			root = _insert(root, item);
			++count;
		}
		bool remove(const Key& key, T& item) {
			Node<T>* t = NULL;
			root = _remove(root, key, t);
			// key is not found
			if (t == NULL) return false;
			item = t->value();
			--count;
			delete t;
			return true;
		}
		bool removeAny(T& item) {
			// delete minimum value
			if (root == NULL) return false;
			Node<T>* t;
			root = removeMin(root, t);
			item = t->value();
			delete t;
			--count;
			return true;
		}
		bool find(const Key& key, T& item) const {
			return _find(root, key, item);
		}
		unsigned int size() {
			return count;
		}
		void print() const {
			if (root == NULL) {
				cout << "The tree is empty." << endl;
			} else {
				_print(root, 0);
			}
		}
	};
	
	// Assume T and Key are the same
	template <typename T>
	class BSTree : public Tree<T, T>
	{
	private:
		bool _find(Node<T>* subroot, const T& key, T& item) const {
			if (subroot == NULL) {
				return false;
			} else if (key < subroot->value) {
				// check left
				return _find(subroot->left, key, item);
			} else if (key > subroot->value) {
				// check right
				return _find(subroot->left, key, item);
			} else {
				// found
				item = subroot->value;
				return true;
			}
		}
		Node<T>* _insert(Node<T>* subroot, const T& item) {
			if (subroot == NULL) {
				// empty tree: create node
				return new Node<T>(item);
			}
			if (item < subroot->value) {
				// insert on left
				subroot->left = _insert(subroot->left, item);
			} else {
				// insert on right
				subroot->right = _insert(subroot->right, item);
			}
			// return subtree with node inserted
			return subroot;
		}
		// return the subtree after the node with the specified key has been removed
		// remove designated node by replacing it with the node
		//  with the least key value greater than the one being removed,
		//  i.e. the node with the least key value in the right branch
		//  (alternatively, one can replace the node with the greatest key value
		//  less than the one being removed, i.e. in the left branch; however,
		//  duplicated values would result in an invalid tree)
		Node<T>* _remove(Node<T>* subroot , const T& key, Node<T>*& item) {
			if (subroot == NULL) {
				// item is not in tree
				return NULL;
			} else if (key < subroot->value) {
				// check left
				subroot->left = _remove(subroot->left, key, item);
			} else if (key > subroot->value) {
				subroot->right = _remove(subroot->right, key, item);
			} else {
				// found it: remove it
				Node<T>* tmp;
				item = subroot;
				if (subroot->left == NULL) {
					// only a right child (or no child): point to it
					subroot = subroot->right;
				} else if (subroot->right == NULL) {
					// only a left child: point to it
					subroot = subroot->left;
				} else {
					// Node to be removed has both children
					// Replace current node's value with the node storing the
					//  least key value that is greater than the node;
					//  store it in $tmp (passed by reference)
					subroot->right = removeMin(subroot->right, tmp);
					// swap values
					T value = subroot->value;
					subroot->value = tmp->value;
					tmp->value = value;
					// Now, subroot is the replacing node
					//   and $tmp stores the removed node: pass it by reference
					item = tmp;
				}
			}
			return subroot;
		}
		// print tree with the root at the left, branching to the right
		void _print(Node<T>* subroot, int level) const {
			if (subroot == NULL) return;  // empty tree: nothing to print
			_print(subroot->left, level+1);
			// indent to level
			for (int i = 0; i < level; ++i) {
				cout << "- ";
			}
			cout << subroot->value << endl;
			_print(subroot->right, level+1);
		}
		
		// return the subroot with the min key value node removed
		//  also returns pointer to the removed node by reference
		Node<T>* removeMin(Node<T>* subroot, Node<T>*& min) {
			if (subroot->left == NULL) {
				// node with min key has no left-child
				// (if it did, it would not have the min key)
				min = subroot;
				// past the right child back to the parent,
				//  s.t. the parent can assign the left child to it,
				//  effectively removing $min from the tree
				return subroot->right;
			} else {
				// continue traversing down the left branch to find the min
				subroot->left = removeMin(subroot->left, min);
				return subroot;
			}
		}
	public:
		BSTree() {}
		~BSTree() {}
	};
	
	
	typedef int* Coordinates;
	
	inline bool equalCoordinates(Coordinates a, Coordinates b, size_t ndim) {
		bool equal = true;
		for (size_t i = 0; i < ndim; ++i) {
			if (a[i] != b[i]) {
				equal = false;
				break;
			}
		}
		return equal;
	}
	
	// T must support T.keys()
	// the size of T.keys() must match ndim
	template <typename T>
	class KDTree : public Tree<T, Coordinates>
	{
	private:
		size_t ndim;
		bool _find(Node<T>* subroot, const Coordinates& keys, T& item) const {
			return this->_find(subroot, keys, item, 0);
		}
		bool _find(Node<T>* subroot, const Coordinates& keys, T& item, int discrim) const {
			if (subroot == NULL) return false;  // empty tree
			Coordinates currKeys = subroot->value.keys();
			if (equalCoordinates(currKeys, keys, ndim)) {
				item = subroot->value;
				return true;
			}
			if (keys[discrim] < currKeys[discrim]) {
				// search left branch, incrementing the discriminant dimension
				return _find(subroot->left, keys, item, (discrim+1)%ndim);
			} else {
				// search righ branch
				return _find(subroot->right, keys, item, (discrim+1)%ndim);
			}
		}
		Node<T>* _insert(Node<T>* subroot, const T& item) {
			return this->_insert(subroot, item, 0);
		}
		Node<T>* _insert(Node<T>* subroot, const T& item, int discrim) {
			if (subroot == NULL) {
				// empty tree: create node
				return new Node<T>(item);
			}
			Coordinates currKeys = subroot->value.keys();
			Coordinates keys = item.keys();
			if (keys[discrim] < currKeys[discrim]) {
				// insert on left
				subroot->left = _insert(subroot->left, item, (discrim+1)%ndim);
			} else {
				// insert on right
				subroot->right = _insert(subroot->right, item, (discrim+1)%ndim);
			}
			// return subtree with node inserted
			return subroot;
		}
		void _print(Node<T>* subroot, int level) const {
			if (subroot == NULL) return;  // empty tree: nothing to print
			_print(subroot->left, level+1);
			// indent to level
			for (int i = 0; i < level; ++i) {
				cout << "- ";
			}
			// print tuple
			cout << "/";
			for (int i = 0; i < ndim; ++i) {
				cout << subroot->value.keys()[i] << "/";
			}
			cout << endl;
			_print(subroot->right, level+1);
		}
		
		Node<T>* _remove(Node<T>* subroot, const Coordinates& keys, Node<T>*& item) {
			return _remove(subroot, keys, item, 0);
		}
		// Remving node N
		// If N has one child...
		//   We cannot simply assign N's parent to point to N's child, since
		//   the lefts would be changed, create an invalid KD tree.
		//   Therefore, replace N by an appropriate node in the left or right tree.
		// If N has only a right child (or both children)...
		//   First find node N, then replace its record by the record in N's right
		//   subtree with the least value of N's discriminator
		// If the right subtree does not exist, it is not satisfactory to
		//   replace N's record with the recording having the greatest value for
		//   the discriminator in the left subtree: this value might appear more than
		//   once in the left subtree.
		//   This would lead to the equal values for the discriminator in N's left tree.
		//   (Recall that nodes in the left tree must have lower values of the
		//   discriminator.)
		// Solution: first move the left subtree of node N to become the right
		//   subtree. Then, process with normal deletion process, replacing the
		//   the record of N to be deleted with the record containing the
		//   *least* value of the discriminator from what is now N's right subtree
		Node<T>* _remove(Node<T>* subroot, const Coordinates& keys, Node<T>*& item, int discrim) {
			
			if (subroot == NULL) {
				// item is not in tree
				return NULL;
			}
			
			Coordinates currKeys = subroot->value.keys();
			if (equalCoordinates(currKeys, keys, ndim)) {
				// Found N == subroot: proceed to remove it
				
				if (subroot->right == NULL && subroot->left == NULL) {
					// N has no children: N can simply be removed
					// Pass N by reference back to caller
					item = subroot;
					// Set subroot to NULL, s.t. calling function can effectively 
					//   remove N from tree while setting its pointer to subroot
					subroot = NULL;
				} else {
					// N has at least one child
				
					// Always replace N's record by the record in N's right subtree with
					//   the least value of N's discriminator
					
					// One special case: right subtree does not exist
					//   Move to the left subtree to the right subtree
					if (subroot->right == NULL) {
						// only a left child: point to it
						subroot->right = subroot->left;
						subroot->left = NULL;
					}
					// Find the node with the minimum key of the N's discriminator,
					//   in the right subtree
					//   $discrim == N's discrminator, since subroot == N
					Node<T>* min = findMin(subroot->right, discrim, discrim);
					// At this point, N must have a right child, because if it only
					//   had a left child, it has been assigned as the right child
					// Since N has a right child, min cannot be NULL
					// After finding min, remove it from the right subtree
					subroot->right = _remove(subroot->right, min->value.keys(), min, discrim);
					// swap values between N and the min Node
					T value = subroot->value;
					subroot->value = min->value;
					min->value = value;
					// Now, min holds of the record of N
					item = min;
				}
				
			} else if (keys[discrim] < currKeys[discrim]) {
				// search left branch, incrementing the discriminant dimension
				// assign left child to the residual left tree with the target node removed
				subroot->left = _remove(subroot->left, keys, item, (discrim+1)%ndim);
			} else {
				// search righ branch
				// assign right child to the residual right tree with the target node removed
				subroot->right = _remove(subroot->right, keys, item, (discrim+1)%ndim);
			}
			
			return subroot;
		}
		// Finds and returns the node with the mininum key value of the
		//   specified discriminator
		// Auxilary function for removing a node
		Node<T>* findMin(Node<T>* subroot, int discrim, int currDiscrim) {
			if (subroot == NULL) return NULL;
			Node<T> *a, *b;
			Coordinates keys, akeys, bkeys;
			keys = subroot->value.keys();
			a = findMin(subroot->left, discrim, (currDiscrim+1)%ndim);
			if (a != NULL) akeys = a->value.keys();
			if (discrim != currDiscrim) {
				// discriminator for the current level is different from the
				//   the specified discrminator
				// Thus, min could be on either side: check the right side too
				// Otherwise, having checked the left side sufficient
				b = findMin(subroot->right, discrim, (currDiscrim+1)%ndim);
				if (b != NULL) bkeys = b->value.keys();
				if (a == NULL || (b != NULL && bkeys[discrim] < akeys[discrim])) {
					// Right side has a smaller key value
					a = b;
					akeys = bkeys;
				}
			}
			// Now, a has the smallest value in children
			if (a == NULL || (keys[discrim] < akeys[discrim])) {
				return subroot;
			}
			return a;
		}
	public:
		KDTree(size_t nDimensions) : ndim(nDimensions) {}
		~KDTree() {}
	};
	
};


#endif